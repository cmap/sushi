import pandas as pd
import argparse
import os
import sys
import glob
import logging
import json
import gzip
import shutil
import boto3
from botocore.exceptions import NoCredentialsError

logger = logging.getLogger('seq_to_mts')
pert_vehicle = "DMSO"
pert_time_unit = "d"


class Config:
    def __init__(self, config_dict):
        for key, value in config_dict.items():
            setattr(self, key, value)


def build_parser():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--build_path', '-b', help='Build path with SUSHI level 3, 4, and 5 data.', required=True)
    parser.add_argument('--out', '-o', help='Output for project level folders', required=True)
    parser.add_argument("--verbose", '-v', help="Whether to print a bunch of output", action="store_true", default=False)
    parser.add_argument('--build_name', '-n', help='Build name.', required=True)
    parser.add_argument('--days', '-d', help='Day timepoints to drop from output data separated by commas.')
    parser.add_argument('--config', '-c', help='Config file for project.', required=True, default='config.json')
    return parser


def read_build_file(search_pattern, args):
    fstr = os.path.join(args.build_path, search_pattern)
    fmatch = glob.glob(fstr)
    assert (len(fmatch) == 1) , "Too many files found: {}".format(fmatch)
    return pd.read_csv(fmatch[0])


def write_key(df):
    df = df[~df['pert_type'].isin(['trt_poscon', 'ctl_vehicle'])]

    df = df[['pert_iname', 'pert_id', 'pert_plate', 'pert_dose', 'x_project_id']]
    distinct_df = df.drop_duplicates().reset_index(drop=True)  
      
    distinct_df['pert_dose'] = distinct_df['pert_dose'].round(10)
    distinct_df.rename(columns={'pert_dose': 'pert_dose_1'}, inplace=True)
    # df = df.fillna(0)

    exclude_columns = [col for col in distinct_df.columns if "pert_dose_1" in col]
    grouped_df = distinct_df.groupby([col for col in distinct_df.columns if col not in exclude_columns]).nunique().reset_index()

    return grouped_df


def read_config(config_path):
    with open(config_path, 'r') as f:
        config_dict = json.load(f)
    return Config(config_dict)


def create_profile_id_column(df, config):
    # Split the sig_cols string into a list of column names
    sig_cols = config.SIG_COLS.split(',')

    # Create the 'profile_id' by concatenating the specified columns with ':' as delimiter
    df['profile_id'] = df[sig_cols].astype(str).agg(':'.join, axis=1)

    return df


def gzip_file(input_file):
    """Gzip a file in place."""
    with open(input_file, 'rb') as f_in:
        with gzip.open(input_file + '.gz', 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    os.remove(input_file)  # Remove the original file after gzipping


def sync_to_s3(local_dir, s3_bucket, s3_prefix):
    """Sync the local directory to S3."""
    s3 = boto3.client('s3')

    for root, dirs, files in os.walk(local_dir):
        for file in files:
            local_path = os.path.join(root, file)
            relative_path = os.path.relpath(local_path, local_dir)
            s3_path = os.path.join(s3_prefix, relative_path)

            try:
                s3.upload_file(local_path, s3_bucket, s3_path)
                logger.info(f"Uploaded {local_path} to s3://{s3_bucket}/{s3_path}")
            except NoCredentialsError:
                logger.error("AWS credentials not found.")
                raise


def remove_invalid_pert_ids(df):
    """
    Removes rows where 'pert_id' is either NaN or 'NONE'.

    Parameters:
    df (pd.DataFrame): Input dataframe containing a 'pert_id' column.

    Returns:
    pd.DataFrame: Dataframe with invalid 'pert_id' rows removed.
    """
    # Drop rows where 'pert_id' is NaN or 'NONE'
    clean_df = df[~df['pert_id'].isna() & (df['pert_id'] != 'NONE')]

    return clean_df


def main(args):
    if os.path.isdir(args.out):
        pass
    else:
        os.makedirs(args.out)

    try:
        print("Reading in data")
        sample_meta = read_build_file("sample_meta.csv", args)
        level_3 = read_build_file("normalized_counts.csv", args)
        level_4 = read_build_file("l2fc.csv", args)
        level_5 = read_build_file("collapsed_l2fc.csv", args)

    except IndexError as err:
        logger.error(err)
        logger.error("Index Error: No file found. Check --build_path arg")
        raise


    # Define the column renaming dictionary
    column_mapping = {
        "project_code": "screen",
        "prism_cell_set": "culture",
        "sig_id": "profile_id",
        "bio_rep": "replicate",
        "day": "pert_time",
        "pert_name": "pert_iname",
        "l2fc": "LFC",
        "median_l2fc": "LFC"
    }

    # Get the columns for profile_id from SIG_COLS
    config_path = args.build_path + '/' + args.config
    config = read_config(config_path=config_path)

    # Add profile_id to levels 3,4 and 5
    level_3 = create_profile_id_column(level_3, config)
    level_4 = create_profile_id_column(level_4, config)

    # Define the list of datasets
    datasets = [sample_meta, level_3, level_4, level_5]

    # Check for existence of "x_project_id" and "pert_plate" before renaming columns
    for ind, dataset in enumerate(datasets):
        dataset_name = ["sample_meta", "level_3", "level_4", "level_5"][ind]
        missing_columns = [col for col in ["x_project_id", "pert_plate"] if col not in dataset.columns]
        print(f'Renaming columns for {dataset_name}...')
        print(f"Original columns {dataset.columns}")
        if missing_columns:
            missing_cols_str = ", ".join(missing_columns)
            raise ValueError(f"Columns '{missing_cols_str}' not found in the '{dataset_name}' dataset. Cannot proceed.")
        else:
            for old_name, new_name in column_mapping.items():
                if old_name in dataset.columns:
                    dataset.rename(columns={old_name: new_name}, inplace=True)
        print(f"Renamed columns {dataset.columns}")

        # Seq projects are all PR500 for now
        dataset["culture"] = "PR500"

    # Define the pert_time values to drop
    if args.days:
        pert_time_to_drop = [int(day) for day in args.days.split(",")]
        print(pert_time_to_drop)
        for dataset in datasets:
            dataset.drop(dataset[dataset['pert_time'].isin(pert_time_to_drop)].index, inplace=True)


    # Setting columns
    print("Reformatting columns...")
    sample_meta = sample_meta.assign(pert_vehicle=pert_vehicle, pert_time_unit = pert_time_unit)

    level_3 = level_3.assign(pert_vehicle=pert_vehicle, pert_time_unit = pert_time_unit)

    level_4 = level_4.assign(pert_vehicle=pert_vehicle, pert_time_unit = pert_time_unit)

    level_5 = level_5.assign(pert_vehicle=pert_vehicle, pert_time_unit = pert_time_unit)

    # Remove invalid 'pert_id' rows, currently 'NONE' and NaN
    sample_meta = remove_invalid_pert_ids(sample_meta)
    level_3 = remove_invalid_pert_ids(level_3)
    level_4 = remove_invalid_pert_ids(level_4)
    level_5 = remove_invalid_pert_ids(level_5)

    # Adding itime/time and idose
    sample_meta["pert_itime"] = sample_meta["pert_time"].astype(str) + " " + sample_meta["pert_time_unit"]
    sample_meta["pert_idose"] = sample_meta["pert_dose"].astype(str) + " " + sample_meta["pert_dose_unit"]

    level_3["pert_itime"] = level_3["pert_time"].astype(str) + " " + level_3["pert_time_unit"]
    level_3["pert_idose"] = level_3["pert_dose"].astype(str) + " " + level_3["pert_dose_unit"]

    level_4["pert_itime"] = level_4["pert_time"].astype(str) + " " + level_4["pert_time_unit"]
    level_4["pert_idose"] = level_4["pert_dose"].astype(str) + " " + level_4["pert_dose_unit"]

    level_5["pert_time"] = level_5["pert_time"].astype(str) + " " + level_5["pert_time_unit"]
    level_5["pert_idose"] = level_5["pert_dose"].astype(str) + " " + level_5["pert_dose_unit"]

    # Sorting columns to resemble MTS style
    level_4.sort_index(axis=1, inplace=True)
    profile_col = level_4.pop('profile_id')
    level_4.insert(1, profile_col.name, profile_col)
    level_4 = level_4[[col for col in level_4.columns if col != 'LFC'] + ['LFC']]

    level_5.sort_index(axis=1, inplace=True)
    level_5 = level_5[[col for col in level_5.columns if col != 'LFC'] + ['LFC']]

    # Writing out modified dataframes
    project = args.build_name

    print("Creating compound key...")

    # Writing out project key
    project_key = write_key(level_3)

    # Saving modified data
    output_files = {
        project + "_inst_info.txt": level_3,
        project + "_LEVEL3_NORMALIZED_COUNTS.csv": level_3,
        project + "_LEVEL4_LFC.csv": level_4,
        project + "_LEVEL5_LFC.csv": level_5,
        project + "_compound_key.csv": project_key,
    }

    for file_name, df in output_files.items():
        output_path = os.path.join(args.out, file_name)
        if file_name.endswith('.csv'):
            df.to_csv(output_path, index=False)
        elif file_name.endswith('.txt'):
            df.to_csv(output_path, sep='\t', index=False)
        gzip_file(output_path)

    # Check for *_EPS_QC_TABLE.csv file and gzip if it exists
    eps_qc_file_pattern = os.path.join(args.build_path, "*_EPS_QC_TABLE.csv")
    eps_qc_files = glob.glob(eps_qc_file_pattern)
    if eps_qc_files:
        for eps_qc_file in eps_qc_files:
            output_eps_qc_file = os.path.join(args.out, os.path.basename(eps_qc_file))
            shutil.copy(eps_qc_file, output_eps_qc_file)
            gzip_file(output_eps_qc_file)

    # Sync the directory to S3
    s3_bucket = "macchiato.clue.io"
    s3_prefix = f"builds/{args.build_name}/build/"
    sync_to_s3(args.out, s3_bucket, s3_prefix)

if __name__ == "__main__":
    args = build_parser().parse_args(sys.argv[1:])
    level = (logging.DEBUG if args.verbose else logging.INFO)
    logging.basicConfig(level=level)
    logger.info("args:  {}".format(args))

    main(args)
