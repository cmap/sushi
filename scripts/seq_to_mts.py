import pandas as pd
import argparse
import os
import sys
import glob
import logging

logger = logging.getLogger('seq_to_mts')
pert_vehicle = "DMSO"
pert_time_unit = "d"

def build_parser():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--build_path', '-b', help='Build path with SUSHI level 3, 4, and 5 data.', required=True)
    parser.add_argument('--out', '-o', help='Output for project level folders', required=True)
    parser.add_argument("--verbose", '-v', help="Whether to print a bunch of output", action="store_true", default=False)
    parser.add_argument('--build_name', '-n', help='Build name.', required=True)
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

def main(args):
    if os.path.isdir(args.out):
        pass
    else:
        os.makedirs(args.out)

    try:
        fstr = os.path.join(args.build_path, '*l2fc*.csv')
        fmatch = glob.glob(fstr)
        assert (len(fmatch) == 1) , "Too many files found"
        print("Reading in data")
        sample_meta = read_build_file("*sample_meta*.csv", args)
        level_3 = read_build_file("*normalized_counts*.csv", args)
        level_4 = read_build_file("*l2fc*.csv", args)
        level_5 = read_build_file("*collapsed_values*.csv", args)

    except IndexError as err:
        logger.error(err)
        logger.error("Index Error: No file found. Check --build_path arg")
        raise


    # Define the column renaming dictionary
    column_mapping = {
        "project_code": "screen",
        "DepMap_ID": "depmap_id",
        "CCLE_name": "ccle_name",
        "prism_cell_set": "culture",
        "trt_type": "pert_type",
        "sig_id": "profile_id",
        "bio_rep": "replicate",
        "day": "pert_time",
        "treatment": "pert_iname",
        "dose": "pert_dose",
        "dose_unit": "pert_dose_unit",
        "l2fc": "LFC",
        "median_l2fc": "LFC"  # Add the mapping for level_5 dataset
    }


    # Define the list of datasets
    datasets = [sample_meta, level_3, level_4, level_5]

    # Check for existence of "x_project_id" and "pert_plate" before renaming columns
    for ind, dataset in enumerate(datasets):
        dataset_name = ["sample_meta", "level_3", "level_4", "level_5"][ind]
        missing_columns = [col for col in ["x_project_id", "pert_plate"] if col not in dataset.columns]
        if missing_columns:
            missing_cols_str = ", ".join(missing_columns)
            raise ValueError(f"Columns '{missing_cols_str}' not found in the '{dataset_name}' dataset. Cannot proceed.")
        else:
            for old_name, new_name in column_mapping.items():
                if old_name in dataset.columns:
                    dataset.rename(columns={old_name: new_name}, inplace=True)

        # Seq projects are all PR500 for now
        dataset["culture"] = "PR500"

    # Define the pert_time values to drop
    pert_time_to_drop = [0, 6]

    for dataset in datasets:
        dataset.drop(dataset[dataset['pert_time'].isin(pert_time_to_drop)].index, inplace=True)


    # Setting columns
    print("Reformatting columns...")
    sample_meta = sample_meta.assign(pert_vehicle=pert_vehicle, pert_time_unit = pert_time_unit, 
                            pert_id = sample_meta["pert_iname"].str.upper(), prc_id = sample_meta["pert_iname"].str.upper())

    level_3 = level_3.assign(pert_vehicle=pert_vehicle, pert_time_unit = pert_time_unit, 
                                pert_id = level_3["pert_iname"].str.upper(), prc_id = level_3["pert_iname"].str.upper())

    level_4 = level_4.assign(pert_vehicle=pert_vehicle, pert_time_unit = pert_time_unit, 
                                pert_id = level_4["pert_iname"].str.upper(), prc_id = level_4["pert_iname"].str.upper())

    level_5 = level_5.assign(pert_vehicle=pert_vehicle, pert_time_unit = pert_time_unit, 
                                pert_id = level_5["pert_iname"].str.upper())

    # Adding itime/time and idose
    sample_meta["pert_itime"] = sample_meta["pert_time"].astype(str) + " " + sample_meta["pert_time_unit"]
    sample_meta["pert_idose"] = sample_meta["pert_dose"].astype(str) + " " + sample_meta["pert_dose_unit"]

    level_3["pert_itime"] = level_3["pert_time"].astype(str) + " " + level_3["pert_time_unit"]
    level_3["pert_idose"] = level_3["pert_dose"].astype(str) + " " + level_3["pert_dose_unit"]

    level_4["pert_itime"] = level_4["pert_time"].astype(str) + " " + level_4["pert_time_unit"]
    level_4["pert_idose"] = level_4["pert_dose"].astype(str) + " " + level_4["pert_dose_unit"]

    level_5["pert_time"] = level_5["pert_time"].astype(str) + " " + level_5["pert_time_unit"]
    level_5["pert_idose"] = level_5["pert_dose"].astype(str) + " " + level_5["pert_dose_unit"]

    # Define a mapping for renaming values
    trt_type_mapping = {
        'poscon': 'trt_poscon',
        'negcon': 'ctl_vehicle'}

    # Replace values in the 'pert_type' column using the mapping
    sample_meta["pert_type"] = sample_meta["pert_type"].replace(trt_type_mapping)
    level_3["pert_type"] = level_3["pert_type"].replace(trt_type_mapping)

    # Sorting columns to resemble MTS style
    level_4.sort_index(axis=1, inplace=True)
    profile_col = level_4.pop('profile_id')
    # pool_col = level_4.pop('pool_id')
    level_4.insert(1, profile_col.name, profile_col)
    # level_4.insert(2, pool_col.name, pool_col)
    level_4 = level_4[[col for col in level_4.columns if col != 'LFC'] + ['LFC']]

    level_5.sort_index(axis=1, inplace=True)
    level_5 = level_5[[col for col in level_5.columns if col != 'LFC'] + ['LFC']]

    # Writing out modified dataframes
    # project = level_4["screen"].unique()[0]
    project = args.build_name

    print("Creating compound key...")

    # Writing out project key
    project_key = write_key(level_3)


    # Saving modified data
    # Add number of cell lines and unique pert_plate/well combinations
    level_3.to_csv(args.out + project + "_inst_info.txt", sep="\t", index=None)
    level_3.to_csv(args.out + project + "_LEVEL3_NORMALIZED_COUNTS.csv", index=0)
    level_4.to_csv(args.out + project + "_LEVEL4_LFC.csv", index=0)
    level_5.to_csv(args.out + "/" + project + "_LEVEL5_LFC.csv", index=0)
    project_key.to_csv(args.out + project + "_compound_key.csv", index=False)
    return level_3, project_key, level_4, level_5


if __name__ == "__main__":
    args = build_parser().parse_args(sys.argv[1:])
    level = (logging.DEBUG if args.verbose else logging.INFO)
    logging.basicConfig(level=level)
    logger.info("args:  {}".format(args))

    main(args)