import pandas as pd
import argparse
import os
import sys
import glob
import json
import boto3
from botocore.exceptions import NoCredentialsError


class Config:
    def __init__(self, config_dict):
        for key, value in config_dict.items():
            setattr(self, key, value)


def build_parser():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "--s3_bucket", "-s", help="S3 bucket for output.", default="macchiato.clue.io"
    )
    parser.add_argument(
        "--build_path", "-b", help="Path to the build directory.", required=True
    )
    return parser


def read_config(config_path):
    with open(config_path, "r") as f:
        config_dict = json.load(f)
    return Config(config_dict)


def read_build_file(search_pattern, args):
    fstr = os.path.join(args.build_path, search_pattern)
    fmatch = glob.glob(fstr)
    assert len(fmatch) == 1, "Too many files found: {}".format(fmatch)
    return pd.read_csv(fmatch[0])


def generate_key(df):
    # Filter out the positive and vehicle controls
    df = df[~df["pert_type"].isin(["trt_poscon", "ctl_vehicle"])]

    # Select the columns we need
    df = df[["pert_name", "pert_id", "pert_plate", "x_project_id"]]

    # Drop duplicates
    distinct_df = df.drop_duplicates().reset_index(drop=True)
    return distinct_df


# Turn the key_df into a json object
def key_to_json(
    df: pd.DataFrame, unique_by: list = None, output_path: str = None
) -> str:
    """
    Convert a pandas DataFrame to a pretty-printed JSON string (records orientation) and write it to a file.

    Parameters
    ----------
    df : pd.DataFrame
        The input DataFrame containing columns like "pert_name", "pert_id", "pert_plate",
        "x_project_id", "pert_dose", etc.
    unique_by : list, optional
        If provided, duplicate rows based on these columns will be dropped before conversion.
    output_path : str, optional
        If provided, the resulting JSON string will be written to this file. Necessary directories
        will be created if they do not exist.

    Returns
    -------
    str
        A pretty-printed JSON string representing the DataFrame.
    """
    if unique_by is not None:
        df = df.drop_duplicates(subset=unique_by)

    json_str = df.to_json(orient="records", indent=2)

    if output_path:
        # Create the directory if it doesn't exist
        directory = os.path.dirname(output_path)
        if directory and not os.path.exists(directory):
            os.makedirs(directory, exist_ok=True)
        with open(output_path, "w") as f:
            f.write(json_str)

    return json_str


# Sync all contents of current directory to S3
def sync_to_s3(local_dir, s3_bucket, s3_prefix, exclude_pattern=None):
    """Sync the local directory to S3."""
    s3 = boto3.client("s3")

    for root, dirs, files in os.walk(local_dir):
        for file in files:
            if exclude_pattern not in file:
                local_path = os.path.join(root, file)
                relative_path = os.path.relpath(local_path, local_dir)
                s3_path = os.path.join(s3_prefix, relative_path)

                try:
                    s3.upload_file(local_path, s3_bucket, s3_path)
                    print(f"Uploaded {local_path} to s3://{s3_bucket}/{s3_path}")
                except NoCredentialsError:
                    print("AWS credentials not found.")
                    raise
            else:
                print(f"Excluding upload of {file}")


def main(args):
    # Get the build path
    build_path = args.build_path

    # Read the config file
    config_path = args.build_path + "/config.json"
    config = read_config(config_path)

    # Get build name from config
    build_name = config.BUILD_NAME

    # Construct the s3 prefix
    s3_prefix = f"builds/{build_name}/"

    # Read the sample_meta file
    sample_meta_path = os.path.join(build_path, config.SAMPLE_META)
    sample_meta = pd.read_csv(sample_meta_path)

    # Generate the key_df
    key_df = generate_key(sample_meta)

    # Convert the key_df to json and write to the build directory
    key_to_json(key_df, output_path=os.path.join(build_path, "compound_key.json"))

    # Sync the build directory to S3
    sync_to_s3(
        local_dir=build_path,
        s3_bucket=args.s3_bucket,
        s3_prefix=s3_prefix,
        exclude_pattern="raw_counts",
    )

    print("Completed sync to S3.")


if __name__ == "__main__":
    args = build_parser().parse_args(sys.argv[1:])

    main(args)
