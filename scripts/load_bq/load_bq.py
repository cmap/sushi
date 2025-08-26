# load_bq.py
import os
import logging
import argparse
from load_bq_functions import (
    get_bigquery_client,
    table_exists,
    load_csv_to_bigquery,
    init_logger,
)

# Define your dataset name here or pass as an arg if needed
BQ_DATASET_ID = "sushi"


def build_parser():
    parser = argparse.ArgumentParser(
        description="Load matching CSV files into BigQuery tables."
    )
    parser.add_argument(
        "--build_path", "-b", help="Path to the build directory.", required=True
    )
    parser.add_argument(
        "--build_name",
        "-n",
        help="Name of the build to delete rows for before loading.",
        required=True,
    )
    parser.add_argument("--screen", "-s", help="Screen name.", required=True)
    return parser


def main(args):
    init_logger()
    client = get_bigquery_client()

    logging.info(f"Recursively scanning: {args.build_path}")
    for root, _, files in os.walk(args.build_path):
        for fname in files:
            if not fname.endswith(".csv"):
                continue

            table_name = f"sushi_{os.path.splitext(fname)[0]}"
            file_path = os.path.join(root, fname)

            if table_exists(client, BQ_DATASET_ID, table_name):
                logging.info(f"Found table match for: {fname} → {table_name}")
                load_csv_to_bigquery(
                    client,
                    BQ_DATASET_ID,
                    table_name,
                    file_path,
                    build_name=args.build_name,
                    screen=args.screen,
                )
            else:
                logging.warning(f"No matching table for: {fname} → {table_name}")


if __name__ == "__main__":
    parser = build_parser()
    args = parser.parse_args()
    main(args)
