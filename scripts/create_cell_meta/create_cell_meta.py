import requests
import pandas as pd
import os
import urllib.parse
import json
import logging
from typing import Optional, Dict, Any
import argparse
import polars as pl


def fetch_metadata(
    url: str, api_key: str, filter_dict: Optional[Dict[str, Any]] = None
):
    headers = {
        "Accept": "application/json",
        "user_key": api_key,
        "prism_key": "prism_mts",
    }

    params = {}
    if filter_dict:
        # Let the `requests` library handle the encoding.
        # It's safer and cleaner.
        params["filter"] = json.dumps(filter_dict)

    response = requests.get(url, headers=headers, params=params)

    if response.status_code == 200:
        return pl.DataFrame(response.json())
    else:
        logging.error(f"Error: {response.status_code}")
        logging.error(response.text)
        response.raise_for_status()


def save_dataframe(df, filename, build_dir):
    # Ensure the directory exists
    os.makedirs(build_dir, exist_ok=True)

    # Define the file path
    filepath = os.path.join(build_dir, f"{filename}.csv")

    # Write the DataFrame to a CSV file
    df.write_csv(filepath, index=False)
    logging.info(f"Wrote metadata to {filepath}")


def main():
    parser = argparse.ArgumentParser(
        description="Fetch data from API and save to a CSV file."
    )
    parser.add_argument(
        "--screen", "-s", type=str, required=True, help="Screen code to query the API"
    )
    parser.add_argument(
        "--api_key", "-k", type=str, required=True, help="API key for authentication"
    )
    parser.add_argument(
        "--build_dir",
        "-b",
        type=str,
        required=True,
        help="Directory to save the CSV file",
    )

    args = parser.parse_args()
    screen = args.screen
    api_key = args.api_key
    build_dir = args.build_dir

    # Construct endpoint urls
    base_url = "https://api.clue.io/"
    sample_meta_url = urllib.parse.urljoin(base_url, "api/v_seq_metadata")
    cell_sets_url = urllib.parse.urljoin(base_url, "api/cell-db/cell-sets")
    control_barcodes_url = urllib.parse.urljoin(
        base_url, "api/cell-db/control-barcodes"
    )

    # Get cell_sets and control_ladders from sample_meta
    sample_meta_df = fetch_metadata(
        url=sample_meta_url,
        api_key=api_key,
        filter_dict={"where": {"project_code": screen}},
    )
    cell_sets = sample_meta_df["cell_set"].unique().to_list()
    logging.info(f"Detected {len(cell_sets)} cell sets in screen {screen}: {cell_sets}")
    control_ladders = sample_meta_df["control_ladder"].unique().to_list()
    logging.info(
        f"Detected {len(control_ladders)} control ladders in screen {screen}: {control_ladders}"
    )

    # Fetch cell_sets metadatas
    cell_sets_df = fetch_metadata(url=cell_sets_url, api_key=api_key).rename(
        {"davepool_id": "cell_set", "barcode_id": "lua"}
    )
    cell_sets_screen = (
        cell_sets_df.filter(pl.col("cell_set").is_in(cell_sets))
        .select(["cell_set", "pool_id", "lua", "depmap_id"])
        .drop_nulls()
        .unique()
    )

    # Fetch control_barcode metadata
    control_barcodes_df = fetch_metadata(url=control_barcodes_url, api_key=api_key).select(
        ["name", "sequence", "log_dose", "set"]
    ).rename(
        {"name": "cb_name",
         "sequence": "forward_read_barcode",
         "set": "control_ladder"}
    )

    # Fetch cell_lines metadata
    cell_lines_url = urllib.parse.urljoin(base_url, "api/cell-db/cell-lines")
    cell_lines_df = fetch_metadata(url=cell_lines_url, api_key=api_key).select(
        ["depmap_id","dna_sequence","lua"]
    ).drop_nulls().unique()

    # Write out the dataframes
    for df, name in [
        (cell_sets_screen, "cell_set_and_pool_meta.csv"),
        (control_barcodes_df, "CB_meta.csv"),
        (cell_lines_df, "cell_line_meta.csv"),
    ]:
        save_dataframe(df, name, build_dir)


if __name__ == "__main__":
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
    )
    main()