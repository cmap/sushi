import argparse
import requests
import pandas as pd
import os


def fetch_data(screen, api_key):
    url = f"https://api.clue.io/api/v_eps_metadata?filter=%7B%22where%22%3A%7B%22project_code%22%3A%22{screen}%22%7D%7D"
    headers = {"Accept": "application/json", "user_key": api_key}

    response = requests.get(url, headers=headers)

    if response.status_code == 200:
        data = response.json()
        df = pd.DataFrame(data)
        return df
    else:
        print(f"Error: {response.status_code}")
        print(response.text)  # Print the response content for more insight
        response.raise_for_status()


def save_dataframe(df, build_dir):
    # Ensure the directory exists
    os.makedirs(build_dir, exist_ok=True)

    # Define the file path
    file_path = os.path.join(build_dir, "sample_meta.csv")

    # Write the DataFrame to a CSV file
    df.to_csv(file_path, index=False)
    print(f"Sample metadata written to {file_path}")


def rename_columns(df):
    rename_map = {"index1": "index_1", "index2": "index_2", "pert_id": "compound_id", "trt_type": "pert_type"}
    return df.rename(columns=rename_map)


def add_control_barcodes(df, control_barcodes):
    control_barcode_df = df.copy()
    control_barcode_df["control_barcodes"] = control_barcodes
    return control_barcode_df


def replace_pert_type(df):
    df["pert_type"] = df["pert_type"].replace(
        {"trt_poscon": "poscon", "ctl_vehicle": "negcon"}
    )
    return df


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
    parser.add_argument(
        "--control_barcodes",
        "-cb",
        type=str,
        required=True,
        help="Type of control barcode ladder, defaults to h-b",
    )

    args = parser.parse_args()
    screen = args.screen
    api_key = args.api_key
    build_dir = args.build_dir
    control_barcodes = args.control_barcodes

    try:
        df = (
            fetch_data(screen, api_key)
            .pipe(rename_columns)
            .pipe(add_control_barcodes, control_barcodes)
            .pipe(replace_pert_type)
        )
        print("Retrieved following sample_metadata from COMET: ")
        print(df.head())
        save_dataframe(df, build_dir)
    except Exception as e:
        print(f"An error occurred: {e}")


if __name__ == "__main__":
    main()
