import argparse
from data_io.file_io import read_csv_file, write_csv_file
from data_io.api_io import fetch_criteria_from_api
from functions.remove_functions import remove_data_based_on_criteria

def parse_arguments():
    """
    Parses command-line arguments.
    """
    parser = argparse.ArgumentParser(description='Remove data based on specified criteria.')
    parser.add_argument('--data_file', '-d', type=str, help='Path to the CSV file containing the main dataset.', required=True)
    parser.add_argument('--criteria_file', '-c', type=str, help='Path to the CSV file containing removal criteria.')
    parser.add_argument('--project_code', '-p', type=str, help='Project code to fetch criteria via API call. Overrides criteria_file.')
    parser.add_argument('--output_file', '-o', type=str, default='output.csv', help='Path to the output CSV file.', required=True)
    
    return parser.parse_args()

def main():
    args = parse_arguments()
    
    # Read the data DataFrame
    data_df = read_csv_file(args.data_file)
    
    # Determine mode based on arguments provided and get the criteria DataFrame accordingly
    if args.criteria_file:
        criteria_df = read_csv_file(args.criteria_file)
    elif args.project_code:
        raise NotImplementedError("This feature will be enabled when removal critera are migrated to the database.")
        # Placeholder for API URL construction. Adjust with actual logic to construct the API URL from the project_code.
        # api_url = f"http://example.com/api/criteria?project_code={args.project_code}"
        # criteria_df = fetch_criteria_from_api(api_url)
    else:
        raise ValueError("Either --criteria_file or --project_code must be provided.")
    
    # Remove data based on criteria
    data_df, n_rows_removed = remove_data_based_on_criteria(data_df, criteria_df)
    print(f"Rows removed: {n_rows_removed}")
    
    # Write the resulting DataFrame to the output file
    write_csv_file(data_df, args.output_file)
    print(f"Output saved to {args.output_file}")

if __name__ == "__main__":
    main()
