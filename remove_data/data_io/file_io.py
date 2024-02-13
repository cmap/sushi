import pandas as pd

def read_csv_file(filepath):
    """
    Reads a CSV file into a pandas DataFrame.
    
    Parameters:
    - filepath: str, path to the CSV file.
    
    Returns:
    - DataFrame containing the CSV data.
    """
    try:
        df = pd.read_csv(filepath, na_values="NA")
        return df
    except Exception as e:
        raise IOError(f"Error reading {filepath}: {e}")
    
def write_csv_file(df, filepath):
    """
    Writes a pandas DataFrame to a CSV file, with error checking.
    
    Parameters:
    - df: DataFrame, the DataFrame to be written to the CSV file.
    - filepath: str, the path where the CSV file will be saved.
    
    Returns:
    - A success message if the file is written successfully, or raises an IOError on failure.
    """
    try:
        df.to_csv(filepath, index=False)
        return f"File successfully written to {filepath}"
    except Exception as e:
        raise IOError(f"Error writing to {filepath}: {e}")
