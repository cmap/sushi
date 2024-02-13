import pandas as pd

def validate_criteria_df(criteria_df, data_df):
    """
    Validates the criteria DataFrame against the data DataFrame.
    
    Parameters:
    - criteria_df: pandas DataFrame, the criteria for data removal.
    - data_df: pandas DataFrame, the main dataset.
    
    Raises:
    - ValueError: If there are columns in criteria_df that do not exist in data_df or if criteria_df has no non-NA values.
    """
    # Check for non-NA/null values in criteria_df
    if criteria_df.dropna(how='all').empty:
        raise ValueError("criteria_df must contain at least one non-NA/null value.")

    # Check that all columns in criteria_df exist in data_df
    missing_columns = set(criteria_df.columns) - set(data_df.columns)
    if missing_columns:
        raise ValueError(f"Columns not found in data_df: {', '.join(missing_columns)}")

def remove_data_based_on_criteria(data_df, criteria_df):
    """
    Removes data rows from data_df based on criteria specified in criteria_df.
    
    Parameters:
    - data_df: pandas DataFrame, the main dataset from which data will be removed.
    - criteria_df: pandas DataFrame, criteria for data removal.
    
    Returns:
    - pandas DataFrame with data removed based on criteria.
    - int, the number of rows removed.
    
    Raises:
    - ValueError: If no rows are removed based on the criteria.
    """
    validate_criteria_df(criteria_df, data_df)

    initial_row_count = len(data_df)
    for index, criteria_row in criteria_df.iterrows():
        non_empty_criteria = criteria_row.dropna()
        query_parts = [f"{col} == {pd.to_numeric(val) if pd.api.types.is_number(val) else repr(val)}" for col, val in non_empty_criteria.items()]
        query_str = " & ".join(query_parts)
        
        if query_str:
            data_df = data_df.query(f"not ({query_str})", engine='python')

    n_rows_removed = initial_row_count - len(data_df)
    if n_rows_removed == 0:
        raise ValueError("No rows were removed based on the provided criteria.")

    return data_df, n_rows_removed
