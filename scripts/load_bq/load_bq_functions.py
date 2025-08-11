# load_bq_functions.py
import logging
from google.cloud import bigquery
import os
import tempfile
import polars as pl
from typing import Dict, Any


def get_polars_dtype_from_bq_field(field: bigquery.SchemaField) -> pl.DataType:
    """Map BigQuery field types to Polars data types."""
    type_mapping = {
        "STRING": pl.Utf8,
        "INTEGER": pl.Int64,
        "INT64": pl.Int64,
        "FLOAT": pl.Float64,
        "FLOAT64": pl.Float64,
        "BOOLEAN": pl.Boolean,
        "BOOL": pl.Boolean,
        "DATE": pl.Date,
        "DATETIME": pl.Datetime,
        "TIMESTAMP": pl.Datetime,
        "TIME": pl.Time,
        "NUMERIC": pl.Float64,
        "BIGNUMERIC": pl.Float64,
    }
    
    bq_type = field.field_type.upper()
    polars_type = type_mapping.get(bq_type, pl.Utf8)  # Default to string if unknown
    
    if field.mode == "REPEATED":
        # For repeated fields, wrap in List
        return pl.List(polars_type)
    
    return polars_type


def coerce_dataframe_types(df: pl.DataFrame, table: bigquery.Table) -> pl.DataFrame:
    """Coerce DataFrame column types to match BigQuery table schema."""
    schema_map = {field.name: field for field in table.schema}
    
    for col in df.columns:
        if col in schema_map:
            field = schema_map[col]
            target_dtype = get_polars_dtype_from_bq_field(field)
            
            try:
                # Only cast if the current type doesn't match the target
                if df[col].dtype != target_dtype:
                    logging.info(f"Coercing column '{col}' from {df[col].dtype} to {target_dtype}")
                    df = df.with_columns(pl.col(col).cast(target_dtype, strict=False))
            except Exception as e:
                logging.warning(f"Failed to coerce column '{col}' to {target_dtype}: {e}")
                # Continue without coercing this column
                pass
    
    return df


def filter_csv_to_matching_columns(
        file_path: str, table: bigquery.Table, build_name: str, screen: str
) -> str:
    """Creates a temp CSV with columns correctly aligned to the BigQuery table schema."""
    schema_map = {field.name: field for field in table.schema}

    df = pl.read_csv(file_path, null_values="NA", infer_schema_length=10000)

    # 1. Identify missing columns and columns to drop
    columns_to_keep = [field.name for field in table.schema]
    original_columns = set(df.columns)

    missing_cols = [col for col in columns_to_keep if col not in original_columns]
    extra_cols = original_columns - set(columns_to_keep)

    # 2. Add missing columns with null values
    for col in missing_cols:
        # Get the correct Polars dtype from the BQ schema
        field = schema_map.get(col)
        if field:
            dtype = get_polars_dtype_from_bq_field(field)
            df = df.with_columns(pl.lit(None).cast(dtype).alias(col))
        else:
            logging.warning(f"Could not find schema for missing column: {col}")

    # 3. Select columns in the correct order
    df = df.select(columns_to_keep)

    # 4. Add sushi_build and screen column
    if build_name:
        df = df.with_columns(pl.lit(build_name).alias("sushi_build"))
        df = df.with_columns(pl.lit(screen).alias("screen"))

    # 5. Coerce data types and write to file...
    df = coerce_dataframe_types(df, table)

    temp_fd, temp_path = tempfile.mkstemp(suffix=".csv")
    os.close(temp_fd)
    df.write_csv(temp_path)
    logging.info(f"Final columns for upload: {df.columns}")
    return temp_path


def init_logger():
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
    )


def get_bigquery_client():
    return bigquery.Client(project="prism-359612")


def table_exists(client, dataset_id, table_id):
    try:
        client.get_table(f"{dataset_id}.{table_id}")
        return True
    except Exception:
        return False


def delete_rows_for_build(client, dataset_id, table_id, build_name):
    query = f"""
    DELETE FROM `{dataset_id}.{table_id}`
    WHERE sushi_build = @build
    """
    job_config = bigquery.QueryJobConfig(
        query_parameters=[bigquery.ScalarQueryParameter("build", "STRING", build_name)]
    )
    logging.info(
        f"Deleting rows from {dataset_id}.{table_id} for sushi_build='{build_name}'"
    )
    client.query(query, job_config=job_config).result()


def load_csv_to_bigquery(client, dataset_id, table_id, file_path, build_name, screen):
    if build_name:
        delete_rows_for_build(client, dataset_id, table_id, build_name)

    table_ref = client.dataset(dataset_id).table(table_id)
    table = client.get_table(table_ref)  # Get full table object including schemaIn the

    # Filter the CSV to only include matching columns
    filtered_csv = filter_csv_to_matching_columns(file_path, table, build_name, screen)

    job_config = bigquery.LoadJobConfig(
        source_format=bigquery.SourceFormat.CSV,
        skip_leading_rows=1,
        autodetect=False,
        schema=table.schema,
        write_disposition=bigquery.WriteDisposition.WRITE_APPEND,
        ignore_unknown_values=True,
    )

    with open(filtered_csv, "rb") as source_file:
        job = client.load_table_from_file(source_file, table_ref, job_config=job_config)

    job.result()
    logging.info(f"Appended {filtered_csv} to {dataset_id}.{table_id}")
