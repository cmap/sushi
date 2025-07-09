# load_bq_functions.py
import logging
from google.cloud import bigquery
import os
import tempfile
import polars as pl


def filter_csv_to_matching_columns(file_path: str, table: bigquery.Table, build_name: str) -> str:
    """Create a temp CSV with only the columns that exist in the BigQuery table."""
    allowed_columns = {field.name for field in table.schema}

    df = pl.read_csv(file_path, null_values="NA", infer_schema_length=10000)
    filtered_cols = [col for col in df.columns if col in allowed_columns]
    df = df.select(filtered_cols)

    # Add sushi_build column
    if build_name:
        df = df.with_columns(pl.lit(build_name).alias("sushi_build"))

    # Write to temp file
    temp_fd, temp_path = tempfile.mkstemp(suffix=".csv")
    os.close(temp_fd)  # Close file descriptor so Polars can write
    df.write_csv(temp_path)
    logging.info(f"Filtered {len(df.columns)} columns for upload: {filtered_cols}")
    return temp_path


def init_logger():
    logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")


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
        query_parameters=[
            bigquery.ScalarQueryParameter("build", "STRING", build_name)
        ]
    )
    logging.info(f"Deleting rows from {dataset_id}.{table_id} for sushi_build='{build_name}'")
    client.query(query, job_config=job_config).result()


def load_csv_to_bigquery(client, dataset_id, table_id, file_path, build_name):
    if build_name:
        delete_rows_for_build(client, dataset_id, table_id, build_name)

    table_ref = client.dataset(dataset_id).table(table_id)
    table = client.get_table(table_ref)  # Get full table object including schema

    # Filter the CSV to only include matching columns
    filtered_csv = filter_csv_to_matching_columns(file_path, table, build_name)

    job_config = bigquery.LoadJobConfig(
        source_format=bigquery.SourceFormat.CSV,
        skip_leading_rows=1,
        autodetect=True,
        write_disposition=bigquery.WriteDisposition.WRITE_APPEND,
        ignore_unknown_values=True
    )

    with open(filtered_csv, "rb") as source_file:
        job = client.load_table_from_file(source_file, table_ref, job_config=job_config)

    job.result()
    logging.info(f"Appended {filtered_csv} to {dataset_id}.{table_id}")

