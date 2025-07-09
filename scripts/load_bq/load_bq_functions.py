# load_bq_functions.py
import logging
from google.cloud import bigquery
import os

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
    job_config = bigquery.LoadJobConfig(
        source_format=bigquery.SourceFormat.CSV,
        skip_leading_rows=1,
        autodetect=True,
        write_disposition=bigquery.WriteDisposition.WRITE_APPEND
    )

    with open(file_path, "rb") as source_file:
        job = client.load_table_from_file(source_file, table_ref, job_config=job_config)

    job.result()
    logging.info(f"Appended {file_path} to {dataset_id}.{table_id}")

