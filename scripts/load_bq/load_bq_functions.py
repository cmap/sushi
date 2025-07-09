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

def load_csv_to_bigquery(client, dataset_id, table_id, file_path):
    table_ref = client.dataset(dataset_id).table(table_id)
    job_config = bigquery.LoadJobConfig(
        source_format=bigquery.SourceFormat.CSV,
        skip_leading_rows=1,
        autodetect=False,
        write_disposition=bigquery.WriteDisposition.WRITE_APPEND,
    )
    with open(file_path, "rb") as source_file:
        job = client.load_table_from_file(source_file, table_ref, job_config=job_config)

    job.result()
    logging.info(f"Loaded {file_path} into {dataset_id}.{table_id}")
