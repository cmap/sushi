# Remove Data Module for SUSHI Pipeline

### John Davis 02/24

The `remove_data` module is a component of the Sushi pipeline deisgned to selectively remove data based on predefined criteria.

## Features

- **Data filtering**: Remove rows from a dataset based on single or multiple criteria.
- **Flexible criteria input**: Supports criteria defined in CSV files or fetched dynamically via an API, based on user defined inputs.
- **Validation and testing**: Ensures integrity of input data and criteria before processing.

## Getting Started

### Prerequisites

- Python 3.8 or later
- pandas
- requests (for API calls)

### Installation

1. Clone the SUSHI pipeline repository:

    ```sh
    git clone https://github.com/yourusername/sushi.git
    cd sushi/remove_data
    ```

2. Install required Python packages:

    ```sh
    pip install -r requirements.txt
    ```

### Usage

The `remove_data` module can be executed directly from the command line with various options for specifying input data and criteria.

1. **Using CSV Files for Criteria**

    ```sh
    python main.py --data_file path/to/data.csv --criteria_file path/to/criteria.csv --output_file path/to/output.csv
    ```

2. **Using an API Call for Criteria (note: to be implemented when this table is added to the database)**

    ```sh
    python main.py --data_file path/to/data.csv --project_code PROJECT123 --output_file path/to/output.csv
    ```

### File Contents
**data_file**: This is your underlying data file. While it can work with any properly formatted csv, generally it will be the output of the `filter_raw_reads` module pre normalization.

**criteria_file**: This is the file detailing the data you wish to be removed. It can be as specific or general as you want, however **all of the column headers in this file must match headers that exist in your data_file**. Each row is used independently to remove data. For example, if the first row has a column `ccle_name: NCIH358_LUNG` and the second row has columns `ccle_name: TE1_OESOPHAGUS` and `pcr_well: A03`, all NCIH358_LUNG rows will be removed, while only TE1_OESOPHAGUS samples that come from pcr_well A03 will be removed.

**output_file**: Where to save the resultant dataframe with the appropriate data removed.

**project_code**: This is the stand-in identifier for querying the API to build the criteria_df. This may change upon implementation.

### Docker Support

The module is designed to be containerized with Docker

- Build the Docker image:

    ```sh
    docker build -t remove_data .
    ```

- Run the module inside a Docker container:

    ```sh
    docker run remove_data [COMMAND]  # Replace [COMMAND] with the appropriate command line arguments.
    ```

## Development

### Running Tests

Ensure the functionality of the module by running unit tests:

```sh
python -m unittest discover tests
```
