import requests
import pandas as pd

def fetch_criteria_from_api(api_url):
    """Fetches removal criteria from an API.
    Placehold to be replaced when feature is implemented."""
    response = requests.get(api_url)
    if response.status_code == 200:
        data = response.json()
        # Assuming the API returns data in a format that pandas can directly convert to DataFrame
        return pd.DataFrame(data)
    else:
        raise Exception(f"Failed to fetch data from API, status code: {response.status_code}")
