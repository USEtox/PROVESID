import requests

class ClassyFireAPI:
    """
    Class to interact with the ClassyFire API. The class is converted from the original Ruby code
    provided at https://bitbucket.org/wishartlab/classyfire_api/src/master/lib/classyfire_api.rb
    It is converted to Python by using the requests library, using copilot to help with the conversion.
    It is then tested by the author (simulkade) to ensure that it works as expected.
    The following is an example of how to use the ClassyFireAPI class:
    ```
    import time
    import json
    from data_extract import ClassyFireAPI
    # SMILES string to be queried
    smiles = "C1=CC(=CC=C1[N+](=O)[O-])Cl"

    # Submit the query
    response = ClassyFireAPI.submit_query("Example Query", smiles)

    # Check if the response is successful
    if response.status_code == 200:
        query_result = response.json()
        query_id = query_result['id']
        print(f"Query submitted successfully. Query ID: {query_id}")

        # Wait for some time to allow the query to be processed
        time.sleep(1)  # Adjust the sleep time as needed

        # Retrieve the classification results
        result_response = ClassyFireAPI.get_query(query_id, format="json")

        if result_response.status_code == 200:
            classification_results = result_response.json()
            print("Classification Results:")
            print(json.dumps(classification_results, indent=2))
        else:
            print(f"Failed to retrieve classification results. Status code: {result_response.status_code}")
    else:
        print(f"Failed to submit query. Status code: {response.status_code}")
    ```
    """
    URL = 'http://classyfire.wishartlab.com'

    @staticmethod
    def submit_query(label, input, type='STRUCTURE'):
        try:
            response = requests.post(
                f"{ClassyFireAPI.URL}/queries",
                json={"label": label, "query_input": input, "query_type": type},
                headers={"Accept": "application/json", "Content-Type": "application/json"}
            )
            response.raise_for_status()
        except requests.exceptions.HTTPError as e:
            return e.response
        except requests.exceptions.RequestException as e:
            return e.response
        return response
    
    @staticmethod
    def query_status(query_id):
        """
        Retrieves the status of a query.
        :param query_id: The ID of the query.
        :return: response
        """
        try:
            response = requests.get(
                f"{ClassyFireAPI.URL}/queries/{query_id}/status.json",
                headers={"Accept": "application/json", "Content-Type": "application/json"}
            )
            response.raise_for_status()
            return response
        except requests.exceptions.RequestException as e:
            return None

    @staticmethod
    def get_query(query_id, format="json"):
        try:
            if format == "json":
                response = requests.get(
                    f"{ClassyFireAPI.URL}/queries/{query_id}.json",
                    headers={"Accept": "application/json"}
                )
            elif format == "sdf":
                response = requests.get(
                    f"{ClassyFireAPI.URL}/queries/{query_id}.sdf",
                    headers={"Accept": "chemical/x-mdl-sdfile"}
                )
            elif format == "csv":
                response = requests.get(
                    f"{ClassyFireAPI.URL}/queries/{query_id}.csv",
                    headers={"Accept": "text/csv"}
                )
            else:
                raise ValueError("Invalid format. Must be one of json, sdf, or csv.")
            response.raise_for_status()
        except requests.exceptions.HTTPError as e:
            return e.response
        except requests.exceptions.RequestException as e:
            return e.response
        return response
