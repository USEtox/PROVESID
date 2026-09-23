# ClassyFire

!!! warning "The ClassyFire service no longer classifies new structures"
    ClassyFire (Wishart lab) has classified nothing new since February 2023.
    On 2026-09-23 a new submission was refused with HTTP 500. Results of
    queries submitted before then are still served. For ChEBI chemical classes
    computed offline, use Chebifier; see
    [Chebifier taxonomy](../../guide/chebifier.md).

`ClassyFireAPI` is kept in the package in case the service returns. Its
methods are static, return the raw `requests.Response`, and do not pace or
retry. A rewrite is planned. The workflow is: submit, poll, fetch.

```python
from provesid import ClassyFireAPI

response = ClassyFireAPI.submit_query("aspirin", "CC(=O)OC1=CC=CC=C1C(=O)O")
query_id = response.json()["id"]                 # when the service accepts it

ClassyFireAPI.query_status(query_id).text        # "Done" when finished
result = ClassyFireAPI.get_query(query_id, format="json").json()
result["entities"][0]["kingdom"]["name"], result["entities"][0]["class"]["name"]
```

An existing query can still be fetched by its ID:

```python
ClassyFireAPI.get_query(1).json()["classification_status"]   # 'Done'
```

`format` may also be `"sdf"` or `"csv"`. A connection failure comes back as
`None`, not as an exception. The [API reference](../../api/classyfire.md)
lists every method.
