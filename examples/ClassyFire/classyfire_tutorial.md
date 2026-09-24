# ClassyFire

!!! warning "The ClassyFire service no longer classifies new structures"
    ClassyFire (Wishart lab) has classified nothing new since February 2023.
    On 2026-09-23 a new submission was refused with HTTP 500. Results of
    queries submitted before then are still served. For ChEBI chemical classes
    computed offline, use Chebifier; see
    [Chebifier taxonomy](../../guide/chebifier.md).

`ClassyFireAPI` is kept in the package in case the service returns. Its
methods are static. The three raw calls return the `requests.Response` and
do not pace or retry; a rewrite is planned. `get_classification` returns the
parsed result, through the shared transport. The workflow is: submit, poll,
fetch.

```python
from provesid import ClassyFireAPI

response = ClassyFireAPI.submit_query("aspirin", "CC(=O)OC1=CC=CC=C1C(=O)O")
query_id = response.json()["id"]                 # when the service accepts it

ClassyFireAPI.query_status(query_id).text        # "Done" when finished
result = ClassyFireAPI.get_classification(query_id)
result["entities"][0]["kingdom"]["name"], result["entities"][0]["class"]["name"]
```

An existing query can still be fetched by its ID:

```python
result = ClassyFireAPI.get_classification(1)   # 66 pages of 10
result["classification_status"], len(result["entities"])
```

The server stops sending after about 108 KB but still answers HTTP 200, so
a query larger than that never arrives whole. `get_query(1)` gets a response
cut off mid-body and returns `None`, with a WARNING in the log.
`get_classification` asks for pages of `ClassyFireAPI.PAGE_SIZE` (10)
entities, which arrive complete, and joins them. It raises `ClassyFireError`
when a page fails. The server also answers HTTP 429 when asked too often;
the transport waits and retries.

`get_query` takes `page` and `per_page` too, and its `format` may also be
`"sdf"` or `"csv"`. The raw calls return `None` on a connection failure, not
an exception. The [API reference](../../api/classyfire.md) lists every
method.
