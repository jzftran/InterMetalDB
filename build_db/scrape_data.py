"""Gets information from RCSB PDB about the number of all structures,
writes this information to the database."""
from django.conf import settings
from tools import *
env_django()
from fill_models import *
import datetime
import requests

from functools import reduce
import operator
from django.db.models import Q


def fetch_PDB_codes(date):

    date = date.strftime('%Y-%m-%dT%H:%M:%S')

    """Fetchs PDB codes for all structures with metal ion in them.
    If the response returned has an 500 error code, an error will
    be comunicated. """
    string = '''
    {
    "query": {
        "type": "group",
        "logical_operator": "and",
        "nodes": [
        {
            "type": "terminal",
            "service": "text",
            "parameters": {
            "operator": "less",
            "value":"'''+str(date)+'''",
            "attribute": "rcsb_accession_info.initial_release_date"
          }
        }
      ]
    },
    "request_options": {

        "return_counts": true

        "return_all_hits": true

    },
    "return_type": "entry"
    }'''

    query_string= json.loads(string)

    url = 'https://search.rcsb.org/rcsbsearch/v2/query'


    response = requests.post(url,json=query_string)

    if response.status_code == 200: #request has succeeded
        jsonresponse = response.json()

        return jsonresponse['total_count']

        # return [jsonresponse['result_set'][code]['identifier'] for code in range(jsonresponse['total_count'])]


    raise Exception("No fetched codes from RCSB.")




with_metals_no_interface = Pdb.objects.filter(has_metal_interface=False).count()
with_metal_interface = Pdb.objects.filter(has_metal_interface=True).count()
database_count = Pdb.objects.all().count()
newest_deposition_date = Pdb.objects.latest("deposition_date").deposition_date







pdb_count = fetch_PDB_codes(newest_deposition_date)

pdb_count = len(fetch_PDB_codes(newest_deposition_date))

no_metal = int(pdb_count) - int(database_count)
update_date = datetime.date.today()
#save results to the database
add_metalPDB_Data(pdb_count,no_metal,newest_deposition_date, update_date)


print("Categorizing DB succeeded.")
make_log(f"DB categorization has succeded at:{datetime.now()}")
