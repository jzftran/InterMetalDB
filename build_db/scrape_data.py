"""Gets information from RCSB PDB about the number of all structures,
writes this information to the database."""
from django.conf import settings
from tools import *
env_django()
from fill_models import *
import datetime
import requests

with_metals_no_interface = Pdb.objects.filter(has_metal_interface=False).count()
with_metal_interface = Pdb.objects.filter(has_metal_interface=True).count()
database_count = Pdb.objects.all().count()
newest_deposition_date = Pdb.objects.latest("deposition_date").deposition_date

r = requests.post("https://www.rcsb.org/pdb/search/smartRow.do",
headers={"Content-Type": "application/x-www-form-urlencoded"},
data={"smartSearchSubtype": "ReleaseDateQuery", "r": "0", "target": "Current",
"pdbx_audit_revision_history.revision_date.comparator": "between",
"pdbx_audit_revision_history.revision_date.max": newest_deposition_date,
})

pdb_count = int(r.text.split("'")[1].split()[0])
no_metal = int(pdb_count) - int(database_count)
update_date = datetime.date.today()
#save results to the database
add_metalPDB_Data(pdb_count,no_metal,newest_deposition_date, update_date)
