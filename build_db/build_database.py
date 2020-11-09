"""Build and update database."""

from tools import *
env_django()
from tqdm import tqdm
import traceback
from django.db import transaction
from limb.models import Pdb
from process_pdb import *
from datetime import datetime

make_log(f"Build database has started at: {datetime.now()}\nFollowing codes could not be processed:\n")
def main():
    # Iterates over list of metals, fetching all PDB codes that have specific metal
    for metal_element in metal_list:
        PDB_codes = fetch_PDB_codes(metal_element)
        print(f"{len(PDB_codes)} PDB codes with {metal_element} found")

        # Gets PDB codes alredy in database
        codes_in_db = Pdb.objects.all().values_list("id", flat=True) #calls Pdb table defined by models, retrieves id, and creates a list of codes already in db
        new_codes = []
        #checks which codes are new
        for pdb_id in PDB_codes:
            if pdb_id not in codes_in_db:
                new_codes.append(pdb_id)

        print(f"{len(new_codes)} of these are new")

        #Part bellow process pdb codes with "process_pdb"
        not_checked = {} #creates set of PDB codes that could not be processed, such codes will be wrtitten into log
        for pdb_id in tqdm(new_codes):
            try:
                with transaction.atomic(): process_pdb(pdb_id)
            except Exception as e: not_checked[pdb_id] = traceback.format_exc()  #adds PDB ID to log
        print("PDBs that could not be processed:\n")
        for pdb_id in not_checked:
            make_log(pdb_id)
            print(pdb_id)

main()
make_log(f"Build database has ended at: {datetime.now()}")
