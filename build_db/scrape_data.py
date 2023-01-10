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
    #date = date.strftime('%Y-%m-%dT%H:%M:%S')
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
    },
    "return_type": "entry"
    }'''

    query_string= json.loads(string)
    url = 'https://search.rcsb.org/rcsbsearch/v2/query'

    response = requests.post(url,json=query_string)

    if response.status_code == 200: #request has succeeded
        jsonresponse = response.json()
        return jsonresponse['total_count']

    raise Exception("No fetched codes from RCSB.")




with_metals_no_interface = Pdb.objects.filter(has_metal_interface=False).count()
with_metal_interface = Pdb.objects.filter(has_metal_interface=True).count()
database_count = Pdb.objects.all().count()
newest_deposition_date = Pdb.objects.latest("deposition_date").deposition_date






pdb_count = fetch_PDB_codes(newest_deposition_date)
no_metal = int(pdb_count) - int(database_count)
update_date = datetime.date.today()
#save results to the database
add_metalPDB_Data(pdb_count,no_metal,newest_deposition_date, update_date)

# not really efficient

def get_DB_composition(**filter_param):
    """Returns querysets of complexes with certain compostion. Protein-protein, protein-DNA etc. """
    #  https://pdb101.rcsb.org/learn/guide-to-understanding-pdb-data/primary-sequences-and-the-pdb-format
    # So far atomium does not contain type identification, thus type finding is more complicated
    amino_acids= ['ALA.', 'CYS.', 'ASP.', 'GLU.', 'PHE.', 'GLY.', 'HIS.', 'ILE.', 'LYS.', 'LEU.', 'MET.', 'ASN.', 'PRO.', 'GLN.', 'ARG.', 'SER.', 'THR.', 'VAL.', 'TRP.', 'TYR.']
    deoxyribonucleotides = ['DA.', 'DC.', 'DG.', 'DT.', 'DI.']

    ribonucleotides = [r'(?<!.)A\.', r'(?<!.)C\.', r'(?<!.)G\.', r'(?<!.)U\.', r'(?<!.)I\.', r'^A\.', r'^C\.', r'^G\.', r'^U\.', r'^I\.', r'\.A\.', r'\.C\.', r'\.G\.', r'\.U\.', r'\.I\.']
    #DNA Protein RNA complexes
    DPR = MetalSite.objects.filter(reduce(operator.or_, (Q(aminoacid_residue_names__regex=x) for x in ribonucleotides)), **filter_param, ready_for_presentation=True).filter(reduce(operator.or_, (Q(aminoacid_residue_names__icontains=x) for x in deoxyribonucleotides))).filter(reduce(operator.or_, (Q(aminoacid_residue_names__icontains=x) for x in amino_acids)))

    # DNA Protein complexes
    DP = MetalSite.objects.filter(reduce(operator.or_, (Q(aminoacid_residue_names__icontains=x) for x in deoxyribonucleotides)), **filter_param, ready_for_presentation=True).filter(reduce(operator.or_, (Q(aminoacid_residue_names__icontains=x) for x in amino_acids)))
    DP = DP.difference(DPR)

    # DNA RNA
    DR = MetalSite.objects.filter(reduce(operator.or_, (Q(aminoacid_residue_names__regex=x) for x in ribonucleotides)), **filter_param, ready_for_presentation=True).filter(reduce(operator.or_, (Q(aminoacid_residue_names__icontains=x) for x in deoxyribonucleotides)))
    DR = DR.difference(DPR)

    # Protein RNA
    PR = MetalSite.objects.filter(reduce(operator.or_, (Q(aminoacid_residue_names__regex=x) for x in ribonucleotides)), **filter_param, ready_for_presentation=True).filter(reduce(operator.or_, (Q(aminoacid_residue_names__icontains=x) for x in amino_acids)))
    PR = PR.difference(DPR)


    # Protein protein
    PP = MetalSite.objects.filter(reduce(operator.or_, (Q(aminoacid_residue_names__icontains=x) for x in amino_acids)), **filter_param, ready_for_presentation=True)
    PP = PP.difference(DP, DR, PR, DPR)


    # DNA DNA
    DD = MetalSite.objects.filter(reduce(operator.or_, (Q(aminoacid_residue_names__icontains=x) for x in deoxyribonucleotides)), **filter_param, ready_for_presentation=True)
    DD = DD.difference(DP, DR, PR, DPR)


    # RNA RNA
    RR = MetalSite.objects.filter(reduce(operator.or_, (Q(aminoacid_residue_names__regex=x) for x in ribonucleotides)), **filter_param, ready_for_presentation=True)
    RR = RR.difference(DP, DR, PR, DPR)

    DP =len(DP)
    DR = len(DR)
    PR = len(PR)
    PP = len(PP)
    DD= len(DD)
    RR= len(RR)
    DPR = len(DPR)
    #other_complexes
    other_complexes = MetalSite.objects.filter(**filter_param, ready_for_presentation=True).count()

    return [DPR, DP, DR, PR, PP, DD, RR, other_complexes-DPR-DP-DR-PR-PP-DD-RR]


elements = ["MG","PT4","IR3","BS3","4TI","BA","AU3","FE2","MN3","CU1","ZN","NA","K","CA","RB","SR","CS","SC","V","CR","MN","FE","CO","NI","CU","LI","Y","ZR", "MO","TC","RU","RH","PD","AG","CD","W","RE","OS","IR","PT","AU","HG","AL","GA", "IN","SB","TL","PB", "CE", "U","TH","NP","PU","AM","CF","NO"]
elements = sorted(elements)
lanthanides = ['LA',"SM","YB","PR","LU","TB", "GD","GD3","HO","HO3","EU3","EU",]
lanthanides = sorted(lanthanides)
all_elements = elements + lanthanides



composition_statistics = get_DB_composition(representative=True)
add_DBcomposition(*composition_statistics,representative=True, element = 'All')
composition_statistics = get_DB_composition()
add_DBcomposition(*composition_statistics,representative=None, element = 'All')

for metal in all_elements:
    composition_statistics = get_DB_composition(representative=True, element = metal)
    add_DBcomposition(*composition_statistics,representative=True, element = metal)
    composition_statistics = get_DB_composition(element = metal)
    add_DBcomposition(*composition_statistics,representative=None, element = metal)
print("Categorizing DB succeeded.")
make_log(f"DB categorization has succeded at:{datetime.now()}")
