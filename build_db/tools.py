import math
import json
import requests
import atomium
import re
metal_list = ["LU","MG","PT4","IR3","HO","HO3","GD3","GD","EU3","EU","BS3","4TI","BA","TB","SM","AU3","FE2","MN3","CU1","ZN","NA","K","CA","RB","SR","CS","SC","V","CR","MN","FE","CO","NI","CU","LI","Y","ZR", "MO","TC","RU","RH","PD","AG","CD","LA","W","RE","OS","IR","PT","AU","HG","AL","GA", "IN","SB","TL","PB", "CE","PR", "U", "YB","TH","NP","PU","AM","CF","NO"]
# ZCM curium (III)
# 4TI titanium (IV)
# BS3 bismuth (III)


def env_django():
    """Django setup -- allows to use django in script."""
    import os, sys, django
    sys.path.append(os.path.join(".."))
    os.environ.setdefault("DJANGO_SETTINGS_MODULE", "mysite.settings")
    django.setup()

def make_log(txt):
    """Returns PDB codes, that could not be processed"""
    with open('data/log.txt','a') as f:
        f.write(txt+'\n')

def fetch_PDB_codes(metal_element):

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
            "operator": "exact_match",
            "value":"'''+metal_element+'''",
            "attribute": "rcsb_chem_comp_container_identifiers.comp_id"
          }
        }
      ]
    },
    "request_options": {
        "return_all_hits": true
    },
    "return_type": "entry"
    }'''
    print(metal_element)
    query_string= json.loads(string)
    url = 'http://search.rcsb.org/rcsbsearch/v1/query'

    response = requests.post(url,json=query_string)

    if response.status_code == 200: #request has succeeded
        jsonresponse = response.json()
        return [jsonresponse['result_set'][code]['identifier'] for code in range(jsonresponse['total_count'])]

    raise Exception("No fetched codes from RCSB.")



def sort_mod(pdb):
    """Selects lowest energy model, and returns that model."""
    sorted_list = pdb.assemblies
    sorted_list.sort(key=lambda assembly: math.inf if assembly["delta_energy"] is None else assembly["delta_energy"]) #
    if sorted_list:
        model = pdb.generate_assembly(sorted_list[0]["id"])
        metals = model.atoms(is_metal=True) #gets metalls in model
        while len(metals) == 0: #if no metal in assembly, assembly is removed and next assembly is checked if contains metal
            del sorted_list[0]
            model = pdb.generate_assembly(sorted_list[0]["id"])
            metals = model.atoms(is_metal=True)
        return model, sorted_list[0]["id"]
    else:
        return pdb.model, None



def metal_of_interest_outside_model(model, pdb):
    """Returns metal atoms that are not in assembly while present in input PDB."""
    metals_of_interest_ids=[]
    metals_of_interest = []
    metals_outside = []
    for element in metal_list:
        metals_of_interest = pdb.model.atoms(element=element)
        metals_of_interest_ids = [atom.id for atom in model.atoms(element=element)]
    for metal in metals_of_interest:
        if metal.id not in metals_of_interest_ids:
            metals_outside.append(metal)
    return metals_outside



#def get_atom_liganding_atoms(metal,def_cutoff=3, pseudo_MFS=False):
#    """Atoms coordinating metal ions in proteins rarely create angle less than 45 degrees
#    between other atoms with metal ion as a vertex. Such atom is remved from the list.
#    Because metals are not coordinated by carbon, hydrogen atoms, these are removed in the first place.
#    Starting with the closest atom, angle check is performed for the rest of atoms."""
#    nearby_atoms = []
#    for atom in metal.nearby_atoms(cutoff=def_cutoff, is_metal=False):
#        if pseudo_MFS == False:
#            if atom.element not in "CH":
#                nearby_atoms.append(atom)
#        else:
#            nearby_atoms.append(atom)
#    nearby_atoms = list(delete_duplicated_atoms(nearby_atoms))
#    nearby_atoms.sort(key=lambda atom: atom.distance_to(metal))
#    if pseudo_MFS == False:
#        for atom in nearby_atoms:
#            for atom2 in nearby_atoms[1:]:
#                if atom2 !=atom:
#                    if metal.angle(atom, atom2) <math.pi / 4:
#                        nearby_atoms.remove(atom2)
#    return nearby_atoms


def get_atom_liganding_atoms(metal,def_cutoff=3, pseudo_MFS=False):
    """Atoms coordinating metal ions in proteins rarely create angle less than 45 degrees
    between other atoms with metal ion as a vertex. Such atom is remved from the list.
    Because metals are not coordinated by carbon, hydrogen atoms, these are removed in the first place.
    Starting with the closest atom, angle check is performed for the rest of atoms."""
    nearby_atoms = []
    for atom in metal.nearby_atoms(cutoff=def_cutoff, is_metal=False):
        if pseudo_MFS == False:
            if atom.element not in "CH":
                nearby_atoms.append(atom)
        else:
            nearby_atoms.append(atom)
    nearby_atoms = list(delete_duplicated_atoms(nearby_atoms))
    nearby_atoms.sort(key=lambda atom: atom.distance_to(metal))

    return nearby_atoms





def remove_non_interfacial_metals(metals_residues):
    """Takes dictionary of metal:coordinating residues, and for each dictionary element checks if metal is coordinated
    by at least two chains. Creates new set of chains, and checks if this set is longer than 1. If chains aren't uniqe, returns
    an empty set."""
    interfacial_metals=metals_residues.copy()
    for i,m in metals_residues.items():
        unique_chain_set =set()
        for residue in m:
            unique_chain_set.add(residue.chain)
            unique_chain_set
        if len(unique_chain_set) <= 1:
            interfacial_metals.pop(i)
    return(interfacial_metals)


def remove_non_protein_atoms(coordinating_atoms):
    """Takes dictionary of metals and liganding atoms. Generates new empty dictionary
    that is filled only with atoms that are only from protein.
    Returns new dictionary"""
    protein_atoms =dict()
    for m,i in coordinating_atoms.items():
        new_set= set()
        for x in m.nearby_hets(cutoff=3,residues=True, ligands=False) : #removes ligands from coordinating
            atom_list = list()
            if re.search(r'\bResidue [A\s|G\s|C\s|U\s|DA\s|DT\s|DG\s|DC\s]\b',str(x)) == None:
                new_set.add(x)
            for k in i:
                if k.het in new_set:
                    atom_list.append(k)
                    protein_atoms[m] = atom_list
    return protein_atoms


class metal_sites:
    """Creates class of residues coordinating metal."""
    def __init__(self, MOI, metals_residues):
        self.sites = sites = [{"metals": {m: l}} for m, l in metals_residues.items()]
        self.sites = [site for site in sites if MOI in [
        atom.element for atom in site["metals"].keys()
        ]]
        self.sites.sort(key=lambda k: min(atom.id for atom in k["metals"].keys())) #sorts sites



def site_residues(site):
    """Takes a dict with metal/liganding atom information, and gets all the
    residues assoicated with those liganding atoms."""
    values = site["metals"].values()
    site_residues = {atom.het for atoms in values for atom in atoms}
    return site_residues



def site_chains(site):
    """Takes a dictionary with residues, and produces chains to which these residues belong. """
    all_chains = set()
    for residue in site["residues"]:
        if isinstance(residue, atomium.Residue):
            all_chains.add(residue.chain.id)
    unique_chains = set()
    for residue in site["residues"]:
        if isinstance(residue, atomium.Residue):
            if residue.chain.id in all_chains:
                unique_chains.add(residue.chain)
    return unique_chains




def delete_duplicated_atoms(atoms):
    """Sometimes PDB files can have multiple representations of the same atom.
    In order to remove this duplicated atoms given set of atoms is transformed into a different elements set
    list, then every atom is compared with every atom in the mono-element list. From each elements list unique atoms are
    returned, and combined in unique atoms set. Uniqueness is estimated by checking if atoms are closer than 1 Angstrom.
    """

    unique_atoms = set()
    elements = set([m.element for m in atoms])
    for element in elements:
        relevant_atoms = set()
        for atom in atoms:
            if atom.element == element:
                relevant_atoms.add(atom)

        relevant_atoms = list(relevant_atoms)
        for atom in relevant_atoms:
            for atom2 in relevant_atoms[1:]:
                if atom2 !=atom:
                    if atom.distance_to(atom2) <1:
                        relevant_atoms.remove(atom2)
        unique_atoms.update(relevant_atoms)
    return unique_atoms

#from limb.models import Chains
#def chains_to_fasta():
#    lines =[]
#    for chain in Chain.objects.all():
#        lines.append(">"+str(chain.id))
#        lines.append(chain.sequence)
#    return lines
#
#fasta = chains_to_fasta()
