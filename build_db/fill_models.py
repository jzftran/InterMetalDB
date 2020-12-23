"""Instructions and functions that populate models defined in models.py"""

from limb.models import *

import atomium
def add_pdb_record(pdb, assembly_id):
    """Creates a Pdb record in the database. The input for this function are atomium file
    and assembly ID generated with Atomium."""

    return Pdb.objects.create(
     id=pdb.code,
     rvalue=pdb.rvalue,
     classification=pdb.classification,
     deposition_date=pdb.deposition_date,
     organism=pdb.source_organism,
     expression_system=pdb.expression_system,
     technique=pdb.technique,
     keywords=", ".join(pdb.keywords) if pdb.keywords else "",
     title=pdb.title,
     resolution=pdb.resolution,
     assembly=assembly_id,
     )



def add_metal_record(atom, pdb_record, site_record=None, omission=None):
    """Adds metal record to the database, uses pdb record as a primary
    key. If there is no reason for omission, site value is added as another
    foreign key."""

    residue = atom.het
    x, y, z = atom.location
    return Metal.objects.create(
     atomium_id=atom.id,
     element=atom.element,
     name=atom.name, x=x, y=y, z=z,
     residue_number=residue.id.split(".")[1],
     chain_id=atom.chain.id,
     residue_name=residue.name,
     pdb=pdb_record,
     site=site_record,
     omission_reason=omission
    )

def add_chain_record(chain, pdb_record, site_record, index):
    """Creates a chains record for each metals site.
    Need to be used for sequence search."""

    return Chain.objects.create(
    id = f"{site_record.id}-{chain.id}-{index}",
    chain_id = chain.id,
    sequence = chain.sequence,
    pdb = pdb_record,
    site = site_record,
    )


def create_site_record(site_dict, metal_element, pdb_record, index, all_site_dict, p_MFS,presentation_ready):
    """Creates a record of a metal site, and metal itself.
    For metal sites chain record is added as well."""

    # Create site record itself
    residue_names = set()
    for residue in site_dict["residues"]:
        residue_names.add(residue.name+'.')
    residue_names = sorted(list(residue_names))

    residue_ids = set()
    for residue in site_dict["residues"]:
        residue_ids.add(residue.name + ' '+ residue.id+' ')
    residue_ids = sorted(list(residue_ids))

    all_residue_names = set()
    for residue in all_site_dict["residues"]:
        all_residue_names.add(residue.name+'.')
    all_residue_names = sorted(list(all_residue_names))

    all_residue_ids = set()
    for residue in all_site_dict["residues"]:
        all_residue_ids.add(residue.name + ' '+ residue.id+' ')
    all_residue_ids = sorted(list(all_residue_ids))




    site_record = MetalSite.objects.create(
         id=f"{pdb_record.id}-{metal_element}-{index}",
         element = metal_element,
         ion = list(site_dict['metals'].keys())[0].het.name, # metal residue name, contains information about oxidation
         binding_family=create_site_family(site_dict["residues"]),
         pseudo_MFS = create_site_family(p_MFS["residues"]),
         pdb_id=pdb_record,
         aminoacid_residue_names="".join(residue_names),
         all_residue_names ="".join(all_residue_names),
         residue_ids ="".join(residue_ids),
         all_residue_ids ="".join(all_residue_ids),
         is_homomer=check_if_homomer(site_dict["chains"]),
         number_of_coordinating_residues = len(all_site_dict["residues"]),
         number_of_coordinating_aminoacid_residues = len(site_dict["residues"]),
         number_of_coordinating_chains = len(site_dict['chains']),
         number_of_nonprotein_residues = len(all_site_dict["residues"]) - len(site_dict["residues"]),
         ready_for_presentation = presentation_ready,

        )

    # Create metals
    metals_keys = site_dict["metals"].keys()

    metals_keys=sorted(metals_keys, key= lambda metal: metal.id)
    metals_dict ={}
    for metal in metals_keys:
        metals_dict[metal.id] = add_metal_record(metal, pdb_record, site_record)



#    #create atom and bonds records

    for metal, atoms in sorted(site_dict["metals"].items(), key=lambda a: a[0].id):
        for atom in sorted(atoms, key=lambda a: a.id):
            atoms_dict = {}
            atoms_dict[atom] =create_atom_record(atom)
            Bonds.objects.create(
             metal=metals_dict[metal.id], atom=atoms_dict[atom]
            )
#
    for residue in sorted(site_dict["residues"], key=lambda residue: residue.id):
        Residues.objects.create(
            name = residue.name,
            pdb = pdb_record,
            site =site_record
             )


    # create chain records
    chains_list = site_dict["chains"]

    unique_chains =set()
    for chain in chains_list:
        if not chain.id in {chain.id for chain in unique_chains}:
            unique_chains.add(chain)
    index_chain = 0
    for chain in unique_chains:

        add_chain_record(chain, pdb_record, site_record, index_chain)
        index_chain = index_chain+1

def create_site_family(residues):
    """Generates a family string (H3, C2H2 etc.) from  a list of residues."""
    codes = list()
    for residue in residues:
        if isinstance(residue, atomium.Residue):
            codes.append(atomium.data.CODES.get(residue.name, "X")) #dunno
    output=''
    for code in sorted(set(codes)):
        to_join = (f"{code}{codes.count(code)}")
        output=output+to_join
    return output


def check_if_homomer(chains):
    """Generates an empty set, then adds sequences to set.
    Then checks if set has one unique sequence, if so, metal site is homomeric."""
    chain_set = set()
    for chain in chains:
        sequence = chain.sequence
        chain_set.add(sequence)
    if len(chain_set) == 1:
        return True
    else:
        return False

def add_has_metal_interface(pdb):
    """Changes the value of has_metal_interface field to True if PDB contains interfacial
    metal binding site."""
    pdb_record_to_update = Pdb.objects.get(id=pdb.code)
    pdb_record_to_update.has_metal_interface = True
    pdb_record_to_update.save()

def add_ready_for_presentation(pdb):
    """Changes the value of ready_for_presentation field to True if PDB is ready for presentation, i.e. contains
    metal site that is ready for presentation.
    PDB is not ready for presentation if interfacial metal binding site is created by a metal that is self.assert_(
    part of another molecule. Such feature will be implemented in the future"""
    pdb_record_to_update = Pdb.objects.get(id=pdb.code)
    pdb_record_to_update.ready_for_presentation = True
    pdb_record_to_update.save()

def create_atom_record(atom):
    """Creates an Atom record from the information provided."""

    x, y, z = atom.location
    return Atom.objects.create(
     atomium_id=atom.id, name=atom.name, x=x, y=y, z=z,
     element=atom.element
    )


def add_metalPDB_Data(metals_in_pdb, no_metal, scraping_date, update_date):
    """Gets meta data about this database and RCSB database. Used in
    statistics view generation."""

    return MetaDataPDB.objects.create(
    metals_in_pdb = metals_in_pdb,
    no_metal = no_metal,
    scraping_date = scraping_date,
    update_date = update_date
    )


def add_DBcomposition(DNA_protein_RNA, DNA_protein, DNA_RNA, protein_RNA,protein_protein,DNA_DNA,RNA_RNA,other_complexes,representative,element):
    """Gets meta data about this database and RCSB database. Used in
    statistics view generation."""

    return DBcomposition.objects.create(
    DNA_protein_RNA = DNA_protein_RNA,
    DNA_protein = DNA_protein,
    DNA_RNA = DNA_RNA,
    protein_RNA = protein_RNA,
    protein_protein = protein_protein,
    DNA_DNA = DNA_DNA,
    RNA_RNA = RNA_RNA,
    other_complexes = other_complexes,
    representative = representative,
    element = element
    )
