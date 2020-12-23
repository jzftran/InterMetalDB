"""This script is called by build_database. It processes pdb codes
calling functions from tools. For filling the models defined in
models.py functions from fill_models are called"""
import atomium
from django.conf import settings
from fill_models import *



from tools import *


def process_pdb(code):

    """Builds all the relevant objects for any given PDB code."""
code = '1i94'
    # Fetch PDB, get lowest energy model, optimise distances in model
    pdb = atomium.fetch(code)
    model, assembly_id = sort_mod(pdb)
    model.optimise_distances() #organises atoms into grids, making all neighbour searching faster
    #add PDB record to database
    pdb_record = add_pdb_record(pdb, assembly_id) #creates PDB record in a database


    # Save any metal of interests not in biological assembly but present in assymetric unit
    #changes
    for metal in sorted(metal_of_interest_outside_model(model, pdb), key=lambda metal: metal.id):
        add_metal_record(metal, pdb_record, omission="Metal is in an asymmetric unit but not in a biological assembly.")

    # Get all metals and remove duplicates
    metals = delete_duplicated_atoms(model.atoms(is_metal=True)) #removes duplicates metals, that may have been generated during symmetry operation

    #get liganding atoms for every metal
    metals_with_liganding_atoms = dict()
    pseudo_MFS =dict()
    for metal in metals:
        metals_with_liganding_atoms[metal]=get_atom_liganding_atoms(metal)

    # Remove non protein atoms
    only_protein_coordinating_residues = remove_non_protein_atoms(metals_with_liganding_atoms)
    only_protein_coordinating_residues
    interfacial_metals_only_protein=remove_non_interfacial_metals(only_protein_coordinating_residues)

interfacial_metals_only_protein
    for metal in interfacial_metals_only_protein:
        pseudo_MFS[metal] = get_atom_liganding_atoms(metal, def_cutoff=5, pseudo_MFS=True)

    # If PDB file contains interfacial metal it is noted in database
    if interfacial_metals_only_protein:
        add_has_metal_interface(pdb)

    for MOI in metal_list:
        # Uses metals sets and creates sets of residues and chains for all binding moieties
        sites_all_residues = metal_sites(MOI, remove_non_interfacial_metals(metals_with_liganding_atoms)).sites
        for site_all_residues in sites_all_residues:
            site_all_residues["residues"] = site_residues(site_all_residues)
            site_all_residues["chains"]= site_chains(site_all_residues)

        # Gets list of binding sites from coordinating_residues dictionary for only protein ligands
        sites_only_protein = metal_sites(MOI, interfacial_metals_only_protein).sites

        # Creates sets of residues and chains
        for site_only_protein in sites_only_protein:
            site_only_protein["residues"] = site_residues(site_only_protein)
            site_only_protein["chains"]= site_chains(site_only_protein)

        sites_only_protein_pseudo_MFS = metal_sites(MOI, pseudo_MFS).sites

        for site_only_protein_pseudo_MFS in sites_only_protein_pseudo_MFS:
            site_only_protein_pseudo_MFS["residues"] =site_residues(site_only_protein_pseudo_MFS)

        #gets two lists, one of metals that are coordinated only by protein ligands, second one are metals that are
        #coordinated by all ligands, including nonprotein ligands
        #Saves sites to database
        index = 1

        for site, site_all, p_MFS in zip(sites_only_protein, sites_all_residues, sites_only_protein_pseudo_MFS):
            #some sites are not ready for presentation, because they come from molecules like iron-cluster sites
            #this feature will be added in the future
            if list(site['metals'].keys())[0].het.name not in metal_list:
                presentation_ready = False
            if list(site['metals'].keys())[0].het.name  in metal_list:
                presentation_ready = True
                #PROBLEM if pdb contains many sites that are considered good will perform this many times... Need to change this in the future
                add_ready_for_presentation(pdb)
            site_record = create_site_record(site, MOI, pdb_record, index, site_all, p_MFS,presentation_ready)
            index = index +1
