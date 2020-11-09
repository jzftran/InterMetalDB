"""Contains instructions for statistics generation. This should be ommited in the future
and be implemented in model instead of file that is imported by views. """
from django.shortcuts import render
from django.db.models import Count, F
from limb.models import *
import datetime
import requests
import numpy as np
from collections import Counter




# %%
try:
    def get_distance(a, b):
        """Returns distances between two given atoms."""
        return np.sqrt((a.x - b.x)**2+(a.y - b.y)**2+(a.z - b.z)**2)
    STEP = 0.01

    def get_distances_for_sites(representative=None):
        sites = MetalSite.objects.annotate(metals=Count('metal'))
        sites = sites.annotate(resolution=F("pdb_id__resolution"))
        #presentation_ready
        sites = sites.filter(ready_for_presentation=True)
        if representative!= None:
            #presentation_ready
            sites = sites.filter(representative=representative,ready_for_presentation=True)


        atoms = []
        for site in sites:
            #site resolution better than 3 and mononuclear site, if there is no resolution 0 is assumed
            if (int(site.resolution or 0)<3 and site.metals == 1) ==True:
                metal = Metal.objects.get(site=site)
                if metal.name == metal.residue_name:
                    bonds = Bonds.objects.filter(metal=metal)
                    for bond in bonds.all():
                        atoms.append({
                        "metal": metal.name, "name": bond.atom.name, "element": bond.atom.element,
                        "distance": get_distance(bond.atom, metal)
                        })
        distance = atoms
        dist_data = []
        for element in "SON":
            distances = [a["distance"] for a in distance if a["element"] == element]
            mean = sum(distances) / len(distances)
            sd = np.std(distances)
            dist_data.append([element, mean, sd])
        return distance

    distance =  get_distances_for_sites(), get_distances_for_sites(representative = True)
    """distance =  ([{'metal': 'ZN',
      'name': 'OD2',
      'element': 'O',
      'distance': 2.4495620016647868},
     {'metal': 'ZN',
      'name': 'OE1',
      'element': 'O',
      'distance': 2.7788965435942363}],
      [{'metal': 'ZN',
        'name': 'OD2',
        'element': 'O',
        'distance': 2.4495620016647868},
       {'metal': 'ZN',
        'name': 'OE1',
        'element': 'O',
        'distance': 2.7788965435942363},
        {'metal': 'ZN',
         'name': 'OE1',
         'element': 'O',
         'distance': 2.7788965435942363},
         {'metal': 'ZN',
          'name': 'OE1',
          'element': 'O',
          'distance': 2.7788965435942363}])"""


    def get_residues(ion = None, representative = None):
        """Returns residues that are coordinating metals"""
        if ion and representative == None:
            #, ready_for_presentation=True
            MetalSites_with_ion = MetalSite.objects.filter(ion = ion, ready_for_presentation=True).values_list('id', flat=True)
            residues = Residues.objects.filter(site_id__in=MetalSites_with_ion)


        if ion and representative !=None:
            #Django backward relationship lookup
            MetalSites_with_ion = MetalSite.objects.filter(ion = ion, representative = representative, ready_for_presentation=True).values_list('id', flat=True)
            residues = Residues.objects.filter(site_id__in=MetalSites_with_ion)

        if ion == None and representative !=None:
            MetalSites_representative= MetalSite.objects.filter( representative = representative, ready_for_presentation=True).values_list('id', flat=True)
            residues = Residues.objects.filter(site_id__in=MetalSites_representative)

        if ion ==None and representative == None:
            #, ready_for_presentation=True
            #residues = Residues.objects.all()
            MetalSites_representative= MetalSite.objects.filter(ready_for_presentation=True).values_list('id', flat=True)
            residues = Residues.objects.filter(site_id__in=MetalSites_representative)


        residues_dict = Counter([str(residue.name).upper() for residue in residues]).most_common(9)
        residues_data = [n[1] for n in residues_dict]
        residues_data.append(len(residues)-sum(residues_data))

        residues_labels = [n[0] for n in residues_dict]
        residues_labels.append('OTHER')
        return (residues_data, residues_labels)
####
    def get_year(ion = None, representative = None):
        if ion and representative == None:
            #Django backward relationship lookup
            #, ready_for_presentation=True
            PDBs_with_ion = MetalSite.objects.filter(ion = ion, ready_for_presentation=True).values_list('pdb_id', flat=True)
            pdbs = Pdb.objects.filter(id__in=PDBs_with_ion)

        if ion and representative !=None:
            PDBs_with_ion = MetalSite.objects.filter(ion = ion, representative = representative, ready_for_presentation=True).values_list('pdb_id', flat=True)
            pdbs = Pdb.objects.filter(id__in=PDBs_with_ion)

        if ion == None and representative !=None:
            PDBs_with_ion = MetalSite.objects.filter(representative = representative, ready_for_presentation=True).values_list('pdb_id', flat=True)
            pdbs = Pdb.objects.filter(id__in=PDBs_with_ion)

        if ion ==None and representative == None:
            pdbs = Pdb.objects.filter(has_metal_interface=True, ready_for_presentation=True)

        year = Counter([str(protein.deposition_date).upper()[:-6] for protein in pdbs]).most_common()
        year.sort(key=lambda pair: pair[0], reverse=False)
        year_labels = [n[0] for n in year]
        year_data =[n[1] for n in year]

        return (year_data, year_labels)

###
    def get_techniques(ion = None, representative = None):
        """Returns techniqes that were used to acquire structural information of the protein."""
        if ion and representative == None:
            #Django backward relationship lookup
            PDBs_with_ion = MetalSite.objects.filter(ion = ion, ready_for_presentation=True).values_list('pdb_id', flat=True)
            pdbs = Pdb.objects.filter(id__in=PDBs_with_ion)

        if ion and representative !=None:
            PDBs_with_ion = MetalSite.objects.filter(ion = ion, representative = representative, ready_for_presentation=True).values_list('pdb_id', flat=True)
            pdbs = Pdb.objects.filter(id__in=PDBs_with_ion)

        if ion == None and representative !=None:
            PDBs_with_ion = MetalSite.objects.filter(representative = representative, ready_for_presentation=True).values_list('pdb_id', flat=True)
            pdbs = Pdb.objects.filter(id__in=PDBs_with_ion)

        if ion ==None and representative == None:
            pdbs = Pdb.objects.filter(has_metal_interface=True, ready_for_presentation=True)

        techniqes = Counter([str(protein.technique).upper() for protein in pdbs]).most_common()
        techniqes_data = [n[1] for n in techniqes]
        techniqes_labels = [n[0] for n in techniqes]
        return (techniqes_data, techniqes_labels)


    def get_classification(ion = None, representative=None):
        """Returns enzyme classifications in selected PDBs, the rest PDBs is called other."""
        if ion and representative == None:
            #Django backward relationship lookup
            PDBs_with_ion = MetalSite.objects.filter(ion = ion, ready_for_presentation=True).values_list('pdb_id', flat=True)
            pdbs = Pdb.objects.filter(id__in=PDBs_with_ion)

        if ion and representative !=None:
            PDBs_with_ion = MetalSite.objects.filter(ion = ion, representative = representative, ready_for_presentation=True).values_list('pdb_id', flat=True)
            pdbs = Pdb.objects.filter(id__in=PDBs_with_ion)

        if ion == None and representative !=None:
            PDBs_with_ion = MetalSite.objects.filter(representative = representative, ready_for_presentation=True).values_list('pdb_id', flat=True)
            pdbs = Pdb.objects.filter(id__in=PDBs_with_ion)

        if ion ==None and representative == None:
            pdbs = Pdb.objects.filter(has_metal_interface=True, ready_for_presentation=True)


        all_pdbs = Counter([protein.classification for protein in pdbs])
        classification_labels = ['HYDROLASE', 'TRANSFERASE', 'OXIDOREDUCTASE', 'LYASE', 'ISOMERASE', 'LIGASE', 'ISOMERASE']
        classification_data =[]
        for enzyme in classification_labels:
            classification_data.append(all_pdbs[enzyme])
        not_enzymes = len(pdbs)-sum(classification_data)
        classification_data.append(not_enzymes)
        classification_labels.append('OTHER')
        return (classification_data, classification_labels)
###
    ##gets family info
    def get_family_info(ion = None, representative=None):

        if ion and representative == None:
            sites = MetalSite.objects.filter(ion=ion, ready_for_presentation=True)
            all_sites = MetalSite.objects.filter(ready_for_presentation=True)
        if ion and representative !=None:
            sites = MetalSite.objects.filter(ion=ion, representative=representative, ready_for_presentation=True)
            all_sites = MetalSite.objects.filter(representative=representative, ready_for_presentation=True)
        if ion == None and representative !=None:
            sites = MetalSite.objects.filter(representative=representative, ready_for_presentation=True)
            all_sites = MetalSite.objects.filter(representative=representative, ready_for_presentation=True)
        if ion ==None and representative == None:
            sites = MetalSite.objects.filter(ready_for_presentation=True)
            all_sites = MetalSite.objects.filter(ready_for_presentation=True)


        rest_binding_sites = len(all_sites) - len(sites)
        common_codes = Counter([site.binding_family for site in sites]).most_common(10)
        all_codes = Counter([site.binding_family for site in sites])
        top_ten_code_count = sum([n[1] for n in common_codes])
        return (common_codes, top_ten_code_count, len(sites), rest_binding_sites, all_codes)


    #gets information
    def number_of_PDB_interfaces(representative=None):
        database_count = Pdb.objects.filter(has_metal_interface=True, ready_for_presentation=True).count()
        scraping_date = MetaDataPDB.objects.latest("scraping_date").scraping_date
        newest_record = MetaDataPDB.objects.latest("scraping_date")
        metals_in_pdb =getattr(newest_record, 'metals_in_pdb')
        no_metal = getattr(newest_record, 'no_metal')
        has_metal_no_interface = Pdb.objects.filter(has_metal_interface=False).count()
        Pdbs_with_representative_sites= Pdb.objects.filter(metalsite__representative=True, ready_for_presentation=True).distinct().count()
        update_date = MetaDataPDB.objects.latest("update_date").update_date
        if representative ==True:
            return  [Pdbs_with_representative_sites, database_count-Pdbs_with_representative_sites], scraping_date, update_date
        else:
            return [database_count, has_metal_no_interface, int(no_metal)], scraping_date, update_date





    def make_histogram_data(j, metal=None, element=None, representative=None):
        distances = [a["distance"] for a in j if not metal  or a["metal"] == metal if not element or a["element"] == element]
        bins = {round(x * STEP, 2): [] for x in range(int(3 * (1 / STEP)))}
        for d in distances:
            for cutoff in reversed(sorted(bins.keys())):
                if d > cutoff:
                    bins[cutoff].append(d)
                    break
        bins = {x: len(y) for x, y in bins.items()}
        return zip(*bins.items())

    def get_organism(ion = None, representative=None):
        if ion and representative == None:
            PDBs_with_ion = MetalSite.objects.filter(ion = ion, ready_for_presentation=True).values_list('pdb_id', flat=True)
            pdbs = Pdb.objects.filter(id__in=PDBs_with_ion, ready_for_presentation=True)
        if ion and representative !=None:
            PDBs_with_ion = MetalSite.objects.filter(ion = ion, representative=representative, ready_for_presentation=True).values_list('pdb_id', flat=True)
            pdbs = Pdb.objects.filter(id__in=PDBs_with_ion, ready_for_presentation=True)
        if ion == None and representative !=None:
            PDBs_with_ion = MetalSite.objects.filter(representative=representative, ready_for_presentation=True).values_list('pdb_id', flat=True)
            pdbs = Pdb.objects.filter(id__in=PDBs_with_ion, ready_for_presentation=True)
        if ion ==None and representative == None:
            pdbs = Pdb.objects.filter(has_metal_interface=True, ready_for_presentation=True)

        most_common_organisms = Counter([' '.join(str(protein.organism).upper().split(' ')[0:2]) for protein in pdbs]).most_common(7)
        organism_labels = [n[0] for n in most_common_organisms]
        organism_data =[n[1] for n in most_common_organisms]
        rest_organisms = len(pdbs)-sum(organism_data)
        organism_data.append(rest_organisms)
        organism_labels.append('OTHER')
        return (organism_data, organism_labels)


except:
  print("An exception occurred")
