"""Contains instructions  to generate imqages for the publication"""

# %% load functions

import sys
import os
import requests
from collections import Counter
from itertools import groupby
sys.path.append(os.path.join(".."))
os.environ.setdefault("DJANGO_SETTINGS_MODULE", "mysite.settings")
os.environ["DJANGO_ALLOW_ASYNC_UNSAFE"] = "true"
import django
django.setup()
from django.db.models import Count, F
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import fisher_exact
from limb.models import *
from django.shortcuts import render
from django.db.models import Count, F
from limb.models import *
import datetime
import requests
import numpy as np
from collections import Counter
import operator
from django.db.models import Q
### opis
from scipy.stats import fisher_exact
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 6})
from matplotlib.ticker import MaxNLocator
import pandas as pd

from tools import elements, lanthanides, metal_list
#####



font = {'family' : 'DejaVu Sans',
        'weight' : 'normal',
        'size'   : 6}
plt.rc('font', **font)


charge_dict= {'MN3': 'Mn³⁺',
'LI': 'Li⁺',
'NA': 'Na⁺',
'MG': 'Mg²⁺',
'K': 'K⁺',
'CA': 'Ca²⁺',
'RB': 'Rb⁺',
'SR': 'Sr²⁺',
'CS': 'Cs⁺',
'BA': 'Ba²⁺',
'V': 'V³⁺',
'CR': 'Cr³⁺',
'MN': 'Mn²⁺',
'FE': 'Fe³⁺',
'FE2': 'Fe²⁺',
'CO': 'Co²⁺',
'NI': 'Ni²⁺',
'CU': 'Cu²⁺',
'ZN': 'Zn²⁺',
'ZR': 'Zr\u2074\u207A',
'RU': 'Ru³⁺',
'RH': 'Rh⁺',
'PD': 'Pd²⁺',
'AG': 'Ag⁺',
'CD': 'Cd²⁺',
'LA': 'La³⁺',
'W': 'W\u2076\u207A',
'OS': 'Os³⁺',
'IR': 'Ir\u2074\u207A',
'IR3': 'Ir\u00B3\u207A',
'PT': 'Pt²⁺',
'PT4': 'Pt\u2074\u207A',
'AU': 'Au⁺',
'HG': 'Hg²⁺',
'AL': 'Al³⁺',
'GA': 'Ga³⁺',
'IN': 'In³⁺',
'SB': 'Sb³⁺',
'TL': 'Tl⁺',
'PB': 'Pb²⁺',
'CE': 'Ce³⁺',
'PR': 'Pr³⁺',
'SM': 'Sm³⁺',
'EU': 'Eu²⁺',
'TB': 'Tb³⁺',
'DY': 'Dy³⁺',
'YB': 'Yb³⁺',
'LU': 'Lu³⁺',
'TH': 'Th\u2074\u207A',
'AM': 'Am³⁺',
'CF': 'Cf³⁺',
'HO': 'Ho',
'HO3':'Ho\u00B3\u207A',
'TB':'Tb\u00B3\u207A',
'GD':'Gd\u00B3\u207A',
'EU3':'Eu\u00B3\u207A',
'EU2':'Eu\u00B2\u207A',
'BS3':'Bs\u00B3\u207A',
'4TI':'Ti\u2074\u207A',
'BA':'Ba\u00B2\u207A',
'TB':'Tb\u00B3\u207A',
'SM':'Sm\u00B3\u207A',
'LU':'Lu\u00B3\u207A',
'AU3':'Au\u00B3\u207A',
'FE2':'Fe\u00B2\u207A',
'CU1':'Cu\u207A',
 'GD':'Gd',
 'RE':'Re',
 'GD3':'GD³⁺',
 'MO': 'Mo'}

def get_distance(a, b):
    """Returns distances between two given atoms."""
    return np.sqrt((a.x - b.x)**2+(a.y - b.y)**2+(a.z - b.z)**2)
STEP = 0.01

def get_distances_for_sites(representative=None):
    sites = MetalSite.objects.annotate(metals=Count('metal'))
    sites = sites.annotate(resolution=F("pdb_id__resolution"))
    #presentation_ready
    sites = sites.objects.filter(ready_for_presentation=True)
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







### get number of sites, element changed to ion, in order to ommit hemes etc.
def get_number_of_sites(element = None, representative=None):

    if element and representative == None:
        sites = MetalSite.objects.filter(ion=element,ready_for_presentation=True)
    if element and representative != None:
        sites = MetalSite.objects.filter(ion=element, representative=representative,ready_for_presentation=True)
    if element == None and representative !=None:
        sites = MetalSite.objects.filter(representative=representative,ready_for_presentation=True)
    if element ==None and representative == None:
        sites = MetalSite.objects.filter(ready_for_presentation=True)
    all_sites = MetalSite.objects.filter(ready_for_presentation=True)
    #all_element_sites = MetalSite.objects.filter(element=name)
    number_of_sites = len(sites)
    return element, number_of_sites


def DB_comp(**filter_param):
    database_composition = DBcomposition.objects.filter(**filter_param).latest('id')
    return [database_composition.DNA_protein_RNA,
    database_composition.DNA_protein,
    database_composition.DNA_RNA ,
    database_composition.protein_RNA ,
    database_composition.protein_protein,
    database_composition.DNA_DNA  ,
    database_composition.RNA_RNA  ,
    database_composition.other_complexes ,
    database_composition.representative,
    database_composition.element ]

    def number_of_PDB_interfaces(**filter_param):
        database_composition = DBcomposition.objects.filter(**filter_param)

        database_count = Pdb.objects.filter(has_metal_interface=True, ready_for_presentation=True).count()
        scraping_date = MetaDataPDB.objects.latest("scraping_date").scraping_date
        newest_record = MetaDataPDB.objects.latest("scraping_date")
        metals_in_pdb =getattr(newest_record, 'metals_in_pdb')
        no_metal = getattr(newest_record, 'no_metal')
        has_metal_no_interface = Pdb.objects.filter(has_metal_interface=False).count()
        Pdbs_with_representative_sites= Pdb.objects.filter(metalsite__representative=True, ready_for_presentation=True).distinct().count()
        update_date = MetaDataPDB.objects.latest("update_date").update_date
        DBcomposition.objects.latest('id').DNA_protein_RNA
        'DNA-protein-RNA', 'DNA-protein', 'DNA-RNA', 'protein-RNA', 'protein-protein', 'DNA-DNA', 'RNA-RNA', 'other complexes'
        if representative ==True:
            return  [Pdbs_with_representative_sites, database_count-Pdbs_with_representative_sites], scraping_date, update_date
        else:
            return [database_count, has_metal_no_interface, int(no_metal)], scraping_date, update_date
### get coordination_numbers

def get_coordination_numbers(element = None, representative=None):

    if element and representative == None:
        sites = MetalSite.objects.filter(element=element,ready_for_presentation=True)
    if element and representative !=None:
        sites = MetalSite.objects.filter(element=element, representative=representative,ready_for_presentation=True)
    if element == None and representative !=None:
        sites = MetalSite.objects.filter(representative=representative,ready_for_presentation=True)
    if element ==None and representative == None:
        sites = MetalSite.objects.filter(ready_for_presentation=True)
    all_sites = MetalSite.objects.filter(ready_for_presentation=True)
    coord_number = Counter([site.number_of_coordinating_aminoacid_residues for site in sites])
    #all_element_sites = MetalSite.objects.filter(element=name)

    number_of_sites = len(sites)
    return element, sorted(coord_number.items())

### get number of coordinating chains
def get_number_of_coordinating_chains(element = None, representative=None):

    if element and representative == None:
        sites = MetalSite.objects.filter(element=element,ready_for_presentation=True)
    if element and representative !=None:
        sites = MetalSite.objects.filter(element=element, representative=representative,ready_for_presentation=True)
    if element == None and representative !=None:
        sites = MetalSite.objects.filter(representative=representative,ready_for_presentation=True)
    if element ==None and representative == None:
        sites = MetalSite.objects.filter(ready_for_presentation=True)
    all_sites = MetalSite.objects.filter(ready_for_presentation=True)
    coord_number = Counter([site.number_of_coordinating_chains for site in sites])

    #all_element_sites = MetalSite.objects.filter(element=name)
    number_of_sites = len(sites)
    return element, sorted(coord_number.items())



### get residues
def get_residues(element = None, representative = None):
    """Returns residues that are coordinating metals"""
    if element and representative == None:
        MetalSites_with_element = MetalSite.objects.filter(element = element,ready_for_presentation=True).values_list('id', flat=True)
        residues = Residues.objects.filter(site_id__in=MetalSites_with_element)
    if element and representative !=None:
        #Django backward relationship lookup
        MetalSites_with_element = MetalSite.objects.filter(element = element, representative = representative,ready_for_presentation=True).values_list('id', flat=True)
        residues = Residues.objects.filter(site_id__in=MetalSites_with_element)
    if element == None and representative !=None:
        MetalSites_representative= MetalSite.objects.filter( representative = representative,ready_for_presentation=True).values_list('id', flat=True)
        residues = Residues.objects.filter(site_id__in=MetalSites_representative)
    if element ==None and representative == None:
        MetalSites_representative= MetalSite.objects.filter(ready_for_presentation=True).values_list('id', flat=True)
        residues = Residues.objects.filter(site_id__in=MetalSites_representative)
    residues_dict = Counter([str(residue.name).upper() for residue in residues]).most_common(6)
    residues_dict
    residues_data = [n[1] for n in residues_dict]
    other = len(residues)-sum(residues_data)
    if other == 0:
        other = np.nan
    residues_data.append(other)
    residues_labels = [n[0] for n in residues_dict]
    residues_labels.append('OTHER')
    return (element, residues_data, residues_labels)







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



    return {'DPR':DPR, 'DP':DP, 'DR':DR, 'PR':PR, 'PP':PP, 'DD':DD, 'RR':RR, 'other_complexes':other_complexes-DPR-DP-DR-PR-PP-DD-RR }


def autolabel(rects):
    for rect in rects:
        width = rect.get_width()
        ax.annotate('{}'.format(width),
                xy=(1.00*rect.get_width(), rect.get_y()+0.0*rect.get_height()),
                 xytext=(2, 0),  # 1 points vertical offset
                 textcoords="offset points", )




# %% Plot pie graph, structures containing metal interfacial_metal
database_count = Pdb.objects.filter(has_metal_interface=True).count()
scraping_date = MetaDataPDB.objects.latest("scraping_date").scraping_date
newest_record = MetaDataPDB.objects.latest("scraping_date")
metals_in_pdb =getattr(newest_record, 'metals_in_pdb')
no_metal = getattr(newest_record, 'no_metal')
has_metal_no_interface = Pdb.objects.filter(has_metal_interface=False).count()

#plot pie graph
plt.pie([database_count, has_metal_no_interface, no_metal], colors=["#472b52", "#FF6347", "gray"], startangle=90, labels=["PDB with interface", "PDB with metal, but no interface","PDB without metal"])
plt.show()






# %% get chain number and plot representative
#tu zmienic na co sie chce
coord_num = pd.DataFrame()
for metal in elements:
    a= get_number_of_coordinating_chains(metal, representative=True)
    dic = dict(a[1])
    dic[7] =0
    if a[1] != []:
        d= pd.DataFrame.from_dict(dic, orient='index', columns = [metal])
        coord_num = pd.concat((coord_num,d), axis=1)
df = coord_num
df = df.fillna(0)
df = df.astype(int)
fig, axes = plt.subplots(nrows=5, ncols=5, figsize=(7, 6.66))
for idx, (col, ax) in enumerate(zip(df.columns, axes.flatten())):
    bars =ax.barh(df.index, df[col])
    autolabel(bars)
    ax.set_title(col.capitalize())
    x_lims = ax.get_xlim()
    x_lims=((x_lims[0], x_lims[1]*1.3))
    ax.set_xlim(x_lims)
    ax.set_xlabel('Chain count')
    ax.set_ylabel('Count')
    ax.xaxis.set_major_locator(MaxNLocator(4,integer=True))
    plt.subplots_adjust(wspace=.6, hspace=.6)
else:
    [ax.set_visible(False) for ax in axes.flatten()[idx+1:]]

fig.tight_layout()
fig.savefig('data/number_of_chains_representative.tiff',dpi = 300)
fig.clf()
# for lanthanides
coord_num = pd.DataFrame()
for metal in lanthanides:
    a= get_number_of_coordinating_chains(metal, representative=True)
    dic = dict(a[1])
    dic[7] =0
    if a[1] != []:
        d= pd.DataFrame.from_dict(dic, orient='index', columns = [metal])
        coord_num = pd.concat((coord_num,d), axis=1)
df = coord_num
df = df.fillna(0)
df = df.astype(int)
fig, axes = plt.subplots(nrows=3, ncols=3, figsize=(3.66, 3.66))
for idx, (col, ax) in enumerate(zip(df.columns, axes.flatten())):
    bars =ax.barh(df.index, df[col])
    autolabel(bars)
    ax.set_title(col.capitalize())
    x_lims = ax.get_xlim()
    x_lims=((x_lims[0], x_lims[1]*1.3))
    ax.set_xlim(x_lims)
    ax.set_xlabel('Chain count')
    ax.set_ylabel('Count')
    ax.xaxis.set_major_locator(MaxNLocator(4,integer=True))
    plt.subplots_adjust(wspace=.6, hspace=.6)
else:
    [ax.set_visible(False) for ax in axes.flatten()[idx+1:]]

fig.tight_layout()
fig.savefig('data/number_of_chains_representative_lanthanides.tiff',dpi = 300)
fig.clf()








# %%  get_coordination_numbers number and plot

coord_num = pd.DataFrame()
for metal in elements:
    a= get_coordination_numbers(metal, representative=True)
    dic = dict(a[1])
    #dic[7] =0
    if a[1] != []:
        d= pd.DataFrame.from_dict(dic, orient='index', columns = [metal])
        coord_num = pd.concat((coord_num,d), axis=1)
df = coord_num
df = df.fillna(0)
df = df.astype(int)
fig, axes = plt.subplots(nrows=5, ncols=5, figsize=(7, 6.66))
for idx, (col, ax) in enumerate(zip(df.columns, axes.flatten())):
    bars =ax.barh(df.index, df[col])
    autolabel(bars)
    ax.set_title(col.capitalize())
    x_lims = ax.get_xlim()
    x_lims=((x_lims[0], x_lims[1]*1.3))
    ax.set_xlim(x_lims)
    ax.set_xlabel('Donors count')
    ax.set_ylabel('Count')
    ax.xaxis.set_major_locator(MaxNLocator(4,integer=True))
    plt.subplots_adjust(wspace=.6, hspace=.6)
else:
    [ax.set_visible(False) for ax in axes.flatten()[idx+1:]]
fig.tight_layout()
fig.savefig('data/coord_num_representative.tiff',dpi = 300)

coord_num = pd.DataFrame()
for metal in lanthanides:
    a= get_coordination_numbers(metal, representative=True)
    dic = dict(a[1])
    #dic[7] =0
    if a[1] != []:
        d= pd.DataFrame.from_dict(dic, orient='index', columns = [metal])
        coord_num = pd.concat((coord_num,d), axis=1)
df = coord_num
df = df.fillna(0)
df = df.astype(int)
fig, axes = plt.subplots(nrows=3, ncols=3, figsize=(3.66, 3.66))
for idx, (col, ax) in enumerate(zip(df.columns, axes.flatten())):
    bars =ax.barh(df.index, df[col])
    autolabel(bars)
    ax.set_title(col.capitalize())
    x_lims = ax.get_xlim()
    x_lims=((x_lims[0], x_lims[1]*1.3))
    ax.set_xlim(x_lims)
    ax.set_xlabel('Donors count')
    ax.set_ylabel('Count')
    ax.xaxis.set_major_locator(MaxNLocator(4,integer=True))
    plt.subplots_adjust(wspace=.6, hspace=.6)
else:
    [ax.set_visible(False) for ax in axes.flatten()[idx+1:]]
fig.tight_layout()
fig.savefig('data/coord_num_representative_lanthanides.tiff',dpi = 300)






# %% Non representative data
"""Non representative data"""


# %% get chain number and plot representative
#tu zmienic na co sie chce
coord_num = pd.DataFrame()
for metal in elements:
    a= get_number_of_coordinating_chains(metal)
    dic = dict(a[1])
    dic[7] =0
    if a[1] != []:
        d= pd.DataFrame.from_dict(dic, orient='index', columns = [metal])
        coord_num = pd.concat((coord_num,d), axis=1)
df = coord_num
df = df.fillna(0)
df = df.astype(int)
fig, axes = plt.subplots(nrows=5, ncols=5, figsize=(7, 6.66))
for idx, (col, ax) in enumerate(zip(df.columns, axes.flatten())):
    bars =ax.barh(df.index, df[col])
    autolabel(bars)
    ax.set_title(col.capitalize())
    x_lims = ax.get_xlim()
    x_lims=((x_lims[0], x_lims[1]*1.4))
    ax.set_xlim(x_lims)
    ax.set_xlabel('Chain count')
    ax.set_ylabel('Count')
    ax.xaxis.set_major_locator(MaxNLocator(4,integer=True))
    plt.subplots_adjust(wspace=.6, hspace=.6)
else:
    [ax.set_visible(False) for ax in axes.flatten()[idx+1:]]

fig.tight_layout()
fig.savefig('data/number_of_chains.tiff',dpi = 300)

# for lanthanides
coord_num = pd.DataFrame()
for metal in lanthanides:
    a= get_number_of_coordinating_chains(metal)
    dic = dict(a[1])
    dic[7] =0
    if a[1] != []:
        d= pd.DataFrame.from_dict(dic, orient='index', columns = [metal])
        coord_num = pd.concat((coord_num,d), axis=1)
df = coord_num
df = df.fillna(0)
df = df.astype(int)
fig, axes = plt.subplots(nrows=3, ncols=3, figsize=(3.66, 3.66))
for idx, (col, ax) in enumerate(zip(df.columns, axes.flatten())):
    bars =ax.barh(df.index, df[col])
    autolabel(bars)
    ax.set_title(col.capitalize())
    x_lims = ax.get_xlim()
    x_lims=((x_lims[0], x_lims[1]*1.3))
    ax.set_xlim(x_lims)
    ax.set_xlabel('Chain count')
    ax.set_ylabel('Count')
    ax.xaxis.set_major_locator(MaxNLocator(4,integer=True))
    plt.subplots_adjust(wspace=.6, hspace=.6)
else:
    [ax.set_visible(False) for ax in axes.flatten()[idx+1:]]

fig.tight_layout()
fig.savefig('data/number_of_chains_lanthanides.tiff',dpi = 300)









# %%  get_coordination_numbers number and plot

coord_num = pd.DataFrame()
for metal in elements:
    a= get_coordination_numbers(metal)
    dic = dict(a[1])
    #dic[7] =0
    if a[1] != []:
        d= pd.DataFrame.from_dict(dic, orient='index', columns = [metal])
        coord_num = pd.concat((coord_num,d), axis=1)
df = coord_num
df = df.fillna(0)
df = df.astype(int)
fig, axes = plt.subplots(nrows=5, ncols=5, figsize=(7, 6.66))
for idx, (col, ax) in enumerate(zip(df.columns, axes.flatten())):
    bars =ax.barh(df.index, df[col])
    autolabel(bars)
    ax.set_title(col.capitalize())
    x_lims = ax.get_xlim()
    x_lims=((x_lims[0], x_lims[1]*1.3))
    ax.set_xlim(x_lims)
    ax.set_xlabel('Donors count')
    ax.set_ylabel('Count')
    ax.xaxis.set_major_locator(MaxNLocator(4,integer=True))
    plt.subplots_adjust(wspace=.6, hspace=.6)
else:
    [ax.set_visible(False) for ax in axes.flatten()[idx+1:]]
fig.tight_layout()
fig.savefig('data/coord_num.tiff',dpi = 300)

coord_num = pd.DataFrame()
for metal in lanthanides:
    a= get_coordination_numbers(metal)
    dic = dict(a[1])
    #dic[7] =0
    if a[1] != []:
        d= pd.DataFrame.from_dict(dic, orient='index', columns = [metal])
        coord_num = pd.concat((coord_num,d), axis=1)
df = coord_num
df = df.fillna(0)
df = df.astype(int)
fig, axes = plt.subplots(nrows=3, ncols=3, figsize=(3.66, 3.66))
for idx, (col, ax) in enumerate(zip(df.columns, axes.flatten())):
    bars =ax.barh(df.index, df[col])
    autolabel(bars)
    ax.set_title(col.capitalize())
    x_lims = ax.get_xlim()
    x_lims=((x_lims[0], x_lims[1]*1.3))
    ax.set_xlim(x_lims)
    ax.set_xlabel('Donors count')
    ax.set_ylabel('Count')
    ax.xaxis.set_major_locator(MaxNLocator(4,integer=True))
    plt.subplots_adjust(wspace=.6, hspace=.6)
else:
    [ax.set_visible(False) for ax in axes.flatten()[idx+1:]]
fig.tight_layout()
fig.savefig('data/coord_num_lanthanides.tiff',dpi = 300)




# %% Get number of sites
#all sites
get_number_of_sites(representative=True)
#representative Ca2+
get_number_of_sites('CA',representative=True)
# Fe2+ , representative
get_number_of_sites('FE2', True)
# Fe3+ , representative
get_number_of_sites('FE', True)

coord_num = pd.DataFrame()
for metal in metal_list:
    a= get_number_of_sites(metal,True), get_number_of_sites(metal)
    if a[1][1] != 0:
        d= pd.DataFrame([[a[0][1]],[a[1][1]]], columns = [charge_dict[a[0][0]]], index=['representative', 'non-representative'])
        coord_num = pd.concat((coord_num,d), axis=1)
        coord_num.T

coord_num.to_csv('data/number_of_sites.csv',encoding="utf-8-sig")


# %% Get number of residues
residues=get_residues(representative=True)
d = pd.DataFrame(residues[1], columns=['Count'],index = residues[2])
print(d)

residues=get_residues()
d = pd.DataFrame(residues[1], columns=['Count'],index = residues[2])
print(d)


# %% Get number of donors

donors =get_coordination_numbers(representative=True)
non_representative_donors =get_coordination_numbers()



for x in donors[1]:
    print(x[1])
[x[1] for x in donors[1]]
d = pd.DataFrame([[x[1] for x in donors[1]], [x[1] for x in non_representative_donors[1]]], columns = [x[0] for x in donors[1]], index = ['Representative','Non-representative'])
print(d)

#get precentage of cysteines
cys_residues=787
(787/sum(residues[1]))*100




# %% Get number of chains. Get number of sites that are coordinated by x number of chains


chains =get_number_of_coordinating_chains(representative=True)
non_representative_chains =get_number_of_coordinating_chains()

for x in chains[1]:
    print(x[1])
[x[1] for x in chains[1]]
d = pd.DataFrame([[x[1] for x in chains[1]], [x[1] for x in non_representative_chains[1]]], columns = [x[0] for x in chains[1]], index = ['Representative','Non-representative'])
print(d)
