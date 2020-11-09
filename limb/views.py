"""Generates web pages, contains charge_dict - dictionary
that will allow propper ion display. This should be
implemented in the model, but is not. Making views bigger..."""

from django.contrib.auth.models import User
from django.shortcuts import render, get_object_or_404
from limb.filters import PdbFilter
from limb.filters import MetalSiteFilter
from limb.models import Pdb
from limb.models import MetalSite
from limb.models import Metal
from django.http import HttpResponse
from django_filters.views import FilterView
from django_tables2.views import SingleTableMixin
from limb.tables import PdbTable
from limb.tables import MetalSiteTable
from limb.models import Chain

from django_tables2.export.views import ExportMixin
import django_tables2 as tables

# statistics generation should be ommited in the future and be included in models
from limb.generate_statistics import *

charge_dict = {'AG': 'Ag⁺',
 'AL': 'Al³⁺',
 'AM': 'Am³⁺',
  'AU': 'Au⁺',
 'AU3': 'Au³⁺',
 'BA': 'Ba²⁺',
 'BS3': 'Bs³⁺',
 'CA': 'Ca²⁺',
 'CD': 'Cd²⁺',
 'CE': 'Ce³⁺',
 'CF': 'Cf³⁺',
 'CO': 'Co²⁺',
 'CR': 'Cr³⁺',
 'CS': 'Cs⁺',
 'CU1': 'Cu⁺',
 'CU': 'Cu²⁺',
 'DY': 'Dy³⁺',
 'EU': 'Eu²⁺',
 'EU3': 'Eu³⁺',
 'FE2': 'Fe²⁺',
 'FE': 'Fe³⁺',
 'GA': 'Ga³⁺',
 'GD': 'Gd',
 'GD3':'GD³⁺',
 'HG': 'Hg²⁺',
 'HO': 'Ho',
 'HO3': 'Ho³⁺',
 'IN': 'In³⁺',
 'IR3': 'Ir³⁺',
 'IR': 'Ir⁴⁺',
 'K': 'K⁺',
 'LA': 'La³⁺',
 'LI': 'Li⁺',
 'LU': 'Lu³⁺',
 'MG': 'Mg²⁺',
 'MN': 'Mn²⁺',
 'MO': 'Mo',
 'MN3': 'Mn³⁺',
 'NA': 'Na⁺',
 'NI': 'Ni²⁺',
 'OS': 'Os³⁺',
 'PB': 'Pb²⁺',
 'PD': 'Pd²⁺',
 'PR': 'Pr³⁺',
 'PT': 'Pt²⁺',
 'PT4': 'Pt⁴⁺',
 'RB': 'Rb⁺',
 'RE': 'Re',
 'RH': 'Rh⁺',
 'RU': 'Ru³⁺',
 'SB': 'Sb³⁺',
 'SM': 'Sm³⁺',
 'SR': 'Sr²⁺',
 'TB': 'Tb³⁺',
 'TH': 'Th⁴⁺',
 '4TI': 'Ti⁴⁺',
 'TL': 'Tl⁺',
 'V': 'V³⁺',
 'W': 'W⁶⁺',
 'YB': 'Yb³⁺',
 'ZN': 'Zn²⁺',
 'ZR': 'Zr⁴⁺'}




def index(request):
    """Defines views for index website. Located in templates/limb"""
    return render(request, 'limb/index.html', {})

def about(request):
    """Defines views for about website. Located in templates/limb"""
    return render(request, 'limb/about.html', {})

def references(request):
    """Defines references for this work. Located in templates/limb"""
    return render(request, 'limb/references.html', {})

def PDB_summary(request, id):
    """Defines views for PDB summary website. Located in templates/limb"""
    pdb = get_object_or_404(Pdb, id=id)
    #metal_sites = MetalSite.objects.get(pdb_id=pdb)
    metal_sites = MetalSite.objects.filter(pdb_id=pdb)
    metal = Metal.objects.filter(pdb = pdb)

    return render(request, 'limb/PDB_summary.html', {'pdb': pdb, 'metal_sites' : metal_sites, 'metal':metal})

def metal_site_summary(request, id):
    """Defines views for metal_site summary website. Located in templates/limb"""
    metal_site = get_object_or_404(MetalSite, id=id)
    metal = Metal.objects.filter(site = metal_site)
    g_id= Chain.objects.filter(site_id = metal_site).values_list('group_id', flat=True).distinct()
    sites_with_the_same_group_id= Chain.objects.filter(group_id =g_id[0])
    metal_sites = MetalSite.objects.filter(id__in=sites_with_the_same_group_id.values_list('site_id', flat=True))
    representative_id=metal_sites.filter(pseudo_MFS=metal_site.pseudo_MFS, representative=True)#.values_list('id', flat=True).distinct()
    non_representative_ids=metal_sites.filter(pseudo_MFS=metal_site.pseudo_MFS, representative=False)#.values_list('id', flat=True).distinct()
    return render(request, 'limb/metal_site_summary.html', { 'metal_site' : metal_site, 'metal':metal, 'representative_id': representative_id, 'non_representative_ids':non_representative_ids})







class FilteredPdbListView(FilterView, ExportMixin, tables.SingleTableView):
    #chaining queries
    queryset = Pdb.objects.exclude(has_metal_interface=False).exclude(ready_for_presentation=False) # filtruje przed nałożeniem kolejnego filtru, tu można zmienić wartość bool
    table_class = PdbTable
    template_name = 'search/search_template.html'
    paginate_by = 20
    filterset_class = PdbFilter


class FilteredMetalSiteListView(FilterView, ExportMixin, tables.SingleTableView):
    queryset = MetalSite.objects.exclude(ready_for_presentation=False)
    table_class = MetalSiteTable
    template_name = 'search/search_template.html'
    paginate_by = 20
    filterset_class = MetalSiteFilter





def generate_element_graphs(request, name=None, representative=None):
    if representative == 'representative':
        representative = True

    metals = Metal.objects.values('residue_name').distinct()

    common_codes , top_ten_element_count, all_element_count, rest_binding_sites, all_codes = get_family_info(name, representative=representative)
    abundance_in_interface = [all_element_count, rest_binding_sites]
    abundance_in_interface_labels = [name, 'other metals']


    if representative == True:
        x, y = make_histogram_data(distance[1], metal = name)
        xs, ys = make_histogram_data(distance[1], metal = name, element="S")
        xn, yn = make_histogram_data(distance[1], metal = name, element="N")
        xo, yo = make_histogram_data(distance[1], metal = name, element="O")
        data, scraping_date, update_date = number_of_PDB_interfaces(representative=True)
        labels = ['Representative PDB with metal interface', 'Non-representative PDB with metal interface']


    if representative == None:
        x, y = make_histogram_data(distance[0], metal = name)
        xs, ys = make_histogram_data(distance[0], metal = name, element="S")
        xn, yn = make_histogram_data(distance[0], metal = name, element="N")
        xo, yo = make_histogram_data(distance[0], metal = name, element="O")
        data, scraping_date, update_date = number_of_PDB_interfaces()
        labels = ['PDB with metal interface','PDB with metal, but no metal interface', 'PDB without metal']

    labels_histo = list(x)
    bond_length = list(y)

    labels_s = list(xs)
    labels_n = list(xn)
    labels_o = list(xo)

    bond_length_s= list(ys)
    bond_length_n = list(yn)
    bond_length_o = list(yo)

    family_names = [n[0] for n in common_codes]
    family_names.append('OTHER')

    #gets the number of rest binding_sites
    other = all_element_count - top_ten_element_count
    family_count =[n[1] for n in common_codes]
    family_count.append(other)

    residues_data , residues_labels = get_residues(name,representative=representative)
    techniqes_data, techniqes_labels = get_techniques(name,representative=representative)
    classification_data, classification_labels = get_classification(name,representative=representative)
    organism_data, organism_labels = get_organism(name,representative=representative)
    year_data, year_labels = get_year(name,representative=representative)

    if name:
        charge = charge_dict[name]
        return render(request, 'limb/statistics_element.html', {

            'data':data,
            'labels':labels,
            'labels_histo': labels_histo,
            'bond_length':bond_length,
            'bond_length_s':bond_length_s,
            'bond_length_n':bond_length_n,
            'bond_length_o':bond_length_o,
            'labels_s' : labels_s,
            'labels_n' : labels_n,
            'labels_o' : labels_o,
            'family_names' : family_names,
            'family_count' : family_count,
            'abundance_in_interface' : abundance_in_interface,
            'abundance_in_interface_labels': abundance_in_interface_labels,
            'residues_data' : residues_data,
            'residues_labels': residues_labels,
            'techniqes_data': techniqes_data,
            'techniqes_labels': techniqes_labels,
            'classification_data': classification_data,
            'classification_labels': classification_labels,
            'organism_data': organism_data,
            'organism_labels': organism_labels,
            'name':name.capitalize(),
            'charge': charge,
            'year_data': year_data,
            'year_labels': year_labels,
            'all_codes': all_codes,
            'representative':representative,
            'update_date': update_date,
            'metals':metals,
            'charge_dict': charge_dict,
            'representative':representative,
            })
    else:
        charge = None
        return render(request, 'limb/statistics.html', {

            'data':data,
            'labels':labels,
            'labels_histo': labels_histo,
            'bond_length':bond_length,
            'bond_length_s':bond_length_s,
            'bond_length_n':bond_length_n,
            'bond_length_o':bond_length_o,
            'labels_s' : labels_s,
            'labels_n' : labels_n,
            'labels_o' : labels_o,
            'family_names' : family_names,
            'family_count' : family_count,
            'abundance_in_interface' : abundance_in_interface,
            'abundance_in_interface_labels': abundance_in_interface_labels,
            'residues_data' : residues_data,
            'residues_labels': residues_labels,
            'techniqes_data': techniqes_data,
            'techniqes_labels': techniqes_labels,
            'classification_data': classification_data,
            'classification_labels': classification_labels,
            'organism_data': organism_data,
            'organism_labels': organism_labels,
            'name':name,
            'charge': charge,
            'metals':metals,
            'year_data': year_data,
            'year_labels': year_labels,
            'charge_dict': charge_dict,
            'all_codes': all_codes,
            'representative':representative,
            'update_date': update_date
            })
