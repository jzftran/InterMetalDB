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


# export statistics
import csv

from django.http import StreamingHttpResponse

class Echo:
    """An object that implements just the write method of the file-like
    interface.
    """
    def write(self, value):
        """Write the value by returning it, instead of storing in a buffer."""
        return value


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
        DB_composition = DB_comp(representative =True)
        if name:
            DB_composition = DB_comp(representative =True, element=name)
        else:
            DB_composition = DB_comp(representative =True, element='All')
    if representative == None:
        x, y = make_histogram_data(distance[0], metal = name)
        xs, ys = make_histogram_data(distance[0], metal = name, element="S")
        xn, yn = make_histogram_data(distance[0], metal = name, element="N")
        xo, yo = make_histogram_data(distance[0], metal = name, element="O")
        data, scraping_date, update_date = number_of_PDB_interfaces()
        labels = ['PDB with metal interface','PDB with metal, but no metal interface', 'PDB without metal']
        DB_composition = DB_comp()
        if name:
            DB_composition = DB_comp(element=name)
        else:
            DB_composition = DB_comp(element='All')

    DB_composition_labels = ['DNA-protein-RNA', 'DNA-protein', 'DNA-RNA', 'protein-RNA', 'protein-protein', 'DNA-DNA', 'RNA-RNA', 'other complexes']

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
            'DB_composition':DB_composition,
            'DB_composition_labels':DB_composition_labels
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
            'update_date': update_date,
            'DB_composition':DB_composition,
            'DB_composition_labels':DB_composition_labels
            })


def some_streaming_csv_view(request, name=None, representative=None):
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
        DB_composition = DB_comp(representative =True)
        if name:
            DB_composition = DB_comp(representative =True, element=name)
        else:
            DB_composition = DB_comp(representative =True, element='All')
    if representative == None:
        x, y = make_histogram_data(distance[0], metal = name)
        xs, ys = make_histogram_data(distance[0], metal = name, element="S")
        xn, yn = make_histogram_data(distance[0], metal = name, element="N")
        xo, yo = make_histogram_data(distance[0], metal = name, element="O")
        data, scraping_date, update_date = number_of_PDB_interfaces()
        labels = ['PDB with metal interface','PDB with metal, but no metal interface', 'PDB without metal']
        DB_composition = DB_comp()
        if name:
            DB_composition = DB_comp(element=name)
        else:
            DB_composition = DB_comp(element='All')

    DB_composition_labels = ['DNA-protein-RNA', 'DNA-protein', 'DNA-RNA', 'protein-RNA', 'protein-protein', 'DNA-DNA', 'RNA-RNA', 'other complexes']

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

    """A view that streams a large CSV file."""
    # Generate a sequence of rows. The range is based on the maximum number of
    # rows that can be handled by a single sheet in most spreadsheet
    # applications.
    if name:
        rows = (
        [f"PDBs containing {name} and other metals"],
        abundance_in_interface_labels,
        abundance_in_interface,

        [f"Structures containing {name} at interface"],
        year_labels,
        year_data,

        [f"Bond lengts for {name} and liganding atom pairs"],
        ["All bonds"],
        labels_histo,
        bond_length,

        ["S bond lengths"],
        labels_s,
        bond_length_s,
        ["N bond lengths"],
        labels_n,
        bond_length_n,
        ["O bond lengths"],
        labels_o,
        bond_length_o,
        [f"{name} bound residues in metal-involved interfaces"],
        family_names,
        family_count,

        [f"Residues coordinating {name} at interfaces"],
        residues_labels,
        residues_data,

        [f"Classification of PDBs containing {name}-involved interface"],
        classification_labels,
        classification_data,

        [f"Gene source for structures containing {name}-involved interface"],
        organism_labels,
        organism_data,

        ["Structure acquisition methods"],
        techniqes_labels,
        techniqes_data,

        ["Complex type in metal binding"],
        DB_composition_labels,
        DB_composition
        )
    else:
        rows = (
        ["Interfaces in RCSB PDB"],
        labels, data,
        ["Structures containing metals at interface"],
        year_labels,
        year_data,
        ["Bond lengts for metals and liganding atom pairs"],
        ["All bonds"],
        labels_histo,
        bond_length,
        ["S bond lengths"],
        labels_s,
        bond_length_s,
        ["N bond lengths"],
        labels_n,
        bond_length_n,
        ["O bond lengths"],
        labels_o,
        bond_length_o,
        ["Bound residues in metal-involved interfaces"],
        family_names,
        family_count,

        ["Residues coordinating metals at interfaces"],
        residues_labels,
        residues_data,

        ["Classification of PDBs containing metal-involved interface"],
        classification_labels,
        classification_data,

        ["Gene source for structures containing metal-involved interface"],
        organism_labels,
        organism_data,

        ["Structure acquisition methods"],
        techniqes_labels,
        techniqes_data,

        ["Complex type in metal binding"],
        DB_composition_labels,
        DB_composition
        )
    pseudo_buffer = Echo()
    writer = csv.writer(pseudo_buffer)
    response = StreamingHttpResponse((writer.writerow(row) for row in rows), content_type="text/csv")
    response['Content-Disposition'] = 'attachment; filename="statistics.csv"'
    return response
