""" Contains filters for  search page. """
from django.db import models
from django.contrib.auth.models import User
from .models import Pdb
from .models import MetalSite
import django_filters
from django_filters.filterset import FilterSet
from django_filters.filters import ModelChoiceFilter
from django.forms.widgets import TextInput






class MyRangeWidget(django_filters.widgets.RangeWidget):

    def __init__(self, from_attrs=None, to_attrs=None, attrs=None):
        super(MyRangeWidget, self).__init__(attrs)
        if from_attrs:
            self.widgets[0].attrs.update(from_attrs)
        if to_attrs:
            self.widgets[1].attrs.update(to_attrs)









class PdbFilter(django_filters.FilterSet):
    """This class defines how to filter Pdb model."""
    id = django_filters.CharFilter(lookup_expr='iexact', label='PDB ID', widget=TextInput(attrs={'placeholder': 'e.g. 1AYM'}))
    title = django_filters.CharFilter(lookup_expr='icontains', label='PDB title', widget=TextInput(attrs={'placeholder': 'e.g. human insulin'}))
    organism = django_filters.CharFilter(lookup_expr='icontains', label ='Gene source organism', widget=TextInput(attrs={'placeholder': 'e.g. Helicobacter pylori'}))
    technique =django_filters.CharFilter(lookup_expr='icontains', label ='Technique contains', widget=TextInput(attrs={'placeholder': 'e.g. X-RAY'}))
    classification = django_filters.CharFilter(lookup_expr='icontains',label= 'PDB classification', widget=TextInput(attrs={'placeholder': 'e.g. hormone'}))
    keywords = django_filters.CharFilter(lookup_expr='icontains', label = 'keywords', widget=TextInput(attrs={'placeholder': 'e.g. metal transport'}))
    resolution = django_filters.RangeFilter(
        label='Resolution [Å]',
        widget=MyRangeWidget(
            from_attrs={'placeholder':'from'},
            to_attrs={'placeholder':'to'},
        )
    )
    deposition_date = django_filters.DateFromToRangeFilter(
        label='Deposition date',
        widget=MyRangeWidget(
            from_attrs={'placeholder':'from (YYYY-MM-DD)'},
            to_attrs={'placeholder':'to (YYYY-MM-DD)'},
        )
    )

    #deposition_date = django_filters.DateFromToRangeFilter(field_name='aa',widget=RangeWidget(attrs={'placeholder': 'YYYY/MM/DD'}))
    #is_homomer = django_filters.BooleanFilter()

    class Meta:
        model = Pdb
        fields = ['id', 'title','organism','technique', 'classification', 'keywords','resolution','deposition_date']



class MetalSiteFilter(django_filters.FilterSet):
    """This class defines how to filter MetalSite model."""
    id = django_filters.CharFilter(lookup_expr='icontains', label='Metal site id', widget=TextInput(attrs={'placeholder': 'e.g. 1AYM-ZN-2 or 1AYM'}))
    element = django_filters.CharFilter(lookup_expr='icontains', label='Metal element', widget=TextInput(attrs={'placeholder': 'e.g. Cd'}))
    number_of_coordinating_chains = django_filters.CharFilter(lookup_expr='iexact', label='Number of bound chains', widget=TextInput(attrs={'placeholder': 'e.g. 3'}))

    number_of_coordinating_aminoacid_residues = django_filters.CharFilter(lookup_expr='iexact', label='Number of bound amino acids or nucleotide residues', widget=TextInput(attrs={'placeholder': 'e.g. 4'}))
    number_of_coordinating_residues = django_filters.CharFilter(lookup_expr='iexact', label='Number of all bound ligands', widget=TextInput(attrs={'placeholder': 'e.g. 4'}))
    binding_family = django_filters.CharFilter(lookup_expr='iexact', label='Bound amino acid or nucleotide residues (is exact)', widget=TextInput(attrs={'placeholder': 'e.g. C3H1, E2H2'}))
    aminoacid_residue_names = django_filters.CharFilter(lookup_expr='icontains', label='Amino acid or nucleotide residue names', widget=TextInput(attrs={'placeholder': 'e.g. CYS.HIS, CYS or ASP'}))
    all_residue_names = django_filters.CharFilter(lookup_expr='icontains', label='All bound ligands names', widget=TextInput(attrs={'placeholder': 'e.g. CL.HIS, DG, HIS or CL'}))
    #is_homomer = django_filters.BooleanFilter(lookup_expr='isnull')
    class Meta:
        model = MetalSite
        fields = ['id', 'element', 'number_of_coordinating_chains','number_of_coordinating_aminoacid_residues', 'number_of_coordinating_residues', 'binding_family', 'binding_family', 'aminoacid_residue_names', 'all_residue_names']
