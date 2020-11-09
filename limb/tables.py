import django_tables2 as tables
from .models import Pdb
from .models import MetalSite
from django_tables2.utils import A
from django_tables2 import LinkColumn
from tqdm import tqdm
import traceback
from django.db import transaction
from limb.models import Pdb
from limb.models import MetalSite
from django.conf import settings
import django_filters

charge_dict = {'MN3': 'Mn³⁺',
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
  'GD3':'GD³⁺'}


#class TruncatedTextColumn(tables.Column):
#    '''A Column to limit to 100 characters and add an ellipsis'''
#    def render(self, value):
#        if len(value) > 102:
#            return value[0:99] + '...'
#        return str(value)

class PdbTable(tables.Table):
    # text=lambda record: record.id pozwala na zamianę tekstu na określony id
    id = tables.LinkColumn('PDB_summary',text=lambda record: record.id, args=[A('id')],
    attrs={
        "td": {
            "data-length": lambda value: len(value)}
        })#, attrs={"td": {"class": "text-overflow: ellipsis"}})
    class Meta:
        model = Pdb
        template_name = 'django_tables2/bootstrap-responsive.html'
        #id = tables.LinkColumn('PDB', text='static text', args=[A('pk')])

        fields = ('id','title','classification','keywords','deposition_date',
        'resolution','rvalue','organism','expression_system','technique','assembly',
        )
        attrs = {"class": "table-condensed responsive-table striped "}


class MetalSiteTable(tables.Table):


    id = tables.LinkColumn('metal_site_summary',text=lambda record: record.id, args=[A('id')])
    #id = tables.Column('metal_site_summary', record.id, args=[A('id')])

    class Meta:

        model = MetalSite
        template_name = 'django_tables2/bootstrap.html'
        fields = ('id','element','ion','binding_family','aminoacid_residue_names','all_residue_names','is_homomer','number_of_coordinating_aminoacid_residues','number_of_coordinating_residues','number_of_coordinating_chains','representative')
        attrs = {"class": "responsive-table striped"}

    # mofify ion display to follow
    def render_ion(self,value):
        return charge_dict[value]
