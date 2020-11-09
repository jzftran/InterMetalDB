"""Defines urls and paths for website"""

from django.conf.urls import url
from limb import views
from django_filters.views import FilterView
from limb.filters import PdbFilter
from limb.filters import MetalSiteFilter
from limb.models import Pdb
from limb.models import MetalSite
from django.urls import path
from limb.views import FilteredPdbListView
from limb.views import FilteredMetalSiteListView




urlpatterns = [
    path('', views.index, name='index'),
    url(r'^search/$', FilteredPdbListView.as_view(), name="search"),
    url(r'^search_metal_site/$', FilteredMetalSiteListView.as_view(), name="searchMetalSite"),
    path('PDB_summary/<str:id>/', views.PDB_summary, name='PDB_summary'),
    path('metal_site_summary/<str:id>/', views.metal_site_summary, name='metal_site_summary'),
    path('statistics/', views.generate_element_graphs, name='statistics'),
    path('statistics/representative', views.generate_element_graphs, {'representative':'representative'}, name='representative'),
    url(r'^statistics/(?P<name>\w+)/(?P<representative>\w+)/$',views.generate_element_graphs, {'representative':'representative'}, name='representative_element'),
    path('statistics/<str:name>/', views.generate_element_graphs, name='generate_element_graphs'),
    path('about', views.about, name='about'),
    path('references', views.references, name='references'),
]
