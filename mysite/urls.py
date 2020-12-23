"""mysite URL Configuration

The `urlpatterns` list routes URLs to views. For more information please see:
    https://docs.djangoproject.com/en/2.2/topics/http/urls/
Examples:
Function views
    1. Add an import:  from my_app import views
    2. Add a URL to urlpatterns:  path('', views.home, name='home')
Class-based views
    1. Add an import:  from other_app.views import Home
    2. Add a URL to urlpatterns:  path('', Home.as_view(), name='home')
Including another URLconf
    1. Import the include() function: from django.urls import include, path
    2. Add a URL to urlpatterns:  path('blog/', include('blog.urls'))
"""
from django.contrib import admin
from django.urls import path, include
#from limb import views # imports views from limb

from django.contrib.sitemaps.views import sitemap
from limb.sitemaps import StaticViewSitemap
from django.contrib.sitemaps import GenericSitemap
# for sitemap
from limb.models import Pdb
from limb. models import MetalSite


Pdb_dict = {
    'queryset': Pdb.objects.filter(ready_for_presentation=True),
}

MetalSite_dict = {
    'queryset': MetalSite.objects.filter(ready_for_presentation=True),
}

sitemaps = {

    'static' : StaticViewSitemap,
    'pdb': GenericSitemap(Pdb_dict, priority=0.5),
    'metalsite': GenericSitemap(MetalSite_dict, priority=0.5),
}

urlpatterns = [
    path('',include('limb.urls')),
    path('sitemap.xml', sitemap, {'sitemaps': sitemaps}, name='django.contrib.sitemaps.views.sitemap'),

]
