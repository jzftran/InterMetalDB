
from django.contrib.sitemaps import Sitemap
from django.shortcuts import reverse
from .models import Pdb
class StaticViewSitemap(Sitemap):
    def items(self):
        return ['index',
        'search',
        'searchMetalSite',
        'statistics',
        'references',
        ]
    def location(self, item):
        return reverse(item)
