"""Clusters sequences and metal sites itself.
Should be run after each database update.
Can be run as an independent script."""

from tools import *
env_django()
import django; django.setup()
from django.db import transaction
from limb.models import Pdb
from limb.models import MetalSite
from limb.models import Chain
from django.conf import settings
from fill_models import *
from django.db.models import F
import subprocess
import glob, os, shutil, sys
from tqdm import tqdm
from datetime import datetime


# makes a log
make_log(f"Clustering has started at:{datetime.now()}")
def chains_to_fasta():
    """Gets chains sequence from the database.
    Returns lines that could be written to a file.
    Might not work well with a large dataset!"""
    lines =[]
    for chain in Chain.objects.all():
        lines.append(">"+str(chain.id))
        sequence = chain.sequence
        while sequence:
            lines.append(sequence[:80])
            sequence = sequence[80:]
    return "\n".join(lines)
fasta = chains_to_fasta()

with open('chains.fasta','w') as f:
    f.write(fasta)
    f.close()

#describes path to mmseqs

#mmseqs = r"/Users/jozeftran/mmseqs/bin/mmseqs"
#mmseqs = "~/mmseqs/bin/mmseqs"
#%%
import subprocess
mmseqs = r"C:\Users\jozef\Documents\MMseqs2-Windows-Unified\mmseqs\bin\mmseqs.exe"
run = subprocess.call(mmseqs + " createdb chains.fasta DB",shell=True, stdout=subprocess.PIPE)
print(run)
#%%
#identity treshold 0.5, coverage 0.8, coverage mode 0 (query and targer)
subprocess.call(mmseqs + " cluster DB DB_clu tmp --min-seq-id 0.5 -c 0.8 --cov-mode 0", shell=True)
#output as tsv
#%%
subprocess.call(mmseqs + " createtsv DB DB DB_clu results.tsv", shell=True, stdout=subprocess.PIPE)
clusters={}
#%%

#cluster chains in DB
with open('results.tsv') as f:
    print('Saving representative chains...')
    for line in tqdm(f):
        cluster_representative, cluster_member = line.strip().split('\t')
        chain = Chain.objects.get(id=cluster_member)
        chain.group_id = cluster_representative
        chain.save()

# removes clustering from metalsites
chains = Chain.objects.all()
metal_sites = MetalSite.objects.all()
metal_sites.update(representative = False)



group_ids = Chain.objects.all().values_list('group_id', flat=True).distinct()
print('Saving representative metal-sites')
for group_id in tqdm(group_ids):
    grouped_sites = Chain.objects.filter(group_id=group_id).values_list('site', flat=True).distinct()
    pseudo_MFS= MetalSite.objects.filter(id__in=grouped_sites, representative = None).values_list('pseudo_MFS',flat=True).distinct()
    #annotates metalsites based on pdb structure resolution, if no resolution information is present such site wont be representative
    annotated_MetalSites = MetalSite.objects.filter(id__in=grouped_sites).annotate(resolution=F("pdb_id__resolution")).order_by(F('resolution').asc(nulls_last=True))

    fiels_values= set()
    grouped = dict()
    for site in  grouped_sites:
        obj = MetalSite.objects.get(id=site)
        field_value = getattr(obj, 'pseudo_MFS')
        fiels_values.add(field_value)
    for binding_family in fiels_values:
        grouped[binding_family] =[]
        grouped[binding_family] += (annotated_MetalSites.filter( pseudo_MFS=binding_family))
    for key,value in grouped.items():
        value[0].representative = True
        value[0].save()



make_log(f"Clustering has ended at:{datetime.now()}")
#remove temporary files
try:
    for f in glob.glob("DB*"):
        os.remove(f)
    os.remove('results.tsv')
    shutil.rmtree('tmp')
except:
    pass
