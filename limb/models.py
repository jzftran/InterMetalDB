from django.db import models
from django.shortcuts import reverse

class Pdb(models.Model):

    """General information about PDB contained in database."""


    class Meta:
        db_table = "PDB_codes_db" #name of the table

    id = models.CharField(primary_key=True, max_length=32)
    title = models.CharField(max_length=1024)
    classification = models.CharField(null=True, blank=True, max_length=1024)
    keywords = models.CharField(null=True, blank=True, max_length=2048)
    deposition_date = models.DateField(null=True, blank=True)
    resolution = models.FloatField(null=True, blank=True)
    rvalue = models.FloatField(null=True, blank=True)
    organism = models.CharField(null=True, blank=True, max_length=1024)
    expression_system = models.CharField(null=True, blank=True, max_length=1024)
    technique = models.CharField(null=True, blank=True, max_length=1024)
    assembly = models.IntegerField(null=True, blank=True)
    has_metal_interface = models.BooleanField(default = False, blank=False)
    ready_for_presentation = models.BooleanField(default = False, blank=False) #data for the database has been downloaded for chemical identifier, so it is not a full picture of whole RCSB PDB database, even though some of the record contain not only standalone ions, this ensures that in the dabase are presented only standalone ions

    def get_absolute_url(self):
        return reverse('PDB_summary', kwargs={'id':self.id})


class MetalSite(models.Model):

    """Metal binding site record."""

    class Meta:
        db_table = "metal_sites"

    id = models.CharField(primary_key=True, max_length=128)
    element = models.CharField(max_length=64, null =True)
    ion  = models.CharField(max_length=64, null =True)
    binding_family = models.CharField(max_length=128)
    pseudo_MFS = models.CharField(max_length=128, null =True)
    aminoacid_residue_names = models.CharField(max_length=512,verbose_name='amino acids or nucleotide residues names')
    all_residue_names = models.CharField(max_length=512, null=True, verbose_name='all bound residues')
    pdb_id = models.ForeignKey(Pdb, on_delete=models.CASCADE) # "on_delete=models.CASCADE" removes remaining stuff after itself. ForeignKey(Pdb links to Pdb table
    is_homomer = models.BooleanField(null=True)
    residue_ids =models.CharField(max_length=512, null=True)
    all_residue_ids =models.CharField(max_length=512, null=True)
    number_of_coordinating_aminoacid_residues = models.IntegerField(null=True, blank=True, verbose_name= 'number of bound amino acids or nucleotide residues')
    number_of_coordinating_chains = models.IntegerField(null=True, blank=True)
    number_of_coordinating_residues = models.IntegerField(null=True, blank=True, verbose_name= 'number of all bound residues')
    number_of_nonprotein_residues = models.IntegerField(null=True, blank=True)
    representative = models.BooleanField(default = False, blank = True, null =True)
    ready_for_presentation = models.BooleanField(default = False, blank=False) #data for the database has been downloaded for chemical identifier, so it is not a full picture of whole RCSB PDB database, even though some of the record contain not only standalone ions, this ensures that in the dabase are presented only standalone ions

    def get_absolute_url(self):
        return reverse('metal_site_summary', kwargs={'id':self.id})

class Metal(models.Model):
    """A metal atom model."""

    class Meta:
        db_table = "metals"

    atomium_id = models.IntegerField()
    name = models.CharField(max_length=32)
    x = models.FloatField()
    y = models.FloatField()
    z = models.FloatField()
    element = models.CharField(max_length=8)
    residue_name = models.CharField(max_length=32)
    residue_number = models.IntegerField()
    chain_id = models.CharField(max_length=128)
    omission_reason = models.TextField(blank=True, null=True)
    pdb = models.ForeignKey(Pdb, on_delete=models.CASCADE)
    site = models.ForeignKey(MetalSite, on_delete=models.CASCADE, blank=True, null=True)


class Chain(models.Model):
    """Chains coordinating metal."""

    class Meta:
        db_table = "chains"

    id = models.CharField(primary_key=True, max_length=128)
    chain_id = models.CharField(max_length=128,blank=True, null=True)
    sequence = models.TextField()
    pdb = models.ForeignKey(Pdb, on_delete=models.CASCADE, blank=True, null=True)
    group_id = models.CharField(max_length=128,blank=True, null=True)
    site = models.ForeignKey(MetalSite, on_delete=models.CASCADE,blank=True, null=True)

class Atom(models.Model):
    """Atoms coordinating metals."""
    class Meta:
        db_table = "atoms"
    atomium_id = models.IntegerField()
    name = models.CharField(max_length=32)
    x = models.FloatField()
    y = models.FloatField()
    z = models.FloatField()
    element = models.CharField(max_length=8)



class Residues(models.Model):
    """Residues coordinating metals."""
    class Meta:
        db_table = "residues"
    name = models.CharField(max_length=32)
    pdb = models.ForeignKey(Pdb, on_delete=models.CASCADE, blank=True, null=True)
    site = models.ForeignKey(MetalSite, on_delete=models.CASCADE,blank=True, null=True)





class Bonds(models.Model):
    """Bonds between atoms and metals."""
    class Meta:
        db_table =  "bonds"
    metal = models.ForeignKey(Metal, on_delete=models.CASCADE)
    atom = models.ForeignKey(Atom, on_delete=models.CASCADE)

class MetaDataPDB(models.Model):
    """Data scraped from RCSB."""
    class Meta:
        db_table =  "MetaDataPDB"

    metals_in_pdb = models.IntegerField()
    no_metal = models.IntegerField()
    scraping_date = models.DateField(null=True, blank=True)
    update_date = models.DateField(null=True, blank=True)


class DBcomposition(models.Model):
    """DB composition"""
    class Meta:
        db_table =  "DBcomposition"

    DNA_protein_RNA = models.IntegerField()
    DNA_protein = models.IntegerField()
    DNA_RNA = models.IntegerField()
    protein_RNA = models.IntegerField()
    protein_protein = models.IntegerField()
    DNA_DNA = models.IntegerField()
    RNA_RNA = models.IntegerField()
    other_complexes = models.IntegerField()
    representative = models.BooleanField(null=True, blank=True)
    element = models.CharField(max_length=64, null =True)
