# Generated by Django 3.0.6 on 2020-10-13 09:26

from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    initial = True

    dependencies = [
    ]

    operations = [
        migrations.CreateModel(
            name='Atom',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('atomium_id', models.IntegerField()),
                ('name', models.CharField(max_length=32)),
                ('x', models.FloatField()),
                ('y', models.FloatField()),
                ('z', models.FloatField()),
                ('element', models.CharField(max_length=8)),
            ],
            options={
                'db_table': 'atoms',
            },
        ),
        migrations.CreateModel(
            name='MetaDataPDB',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('metals_in_pdb', models.IntegerField()),
                ('no_metal', models.IntegerField()),
                ('scraping_date', models.DateField(blank=True, null=True)),
                ('update_date', models.DateField(blank=True, null=True)),
            ],
            options={
                'db_table': 'MetaDataPDB',
            },
        ),
        migrations.CreateModel(
            name='MetalSite',
            fields=[
                ('id', models.CharField(max_length=128, primary_key=True, serialize=False)),
                ('element', models.CharField(max_length=64, null=True)),
                ('binding_family', models.CharField(max_length=128)),
                ('pseudo_MFS', models.CharField(max_length=128, null=True)),
                ('aminoacid_residue_names', models.CharField(max_length=512, verbose_name='amino acids or nucleotide residues names')),
                ('all_residue_names', models.CharField(max_length=512, null=True, verbose_name='all bound residues')),
                ('is_homomer', models.BooleanField(null=True)),
                ('residue_ids', models.CharField(max_length=512, null=True)),
                ('all_residue_ids', models.CharField(max_length=512, null=True)),
                ('number_of_coordinating_aminoacid_residues', models.IntegerField(blank=True, null=True, verbose_name='number of bound amino acids or nucleotide residues')),
                ('number_of_coordinating_chains', models.IntegerField(blank=True, null=True)),
                ('number_of_coordinating_residues', models.IntegerField(blank=True, null=True, verbose_name='number of all bound residues')),
                ('number_of_nonprotein_residues', models.IntegerField(blank=True, null=True)),
                ('representative', models.BooleanField(default=False)),
            ],
            options={
                'db_table': 'metal_sites',
            },
        ),
        migrations.CreateModel(
            name='Pdb',
            fields=[
                ('id', models.CharField(max_length=32, primary_key=True, serialize=False)),
                ('title', models.CharField(max_length=1024)),
                ('classification', models.CharField(blank=True, max_length=1024, null=True)),
                ('keywords', models.CharField(blank=True, max_length=2048, null=True)),
                ('deposition_date', models.DateField(blank=True, null=True)),
                ('resolution', models.FloatField(blank=True, null=True)),
                ('rvalue', models.FloatField(blank=True, null=True)),
                ('organism', models.CharField(blank=True, max_length=1024, null=True)),
                ('expression_system', models.CharField(blank=True, max_length=1024, null=True)),
                ('technique', models.CharField(blank=True, max_length=1024, null=True)),
                ('assembly', models.IntegerField(blank=True, null=True)),
                ('has_metal_interface', models.BooleanField(default=False)),
            ],
            options={
                'db_table': 'PDB_codes_db',
            },
        ),
        migrations.CreateModel(
            name='Residues',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('name', models.CharField(max_length=32)),
                ('pdb', models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.CASCADE, to='limb.Pdb')),
                ('site', models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.CASCADE, to='limb.MetalSite')),
            ],
            options={
                'db_table': 'residues',
            },
        ),
        migrations.AddField(
            model_name='metalsite',
            name='pdb_id',
            field=models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='limb.Pdb'),
        ),
        migrations.CreateModel(
            name='Metal',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('atomium_id', models.IntegerField()),
                ('name', models.CharField(max_length=32)),
                ('x', models.FloatField()),
                ('y', models.FloatField()),
                ('z', models.FloatField()),
                ('element', models.CharField(max_length=8)),
                ('residue_name', models.CharField(max_length=32)),
                ('residue_number', models.IntegerField()),
                ('chain_id', models.CharField(max_length=128)),
                ('omission_reason', models.TextField(blank=True, null=True)),
                ('pdb', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='limb.Pdb')),
                ('site', models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.CASCADE, to='limb.MetalSite')),
            ],
            options={
                'db_table': 'metals',
            },
        ),
        migrations.CreateModel(
            name='Chain',
            fields=[
                ('id', models.CharField(max_length=128, primary_key=True, serialize=False)),
                ('chain_id', models.CharField(blank=True, max_length=128, null=True)),
                ('sequence', models.TextField()),
                ('group_id', models.CharField(blank=True, max_length=128, null=True)),
                ('pdb', models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.CASCADE, to='limb.Pdb')),
                ('site', models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.CASCADE, to='limb.MetalSite')),
            ],
            options={
                'db_table': 'chains',
            },
        ),
        migrations.CreateModel(
            name='Bonds',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('atom', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='limb.Atom')),
                ('metal', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='limb.Metal')),
            ],
            options={
                'db_table': 'bonds',
            },
        ),
    ]
