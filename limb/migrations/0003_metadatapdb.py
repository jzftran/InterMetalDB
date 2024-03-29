# Generated by Django 3.0.6 on 2020-10-13 09:31

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('limb', '0002_delete_metadatapdb'),
    ]

    operations = [
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
    ]
