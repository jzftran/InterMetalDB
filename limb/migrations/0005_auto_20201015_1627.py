# Generated by Django 3.0.6 on 2020-10-15 14:27

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('limb', '0004_metalsite_ion'),
    ]

    operations = [
        migrations.AddField(
            model_name='metalsite',
            name='ready_for_presentation',
            field=models.BooleanField(default=False),
        ),
        migrations.AddField(
            model_name='pdb',
            name='ready_for_presentation',
            field=models.BooleanField(default=False),
        ),
    ]
