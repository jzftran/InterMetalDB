# Generated by Django 3.0.6 on 2020-10-15 11:15

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('limb', '0003_metadatapdb'),
    ]

    operations = [
        migrations.AddField(
            model_name='metalsite',
            name='ion',
            field=models.CharField(max_length=64, null=True),
        ),
    ]