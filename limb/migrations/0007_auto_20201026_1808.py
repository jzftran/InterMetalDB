# Generated by Django 3.0.6 on 2020-10-26 17:08

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('limb', '0006_auto_20201026_1805'),
    ]

    operations = [
        migrations.AlterField(
            model_name='metalsite',
            name='representative',
            field=models.BooleanField(blank=True, default=False, null=True),
        ),
    ]
