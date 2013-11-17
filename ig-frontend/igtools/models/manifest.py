__author__ = 'Kos'
from django.db import models


class Manifest(models.Model):
    tool_name = models.CharField(max_length=256, null=True, blank=True)
    tool_path = models.CharField(max_length=256, null=True, blank=True)
    manifest = models.CharField(max_length=4096, null=True, blank=True)

    class Meta:
        app_label = 'igtools'