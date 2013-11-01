from django.db import models
from django.forms.models import model_to_dict
from django import forms

import json
import logging
log = logging.getLogger('all')


class StorageItem(models.Model):
    file_id = models.CharField(max_length=256, null=True, blank=True)
    comment = models.CharField(max_length=256, null=True, blank=True)
    path = models.CharField(max_length=256, null=True, blank=True)

    def __str__(self):
        return json.dumps(model_to_dict(self, fields=[], exclude=[]))


class StorageItemForm(forms.Form):
    file_id = forms.CharField(widget=forms.TextInput(attrs={'class': 'form-control', 'type': 'text'}), label='File ID', required=True)
    comment = forms.CharField(widget=forms.TextInput(attrs={'class': 'form-control', 'type': 'text'}), label='Comment', required=True)
    path = forms.CharField(widget=forms.TextInput(attrs={'class': 'form-control', 'type': 'text'}), label='Path', required=True)

