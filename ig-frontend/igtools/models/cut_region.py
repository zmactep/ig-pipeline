__author__ = 'Kos'

from django.db import models
from django import forms
from django.forms.models import model_to_dict
from igtools.models.general import CustomModelChoiceField
from igstorage.models import StorageItem
import json


import logging
log = logging.getLogger('all')


class CutRegion(models.Model):
    backend_id = models.IntegerField(null=True)
    name = models.CharField(max_length=256, null=True, blank=True)

    # files and paths
    input = models.CharField(max_length=256, null=True, blank=True)
    kabat = models.CharField(max_length=256, null=True, blank=True)

    # params
    reg_num = models.IntegerField(null=True, blank=True)

    # metadata
    group = models.CharField(max_length=256)
    comment = models.CharField(max_length=256)

    def __str__(self):
        return json.dumps(model_to_dict(self, fields=[], exclude=[]))

    def read_params(self, params_map):
        if 'input' in params_map:
            self.input = params_map['input'].path

        if 'kabat' in params_map:
            self.kabat = params_map['kabat'].path

        if 'reg_num' in params_map:
            self.reg_num = params_map['reg_num']

        if 'group' in params_map:
            self.group = params_map['group']

        if 'comment' in params_map:
            self.comment = params_map['comment']

    def get_backend_request(self):

        request = {
                  "commands":[
                      {"executable": "ig-snooper/ig_snooper_utils/cut_region.py",
                        "input": {
                           "params": [
                                {"name": "input", "value": self.input},
                                {"name": "kabat", "value": self.kabat},
                                {"name": "reg_num", "value": str(self.reg_num)}
                           ],
                           "comment": self.comment,
                           "group": self.group
                       }
                    }
                  ]}

        return json.dumps(request)

    class Meta:
        app_label = 'igtools'


class CutRegionForm(forms.Form):
    name            = forms.CharField(widget=forms.TextInput(attrs={'class': 'label', 'type': 'text', 'style': 'display: none;'}), initial='CutRegion', label='CutRegion')
    input           = CustomModelChoiceField(queryset=StorageItem.objects.using('ig').filter(path__endswith='fasta').order_by('file_id'), label='Input FASTA file', required=False)
    kabat           = CustomModelChoiceField(queryset=StorageItem.objects.using('ig').filter(path__endswith='kabat').order_by('file_id'), label='Input KABAT file', required=False)
    reg_num         = forms.IntegerField(widget=forms.TextInput(attrs={'class': 'form-control', 'type': 'text'}), initial=6, label='0-based index of region to cut')
    group           = forms.CharField(widget=forms.TextInput(attrs={'class': 'form-control', 'type': 'text'}), initial='testrun', label='Group')
    comment         = forms.CharField(widget=forms.TextInput(attrs={'class': 'form-control', 'type': 'text'}), initial='Some useful comment', label='Your comment')

    def clean(self):
        try:
            if not ('input' in self.cleaned_data and self.cleaned_data['input']):
                raise forms.ValidationError('fasta parameters missing')
            if not ('kabat' in self.cleaned_data and self.cleaned_data['kabat']):
                raise forms.ValidationError('kabat parameters missing')
            if not ('reg_num' in self.cleaned_data):
                raise forms.ValidationError('reg_num parameters missing')
            if not ('group' in self.cleaned_data):
                raise forms.ValidationError('group parameters missing')

        except forms.ValidationError as e:
            log.debug('Error: %s' % e.messages)
            raise
        return self.cleaned_data