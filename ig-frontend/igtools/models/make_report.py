__author__ = 'Kos'

from django.db import models
from django import forms
from django.forms.models import model_to_dict
from igtools.models.general import CustomModelChoiceField
from igstorage.models import StorageItem
import json


import logging
log = logging.getLogger('all')


class Report(models.Model):
    backend_id = models.IntegerField(null=True)
    name = models.CharField(max_length=256, null=True, blank=True)

    # files and paths
    ldir = models.CharField(max_length=256, null=True, blank=True)
    hdir = models.CharField(max_length=256, null=True, blank=True)

    # params
    fix_suffix = models.CharField(max_length=256)

    # metadata
    group = models.CharField(max_length=256)
    comment = models.CharField(max_length=256)

    def __str__(self):
        return json.dumps(model_to_dict(self, fields=[], exclude=[]))

    def read_params(self, params_map):
        if 'ldir' in params_map:
            self.ldir = params_map['ldir'].path

        if 'hdir' in params_map:
            self.hdir = params_map['hdir'].path

        if 'group' in params_map:
            self.group = params_map['group']

        if 'comment' in params_map:
            self.comment = params_map['comment']

    def get_backend_request(self):

        request = {
                  "commands":[
                      {"executable": "ig-simplecluster/make_report.py",
                        "input": {
                           "params": [
                                {"name": "ldir", "value": self.ldir},
                                {"name": "hdir", "value": self.hdir}
                           ],
                           "comment": self.comment,
                           "group": self.group
                       }
                    }
                  ]}

        return json.dumps(request)

    class Meta:
        app_label = 'igtools'


class ReportForm(forms.Form):
    name            = forms.CharField(widget=forms.TextInput(attrs={'class': 'label', 'type': 'text', 'style': 'display: none;'}), initial='Report', label='Report')
    ldir            = CustomModelChoiceField(queryset=StorageItem.objects.using('ig').filter(path__endswith='/').order_by('file_id'), label='light chain dir', required=True)
    hdir            = CustomModelChoiceField(queryset=StorageItem.objects.using('ig').filter(path__endswith='/').order_by('file_id'), label='heavy chain dir', required=True)
    group           = forms.CharField(widget=forms.TextInput(attrs={'class': 'form-control', 'type': 'text'}), initial='testrun', label='Group')
    comment         = forms.CharField(widget=forms.TextInput(attrs={'class': 'form-control', 'type': 'text'}), initial='Some useful comment', label='Your comment')

    def clean(self):
        try:
            if not ('ldir' in self.cleaned_data and self.cleaned_data['ldir']):
                raise forms.ValidationError('ldir parameters missing')
            if not ('hdir' in self.cleaned_data and self.cleaned_data['hdir']):
                raise forms.ValidationError('hdir parameters missing')
            if not ('group' in self.cleaned_data):
                raise forms.ValidationError('group parameters missing')

        except forms.ValidationError as e:
            log.debug('Error: %s' % e.messages)
            raise
        return self.cleaned_data