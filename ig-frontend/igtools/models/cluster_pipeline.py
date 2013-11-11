__author__ = 'Kos'

from django.db import models
from django import forms
from django.forms.models import model_to_dict
from igtools.models.general import CustomModelChoiceField
from igstorage.models import StorageItem
import json


import logging
log = logging.getLogger('all')


class ClusterPipeline(models.Model):
    backend_id = models.IntegerField(null=True)
    name = models.CharField(max_length=256, null=True, blank=True)

    # files and paths
    lfasta = models.CharField(max_length=256, null=True, blank=True)
    hfasta = models.CharField(max_length=256, null=True, blank=True)

    lkabat = models.CharField(max_length=256, null=True, blank=True)
    hkabat = models.CharField(max_length=256, null=True, blank=True)

    # metadata
    group = models.CharField(max_length=256)
    comment = models.CharField(max_length=256)

    def __str__(self):
        return json.dumps(model_to_dict(self, fields=[], exclude=[]))

    def read_params(self, params_map):
        if 'lfasta' in params_map:
            self.lfasta = params_map['lfasta'].path

        if 'hfasta' in params_map:
            self.hfasta = params_map['hfasta'].path

        if 'lkabat' in params_map:
            self.lkabat = params_map['lkabat'].path

        if 'hkabat' in params_map:
            self.hkabat = params_map['hkabat'].path

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
                                {"name": "input", "value": self.lfasta},
                                {"name": "kabat", "value": self.lkabat},
                                {"name": "reg_num", "value": "6"}
                           ],
                           "comment": self.comment,
                           "group": self.group
                       }
                      },
                      {"executable": "ig-snooper/ig_snooper_utils/cut_region.py",
                        "input": {
                           "params": [
                                {"name": "input", "value": self.hfasta},
                                {"name": "kabat", "value": self.hkabat},
                                {"name": "reg_num", "value": "6"}
                           ],
                           "comment": self.comment,
                           "group": self.group
                       }
                      },
                      {"executable": "ig-simplecluster/clusterize.py",
                        "input": {
                           "params":[
                               {"name": "src", "value": "{0}/output.fasta"}
                           ],
                           "comment": self.comment,
                           "group": self.group
                        }
                       },
                       {"executable": "ig-simplecluster/clusterize.py",
                        "input": {
                           "params":[
                               {"name": "src", "value": "{1}/output.fasta"}
                           ],
                           "comment": self.comment,
                           "group": self.group
                        }
                      },
                      {"executable": "ig-simplecluster/make_report.py",
                        "input": {
                           "params": [
                                {"name": "ldir", "value": "{2}"},
                                {"name": "hdir", "value": "{3}"}
                           ],
                           "comment": self.comment,
                           "group": self.group
                       }
                      }
                  ]}

        return json.dumps(request)

    class Meta:
        app_label = 'igtools'


class ClusterPipelineForm(forms.Form):
    name            = forms.CharField(widget=forms.TextInput(attrs={'class': 'label', 'type': 'text', 'style': 'display: none;'}), initial='ClusterPipeline', label='ClusterPipeline')
    lfasta          = CustomModelChoiceField(queryset=StorageItem.objects.using('ig').filter(path__endswith='fasta').order_by('file_id'), label='light chain fasta', required=True)
    hfasta          = CustomModelChoiceField(queryset=StorageItem.objects.using('ig').filter(path__endswith='fasta').order_by('file_id'), label='heavy chain fasta', required=True)
    lkabat          = CustomModelChoiceField(queryset=StorageItem.objects.using('ig').filter(path__endswith='kabat').order_by('file_id'), label='light chain kabat', required=True)
    hkabat          = CustomModelChoiceField(queryset=StorageItem.objects.using('ig').filter(path__endswith='kabat').order_by('file_id'), label='heavy chain kabat', required=True)
    group           = forms.CharField(widget=forms.TextInput(attrs={'class': 'form-control', 'type': 'text'}), initial='testrun', label='Group')
    comment         = forms.CharField(widget=forms.TextInput(attrs={'class': 'form-control', 'type': 'text'}), initial='Some useful comment', label='Your comment')

    def clean(self):
        try:
            if not ('lfasta' in self.cleaned_data and self.cleaned_data['lfasta']):
                raise forms.ValidationError('lfasta parameters missing')
            if not ('hfasta' in self.cleaned_data and self.cleaned_data['hfasta']):
                raise forms.ValidationError('hfasta parameters missing')
            if not ('lkabat' in self.cleaned_data and self.cleaned_data['lkabat']):
                raise forms.ValidationError('lkabat parameters missing')
            if not ('hkabat' in self.cleaned_data and self.cleaned_data['hkabat']):
                raise forms.ValidationError('hkabat parameters missing')
            if not ('group' in self.cleaned_data):
                raise forms.ValidationError('group parameters missing')

        except forms.ValidationError as e:
            log.debug('Error: %s' % e.messages)
            raise
        return self.cleaned_data