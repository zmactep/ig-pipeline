from django.db import models
from django import forms
from django.forms.models import model_to_dict
from igtools.models.general import CustomModelChoiceField
from igstorage.models import StorageItem
import json


import logging
log = logging.getLogger('all')


class SimpleCluster(models.Model):
    backend_id = models.IntegerField(null=True)
    name = models.CharField(max_length=256, null=True, blank=True)

    # files and paths
    src = models.CharField(max_length=256, null=True, blank=True)

    # params
    sim = models.CharField(max_length=256, null=True, blank=True)
    shortest_cons = models.CharField(max_length=256, null=True, blank=True)
    skip_first = models.IntegerField(null=True, blank=True)
    min_len = models.IntegerField(null=True, blank=True)
    use_prct = models.IntegerField(null=True, blank=True)

    # metadata
    group = models.CharField(max_length=256)
    comment = models.CharField(max_length=256)

    def __str__(self):
        return json.dumps(model_to_dict(self, fields=[], exclude=[]))

    def read_params(self, params_map):
        if 'src' in params_map:
            self.src = params_map['src'].path

        if 'sim' in params_map:
            self.sim = params_map['sim']

        if 'skip_first' in params_map:
            self.skip_first = params_map['skip_first']

        if 'min_len' in params_map:
            self.min_len = params_map['min_len']

        if 'use_prct' in params_map:
            self.use_prct = params_map['use_prct']

        if 'shortest_cons' in params_map:
            self.shortest_cons = params_map['shortest_cons']

        if 'group' in params_map:
            self.group = params_map['group']

        if 'comment' in params_map:
            self.comment = params_map['comment']

    def get_backend_request(self):
        params = [{"name": "src", "value": self.src}]
        if self.sim:
            params += [{"name": "sim", "value": self.sim}]

        if self.skip_first:
            params += [{"name": "skip-first", "value": str(self.skip_first)}]

        if self.min_len:
            params += [{"name": "min-len", "value": str(self.min_len)}]

        if self.use_prct:
            params += [{"name": "use-prct", "value": str(self.use_prct)}]

        if self.shortest_cons:
            params += [{"name": "shortest-cons", "value": str(self.shortest_cons)}]

        request = {"executable": "ig-simplecluster/clusterize.py",
                        "input": {
                           "params": params,
                           "comment": self.comment,
                           "group": self.group
                       }
                    }

        return request

    class Meta:
        app_label = 'igtools'


class SimpleClusterForm(forms.Form):
    name            = forms.CharField(widget=forms.TextInput(attrs={'class': 'label', 'type': 'text', 'style': 'display: none;'}), initial='SimpleCluster', label='SimpleCluster')
    src             = CustomModelChoiceField(queryset=StorageItem.objects.using('ig').filter(path__endswith='fasta').order_by('file_id'), label='Input FASTA file', required=False)
    sim             = forms.CharField(widget=forms.TextInput(attrs={'class': 'form-control', 'type': 'text'}), label='Sim', required=False)
    shortest_cons   = forms.CharField(widget=forms.TextInput(attrs={'class': 'form-control', 'type': 'text'}), label='Shortest Cons', required=False)
    skip_first      = forms.IntegerField(widget=forms.TextInput(attrs={'class': 'form-control', 'type': 'text'}), label='Skip First', required=False)
    min_len         = forms.IntegerField(widget=forms.TextInput(attrs={'class': 'form-control', 'type': 'text'}), label='Minimal Length', required=False)
    use_prct        = forms.IntegerField(widget=forms.TextInput(attrs={'class': 'form-control', 'type': 'text'}), label='Use %', required=False)
    group           = forms.CharField(widget=forms.TextInput(attrs={'class': 'form-control', 'type': 'text'}), initial='testrun', label='Group')
    comment         = forms.CharField(widget=forms.TextInput(attrs={'class': 'form-control', 'type': 'text'}), initial='Some useful comment', label='Your comment', required=False)

    def clean(self):
        try:
            if not ('src' in self.cleaned_data and self.cleaned_data['src']):
                raise forms.ValidationError('src parameters missing')

        except forms.ValidationError as e:
            log.debug('Error: %s' % e.messages)
            raise
        return self.cleaned_data