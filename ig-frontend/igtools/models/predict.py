from django.db import models
from django import forms
from django.forms.models import model_to_dict
from igtools.models.general import CustomModelChoiceField
from igstorage.models import StorageItem
import json


import logging
log = logging.getLogger('all')


class Predict(models.Model):
    backend_id = models.IntegerField(null=True)
    name = models.CharField(max_length=256, null=True, blank=True)

    # files and paths
    fasta = models.CharField(max_length=256, null=True, blank=True)
    model_path = models.CharField(max_length=256, null=True, blank=True)

    # algo params
    ml_window_size = models.IntegerField(null=True, blank=True)
    avg_window_size = models.IntegerField(null=True, blank=True)
    merge_threshold = models.IntegerField(null=True, blank=True)

    # metadata
    group = models.CharField(max_length=256)
    comment = models.CharField(max_length=256)

    def __str__(self):
        return json.dumps(model_to_dict(self, fields=[], exclude=[]))

    def get_backend_request(self):
        request = {
                  "commands":[
                      {"executable": "ig-snooper/predict.py",
                        "input": {
                           "params": [
                                {"name": "fasta", "value": self.fasta},
                                {"name": "model_path", "value": self.model_path},
                                {"name": "ml_window_size", "value": str(self.ml_window_size)},
                                {"name": "merge_threshold", "value": str(self.merge_threshold)},
                                {"name": "avg_window_size", "value": str(self.avg_window_size)}
                           ],
                           "comment": self.comment,
                           "group": self.group
                       }
                    }
                  ]}

        return json.dumps(request)

    class Meta:
        app_label = 'igtools'


class PredictForm(forms.Form):
    name            = forms.CharField(widget=forms.TextInput(attrs={'class': 'label', 'type': 'text', 'style': 'display: none;'}), initial='Predict', label='Predict')
    fasta           = CustomModelChoiceField(queryset=StorageItem.objects.using('ig').filter(path__endswith='fasta').order_by('file_id'), label='Input FASTA file', required=False)
    model_path      = CustomModelChoiceField(queryset=StorageItem.objects.using('ig').filter(path__endswith='model').order_by('file_id'), label='Model file', required=False)
    group           = forms.CharField(widget=forms.TextInput(attrs={'class': 'form-control', 'type': 'text'}), initial='testrun', label='Group')
    comment         = forms.CharField(widget=forms.TextInput(attrs={'class': 'form-control', 'type': 'text'}), initial='Some useful comment', label='Your comment')
    ml_window_size  = forms.IntegerField(widget=forms.TextInput(attrs={'class': 'form-control', 'type': 'text'}), initial=13, label='Machine learning window size')
    merge_threshold = forms.IntegerField(widget=forms.TextInput(attrs={'class': 'form-control', 'type': 'text'}), initial=10, label='Merge threshold window size')
    avg_window_size = forms.IntegerField(widget=forms.TextInput(attrs={'class': 'form-control', 'type': 'text'}), initial=10, label='Average window size')

    def clean(self):
        try:
            if not ('fasta' in self.cleaned_data and self.cleaned_data['fasta']):
                raise forms.ValidationError('fasta parameters missing')
            if not ('model_path' in self.cleaned_data and self.cleaned_data['model_path']):
                raise forms.ValidationError('model_path parameters missing')
            if not ('ml_window_size' in self.cleaned_data):
                raise forms.ValidationError('ml_window_size parameters missing')
            if not ('merge_threshold' in self.cleaned_data):
                raise forms.ValidationError('merge_threshold parameters missing')
            if not ('avg_window_size' in self.cleaned_data):
                raise forms.ValidationError('avg_window_size parameters missing')
            if not ('group' in self.cleaned_data):
                raise forms.ValidationError('group parameters missing')

        except forms.ValidationError as e:
            log.debug('Error: %s' % e.messages)
            raise
        return self.cleaned_data