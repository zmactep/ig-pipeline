from django.db import models
from django import forms
from django.forms.models import model_to_dict
from igtools.forms import CachedModelChoiceField
from igstorage.models import StorageItem
import json
import os

import logging
log = logging.getLogger('all')


class Train(models.Model):
    backend_id = models.IntegerField(null=True)
    name = models.CharField(max_length=256, null=True, blank=True)

    # files and paths
    fasta = models.CharField(max_length=256, null=True, blank=True)
    kabat = models.CharField(max_length=256, null=True, blank=True)

    # algo params
    ml_window_size = models.IntegerField(null=True, blank=True)

    # metadata
    group = models.CharField(max_length=256)
    comment = models.CharField(max_length=256)

    def read_params(self, form, index):
        params_map = form.data
        if str(index) + '-fasta' in params_map:
            self.fasta = params_map[str(index) + '-fasta']

        if str(index) + '-kabat' in params_map:
            self.kabat = params_map[str(index) + '-kabat']

        if str(index) + '-ml_window_size' in params_map:
            self.ml_window_size = params_map[str(index) + '-ml_window_size']

        if str(index) + '-group' in params_map:
            self.group = params_map[str(index) + '-group']

        if str(index) + '-comment' in params_map:
            self.comment = params_map[str(index) + '-comment']

    def __str__(self):
        return json.dumps(model_to_dict(self, fields=[], exclude=[]))

    def get_backend_request(self):
        request = {"executable": "ig-snooper/train.py",
                        "input": {
                           "params": [
                                {"name": "fasta", "value": self.fasta},
                                {"name": "kabat", "value": self.kabat},
                                {"name": "ml_window_size", "value": str(self.ml_window_size)}
                           ],
                           "comment": self.comment,
                           "group": self.group
                       }
                    }

        return request

    class Meta:
        app_label = 'igtools'


class TrainForm(forms.Form):
    name            = forms.CharField(widget=forms.TextInput(attrs={'class': 'form-control', 'type': 'text', 'style': 'display: none;'}), initial='Train', label='Train')
    fasta           = CachedModelChoiceField(widget=forms.Select(attrs={'class': 'form-control', 'type': 'text'}),
                    objects=lambda: {item.path: str(item.file_id + ' - ' + item.group + ' - ' + item.run) for item in StorageItem.objects.using('ig').filter(path__endswith='fasta').order_by('file_id')},
                    label='Входной FASTA файл', empty_label=None, required=True)
    kabat           = CachedModelChoiceField(widget=forms.Select(attrs={'class': 'form-control', 'type': 'text'}),
                    objects=lambda: {item.path: str(item.file_id + ' - ' + item.group + ' - ' + item.run) for item in StorageItem.objects.using('ig').filter(path__endswith='kabat').order_by('file_id')},
                    label='Входной KABAT файл', empty_label=None, required=True)
    group           = forms.CharField(widget=forms.TextInput(attrs={'class': 'form-control', 'type': 'text',
                    'data-validation': 'length', 'data-validation-length': 'min5',
                    'data-validation-error-msg': 'Введите, пожалуйста, название группы не короче 5 символов'}),
                    initial='testrun', label='Группа')
    comment         = forms.CharField(widget=forms.TextInput(attrs={'class': 'form-control', 'type': 'text',
                    'data-validation': 'length', 'data-validation-length': 'min5',
                    'data-validation-error-msg': 'Введите, пожалуйста, комментарий не короче 5 символов'}),
                    initial='comment', label='Комментарий')
    ml_window_size  = forms.IntegerField(widget=forms.TextInput(attrs={'class': 'form-control', 'type': 'text',
                    'data-validation': 'number',
                    'data-validation-error-msg': 'Введите, пожалуйста, размер, на который будут разделены риды'}),
                    initial=13, label='Размер окна для ML')

    def set_additional_files(self, fasta, kabat, model):
        new_fasta = {file: str(os.path.basename(file) + ' - pipeline from stage ' + stage + ' - 0') for file, stage in fasta}
        fasta_in_db = {item.path: str(item.file_id + ' - ' + item.group + ' - ' + item.run) for item in StorageItem.objects.using('ig').filter(path__endswith='fasta').order_by('file_id')}
        fasta_in_db.update(new_fasta)
        self.fields['fasta'].objects = fasta_in_db

        new_kabat = {file: str(os.path.basename(file) + ' - pipeline from stage ' + stage + ' - 0') for file, stage in kabat}
        kabat_in_db = {item.path: str(item.file_id + ' - ' + item.group + ' - ' + item.run) for item in StorageItem.objects.using('ig').filter(path__endswith='kabat').order_by('file_id')}
        kabat_in_db.update(new_kabat)
        self.fields['kabat'].objects = kabat_in_db

    def clean(self):
        try:
            if not ('fasta' in self.cleaned_data and self.cleaned_data['fasta']):
                raise forms.ValidationError('fasta parameters missing')
            if not ('kabat' in self.cleaned_data and self.cleaned_data['kabat']):
                raise forms.ValidationError('kabat parameters missing')
            if not ('ml_window_size' in self.cleaned_data):
                raise forms.ValidationError('ml_window_size parameters missing')
            if not ('group' in self.cleaned_data):
                raise forms.ValidationError('group parameters missing')

        except forms.ValidationError as e:
            log.debug('Error in Train: %s' % e.messages)
            raise
        return self.cleaned_data