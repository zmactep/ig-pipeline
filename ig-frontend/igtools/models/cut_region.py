from igtools.fields import CachedModelChoiceField

__author__ = 'Kos'

from django.db import models
from django import forms
from django.forms.models import model_to_dict
from igstorage.models import StorageItem
import json
import os


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

    def read_params(self, form, index):
        params_map = form.data
        if str(index) + '-input' in params_map:
            self.input = params_map[str(index) + '-input']

        if str(index) + '-kabat' in params_map:
            self.kabat = params_map[str(index) + '-kabat']

        if str(index) + '-reg_num' in params_map:
            self.reg_num = params_map[str(index) + '-reg_num']

        if str(index) + '-group' in params_map:
            self.group = params_map[str(index) + '-group']

        if str(index) + '-comment' in params_map:
            self.comment = params_map[str(index) + '-comment']

    def get_backend_request(self):

        request = {"executable": "ig-snooper/ig_snooper_utils/cut_region.py",
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

        return request

    class Meta:
        app_label = 'igtools'


class CutRegionForm(forms.Form):
    name            = forms.CharField(widget=forms.TextInput(attrs={'class': 'label', 'type': 'text', 'style': 'display: none;'}), initial='CutRegion', label='CutRegion')
    input           = CachedModelChoiceField(widget=forms.Select(attrs={'class': 'form-control', 'type': 'text'}),
                    objects=lambda: {item.path: str(item.file_id + ' - ' + item.group + ' - ' + item.run) for item in StorageItem.objects.using('ig').filter(path__endswith='fasta').order_by('file_id')},
                    label='Входной FASTA файл', empty_label=None, required=True)
    kabat           = CachedModelChoiceField(widget=forms.Select(attrs={'class': 'form-control', 'type': 'text'}),
                    objects=lambda: {item.path: str(item.file_id + ' - ' + item.group + ' - ' + item.run) for item in StorageItem.objects.using('ig').filter(path__endswith='kabat').order_by('file_id')},
                    label='Входной KABAT файл', empty_label=None, required=True)
    reg_num         = forms.IntegerField(widget=forms.TextInput(attrs={'class': 'form-control', 'type': 'text',
                    'data-validation': 'number', 'data-validation-allowing': 'range[0;6]',
                    'data-validation-error-msg': 'Введите, пожалуйста, размер, индекс региона (0-6), который необходимо вырезать. Нумерация с нуля'}),
                    initial=6, label='Номер региона (с нуля)')
    group           = forms.CharField(widget=forms.TextInput(attrs={'class': 'form-control', 'type': 'text',
                    'data-validation': 'length', 'data-validation-length': 'min5',
                    'data-validation-error-msg': 'Введите, пожалуйста, название группы не короче 5 символов'}),
                    initial='testrun', label='Группа')
    comment         = forms.CharField(widget=forms.TextInput(attrs={'class': 'form-control', 'type': 'text',
                    'data-validation': 'length', 'data-validation-length': 'min5',
                    'data-validation-error-msg': 'Введите, пожалуйста, комментарий не короче 5 символов'}),
                    initial='comment', label='Комментарий')

    def set_additional_files(self, fasta, kabat, model):
        new_fasta = {file: str(os.path.basename(file) + ' - pipeline from stage ' + stage + ' - 0') for file, stage in fasta}
        fasta_in_db = {item.path: str(item.file_id + ' - ' + item.group + ' - ' + item.run) for item in StorageItem.objects.using('ig').filter(path__endswith='fasta').order_by('file_id')}
        fasta_in_db.update(new_fasta)
        self.fields['input'].objects = fasta_in_db

        new_kabat = {file: str(os.path.basename(file) + ' - pipeline from stage ' + stage + ' - 0') for file, stage in kabat}
        kabat_in_db = {item.path: str(item.file_id + ' - ' + item.group + ' - ' + item.run) for item in StorageItem.objects.using('ig').filter(path__endswith='kabat').order_by('file_id')}
        kabat_in_db.update(new_kabat)
        self.fields['kabat'].objects = kabat_in_db

    def clean(self):
        try:
            if not ('input' in self.cleaned_data and self.cleaned_data['input']):
                raise forms.ValidationError('input fasta parameter missing')
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