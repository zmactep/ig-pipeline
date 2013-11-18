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

    def read_params(self, form, index):
        params_map = form.data
        if str(index) + '-ldir' in params_map:
            self.ldir = params_map[str(index) + '-ldir']

        if str(index) + '-hdir' in params_map:
            self.hdir = params_map[str(index) + '-hdir']

        if str(index) + '-group' in params_map:
            self.group = params_map[str(index) + '-group']

        if str(index) + '-comment' in params_map:
            self.comment = params_map[str(index) + '-comment']

    def get_backend_request(self):

        request = {"executable": "ig-simplecluster/make_report.py",
                        "input": {
                           "params": [
                                {"name": "ldir", "value": self.ldir},
                                {"name": "hdir", "value": self.hdir}
                           ],
                           "comment": self.comment,
                           "group": self.group
                       }
                    }

        return request

    class Meta:
        app_label = 'igtools'


class ReportForm(forms.Form):
    name            = forms.CharField(widget=forms.TextInput(attrs={'class': 'label', 'type': 'text', 'style': 'display: none;'}), initial='Report', label='Report')
    ldir             = CachedModelChoiceField(widget=forms.Select(attrs={'class': 'form-control', 'type': 'text'}),
                    label='Директория легкой цепи', empty_label=None, required=True)
    hdir             = CachedModelChoiceField(widget=forms.Select(attrs={'class': 'form-control', 'type': 'text'}),
                    label='Директория тяжелой цепи', empty_label=None, required=True)
    group           = forms.CharField(widget=forms.TextInput(attrs={'class': 'form-control', 'type': 'text',
                    'data-validation': 'length', 'data-validation-length': 'min5',
                    'data-validation-error-msg': 'Введите, пожалуйста, название группы не короче 5 символов'}),
                    initial='testrun', label='Группа', required=True)
    comment         = forms.CharField(widget=forms.TextInput(attrs={'class': 'form-control', 'type': 'text',
                    'data-validation': 'length', 'data-validation-length': 'min5',
                    'data-validation-error-msg': 'Введите, пожалуйста, комментарий не короче 5 символов'}),
                    initial='comment', label='Комментарий', required=True)

    def set_additional_files(self, pipelined_files_map):
        dir = pipelined_files_map['dir']
        new_dir = {file: str(os.path.basename(file) + ' - pipeline from stage ' + stage + ' - 0') for file, stage in dir}
        dir_in_db = {item.path: str(item.file_id + ' - ' + item.group + ' - ' + item.run) for item in StorageItem.objects.using('ig').filter(path__endswith='/').order_by('file_id')}
        dir_in_db.update(new_dir)
        self.fields['ldir'].objects = new_dir
        self.fields['hdir'].objects = new_dir

    def clean(self):
        try:
            if not ('ldir' in self.cleaned_data and self.cleaned_data['ldir']):
                raise forms.ValidationError('ldir parameters missing')
            if not ('hdir' in self.cleaned_data and self.cleaned_data['hdir']):
                raise forms.ValidationError('hdir parameters missing')
            if not ('group' in self.cleaned_data):
                raise forms.ValidationError('group parameters missing')
            if not ('comment' in self.cleaned_data):
                raise forms.ValidationError('comment parameters missing')

        except forms.ValidationError as e:
            log.debug('Error: %s' % e.messages)
            raise
        return self.cleaned_data