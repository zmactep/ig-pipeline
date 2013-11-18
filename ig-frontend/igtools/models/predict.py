from django.db import models
from django import forms
from django.forms.models import model_to_dict
from igstorage.models import StorageItem
from igtools.forms import CachedModelChoiceField
import json
import os


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

    def read_params(self, form, index):
        params_map = form.data
        if str(index) + '-fasta' in params_map:
            self.fasta = params_map[str(index) + '-fasta']

        if str(index) + '-model_path' in params_map:
            self.model_path = params_map[str(index) + '-model_path']

        if str(index) + '-merge_threshold' in params_map:
            self.merge_threshold = params_map[str(index) + '-merge_threshold']

        if str(index) + '-avg_window_size' in params_map:
            self.avg_window_size = params_map[str(index) + '-avg_window_size']

        if str(index) + '-ml_window_size' in params_map:
            self.ml_window_size = params_map[str(index) + '-ml_window_size']

        if str(index) + '-group' in params_map:
            self.group = params_map[str(index) + '-group']

        if str(index) + '-comment' in params_map:
            self.comment = params_map[str(index) + '-comment']

    def __str__(self):
        return json.dumps(model_to_dict(self, fields=[], exclude=[]))

    def get_backend_request(self):
        request = {"executable": "ig-snooper/predict.py",
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

        return request

    class Meta:
        app_label = 'igtools'


class PredictForm(forms.Form):
    name            = forms.CharField(widget=forms.TextInput(attrs={'class': 'label', 'type': 'text', 'style': 'display: none;'}), initial='Predict', label='Predict')
    fasta           = CachedModelChoiceField(widget=forms.Select(attrs={'class': 'form-control', 'type': 'text'}),
                        label='Входной FASTA файл', empty_label=None, required=True)
    model_path      = CachedModelChoiceField(widget=forms.Select(attrs={'class': 'form-control', 'type': 'text'}),
                        label='Файл с моделью', empty_label=None, required=True)
    group           = forms.CharField(widget=forms.TextInput(attrs={'class': 'form-control', 'type': 'text',
                        'data-validation': 'length', 'data-validation-length': 'min5',
                        'data-validation-error-msg': 'Введите, пожалуйста, название группы не короче 5 символов'}),
                        initial='testrun', label='Группа', required=True)
    comment         = forms.CharField(widget=forms.TextInput(attrs={'class': 'form-control', 'type': 'text',
                        'data-validation': 'length', 'data-validation-length': 'min5',
                        'data-validation-error-msg': 'Введите, пожалуйста, комментарий не короче 5 символов'}),
                        initial='comment', label='Комментарий', required=True)
    ml_window_size  = forms.IntegerField(widget=forms.TextInput(attrs={'class': 'form-control', 'type': 'text',
                    'data-validation': 'number',
                    'data-validation-error-msg': 'Введите, пожалуйста, размер, на который будут разделены риды'}),
                    initial=13, label='Размер окна для ML', required=True)
    merge_threshold = forms.IntegerField(widget=forms.TextInput(attrs={'class': 'form-control', 'type': 'text',
                    'data-validation': 'number',
                    'data-validation-error-msg': 'Введите, пожалуйста, порог скленивания коротких участков'}),
                    initial=10, label='Порог скленивания коротких участков', required=True)
    avg_window_size = forms.IntegerField(widget=forms.TextInput(attrs={'class': 'form-control', 'type': 'text',
                    'data-validation': 'number',
                    'data-validation-error-msg': 'Введите, пожалуйста, размер усредняющего окна'}),
                    initial=10, label='Размер усредняющего окна', required=True)

    def set_additional_files(self, pipelined_files_map):
        fasta = pipelined_files_map['fasta']
        new_fasta = {file: str(os.path.basename(file) + ' - pipeline from stage ' + stage + ' - 0') for file, stage in fasta}
        fasta_in_db = {item.path: str(item.file_id + ' - ' + item.group + ' - ' + item.run) for item in StorageItem.objects.using('ig').filter(path__endswith='fasta').order_by('file_id')}
        fasta_in_db.update(new_fasta)
        self.fields['fasta'].objects = fasta_in_db

        model = pipelined_files_map['model']
        new_model = {file: str(os.path.basename(file) + ' - pipeline from stage ' + stage + ' - 0') for file, stage in model}
        model_in_db = {item.path: str(item.file_id + ' - ' + item.group + ' - ' + item.run) for item in StorageItem.objects.using('ig').filter(path__endswith='model').order_by('file_id')}
        model_in_db.update(new_model)
        self.fields['model_path'].objects = model_in_db

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
            if not ('comment' in self.cleaned_data):
                raise forms.ValidationError('group parameters missing')

        except forms.ValidationError as e:
            log.debug('Error in Predict: %s' % e.messages)
            raise
        return self.cleaned_data