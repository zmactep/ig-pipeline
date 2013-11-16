from django.db import models
from django import forms
from django.forms.models import model_to_dict
from igtools.models.general import CustomModelChoiceField
from igstorage.models import StorageItem
import json


import logging
log = logging.getLogger('all')


class Train(models.Model):
    backend_id = models.IntegerField(null=True)
    name = models.CharField(max_length=256, null=True, blank=True)

    # files and paths
    fasta = models.CharField(max_length=256, null=True, blank=True)
    kabat = models.CharField(max_length=256, null=True, blank=True)
    model_name = models.CharField(max_length=256, null=True, blank=True)

    # algo params
    ml_window_size = models.IntegerField(null=True, blank=True)

    # metadata
    group = models.CharField(max_length=256)
    comment = models.CharField(max_length=256)

    def read_params(self, params_map):
        if 'fasta' in params_map:
            self.fasta = params_map['fasta'].path

        if 'kabat' in params_map:
            self.kabat = params_map['kabat'].path

        if 'model_name' in params_map:
            self.model_name = params_map['model_name']

        if 'ml_window_size' in params_map:
            self.ml_window_size = params_map['ml_window_size']

        if 'group' in params_map:
            self.group = params_map['group']

        if 'comment' in params_map:
            self.comment = params_map['comment']

    def __str__(self):
        return json.dumps(model_to_dict(self, fields=[], exclude=[]))

    def get_backend_request(self):
        request = {
                  "commands":[
                      {"executable": "ig-snooper/train.py",
                        "input": {
                           "params": [
                                {"name": "fasta", "value": self.fasta},
                                {"name": "kabat", "value": self.kabat},
                                {"name": "model_name", "value": self.model_name},
                                {"name": "ml_window_size", "value": str(self.ml_window_size)}
                           ],
                           "comment": self.comment,
                           "group": self.group
                       }
                    }
                  ]}

        return json.dumps(request)

    class Meta:
        app_label = 'igtools'


class TrainForm(forms.Form):
    name            = forms.CharField(widget=forms.TextInput(attrs={'class': 'form-control', 'type': 'text', 'style': 'display: none;'}), initial='Train', label='Train')
    fasta           = CustomModelChoiceField(widget=forms.Select(attrs={'class': 'form-control', 'type': 'text'}),
                    queryset=StorageItem.objects.using('ig').filter(path__endswith='fasta').order_by('file_id'),
                    label='Входной FASTA файл', empty_label=None, required=True)
    kabat           = CustomModelChoiceField(widget=forms.Select(attrs={'class': 'form-control', 'type': 'text'}),
                    queryset=StorageItem.objects.using('ig').filter(path__endswith='kabat').order_by('file_id'),
                    label='Входной KABAT файл', empty_label=None, required=True)
    model_name      = forms.CharField(widget=forms.TextInput(attrs={'class': 'form-control', 'type': 'text',
                    'data-validation': 'length', 'data-validation-length': 'min7',
                    'data-validation-error-msg': 'Введите, пожалуйста, имя модели. Оно должно оканчиваться на .model'}),
                    label='Имя файла модели', required=True, initial='model.model')
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

    def clean(self):
        try:
            if not ('fasta' in self.cleaned_data and self.cleaned_data['fasta']):
                raise forms.ValidationError('fasta parameters missing')
            if not ('kabat' in self.cleaned_data and self.cleaned_data['kabat']):
                raise forms.ValidationError('kabat parameters missing')
            if not ('model_name' in self.cleaned_data and self.cleaned_data['model_name']):
                raise forms.ValidationError('model_name parameters missing')
            if not ('ml_window_size' in self.cleaned_data):
                raise forms.ValidationError('ml_window_size parameters missing')
            if not ('group' in self.cleaned_data):
                raise forms.ValidationError('group parameters missing')

        except forms.ValidationError as e:
            log.debug('Error: %s' % e.messages)
            raise
        return self.cleaned_data