from django.db import models
from django import forms
from django.forms.models import model_to_dict
from igstorage.models import StorageItem
import json
import os


import logging
from igtools.fields import CachedModelChoiceField

log = logging.getLogger('all')


class SimpleCluster(models.Model):
    backend_id = models.IntegerField(null=True)
    name = models.CharField(max_length=256, null=True, blank=True)

    # files and paths
    src = models.CharField(max_length=256, null=True, blank=True)

    # params
    sim = models.BooleanField()
    shortest_cons = models.BooleanField()
    skip_first = models.IntegerField(null=True, blank=True)
    min_len = models.IntegerField(null=True, blank=True)
    use_prct = models.IntegerField(null=True, blank=True)

    # metadata
    group = models.CharField(max_length=256)
    comment = models.CharField(max_length=256)

    def __str__(self):
        return json.dumps(model_to_dict(self, fields=[], exclude=[]))

    def read_params(self, form, index):
        params_map = form.data
        if str(index) + '-src' in params_map:
            self.src = params_map[str(index) + '-src']

        if str(index) + '-sim' in params_map:
            self.sim = params_map[str(index) + '-sim']

        if str(index) + '-skip_first' in params_map and params_map[str(index) + '-skip_first']:
            self.skip_first = params_map[str(index) + '-skip_first']

        if str(index) + '-min_len' in params_map:
            self.min_len = params_map[str(index) + '-min_len']

        if str(index) + '-use_prct' in params_map:
            self.use_prct = params_map[str(index) + '-use_prct']

        if str(index) + '-shortest_cons' in params_map:
            self.shortest_cons = params_map[str(index) + '-shortest_cons']

        if str(index) + '-group' in params_map:
            self.group = params_map[str(index) + '-group']

        if str(index) + '-comment' in params_map:
            self.comment = params_map[str(index) + '-comment']

    def get_backend_request(self):
        params = [{"name": "src", "value": self.src}]
        if self.sim: # checkbox is converted to "name" without params if true. Backend will convert it to "--name  "
            params += [{"name": "sim", "value": ""}]

        if self.shortest_cons:
            params += [{"name": "shortest-cons", "value": ""}]

        if self.skip_first:
            params += [{"name": "skip-first", "value": str(self.skip_first)}]

        if self.min_len:
            params += [{"name": "min-len", "value": str(self.min_len)}]

        if self.use_prct:
            params += [{"name": "use-prct", "value": str(self.use_prct)}]

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
    src             = CachedModelChoiceField(widget=forms.Select(attrs={'class': 'form-control', 'type': 'text'}),
                    label='Входной FASTA файл', empty_label=None, required=True)
    sim             = forms.BooleanField(widget=forms.CheckboxInput(attrs={'class': 'form-control', 'type': 'checkbox'}), label='Поиск дубликатов вместо кластеризации (sim)', required=False)
    shortest_cons   = forms.BooleanField(widget=forms.CheckboxInput(attrs={'class': 'form-control', 'type': 'checkbox'}), label='Shortest Cons', required=False) # нужно если sim = true
    skip_first      = forms.IntegerField(widget=forms.TextInput(attrs={'class': 'form-control', 'type': 'text'}), label='Skip First', required=False) # нужно если sim = true
    min_len         = forms.IntegerField(widget=forms.TextInput(attrs={'class': 'form-control', 'type': 'text',
                    'data-validation': 'number', 'data-validation-allowing': 'range[0;500]',
                    'data-validation-error-msg': 'Введите, пожалуйста, размер, минимальную длину в диапазоне 0-500'}),
                    initial=100, label='Minimal Length', required=True)
    use_prct        = forms.IntegerField(widget=forms.TextInput(attrs={'class': 'form-control', 'type': 'text',
                    'data-validation': 'number', 'data-validation-allowing': 'range[0;100]',
                    'data-validation-error-msg': 'Введите, пожалуйста, число в диапазоне 0-100'}), label='Use %', initial=100, required=False)  # нужно если sim = true

    group           = forms.CharField(widget=forms.TextInput(attrs={'class': 'form-control', 'type': 'text',
                    'data-validation': 'length', 'data-validation-length': 'min5',
                    'data-validation-error-msg': 'Введите, пожалуйста, название группы не короче 5 символов'}),
                    initial='testrun', label='Группа', required=True)
    comment         = forms.CharField(widget=forms.TextInput(attrs={'class': 'form-control', 'type': 'text',
                    'data-validation': 'length', 'data-validation-length': 'min5',
                    'data-validation-error-msg': 'Введите, пожалуйста, комментарий не короче 5 символов'}),
                    initial='comment', label='Комментарий', required=True)

    def set_additional_files(self, pipelined_files_map):
        fasta = pipelined_files_map['fasta']
        new_fasta = {file: str(os.path.basename(file) + ' - pipeline from stage ' + stage + ' - 0') for file, stage in fasta}
        fasta_in_db = {item.path: str(item.file_id + ' - ' + item.group + ' - ' + item.run) for item in StorageItem.objects.using('ig').filter(path__endswith='fasta').order_by('file_id')}
        fasta_in_db.update(new_fasta)
        self.fields['src'].objects = fasta_in_db


    def clean(self):
        try:
            if not ('src' in self.cleaned_data and self.cleaned_data['src']):
                raise forms.ValidationError('src parameters missing')
            if not ('group' in self.cleaned_data):
                raise forms.ValidationError('group parameters missing')
            if not ('comment' in self.cleaned_data):
                raise forms.ValidationError('comment parameters missing')
            if not ('min_len' in self.cleaned_data):
                raise forms.ValidationError('min_len parameters missing')
        except forms.ValidationError as e:
            log.debug('Error: %s' % e.messages)
            raise
        return self.cleaned_data