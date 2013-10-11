from django.db import models
from django import forms
from django.forms.models import model_to_dict

import json
import logging
log = logging.getLogger('all')

root_path = '/Users/Kos/Dropbox/Biocad/ig-pipeline/tmp2/'

class TaskRequest(models.Model):
    FIND_PATTERNS = 1
    GENERATE_MODEL = 2
    MODEL_LIST = 3
    TASK_CHOICES = (
        (FIND_PATTERNS, 'find patterns'),
        (GENERATE_MODEL, 'generate model'),
        (MODEL_LIST, 'model list')
    )

    task = models.IntegerField(choices=TASK_CHOICES, default=MODEL_LIST, null=False)

    # files and paths
    input_file = models.FilePathField(null=True)
    model_path = models.FilePathField(null=True)
    out_dir = models.CharField(max_length=256, null=True)

    # algo params
    algo = models.CharField(max_length=256, null=True)
    algo_params = models.CharField(max_length=256, null=True)
    avg_window_size = models.IntegerField(null=True)
    ml_window_size = models.IntegerField(null=True)
    merge_threshold = models.IntegerField(null=True)

    # metadata
    group = models.CharField(max_length=256, null=True)
    comment = models.CharField(max_length=256, null=True)

    def __str__(self):
        return json.dumps(model_to_dict(self, fields=[], exclude=[]))


class TaskRequestForm(forms.Form):
    server          = forms.CharField(widget=forms.TextInput(attrs={'class': 'form-control', 'type': 'text', 'id': 'server_id'}), initial='127.0.0.1')
    port            = forms.IntegerField(widget=forms.TextInput(attrs={'class': 'form-control', 'type': 'text', 'id': 'port_id'}), initial=8080)
    task            = forms.ChoiceField(widget=forms.Select(attrs={'class': 'form-control', 'type': 'text', 'id': 'task_id'}),
                             choices=([('1', 'find patterns'), ('2', 'generate model'), ('3', 'model list'), ]),
                             initial='3',
                             required=True)
    input_file      = forms.FilePathField(widget=forms.TextInput(attrs={'class': 'form-control', 'type': 'text', 'id': 'input_file_id', 'style': 'display: none;'}), path=root_path)
    model_path      = forms.FilePathField(widget=forms.TextInput(attrs={'class': 'form-control', 'type': 'text', 'id': 'model_path_id', 'style': 'display: none;'}), path=root_path)
    out_dir         = forms.CharField(widget=forms.TextInput(attrs={'class': 'form-control', 'type': 'text', 'id': 'out_dir_id', 'style': 'display: none;'}), initial='output')
    group           = forms.CharField(widget=forms.TextInput(attrs={'class': 'form-control', 'type': 'text', 'id': 'group_id'}), initial='currently ignored')
    comment         = forms.CharField(widget=forms.TextInput(attrs={'class': 'form-control', 'type': 'text', 'id': 'comment_id', 'style': 'display: none;'}))
    algo            = forms.CharField(widget=forms.TextInput(attrs={'class': 'form-control', 'type': 'text', 'id': 'algo_id', 'style': 'display: none;'}), initial='random forest')
    algo_params     = forms.CharField(widget=forms.TextInput(attrs={'class': 'form-control', 'type': 'text', 'id': 'algo_params_id', 'style': 'display: none;'}), initial='-l 10 -S 0')
    avg_window_size = forms.IntegerField(widget=forms.TextInput(attrs={'class': 'form-control', 'type': 'text', 'id': 'avg_window_size_id', 'style': 'display: none;'}), initial=1)
    ml_window_size  = forms.IntegerField(widget=forms.TextInput(attrs={'class': 'form-control', 'type': 'text', 'id': 'ml_window_size_id', 'style': 'display: none;'}), initial=13)
    merge_threshold = forms.IntegerField(widget=forms.TextInput(attrs={'class': 'form-control', 'type': 'text', 'id': 'merge_threshold_id', 'style': 'display: none;'}), initial=7)
