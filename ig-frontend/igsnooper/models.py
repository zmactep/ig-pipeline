from django.db import models
from django import forms
from django.forms.models import model_to_dict

import json
import logging
log = logging.getLogger('all')


class TaskRequest(models.Model):
    FIND_PATTERNS = 1
    GENERATE_MODEL = 2
    MODEL_LIST = 3
    TASK_CHOICES = (
        (FIND_PATTERNS, 'find patterns'),
        (GENERATE_MODEL, 'generate model'),
        (MODEL_LIST, 'model list')
    )

    backend_id = models.IntegerField(null=True)
    task = models.IntegerField(choices=TASK_CHOICES, default=MODEL_LIST, null=False)

    # files and paths
    input_file_fasta = models.CharField(max_length=256, null=True, blank=True)
    input_file_kabat = models.CharField(max_length=256, null=True, blank=True)
    model_name = models.CharField(max_length=256, null=True, blank=True)
    model_path = models.CharField(max_length=256, null=True, blank=True)
    out_dir = models.CharField(max_length=256, null=True, blank=True)

    # algo params
    algo = models.CharField(max_length=256, null=True, blank=True)
    algo_params = models.CharField(max_length=256, null=True, blank=True)
    avg_window_size = models.IntegerField(null=True, blank=True)
    ml_window_size = models.IntegerField(null=True, blank=True)
    merge_threshold = models.IntegerField(null=True, blank=True)

    # metadata
    group = models.CharField(max_length=256)
    comment = models.CharField(max_length=256)

    def __str__(self):
        return json.dumps(model_to_dict(self, fields=[], exclude=[]))

    def get_backend_request(self):
        request = None
        if int(self.FIND_PATTERNS) == int(self.task):
            request = {"task": "1",
                       "input": {
                           "files": [self.input_file_fasta, self.input_file_kabat],
                           "params": {
                               "mlWindowsize": str(self.ml_window_size),
                               "avgWidowsize": str(self.avg_window_size),
                               "mergeThreshold": str(self.merge_threshold),
                               "modelPath": self.model_path
                           }
                       },
                       "output": {
                           "outdir": self.out_dir
                       }
                    }

        if int(self.GENERATE_MODEL) == int(self.task):
            request = {"task": "2",
                       "input": {
                           "files": [self.input_file_fasta, self.input_file_kabat],
                           "params": {
                               "mlWindowsize": str(self.ml_window_size),
                               "algo": self.algo,
                               "modelName": self.model_name,
                               "algoParams": self.algo_params
                           },
                           "comment": self.comment,
                           "group": self.group
                       },
                       "output": {
                           "outdir": self.out_dir
                       }
                    }

        if int(self.MODEL_LIST) == int(self.task):
            request = {"task": "3", "input": {"group": self.group}}

        return json.dumps(request)


class TaskRequestForm(forms.Form):
    task            = forms.ChoiceField(widget=forms.Select(attrs={'class': 'form-control', 'type': 'text', 'id': 'task_id'}), label='Choose action',
                             choices=([('1', 'find patterns'), ('2', 'generate model'), ('3', 'model list'), ]),
                             initial='3',
                             required=True)
    input_file_fasta= forms.CharField(widget=forms.TextInput(attrs={'class': 'form-control', 'type': 'text', 'id': 'input_file_fasta_id', 'style': 'display: none;'}), label='Input FASTA file', required=False, initial='/Users/Kos/Dropbox/Biocad/ig-pipeline/data/train/VDJH_train.fasta')
    input_file_kabat= forms.CharField(widget=forms.TextInput(attrs={'class': 'form-control', 'type': 'text', 'id': 'input_file_kabat_id', 'style': 'display: none;'}), label='Input KABAT file', required=False, initial='/Users/Kos/Dropbox/Biocad/ig-pipeline/data/train/VDJH_train.kabat')
    model_name      = forms.CharField(widget=forms.TextInput(attrs={'class': 'form-control', 'type': 'text', 'id': 'model_name_id', 'style': 'display: none;'}), label='Model file name', required=False, initial='model.model')
    model_path      = forms.CharField(widget=forms.TextInput(attrs={'class': 'form-control', 'type': 'text', 'id': 'model_path_id', 'style': 'display: none;'}), label='Model file path', required=False, initial='task1/model.model')
    out_dir         = forms.CharField(widget=forms.TextInput(attrs={'class': 'form-control', 'type': 'text', 'id': 'out_dir_id', 'style': 'display: none;'}), label='Output dir', required=False, initial='task1/')
    group           = forms.CharField(widget=forms.TextInput(attrs={'class': 'form-control', 'type': 'text', 'id': 'group_id'}), initial='currently ignored', label='Group')
    comment         = forms.CharField(widget=forms.TextInput(attrs={'class': 'form-control', 'type': 'text', 'id': 'comment_id', 'style': 'display: none;'}), initial='Some useful comment', label='Your comment')
    algo            = forms.CharField(widget=forms.TextInput(attrs={'class': 'form-control', 'type': 'text', 'id': 'algo_id', 'style': 'display: none;'}), initial='random forest', label='Machine learning algorithm')
    algo_params     = forms.CharField(widget=forms.TextInput(attrs={'class': 'form-control', 'type': 'text', 'id': 'algo_params_id', 'style': 'display: none;'}), initial='-l 10 -S 0', label='Parameters for algorithm')
    avg_window_size = forms.IntegerField(widget=forms.TextInput(attrs={'class': 'form-control', 'type': 'text', 'id': 'avg_window_size_id', 'style': 'display: none;'}), initial=1, label='Averaging window size')
    ml_window_size  = forms.IntegerField(widget=forms.TextInput(attrs={'class': 'form-control', 'type': 'text', 'id': 'ml_window_size_id', 'style': 'display: none;'}), initial=13, label='Machine learning window size')
    merge_threshold = forms.IntegerField(widget=forms.TextInput(attrs={'class': 'form-control', 'type': 'text', 'id': 'merge_threshold_id', 'style': 'display: none;'}), initial=7, label='Merge threshold')

    def clean(self):
        try:
            if int(TaskRequest.FIND_PATTERNS) == int(self.cleaned_data['task']):
                if not ('input_file_fasta' in self.cleaned_data and self.cleaned_data['input_file_fasta']):
                    raise forms.ValidationError('input_file_fasta parameters missing')
                if not ('ml_window_size' in self.cleaned_data):
                    raise forms.ValidationError('ml_window_size parameters missing')
                if not ('avg_window_size' in self.cleaned_data):
                    raise forms.ValidationError('avg_window_size parameters missing')
                if not ('merge_threshold' in self.cleaned_data):
                    raise forms.ValidationError('merge_threshold parameters missing')
                if not ('model_path' in self.cleaned_data and self.cleaned_data['model_path']):
                    raise forms.ValidationError('model_path parameters missing')
                if not ('out_dir' in self.cleaned_data and self.cleaned_data['out_dir']):
                    raise forms.ValidationError('out_dir parameters missing')

            if int(TaskRequest.GENERATE_MODEL) == int(self.cleaned_data['task']):
                if not ('input_file_fasta' in self.cleaned_data and self.cleaned_data['input_file_fasta']):
                    raise forms.ValidationError('input_file_fasta parameters missing')
                if not ('input_file_kabat' in self.cleaned_data and self.cleaned_data['input_file_kabat']):
                    raise forms.ValidationError('input_file_kabat parameters missing')
                if not ('ml_window_size' in self.cleaned_data):
                    raise forms.ValidationError('ml_window_size parameters missing')
                if not ('algo' in self.cleaned_data and self.cleaned_data['algo']):
                    raise forms.ValidationError('algo parameters missing')
                if not ('algo_params' in self.cleaned_data and self.cleaned_data['algo_params']):
                    raise forms.ValidationError('algo_params parameters missing')
                if not ('out_dir' in self.cleaned_data and self.cleaned_data['out_dir']):
                    raise forms.ValidationError('out_dir parameters missing')
                if not ('model_name' in self.cleaned_data and self.cleaned_data['model_name']):
                    raise forms.ValidationError('model_name parameters missing')

        except forms.ValidationError as e:
            log.debug('Error: %s' % e.messages)
            raise
        return self.cleaned_data
