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

    task = models.IntegerField(choices=TASK_CHOICES, default=MODEL_LIST)

    def __str__(self):
        return json.dumps(model_to_dict(self, fields=[], exclude=[]))


class TaskRequestForm(forms.Form):
    server = forms.CharField(widget=forms.TextInput())
    port = forms.IntegerField()
    task = forms.ChoiceField(widget = forms.Select(),
                             choices = ([('1','find patterns'), ('2','generate model'),('3','model list'), ]),
                             initial ='3',
                             required = True)