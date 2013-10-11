from django.db import models
from django.core import serializers
from django.forms import ModelForm

import json

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
        data = self.get_queryset()[:1]
        return json.dumps(serializers.serialize('json', data))

class TaskRequestForm(ModelForm):
    class Meta:
        model = TaskRequest
        fields = ['task']