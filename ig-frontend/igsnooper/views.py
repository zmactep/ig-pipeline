from django.http import HttpResponse
from igsnooper.models import TaskRequestForm
from igsnooper.models import TaskRequest
from django.shortcuts import render
from django.views import generic
from django.views.decorators.csrf import csrf_exempt

import urllib.parse
import urllib.request
import logging
import json
import sys

log = logging.getLogger('all')
backend_uri = 'http://127.0.0.1:8080'


class TaskRequestView(generic.ListView):
    template_name = 'igsnooper/task_request.html'
    paginate_by = 25
    context_object_name = 'tasks'

    def get_queryset(self):
        return TaskRequest.objects.all()


def create(request):
    if request.method == 'POST':  # If the form has been submitted...
        form = TaskRequestForm(request.POST, request.FILES)
        if form.is_valid():
            data = form.cleaned_data
            task_request = None
            if int(data['task']) == int(TaskRequest.FIND_PATTERNS):
                task_request = TaskRequest(task=data['task'], input_file_fasta=data['input_file_fasta'],
                                           input_file_kabat=data['input_file_kabat'],
                                           model_path=data['model_path'], ml_window_size=data['ml_window_size'],
                                           avg_window_size=data['avg_window_size'], out_dir=data['out_dir'],
                                           merge_threshold=data['merge_threshold'])

            if int(data['task']) == int(TaskRequest.GENERATE_MODEL):
                task_request = TaskRequest(task=data['task'], input_file_fasta=data['input_file_fasta'],
                                           input_file_kabat=data['input_file_kabat'], algo=data['algo'],
                                           algo_params=data['algo_params'], out_dir=data['out_dir'],
                                           ml_window_size=data['ml_window_size'])
            if int(data['task']) == int(TaskRequest.MODEL_LIST):
                task_request = TaskRequest(task=data['task'], group=data['group'])

            response = ask_server(task_request.get_backend_request())
            task_request.backend_id = json.loads(response)['id']
            task_request.save()
            return HttpResponse(response, content_type="application/json")
    else:
        form = TaskRequestForm()  # An unbound form

    return render(request, 'igsnooper/send_request.html', dictionary={'form': form})


@csrf_exempt
def ask_backend(request):
    back_id = request.POST.get('id')
    return HttpResponse(json.dumps(ask_server(json.dumps({'result_for': back_id}))), content_type="application/json")


def ask_server(query):
    try:
        log.debug("Request to ig-backend: %s" % query)
        req = urllib.request.Request(backend_uri, query.encode('utf-8'), method='POST')

        response = urllib.request.urlopen(req)
        result = response.read().decode()
        log.debug("Response from ig-backend: %s" % result)
        return result
    except:
        _, value, _ = sys.exc_info()
        log.debug('%s', value)
        return {'status': 'fail', 'text': "Server error"}