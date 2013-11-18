from django.http import HttpResponse
from igtools.models.predict import PredictForm, Predict
from igtools.models.train import TrainForm, Train
from igtools.models.cut_region import CutRegionForm, CutRegion
from igtools.models.make_report import ReportForm, Report
from igtools.models.simple_cluster import SimpleClusterForm, SimpleCluster
from igtools.models.manifest import Manifest
from django.shortcuts import render
from django.views import generic
from django.views.decorators.csrf import csrf_exempt
from django.core.paginator import Paginator, EmptyPage, PageNotAnInteger
from django.shortcuts import render_to_response
from itertools import chain

import urllib.parse
import urllib.request
import linecache
import logging
import json
import sys
import os

log = logging.getLogger('all')
backend_uri = 'http://127.0.0.1:8080'


class TaskRequestView(generic.ListView):
    template_name = 'igtools/results.html'
    paginate_by = 25
    context_object_name = 'tasks'

    def get_queryset(self):
        # TODO unhardcode me
        train_list = Train.objects.using('ig').all()
        predict_list = Predict.objects.using('ig').all()
        cluster_list = SimpleCluster.objects.using('ig').all()
        cut_region_list = CutRegion.objects.using('ig').all()
        report_list = Report.objects.using('ig').all()
        return sorted(chain(train_list, predict_list, cluster_list, cut_region_list, report_list), key=lambda record: record.backend_id)


# collect all possible output files from previous forms (by manifest data)
def add_previous_steps_files(forms_in_pipe):
    pipelined_files_map = {'fasta': [], 'kabat': [], 'model': [], 'dir': []}
    for counter, form in zip(range(10), forms_in_pipe):
        tool_name = form.data[str(counter) + '-name']
        tool = Manifest.objects.using('ig').get(tool_name=tool_name)
        for file_info in json.loads(tool.manifest)['output']:
            if file_info['pipelined'] and file_info['type'] in pipelined_files_map:
                pipelined_files_map[file_info['type']] += [('{' + str(counter) + '}/' + file_info['name'], str(counter))]
    return pipelined_files_map


def create(request):
    forms_in_pipe = []
    response = 'Nothing happened yet'
    if request.method == 'POST':  # If the form has been submitted...
        post = request.POST
        #restore all forms from pipeline
        for i in range(10):
            form_name_tag = str(i) + '-name'
            if form_name_tag in post:
                last_form = eval(post[form_name_tag] + 'Form')(post, request.FILES, prefix=str(i))
                last_form.set_additional_files(add_previous_steps_files(forms_in_pipe))
                forms_in_pipe += [last_form]

        if 'btn_delete' in request.POST:
            if len(forms_in_pipe) > 0:
                forms_in_pipe.pop(len(forms_in_pipe) - 1)

        if 'btn_add' in request.POST:
            if len(forms_in_pipe) < 10:
                # manually add files to input lists:
                new_form = eval(post['tools_select'] + 'Form')(prefix=str(len(forms_in_pipe)))
                new_form.set_additional_files(add_previous_steps_files(forms_in_pipe))
                forms_in_pipe += [new_form]

        if 'btn_start' in request.POST:
            list_of_requests = []
            for form, index in zip(forms_in_pipe, range(len(forms_in_pipe))):
                if form.is_valid():
                    params_map = form.cleaned_data
                    req = eval(params_map['name'])()
                    req.read_params(form, index)
                    list_of_requests += [req.get_backend_request()]

            result = ask_server(json.dumps({"commands": list_of_requests}))
            try:
                req.backend_id = json.loads(result)['id']
                req.save(using='ig')
                response = result
            except Exception as e:
                response = 'Bad response from server: %s. Probably, backend is shut down' % e

    tools = Manifest.objects.using('ig').all()
    return render(request, 'igtools/send_request.html', dictionary={'forms': forms_in_pipe, 'tools': tools, 'response': response})


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


lines_per_page = 24


#TODO limit the ability to view any file on filesystem
def listing(request):
    data = []
    params = request.GET
    filename = params.get('file')

    if os.path.exists(filename):
        with open(filename, 'r') as data_source:
            data = [data_source.readline() for i in range(lines_per_page + 1)] # Extra line to create .next in Paginator

    paginator = Paginator(data, lines_per_page)

    return render_to_response('igtools/file_viewer.html', {"data": paginator.page(1), 'file': filename})


def next_listing(request):
    params = request.GET
    page = params.get('page')
    filename = params.get('file')
    data = ['' for i in range((int(page) - 1) * lines_per_page)]

    if os.path.exists(filename):
        data += [linecache.getline(filename, (int(page) - 1) * lines_per_page + i + 1) for i in range(lines_per_page + 1)]

    paginator = Paginator(data, lines_per_page)

    try:
        p = paginator.page(page)
    except PageNotAnInteger:
        p = paginator.page(1)
    except EmptyPage:
        p = paginator.page(paginator.num_pages)

    return render_to_response('igtools/file_viewer_next_page.html', {"data": p, 'file': filename})