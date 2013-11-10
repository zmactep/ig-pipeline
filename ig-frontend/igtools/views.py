from django.http import HttpResponse
from igtools.models.predict import PredictForm, Predict
from igtools.models.train import TrainForm, Train
from igtools.models.simplecluster import SimpleClusterForm, SimpleCluster
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
        return sorted(chain(train_list, predict_list, cluster_list), key=lambda record: record.backend_id)


def create(request):
    if request.method == 'POST':  # If the form has been submitted...
        post = request.POST
        if 'tools_select' in post:
            form = eval(post['tools_select'])
            return render(request, 'igtools/send_request.html', dictionary={'form': form})

        else:
            if post['name']:
                name = post['name']
                form = eval(name + 'Form')(post, request.FILES)
                if form.is_valid():
                    params_map = form.cleaned_data
                    request = eval(name)()
                    request.read_params(params_map)

            response = ask_server(request.get_backend_request())
            request.backend_id = json.loads(response)['id']
            request.save(using='ig')
            return HttpResponse(response, content_type="application/json")
    else:
        form = TrainForm()  # An unbound form

    return render(request, 'igtools/send_request.html', dictionary={'form': form})


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