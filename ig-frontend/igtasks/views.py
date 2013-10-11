from django.http import HttpResponse
from igtasks.models import TaskRequestForm
from igtasks.models import TaskRequest
from django.shortcuts import render

import urllib.parse
import urllib.request
import logging
import json
import sys

log = logging.getLogger('all')


def create(request):
    if request.method == 'POST': # If the form has been submitted...
        form = TaskRequestForm(request.POST)
        if form.is_valid():
            data = form.cleaned_data
            task_request = TaskRequest(task=data['task'])
            task_request.save()
            log.debug("Got request: %s" % task_request)
            return HttpResponse(json.dumps(ask_server(data['server'], data['port'], task_request.__str__())),
                                content_type="application/json")

    else:
        form = TaskRequestForm() # An unbound form

    return render(request, 'send_request.html', dictionary={'form': form})


def send_request(request):
    server = request.POST.get('server')
    port = int(request.POST.get('port', 0))
    rtext = request.POST.get('rtext')

    log.debug("Got request: server = %s; port = %d; rtext = %s" % (server, port, rtext))
    return HttpResponse(json.dumps(ask_server(server, port, rtext)), content_type="application/json")


def ask_server(server, port, query):
    try:
        uri = 'http://' + server + ':' + str(port)
        log.debug("Request to ig-backend @ %s: %s" % (uri, query))
        req = urllib.request.Request(uri, query.encode('utf-8'), method='POST')

        response = urllib.request.urlopen(req)
        result = response.read().decode()
        log.debug("Response from ig-backend: %s" % result)
        return result
    except:
        type, value, tb = sys.exc_info()
        log.debug('%s', value)
        return {'status': 'fail', 'text': "Server error"}