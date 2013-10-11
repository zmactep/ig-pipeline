from django.views.generic.edit import CreateView
from django.http import HttpRequest, HttpResponse
from igtasks.models import TaskRequest
from igtasks.models import TaskRequestForm

import urllib.parse
import urllib.request
import logging
import json
import sys

log = logging.getLogger('all')

class TaskCreateView(CreateView):
    model = TaskRequest
    form_class = TaskRequestForm
    template_name = "send_request.html"

def send_request(request):
    server = request.POST.get('server')
    port = int(request.POST.get('port', 0))
    rtext = request.POST.get('rtext')

    log.debug("Got request: server = %s; port = %d; rtext = %s" % (server, port, rtext))

    if not server or not port or not rtext:
        return HttpResponse(json.dumps({'status': 'fail', 'text': "wrong request"}), content_type="application/json")

    try:
        uri = 'http://' + server + ':' + str(port)
        log.debug("Request to ig-backend @ %s: %s" % (uri, rtext))
        req = urllib.request.Request(uri, rtext.encode('utf-8'), method='POST')

        response = urllib.request.urlopen(req)
        result = response.read().decode()
        log.debug("Response from ig-backend: %s" % result)

        return HttpResponse(json.dumps(result), content_type="application/json")
    except:
        log.debug(sys.exec_info()[0])
        return HttpResponse(json.dumps({'status': 'fail', 'text': "Server error"}), content_type="application/json")
