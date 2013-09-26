# Create your views here.
from django.http import HttpRequest, HttpResponse
from django.shortcuts import render


import socket
import json


def index(request):
    return render(request, "main/sendrequest.html")


def send_request(request):
    print("OK!")
    server = request.POST.get('server')
    port = int(request.POST.get('port', 0))
    rtext = request.POST.get('rtext')

    if not server or not port or not rtext:
        return HttpResponse("", content_type="application/json")

    s = socket.socket()
    s.connect((server, port))

    s.send(rtext.encode())
    result = s.recv(len(rtext))

    s.close()

    return HttpResponse(json.dumps({'status': 'ok', 'text': result.decode()}), content_type="application/json")