# Create your views here.
from django.http import HttpRequest, HttpResponse
from django.shortcuts import render

import socket


def index(request):
    return render(request, "main/sendrequest.html")


def send_request(request):
    server = request.POST['server']
    port = int(request.POST['port'])
    rtext = request.POST['rtext']

    s = socket.socket()
    s.connect((server, port))

    s.send(rtext.encode())
    result = s.recv(len(rtext))

    s.close()

    return render(request, "main/showresponse.html", {"result": result})