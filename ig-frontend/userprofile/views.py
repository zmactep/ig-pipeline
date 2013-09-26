# Create your views here.
from django.http import HttpResponse
from django.shortcuts import render


def profile(request, username):
    return render(request, "userprofile/profile.html", {"viewuser": username})


def userlist(request):
    return render(request, "userprofile/list.html")