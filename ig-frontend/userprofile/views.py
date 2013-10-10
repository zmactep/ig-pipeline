# Create your views here.
from django.shortcuts import render
from django.contrib.auth import get_user_model


def profile(request, username):
    if not get_user_model().objects.get(username=username):
        return render(request, "main/no_user.html", {"username": username})
    return render(request, "userprofile/profile.html", {"viewuser": username})


def userlist(request):
    ulist = get_user_model().objects.all()
    return render(request, "userprofile/list.html", {'userlist': ulist})