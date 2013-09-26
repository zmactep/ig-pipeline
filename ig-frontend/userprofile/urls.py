__author__ = 'mactep'
from django.conf.urls import patterns, url

from userprofile import views

urlpatterns = patterns('',
    url(r'^list/$', views.userlist, name="userlist"),
    url(r'^(?P<username>\w+)$', views.profile, name="profile"),
    )