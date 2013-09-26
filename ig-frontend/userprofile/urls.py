__author__ = 'mactep'
from django.conf.urls import patterns, url

from userprofile import views

urlpatterns = patterns('',
    url(r'^(?P<username>\w+)$', views.profile, name="profile"),
    )