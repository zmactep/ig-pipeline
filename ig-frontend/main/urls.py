__author__ = 'mactep'

from django.conf.urls import patterns, url
from main import views

urlpatterns = patterns('',
    url(r'^about/$', views.about, name="about"),
    url(r'^contact/$', views.contact, name="contact"),
    )