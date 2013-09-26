__author__ = 'mactep'

from django.conf.urls import patterns, url

from main import views

urlpatterns = patterns('',
    url(r'^$', views.index, name="index"),
    url(r'^result/$', views.send_request, name="send_request"),
    url(r'^about/$', views.about, name="about"),
    url(r'^contact/$', views.contact, name="contact"),
    )