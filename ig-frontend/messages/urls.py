__author__ = 'mactep'
from django.conf.urls import patterns, url

from messages import views

urlpatterns = patterns('',
    url(r'^list/$', views.msglist, name="list"),
    url(r'^(?P<username>\w+)$', views.dialog, name="dialog"),
    url(r'^markread/(?P<dialog_id>\d+)$', views.mark_read, name="mark_read"),
    url(r'^send_message/$', views.send_message, name="send_message"),
    )