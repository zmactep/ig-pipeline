from django.conf.urls import patterns, url

from igtasks import views

urlpatterns = patterns('',
    url(r'^$', views.create, name='create_view'),
    url(r'^result/$', views.send_request, name="send_request")
)