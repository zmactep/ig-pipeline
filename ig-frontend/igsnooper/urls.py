from django.conf.urls import patterns, url

from igsnooper import views

urlpatterns = patterns('',
    url(r'^$', views.create, name='create_view'),
    url(r'^result/$', views.ask_backend, name="ask_backend"),
    url(r'^tasks/$', views.TaskRequestView.as_view(), name="task_request")
)