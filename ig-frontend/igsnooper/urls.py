from django.conf.urls import patterns, url

from igsnooper import views

urlpatterns = patterns('',
    url(r'^$', views.create, name='create_view'),
    url(r'^result/$', views.ask_backend, name="ask_backend"),
    url(r'^tasks/$', views.TaskRequestView.as_view(), name="task_request"),
    url(r'^listing/$', views.listing, name="view_results"),
    url(r'^listing/next/$', views.next_listing, name="view_results_next_page")
)