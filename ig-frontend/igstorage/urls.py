from django.conf.urls import patterns, url

from igstorage import views

urlpatterns = patterns('',
    url(r'^storage/$', views.view, name="show_storage_items"),
    url(r'^storage/listdir/$', views.dir_list, name="dir_list"),
)