from django.conf.urls import patterns, url

from igstorage import views

urlpatterns = patterns('',
    url(r'^storage/$', views.StorageItemView.as_view(), name="show_storage_items"),
    url(r'^storage/modify/$', views.modify, name='modify_view')
)