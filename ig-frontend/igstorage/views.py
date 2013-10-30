from django.http import HttpResponse
from django.views import generic
from django.shortcuts import render

import logging
from igstorage.models import StorageItem, StorageItemForm

log = logging.getLogger('all')


class StorageItemView(generic.ListView):
    template_name = 'igstorage/show_storage_items.html'
    paginate_by = 25
    context_object_name = 'items'

    def get_queryset(self):
        return StorageItem.objects.all()


def modify(request):
    if request.method == 'POST':  # If the form has been submitted...
        form = StorageItemForm(request.POST)
        if form.is_valid():
            data = form.cleaned_data
            if 'btn_create' in request.POST:
                try:
                    StorageItem.objects.get(file_id=data['file_id'])
                    return render(request, 'igstorage/modify_storage_items.html', dictionary={'form': form, 'result': 'Already exists'})
                except StorageItem.DoesNotExist:
                    item = StorageItem(file_id=data['file_id'], comment=data['comment'], path=data['path'])
                    item.save()

            elif 'btn_update' in request.POST:
                try:
                    item = StorageItem.objects.get(file_id=data['file_id'])
                    item.comment = data['comment']
                    item.path = data['path']
                    item.save()
                except StorageItem.DoesNotExist:
                    return render(request, 'igstorage/modify_storage_items.html', dictionary={'form': form, 'result': 'Not found'})

            elif 'btn_delete' in request.POST:
                item = StorageItem.objects.get(data['file_id'])
                item.delete()
            return render(request, 'igstorage/modify_storage_items.html', dictionary={'form': form, 'result': 'OK'})
    else:
        form = StorageItemForm()  # An unbound form

    return render(request, 'igstorage/modify_storage_items.html', dictionary={'form': form})