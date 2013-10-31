from django.http import HttpResponse
from django.views import generic
from django.shortcuts import render

import logging
from igstorage.models import StorageItem, StorageItemForm

log = logging.getLogger('all')


def view(request):
    if request.method == 'POST':  # If the form has been submitted...
        if 'btn_delete' in request.POST:
            item_id = request.POST.get('id', 'NA')
            item = StorageItem.objects.using('ig').get(file_id=item_id)
            item.delete(using='ig')
            form = StorageItemForm()
        else:
            form = StorageItemForm(request.POST)

            if form.is_valid():
                data = form.cleaned_data
                if 'btn_create' in request.POST:
                    try:
                        StorageItem.objects.using('ig').get(file_id=data['file_id'])
                        return render(request, 'igstorage/show_storage_items.html', dictionary={'form': form, 'result': 'Already exists'})
                    except StorageItem.DoesNotExist:
                        item = StorageItem(file_id=data['file_id'], comment=data['comment'], path=data['path'])
                        item.save(using='ig')

                    items = StorageItem.objects.using('ig').all()
                    return render(request, 'igstorage/show_storage_items.html', dictionary={'form': form, 'result': 'OK', 'items': items})
    else:
        form = StorageItemForm()  # An unbound form

    items = StorageItem.objects.using('ig').all()
    return render(request, 'igstorage/show_storage_items.html', dictionary={'form': form, 'items': items})
