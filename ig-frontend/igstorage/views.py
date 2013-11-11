from django.shortcuts import render

import logging
from igstorage.models import StorageItem, StorageItemForm

log = logging.getLogger('all')


def view(request):
    if request.method == 'POST':  # If the form has been submitted...
        if 'btn_delete' in request.POST:
            file_id = request.POST.get('file_id', 'NA')
            item = StorageItem.objects.using('ig').get(file_id=file_id)
            item.delete(using='ig')
            form = StorageItemForm()
        elif 'btn_update' in request.POST:
            file_id = request.POST.get('file_id', 'NA')
            new_file_id = request.POST.get('modified_id', 'NA')
            comment = request.POST.get('comment', 'NA')
            group = request.POST.get('group', 'NA')
            run = request.POST.get('run', 'NA')

            item = StorageItem.objects.using('ig').get(file_id=file_id, group=group, run=run)
            path = item.path
            group = item.group
            run = item.run
            item.delete(using='ig')

            new_item = StorageItem(file_id=new_file_id, comment=comment, path=path, run=run, group=group)
            new_item.save(using='ig')

            form = StorageItemForm()
        else:
            form = StorageItemForm(request.POST)

            if form.is_valid():
                data = form.cleaned_data
                if 'btn_create' in request.POST:
                    try:
                        StorageItem.objects.using('ig').get(file_id=data['file_id'])
                        return render(request, 'igstorage/show_storage_items.html', dictionary={'form': form, 'status': 'Already exists'})
                    except StorageItem.DoesNotExist:
                        item = StorageItem(file_id=data['file_id'], comment=data['comment'], path=data['path'], group=data['group'], run=data['run'])
                        item.save(using='ig')

                    items = StorageItem.objects.using('ig').all()
                    return render(request, 'igstorage/show_storage_items.html', dictionary={'form': form, 'status': 'OK', 'items': items})
    else:
        form = StorageItemForm()  # An unbound form

    items = StorageItem.objects.using('ig').all()
    return render(request, 'igstorage/show_storage_items.html', dictionary={'form': form, 'items': items})
