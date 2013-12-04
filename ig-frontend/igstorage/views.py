from django.http import HttpResponse
from django.shortcuts import render

import os
import logging
from django.views.decorators.csrf import csrf_exempt
from igcad import settings
from igstorage.models import StorageItem, StorageItemForm

log = logging.getLogger('all')


def view(request):
    if request.method == 'POST':  # If the form has been submitted...
        if 'btn_delete' in request.POST:
            file_id = request.POST.get('file_id', 'NA')
            group = request.POST.get('group', 'NA')
            run = request.POST.get('run', 'NA')
            item = StorageItem.objects.using('ig').get(file_id=file_id, group=group, run=run)
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
                    for file in scan_dir(data['path']):
                        try:
                            StorageItem.objects.using('ig').get(file_id=get_basename(file), group=data['group'], run=data['run'])
                        except StorageItem.DoesNotExist:
                            item = StorageItem(file_id=get_basename(file), comment=data['comment'], path=file, group=data['group'], run=data['run'])
                            item.save(using='ig')

                    items = StorageItem.objects.using('ig').all().order_by('group', 'run')
                    return render(request, 'igstorage/show_storage_items.html', dictionary={'form': form, 'status': 'OK', 'items': items, 'storage_root': settings.STORAGE_ROOT})
    else:
        form = StorageItemForm()  # An unbound form

    items = StorageItem.objects.using('ig').all().order_by('group', 'run')
    return render(request, 'igstorage/show_storage_items.html', dictionary={'form': form, 'items': items, 'storage_root': settings.STORAGE_ROOT})


@csrf_exempt
def dir_list(request):
    r = ['<ul class="jqueryFileTree" style="display: none;">']
    try:
        r = ['<ul class="jqueryFileTree" style="display: none;">']
        d = request.POST.get('dir')
        for f in os.listdir(d):
            ff=os.path.join(d, f)
            if os.path.isdir(ff):
                r.append('<li class="directory collapsed"><a href="#" rel="%s/">%s</a></li>' % (ff, f))
            else:
                e = os.path.splitext(f)[1][1:] # get .ext and remove dot
                r.append('<li class="file ext_%s"><a href="#" rel="%s">%s</a></li>' % (e, ff, f))
        r.append('</ul>')
    except Exception as e:
        r.append('Could not load directory: %s' % str(e))
    r.append('</ul>')
    return HttpResponse(''.join(r))


def get_basename(name):
    if name.endswith('/'):
        return os.path.basename(name[:-1])
    else:
        return os.path.basename(name)


# TODO rewrite for remote storage
def scan_dir(directory):
    result = [directory]
    if os.path.isdir(directory):
        for path in os.listdir(directory):
            result += [os.path.join(directory, path)]
    return result