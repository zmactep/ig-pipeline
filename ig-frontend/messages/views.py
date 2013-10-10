# Create your views here.
from django.http import HttpResponse
from django.shortcuts import render
from django.contrib.auth.decorators import login_required
from django.contrib.auth import get_user_model
from django.utils import timezone
import json

from messages.models import Dialog, Message


@login_required
def msglist(request):
    dialog_list = map(lambda x: (x, x.get_last_message(), x.get_collocutor(request.user.username)),
                      Dialog.get_dialogs_with(request.user.username))
    dialog_user_list = map(lambda x: x[2], dialog_list)
    return render(request, "messages/list.html", {'dialog_list': dialog_list, 'dialog_user_list': dialog_user_list})


@login_required
def dialog(request, username):
    if not get_user_model().objects.get(username=username):
        return render(request, "main/no_user.html", {"username": username})

    dialogs = Dialog.get_dialog_of(request.user.username, username)
    if not dialogs:
        messages = []
        dialog_id = 0
    else:
        messages = dialogs[0].get_messages()
        dialog_id = dialogs[0].id

    return render(request, "messages/dialog.html", {"username": username, "dialog_id": dialog_id, "messages": messages})


def mark_read(request, dialog_id):
    d = Dialog.objects.get(id=dialog_id)

    if not request.user.is_authenticated() or not dialog_id or \
            not d.in_dialog(request.user.username):
        return HttpResponse(json.dumps({"status": "failed"}), content_type="application/json")

    mid = []
    for msg in d.get_messages():
        if msg.sender.username != request.user.username or d.is_schizophrenia():
            msg.is_read = True
            msg.save()
            mid.append(msg.id)

    return HttpResponse(json.dumps({"status": "ok", "mid": mid}), content_type="application/json")


def send_message(request):
    did = int(request.POST.get('dialog_id'))
    username = request.POST.get('username')
    msg = request.POST.get('msg')

    if not request.user.is_authenticated() or not msg:
        return HttpResponse(json.dumps({"status": "failed"}), content_type="application/json")

    if not did:
        d = Dialog(first_user=request.user, second_user=get_user_model().objects.get(username=username))
        d.save()
    else:
        d = Dialog.objects.get(id=did)

    newmsg = Message(dialog=d, sender=request.user, is_read=False, send_date=timezone.now(), text=msg)
    newmsg.save()

    return HttpResponse(json.dumps({"status": "ok", "date": "Just now" }),
                        content_type="application/json")