__author__ = 'mactep'
from django.contrib import admin
from messages.models import Dialog, Message


class MessageInline(admin.StackedInline):
    model = Message


class DialogAdmin(admin.ModelAdmin):
    fields = ['first_user', 'second_user']
    inlines = [MessageInline]

admin.site.register(Dialog, DialogAdmin)