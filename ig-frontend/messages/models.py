from django.db import models
from django.contrib.auth import get_user_model
# Create your models here.


class Dialog(models.Model):
    first_user = models.ForeignKey(get_user_model(), related_name="first_user")
    second_user = models.ForeignKey(get_user_model(), related_name="second_user")

    def __unicode__(self):
        return self.first_user.username + ' <-> ' + self.second_user.username

    def __str__(self):
        return self.first_user.username + ' <-> ' + self.second_user.username

    @staticmethod
    def get_dialogs_with(username):
        return Dialog.objects.filter(models.Q(first_user__username=username) | models.Q(second_user__username=username))

    @staticmethod
    def get_dialog_of(name1, name2):
        return Dialog.objects.filter(models.Q(first_user__username=name1, second_user__username=name2) |
                                     models.Q(first_user__username=name2, second_user__username=name1))

    def get_last_message(self):
        return Message.objects.filter(dialog=self).latest('send_date')

    def get_collocutor(self, username):
        return self.second_user if self.first_user.username == username else self.first_user

    def get_messages(self):
        return Message.objects.filter(dialog=self)

    def in_dialog(self, username):
        return bool(self.first_user.username == username or self.second_user.username == username)

    def is_schizophrenia(self):
        return self.first_user.username == self.second_user.username


class Message(models.Model):
    dialog = models.ForeignKey(Dialog)
    sender = models.ForeignKey(get_user_model(), related_name='sender')
    is_read = models.BooleanField()
    send_date = models.DateTimeField()
    text = models.TextField()

    def __unicode__(self):
        return str(self.dialog) + "[{}]".format(self.sender.username) + ": " + \
            self.send_date.strftime("%d.%m.%Y %H:%M:%S")

    def __str__(self):
        return str(self.dialog) + "[{}]".format(self.sender.username) + ": " + \
            self.send_date.strftime("%d.%m.%Y %H:%M:%S")