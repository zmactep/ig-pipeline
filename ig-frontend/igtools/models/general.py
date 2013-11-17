from django import forms


class CustomModelChoiceField(forms.ModelChoiceField):
    def label_from_instance(self, obj):
        return "%s - %s - %s" % (obj.file_id, obj.group, obj.run)

    class Meta:
        app_label = 'igtools'
