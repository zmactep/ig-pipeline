from django.contrib import admin
from polls.models import Poll
from polls.models import Choise

ass ChoiseInline(admin.StackedInline):
    model = Choise
    extra = 3

class PollAdmin(admin.ModelAdmin):
    fieldsets = [
        (None,               {'fields': ['question']}),
        ('Date information', {'fields': ['pub_date'], 'classes': ['collapse']}),
    ]
    inlines = [ChoiseInline]

admin.site.register(Poll, PollAdmin)
