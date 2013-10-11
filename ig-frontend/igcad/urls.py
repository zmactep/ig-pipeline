from django.conf.urls import patterns, include, url

# Uncomment the next two lines to enable the admin:
from django.contrib import admin
admin.autodiscover()

from igcad import settings

from registration.forms import RegistrationFormUniqueEmail
from registration.backends.default.views import RegistrationView


urlpatterns = patterns(
    '',

    # Admin
    url(r'^admin/doc/', include('django.contrib.admindocs.urls')),
    url(r'^admin/', include(admin.site.urls)),

    # Auth system
    url(r'^register/$', RegistrationView.as_view(form_class=RegistrationFormUniqueEmail),
        name='registration_register'),
    url(r'', include('registration.backends.default.urls')),

    # Tasks dashboard
    url(r'^', include('igtasks.urls', namespace='igtasks')),

    # Info page content
    url(r'^info/', include('main.urls', namespace='main')),

    # User profile
    url(r'^user/', include('userprofile.urls', namespace='userprofile')),

    # Messages
    (r'^messages/', include('messages.urls', namespace='messages')),

    url(r'^static/', 'django.views.static.serve', {'document_root': settings.STATIC_ROOT}),
)
