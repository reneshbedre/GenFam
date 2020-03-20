from django.conf.urls import url
from . import views
from anot import settings
from django.conf.urls.static import static
import random


app_name = 'family'
# userid = ''.join(random.choice('0123456789ABCDEF') for i in range(16))
# userid = "MMMM"

urlpatterns = [
    url('^$', views.index, name='index'),
    url(r'^results/$', views.results, name='results'),
    url(r'^results/download_enrich/$', views.download_enrich, name='download_enrich'),
    url(r'^results/download_all/$', views.download_all, name='download_all'),
    url(r'^results/download_fig/$', views.download_fig, name='download_fig'),
] + static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT)

