from django.conf.urls import url, include
from django.contrib import admin
from anot import settings
from django.conf.urls.static import static

urlpatterns = [
    url(r'^admin/', admin.site.urls),
    # url(r'^home/', include('accounts.urls')),
    url(r'^', include('accounts.urls')),
    url(r'^family/', include('family.urls')),
] + static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT)
