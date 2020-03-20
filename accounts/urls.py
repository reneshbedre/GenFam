from django.conf.urls import url, include
from . import views
from django.contrib.auth import views as auth_views

app_name = 'accounts'

urlpatterns = [
    url(r'^$', views.index, name='index'),
    url(r'^doc/$', views.doc, name='doc'),
    url(r'^dload/$', views.dload, name='dload'),
    url(r'^news/$', views.news, name='news'),
    url(r'^login/$', auth_views.LoginView, {'template_name': 'accounts/login.html'}, name='login'),
    url(r'^register/$', views.register, name='register'),
    url(r'^register/success/$', views.register_success, name='reg_suc'),
    # url(r'^logout/$', views.logout_page, {'next_page': 'accounts/index.html'}, name='logout'),
    url(r'^logout/$', auth_views.LogoutView, {'next_page': '/home/login'}, name='logout'),
]
