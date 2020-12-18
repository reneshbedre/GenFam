#views.py
# here .form is required for python 3
# form only will work with python 2
# from .form import RegistrationForm, User
# from .models import user
from django.contrib import messages
from django.contrib.auth.decorators import login_required
from django.contrib.auth import logout
from django.views.decorators.csrf import csrf_protect
# from django.shortcuts import render_to_response
from django.http import HttpResponseRedirect, HttpResponse
from django.template import RequestContext
from django.shortcuts import render
import random
# from django.core.urlresolvers import reverse
from django.urls import reverse


# @login_required
def index(request):
    return render(request, 'accounts/index_2.html', {'user': request.user})


def doc(request):
    return render(request, 'accounts/doc.html', {'user': request.user})


def dload(request):
    return render(request, 'accounts/dload.html', {'user': request.user})

def news(request):
    return render(request, 'accounts/news.html', {'user': request.user})

def register(request):
    form = RegistrationForm(request.POST or None)

    if form.is_valid():
        user = User.objects.create_user(username=form.cleaned_data['username'],
                                        password=form.cleaned_data['password1'],
                                        email=form.cleaned_data['email']
                                            )

        context = {'form': RegistrationForm(),
                   'reply': 'You Have Registered Successfully'}
        messages.success(request, "You Have Registered Successfully!")
        return HttpResponseRedirect(reverse("accounts:register"))
    else:
        # context = {'form': RegistrationForm()}
        return render(request, 'accounts/register.html', {'form':form})


def register_success(request):
    return render(request, 'accounts/success.html')

