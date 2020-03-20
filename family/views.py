from .form import SeqForm
from django.shortcuts import render
from django.http import HttpResponseRedirect, HttpResponse
from django.urls import reverse
from django.contrib import messages
import os
from anot.settings import BASE_DIR, MEDIA_ROOT
from wsgiref.util import FileWrapper
from custom_functions import fam_enrich, figures, figures_2
from easy_timezones.utils import is_valid_ip, is_local_ip
import shutil
from itertools import chain
from decimal import Decimal


# @login_required(login_url='/accounts/login/')
def get_ip_address_from_request(request):
    # get IP
    PRIVATE_IPS_PREFIX = ('10.', '172.', '192.', '127.')
    ip_address = ''
    x_forwarded_for = request.META.get('HTTP_X_FORWARDED_FOR', '')
    if x_forwarded_for and ',' not in x_forwarded_for:
        if not x_forwarded_for.startswith(PRIVATE_IPS_PREFIX) and is_valid_ip(x_forwarded_for):
            ip_address = x_forwarded_for.strip()
    else:
        ips = [ip.strip() for ip in x_forwarded_for.split(',')]
        for ip in ips:
            if ip.startswith(PRIVATE_IPS_PREFIX):
                continue
            elif not is_valid_ip(ip):
                continue
            else:
                ip_address = ip
                break
    if not ip_address:
        x_real_ip = request.META.get('HTTP_X_REAL_IP', '')
        if x_real_ip:
            if not x_real_ip.startswith(PRIVATE_IPS_PREFIX) and is_valid_ip(x_real_ip):
                ip_address = x_real_ip.strip()
    if not ip_address:
        remote_addr = request.META.get('REMOTE_ADDR', '')
        if remote_addr:
            if not remote_addr.startswith(PRIVATE_IPS_PREFIX) and is_valid_ip(remote_addr):
                ip_address = remote_addr.strip()
    if not ip_address:
        ip_address = '127.0.0.1'
    ip_address = ip_address.replace('.', '')
    return ip_address


def index(request):
    form = SeqForm(request.POST or None, request.FILES or None)
    if request.session.session_key is None:
        request.session.save()
    session_key = request.session.session_key
    num_visits = request.session.get('num_visits', 0)
    request.session['num_visits'] = num_visits + 1
    ip = get_ip_address_from_request(request)
    user = str(ip) + str(request.session['num_visits'])

    if form.is_valid():
        if not os.path.exists(os.path.join(MEDIA_ROOT, "input_files", str(user))):
            os.makedirs(os.path.join(MEDIA_ROOT, "input_files", str(user)))
        else:
            shutil.rmtree(os.path.join(MEDIA_ROOT, "input_files", str(user)))
            os.makedirs(os.path.join(MEDIA_ROOT, "input_files", str(user)))

        if not os.path.exists(os.path.join(MEDIA_ROOT, "output_files", str(user))):
            os.makedirs(os.path.join(MEDIA_ROOT, "output_files", str(user)))
        else:
            shutil.rmtree(os.path.join(MEDIA_ROOT, "output_files", str(user)))
            os.makedirs(os.path.join(MEDIA_ROOT, "output_files", str(user)))

        _gene_ids = form.cleaned_data['gene_ids']
        _gene_ids = _gene_ids.strip('\n')
        _plant_select = form.cleaned_data['plant_select']
        _id_type = form.cleaned_data['ID_select']
        _stat_sign_test = form.cleaned_data['stat_sign_test']
        _multi_test_corr = form.cleaned_data['multi_test_corr']
        _gene_ids_file = os.path.join(MEDIA_ROOT, "input_files", str(user), "gene_id_file.txt")
        _gene_ids_file = open(_gene_ids_file, "w")
        _count_ids = 0
        _gene_ids_file.write(_gene_ids)
        _gene_ids_file.close()

        # before running analysis check for input fields provided by user
        if _plant_select == "z":
            messages.warning(request, "Error: Select plant species",
                             extra_tags='safe')
            return render(request, 'family/index.html', {'form': form})
        elif _id_type == "z":
            messages.warning(request, "Error: Select gene ID type (If you are not sure about gene ID type, please check "
                                      "allowed IDs for each plant species)",
                             extra_tags='safe')
            return render(request, 'family/index.html', {'form': form})
        elif _stat_sign_test == "z":
            messages.warning(request, "Error: Select statistical significance test",
                             extra_tags='safe')
            return render(request, 'family/index.html', {'form': form})
        elif _multi_test_corr == "z":
            messages.warning(request, "Error: Select multiple testing correction method",
                             extra_tags='safe')
            return render(request, 'family/index.html', {'form': form})

        # common_function_2.fam_enrich(os.path.join(BASE_DIR, "input_files", str(user), "gene_id_file.txt"),
        #                             _plant_select, _stat_sign_test, _multi_test_corr, str(user))
        value = fam_enrich.fam_enrich(os.path.join(MEDIA_ROOT, "input_files", str(user), "gene_id_file.txt"),
                                     _plant_select, _stat_sign_test, _multi_test_corr, str(user), _id_type)

        if value == "okay":
            return HttpResponseRedirect(reverse("family:results"))
        else:
            messages.warning(request, "Analysis stopped due to following reasons:"
                                      "<ul><li>Check if correct plant species is selected</li>"
                                      "<li>Check the number of input IDs "
                                      "(Mappable input IDs must be greater than 5 for reliable results)</li>"
                                      "<li>Check if correct ID type selected for given plant species</li>"
                                      "<li>Check if "
                                      "correct version of plant IDs are selected from Phytozome database</li></ul>",
                             extra_tags='safe')
            return render(request, 'family/index.html', {'form': form})
    else:
        return render(request, 'family/index.html', {'error': 'Error: Provide Input Ids or Upload ID File',
                                                     'form': form})


def results(request):
    ip = get_ip_address_from_request(request)
    user = str(ip) + str(request.session['num_visits'])
    _fam_enrich_out_file = open(os.path.join(MEDIA_ROOT, "output_files", str(user), "fam_enrich_out.txt"), "rU")
    _fam_all_out_file = open(os.path.join(MEDIA_ROOT, "output_files", str(user), "fam_all_out.txt"), "rU")
    _fam_enrich_content_read = _fam_enrich_out_file.readlines()
    _fam_all_content_read = _fam_all_out_file.readlines()
    _fam_enrich_content = []
    _fam_all_content = []
    line_count = 0
    cell_ct_per = 0

    for line in _fam_enrich_content_read:
        line = line.strip().split('\t')
        if not line[0].startswith("Gene"):
            _fam_enrich_content.append([line[0], line[1], '%.2E' % Decimal(float(line[8])), '%.2E' % Decimal(float(line[9]))])

        query_id_count = line[2]
        # get the annotation count. number of genes from user input present in genfam database
        _anot_count = line[3]
        line_count += 1
        # check if cell frequency column is present and count the tests if it is failed due to expected cell count < 5
        if len(line) == 15 and 'Cells' not in line[14] and float(line[14]) > float(20):
            cell_ct_per += 1

    for line in _fam_all_content_read:
        line = line.strip().split('\t')
        if not line[0].startswith("Gene"):
            _fam_all_content.append([line[0], line[1], '%.2E' % Decimal(float(line[8])), '%.2E' % Decimal(float(line[9]))])

        query_id_count = line[2]
        _anot_count = line[3]
        # check if cell frequency column is present and count the tests if it is failed due to expected cell count < 5
        # this is only for chi test; chi test detected using column count
        if len(line) == 15 and 'Cells' not in line[14] and float(line[14]) > float(20):
            cell_ct_per += 1

    # delete first 3 elements which are column names to avoid outputting them on web page
    # sort list based on fdr
    _fam_enrich_content = sorted(_fam_enrich_content, key=lambda x: float(x[3]))
    # flatten 2d to 1d
    _fam_enrich_content = list(chain.from_iterable(_fam_enrich_content))

    # print(_fam_all_content, "xx")
    # _fam_enrich_content = _fam_enrich_content[4:]
    # sort list based on fdr
    _fam_all_content = sorted(_fam_all_content, key=lambda x: float(x[3]))
    # print(_fam_all_content, "yy")
    # _fam_all_content = _fam_all_content[4:]
    # flatten 2d to 1d
    _fam_all_content = list(chain.from_iterable(_fam_all_content))

    # generate figures
    # common_function_3.bar_enrich(os.path.join(BASE_DIR, "output_files", str(user), "fam_enrich_out.txt"), str(user))
    # allow figure if enrichment terms found
    if line_count > 1:
        figures_2.bar_enrich(os.path.join(MEDIA_ROOT, "output_files", str(user), "fam_enrich_out.txt"), str(user))

    _fam_enrich_out_file.close()
    _fam_all_out_file.close()

    context = {'form': SeqForm(),
               'header': 'Success! Data Analysis Completed',
               'results': 'Get Results Table',
               'figures': os.path.join(BASE_DIR, "output_files", str(user), "fam_enrich_out.png"),
               'go_all_content': _fam_all_content,
               'go_enrich_content': _fam_enrich_content,
               'out_file': os.path.join(MEDIA_ROOT, "output_files", "go_enrich_out.txt"),
               'user': user,
               'anot_count': _anot_count,
               'query_id_count': query_id_count,
               'line_count': line_count,
               'chi_cell_count': cell_ct_per,
               #'plant': "xx"
               }
    return render(request, 'family/results.html', context)


def download_enrich(request):
    ip = get_ip_address_from_request(request)
    user = str(ip) + str(request.session['num_visits'])
    path_to_file = os.path.join(MEDIA_ROOT, "output_files", str(user), "fam_enrich_out.txt")
    with open(path_to_file, 'r') as f:
        myfile = f.read()
    response = HttpResponse(myfile, content_type='application/force-download')
    # response['Content-Disposition'] = 'attachment; filename=' + "fam_enrich_out.txt"
    response['Content-Disposition'] = 'attachment; filename=' + os.path.basename(path_to_file)
    return response


def download_all(request):
    ip = get_ip_address_from_request(request)
    user = str(ip) + str(request.session['num_visits'])
    path_to_file = os.path.join(MEDIA_ROOT, "output_files", str(user), "fam_all_out.txt")
    with open(path_to_file, 'r') as f:
        myfile = f.read()
    response = HttpResponse(myfile, content_type='application/force-download')
    response['Content-Disposition'] = 'attachment; filename=' + os.path.basename(path_to_file)
    # response['X-Sendfile'] = smart_str(path_to_file)
    return response


def download_fig(request):
    ip = get_ip_address_from_request(request)
    user = str(ip) + str(request.session['num_visits'])
    path_to_file = os.path.join(MEDIA_ROOT, "output_files", str(user), "fam_enrich_out.png")
    wrapper = FileWrapper(open(path_to_file, "rb"))
    response = HttpResponse(wrapper, content_type='application/force-download')
    response['Content-Disposition'] = 'attachment; filename=' + os.path.basename(path_to_file)
    return response
