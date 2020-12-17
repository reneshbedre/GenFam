#!/usr/bin/python

from __future__ import division
import os
import psycopg2
# for python 3.6
from . import go_enrich
# for python 2.7
import scipy.stats as stats
from anot.settings import BASE_DIR, MEDIA_ROOT


def enrichment_analysis(_uniq_id_count_dict, _get_user_id_count_for_gene_fam, _gene_fam_count_dict, _bg_gene_count,
                        _gene_fam_dict, _stat_sign_test, _multi_test_corr, _uniq_p, _uniq_f, _uniq_c,
                        _get_gene_ids_from_user_dict, _short_fam_dict):
    pvalues = []
    enrichment_result = []
    # get total mapped genes from user list
    mapped_query_ids = sum(_get_user_id_count_for_gene_fam.values())

    for key in _get_user_id_count_for_gene_fam.keys():
        # first row and first column (k)
        # from user list
        # number of genes mapped to family in subset
        gene_in_category = _get_user_id_count_for_gene_fam[key]
        short_fam = _short_fam_dict[key]
        # first row and second column (m-k)  _uniq_id_count_dict (m)
        # n-k (n total mappable deg)
        # from user list
        # total subset - gene mapped to family
        # gene_not_in_category_but_in_sample = _uniq_id_count_dict-gene_in_category
        # calculate stat sign based on only ids mapped and not all query ids as in agrigo
        gene_not_in_category_but_in_sample = mapped_query_ids - gene_in_category
        # second row and first column (n-k)
        # m-k
        # genes that mapped to family are in genome - genes in category from user list
        gene_not_in_catgory_but_in_genome = _gene_fam_count_dict[key]-gene_in_category
        bg_gene_fam_ids = _gene_fam_count_dict[key]
        # second row and second column gene_in_category+gene_not_in_catgory_and_in_genome (n)
        # N-m-n+k
        bg_in_genome = \
            _bg_gene_count-mapped_query_ids-(gene_in_category+gene_not_in_catgory_but_in_genome)+gene_in_category

        # for go terms
        process = _uniq_p[key]
        function = _uniq_f[key]
        comp = _uniq_c[key]

        # gene ids from user list mapped to particular gene family
        gene_ids = _get_gene_ids_from_user_dict[key]

        if _stat_sign_test == "a":
            # fisher exact test
            '''
            if "Glycoside" in key:
                print(key)
                print(gene_in_category, gene_not_in_category_but_in_sample, gene_not_in_catgory_but_in_genome, bg_in_genome)
            '''

            # run analysis if only mappable ids are >= 5
            # if mapped_query_ids >= 5:
            oddsratio, pvalue = stats.fisher_exact([[gene_in_category, gene_not_in_category_but_in_sample],
                                                    [gene_not_in_catgory_but_in_genome, bg_in_genome]], 'greater')
            if int(gene_in_category) > 0:
                '''
                enrichment_result.append([key, _uniq_id_count_dict, gene_in_category, bg_gene_fam_ids, _bg_gene_count,
                                          gene_in_category, gene_not_in_category_but_in_sample,
                                          gene_not_in_catgory_but_in_genome, bg_in_genome, oddsratio, pvalue, process,
                                          function, comp, gene_ids])
                '''
                enrichment_result.append(
                    [key, short_fam, _uniq_id_count_dict, gene_in_category, _bg_gene_count, bg_gene_fam_ids, oddsratio,
                     pvalue, process, function, comp, gene_ids])
                pvalues.append(pvalue)
        # Hypergeometric
        elif _stat_sign_test == "b":
            '''
            if "Homeobox" in key:
                print(key)
                print(gene_in_category, _bg_gene_count, _gene_fam_count_dict[key], mapped_query_ids)
            '''
            # run analysis if only mappable ids are >= 5
            # if mapped_query_ids >= 5:
            pvalue = stats.hypergeom.sf(gene_in_category - 1, _bg_gene_count, _gene_fam_count_dict[key],
                                        mapped_query_ids)
            # pvalue = stats.hypergeom.sf(gene_in_category-1, bg_in_genome,
            #                            gene_not_in_catgory_but_in_genome, gene_not_in_category_but_in_sample)
            oddsratio = "NA"
            if int(gene_in_category) > 0:
                '''
                enrichment_result.append([key, _uniq_id_count_dict, gene_in_category, bg_gene_fam_ids, _bg_gene_count,
                                          gene_in_category, gene_not_in_category_but_in_sample,
                                          gene_not_in_catgory_but_in_genome, bg_in_genome, oddsratio, pvalue, process,
                                          function, comp, gene_ids])
                '''
                enrichment_result.append(
                    [key, short_fam, _uniq_id_count_dict, gene_in_category, _bg_gene_count, bg_gene_fam_ids, oddsratio,
                     pvalue, process, function, comp, gene_ids])
                pvalues.append(pvalue)
        # Binomial
        elif _stat_sign_test == "c":
            # probability from the reference set for particular category
            '''
            if "Homeobox" in key:
                print(key)
                print(_gene_fam_count_dict[key], _bg_gene_count, gene_in_category, mapped_query_ids)
            '''
            # run analysis if only mappable ids are >= 5
            # if mapped_query_ids >= 5:
            exp_pvalue = (_gene_fam_count_dict[key]/_bg_gene_count)
                # pvalue = stats.binom_test(gene_in_category, _uniq_id_count_dict, exp_pvalue, alternative='greater')
                # consider only mapped query ids
            pvalue = stats.binom_test(gene_in_category, mapped_query_ids, exp_pvalue, alternative='greater')
            oddsratio = "NA"
            # print key, exp_pvalue, _gene_fam_count_dict[key], _bg_gene_count, gene_in_category, _uniq_id_count_dict
            if int(gene_in_category) > 0:
                '''
                enrichment_result.append([key, _uniq_id_count_dict, gene_in_category, bg_gene_fam_ids, _bg_gene_count,
                                          gene_in_category, gene_not_in_category_but_in_sample,
                                          gene_not_in_catgory_but_in_genome, bg_in_genome, oddsratio, pvalue, process,
                                          function, comp, gene_ids])
                '''
                enrichment_result.append(
                    [key, short_fam, _uniq_id_count_dict, gene_in_category, _bg_gene_count, bg_gene_fam_ids, oddsratio,
                     pvalue, process, function, comp, gene_ids])
                pvalues.append(pvalue)
        # Chi-Square
        elif _stat_sign_test == "d":
            #   if "Glutathione" in key:
                 # print(key)
                 # print(gene_in_category, gene_not_in_category_but_in_sample, gene_not_in_catgory_but_in_genome,
                 #      bg_in_genome)
            # run analysis if only mappable ids are >= 5
            # if mapped_query_ids >=5:
            chi2, pvalue, dof, exp = stats.chi2_contingency([[gene_in_category, gene_not_in_category_but_in_sample],
                                                             [gene_not_in_catgory_but_in_genome, bg_in_genome]],
                                                            correction=False)
                # count cells where expected frequency < 5
            cell_ct = sum(sum(i<5 for i in exp))
            cell_ct_per = 100* float(cell_ct)/float(4)
                # print(key, exp, cell_ct, cell_ct_per, "count")
            oddsratio = "NA"
            # only report the family where observed count is more than expected gene in category count
            # this is for getting highly enriched genes other under-represented genes will also be found
            if int(gene_in_category) > 0 and gene_in_category >= exp[0][0]:
                '''
                enrichment_result.append([key, _uniq_id_count_dict, gene_in_category, bg_gene_fam_ids, gene_in_category,
                                          _bg_gene_count, gene_not_in_category_but_in_sample,
                                          gene_not_in_catgory_but_in_genome, bg_in_genome, chi2, pvalue, process,
                                          function, comp, gene_ids])
                '''
                enrichment_result.append(
                    [key, short_fam, _uniq_id_count_dict, gene_in_category, _bg_gene_count, bg_gene_fam_ids, oddsratio,
                     pvalue, process, function, comp, gene_ids, cell_ct_per])
                # print key, chi2, pvalue
                pvalues.append(pvalue)

    # FDR Bonferroni
    if _multi_test_corr == "a":
        fdr = go_enrich.correct_pvalues_for_multiple_testing(pvalues, "Bonferroni")
    # FDR Bonferroni-Holm
    elif _multi_test_corr == "b":
        fdr = go_enrich.correct_pvalues_for_multiple_testing(pvalues, "Bonferroni-Holm")
    # FDR Benjamini-Hochberg
    elif _multi_test_corr == "c":
        # print(pvalues)
        fdr = go_enrich.correct_pvalues_for_multiple_testing(pvalues, "Benjamini-Hochberg")
    return enrichment_result, fdr, mapped_query_ids


def get_bg_gene_count(gfam_group):
    try:
        # con = psycopg2.connect(database='family_db_72019', user='', password='', host="127.0.0.1", port="5432")
        # con = psycopg2.connect(database='family_db_72019', user='mlab', password='mlab92', host="127.0.0.1", port="5432")
        con = psycopg2.connect(database='family_db_72019', user='super', password='ps92712268',
                               host='mandadilab-1984.postgres.pythonanywhere-services.com', port='11984')
    except ValueError:
        print("Not able to connect database")
    cur = con.cursor()
    q_loc = "select sum(loc_len) from {}".format(gfam_group)
    cur.execute(q_loc)
    bg_gene_count = cur.fetchone()[0]
    t_loc = "select sum(trn_len) from {}".format(gfam_group)
    cur.execute(t_loc)
    bg_trn_count = cur.fetchone()[0]
    p_loc = "select sum(phyt_id_len) from {}".format(gfam_group)
    cur.execute(p_loc)
    bg_phytid_count = cur.fetchone()[0]
    return bg_gene_count, bg_trn_count, bg_phytid_count


def fam_enrich(_id_file, _plant_select, _stat_sign_test, _multi_test_corr, user, _id_type):
    # input id file obtained from form
    id_file = open(_id_file, "rU")

    if _plant_select == "z":
        # stop the analysis as no plant species selected
        return "stop"

    # create the folder if it doesn't exist.
    try:
        # os.mkdir(os.path.join(MEDIA_ROOT, "output_files"))
        os.mkdir(os.path.join(MEDIA_ROOT, "output_files", user))
    except:
        pass

    # connect to psql db
    try:
       # con = psycopg2.connect(database='plant_db', user='renesh', password='', host="127.0.0.1", port="5432")
       # con = psycopg2.connect(database='family_db', user='renesh', password='', host="127.0.0.1", port="5432")
       # con = psycopg2.connect(database='family_db_62018', user='renesh', password='', host="127.0.0.1", port="5432")
       # con = psycopg2.connect(database='family_db_12019', user='renesh', password='', host="127.0.0.1",
       #                       port="5432")
       # con = psycopg2.connect(database='family_db_32019', user='renesh', password='', host="127.0.0.1", port="5432")
       # con = psycopg2.connect(database='family_db_32019', user='denise', password='', host="127.0.0.1", port="5432")
       # con = psycopg2.connect(database='family_db_72019', user='denise', password='', host="127.0.0.1", port="5432")
       # con = psycopg2.connect(database='family_db_72019', user='mlab', password='mlab92', host="127.0.0.1", port="5432")
       con = psycopg2.connect(database='family_db_72019', user='super', password='ps92712268',
                              host='mandadilab-1984.postgres.pythonanywhere-services.com', port='11984')
    except ValueError:
        print("Not able to connect database")

    cur = con.cursor()
    query1 = ''
    bg_gene_count = 0
    plant_name = ""

    # dicots
    if _plant_select == 'e':
        query1 = "select * from atha_phyt12_gene_fam_group"
        # bg_gene_count = 6650
        # bg_gene_count = 8579
        # bg_trn_count = 11042
        # bg_phytid_count = 11042
        bg_gene_count, bg_trn_count, bg_phytid_count = get_bg_gene_count("atha_phyt12_gene_fam_group")
        plant_name = "Arabidopsis thaliana"
    elif _plant_select == 'ae':
        # query1 = "select * from graimon_phyt12_gene_fam_group"
        # bg_gene_count = 10408
        query1 = "select * from grai_phyt12_gene_fam_group"
        # bg_gene_count = 12134
        # bg_gene_count = 37506
        # bg_trn_count = 24436
        # bg_phytid_count = 24436
        bg_gene_count, bg_trn_count, bg_phytid_count = get_bg_gene_count("grai_phyt12_gene_fam_group")
        # print(bg_gene_count, bg_trn_count, bg_phytid_count)
        plant_name = "Gossypium raimondii"
    elif _plant_select == 'ad':
        # query1 = "select * from ghirs_phyt12_gene_fam_group"
        # bg_gene_count = 18642
        query1 = "select * from ghir_phyt12_gene_fam_group"
        # bg_gene_count = 21453
        # bg_trn_count = 27563
        # bg_phytid_count = 27563
        bg_gene_count, bg_trn_count, bg_phytid_count = get_bg_gene_count("ghir_phyt12_gene_fam_group")
        plant_name = "Gossypium hirsutum"
    elif _plant_select == 'af':
        query1 = "select * from gmax_phyt12_gene_fam_group"
        # bg_gene_count = 14610
        # bg_gene_count = 16984
        # bg_trn_count = 27206
        # bg_phytid_count = 27206
        bg_gene_count, bg_trn_count, bg_phytid_count = get_bg_gene_count("gmax_phyt12_gene_fam_group")
        plant_name = "Glycine max"
    elif _plant_select == 'a':
        # query1 = "select * from ahypo_phyt12_gene_fam_group"
        # bg_gene_count = 6238
        query1 = "select * from ahyp_phyt12_gene_fam_group"
        # bg_gene_count = 7559
        # bg_trn_count = 7570
        # bg_phytid_count = 7570
        bg_gene_count, bg_trn_count, bg_phytid_count = get_bg_gene_count("ahyp_phyt12_gene_fam_group")
        plant_name = "Amaranthus hypochondriacus"
    elif _plant_select == 'c':
        # query1 = "select * from aoccidentale_phyt12_gene_fam_group"
        # bg_gene_count = 11236
        query1 = "select * from aocc_phyt12_gene_fam_group"
        # bg_gene_count = 12908
        # bg_trn_count = 25530
        # bg_phytid_count = 25530
        bg_gene_count, bg_trn_count, bg_phytid_count = get_bg_gene_count("aocc_phyt12_gene_fam_group")
        plant_name = "Anacardium occidentale"
    elif _plant_select == 'h':
        # query1 = "select * from acoerulea_phyt12_gene_fam_group"
        # bg_gene_count = 6668
        query1 = "select * from acoe_phyt12_gene_fam_group"
        # bg_gene_count = 8666
        # bg_trn_count = 12677
        # bg_phytid_count = 12677
        bg_gene_count, bg_trn_count, bg_phytid_count = get_bg_gene_count("acoe_phyt12_gene_fam_group")
        plant_name = "Aquilegia coerulea"
    elif _plant_select == 'f':
        # query1 = "select * from ahalleri_phyt12_gene_fam_group"
        # bg_gene_count = 6131
        query1 = "select * from ahal_phyt12_gene_fam_group"
        # bg_gene_count = 8279
        # bg_trn_count = 8869
        # bg_phytid_count = 8869
        bg_gene_count, bg_trn_count, bg_phytid_count = get_bg_gene_count("ahal_phyt12_gene_fam_group")
        plant_name = "Arabidopsis halleri"
    elif _plant_select == 'g':
        # query1 = "select * from alyrata_phyt12_gene_fam_group"
        # bg_gene_count = 7383
        query1 = "select * from alyr_phyt12_gene_fam_group"
        # bg_gene_count = 9866
        # bg_trn_count = 10473
        # bg_phytid_count = 10473
        bg_gene_count, bg_trn_count, bg_phytid_count = get_bg_gene_count("alyr_phyt12_gene_fam_group")
        plant_name = "Arabidopsis lyrata"
    elif _plant_select == 'l':
        # query1 = "select * from boleraceacapitata_phyt12_gene_fam_group"
        # bg_gene_count = 9395
        query1 = "select * from bole_phyt12_gene_fam_group"
        # bg_gene_count = 11558
        # bg_trn_count = 11558
        # bg_phytid_count = 11558
        bg_gene_count, bg_trn_count, bg_phytid_count = get_bg_gene_count("bole_phyt12_gene_fam_group")
        plant_name = "Brassica oleracea"
    elif _plant_select == 'm':
        # query1 = "select * from brapa_phyt12_gene_fam_group"
        # bg_gene_count = 11019
        query1 = "select * from brap_phyt12_gene_fam_group"
        # bg_gene_count = 14179
        # bg_trn_count = 15247
        # bg_phytid_count = 15247
        bg_gene_count, bg_trn_count, bg_phytid_count = get_bg_gene_count("brap_phyt12_gene_fam_group")
        plant_name = "Brassica rapa"
    elif _plant_select == 'n':
        # query1 = "select * from cgrandiflora_phyt12_gene_fam_group"
        # bg_gene_count = 6408
        query1 = "select * from cgra_phyt12_gene_fam_group"
        # bg_gene_count = 8459
        # bg_trn_count = 9041
        # bg_phytid_count = 9041
        bg_gene_count, bg_trn_count, bg_phytid_count = get_bg_gene_count("cgra_phyt12_gene_fam_group")
        plant_name = "Capsella grandiflora"
    elif _plant_select == 'o':
        # query1 = "select * from crubella_phyt12_gene_fam_group"
        # bg_gene_count = 6802
        query1 = "select * from crub_phyt12_gene_fam_group"
        # bg_gene_count = 10227
        # bg_trn_count = 10911
        # bg_phytid_count = 10911
        bg_gene_count, bg_trn_count, bg_phytid_count = get_bg_gene_count("crub_phyt12_gene_fam_group")
        plant_name = "Capsella rubella"
    elif _plant_select == 'p':
        # query1 = "select * from cpapaya_phyt12_gene_fam_group"
        # bg_gene_count = 5916
        query1 = "select * from cpap_phyt12_gene_fam_group"
        # bg_gene_count = 9093
        # bg_trn_count = 9680
        # bg_phytid_count = 9680
        bg_gene_count, bg_trn_count, bg_phytid_count = get_bg_gene_count("cpap_phyt12_gene_fam_group")
        plant_name = "Carica papaya"
    elif _plant_select == 'q':
        # query1 = "select * from cquinoa_phyt12_gene_fam_group"
        # bg_gene_count = 11089
        query1 = "select * from cqui_phyt12_gene_fam_group"
        # bg_gene_count = 13879
        # bg_trn_count = 13879
        # bg_phytid_count = 13879
        bg_gene_count, bg_trn_count, bg_phytid_count = get_bg_gene_count("cqui_phyt12_gene_fam_group")
        plant_name = "Chenopodium quinoa"
    elif _plant_select == 't':
        # query1 = "select * from cclementina_phyt12_gene_fam_group"
        # bg_gene_count = 6584
        query1 = "select * from ccle_phyt12_gene_fam_group"
        # bg_gene_count = 8198
        # bg_trn_count = 11136
        # bg_phytid_count = 11136
        bg_gene_count, bg_trn_count, bg_phytid_count = get_bg_gene_count("ccle_phyt12_gene_fam_group")
        plant_name = "Citrus clementina"
    elif _plant_select == 'u':
        # query1 = "select * from csinensis_phyt12_gene_fam_group"
        # bg_gene_count = 6671
        query1 = "select * from csin_phyt12_gene_fam_group"
        # bg_gene_count = 8223
        # bg_trn_count = 14959
        # bg_phytid_count = 14959
        bg_gene_count, bg_trn_count, bg_phytid_count = get_bg_gene_count("csin_phyt12_gene_fam_group")
        plant_name = "Citrus sinensis"
    elif _plant_select == 'w':
        # query1 = "select * from dcarota_phyt12_gene_fam_group"
        # bg_gene_count = 8483
        query1 = "select * from dcar_phyt12_gene_fam_group"
        # bg_gene_count = 10393
        # bg_trn_count = 10394
        # bg_phytid_count = 10394
        bg_gene_count, bg_trn_count, bg_phytid_count = get_bg_gene_count("dcar_phyt12_gene_fam_group")
        plant_name = "Daucus carota"
    elif _plant_select == 'y':
        # query1 = "select * from egrandis_phyt12_gene_fam_group"
        # bg_gene_count = 9728
        query1 = "select * from egra_phyt12_gene_fam_group"
        # bg_gene_count = 11935
        # bg_trn_count = 15008
        # bg_phytid_count = 15008
        bg_gene_count, bg_trn_count, bg_phytid_count = get_bg_gene_count("egra_phyt12_gene_fam_group")
        plant_name = "Eucalyptus grandis"
    elif _plant_select == 'ab':
        # query1 = "select * from esalsugineum_phyt12_gene_fam_group"
        # bg_gene_count = 6669
        query1 = "select * from esal_phyt12_gene_fam_group"
        # bg_gene_count = 8896
        # bg_trn_count = 9830
        # bg_phytid_count = 9830
        bg_gene_count, bg_trn_count, bg_phytid_count = get_bg_gene_count("esal_phyt12_gene_fam_group")
        plant_name = "Eutrema salsugineum"
    elif _plant_select == 'ac':
        # query1 = "select * from fvesca_phyt12_gene_fam_group"
        # bg_gene_count = 6593
        query1 = "select * from fves_phyt12_gene_fam_group"
        # bg_gene_count = 8573
        # bg_trn_count = 8573
        # bg_phytid_count = 8573
        bg_gene_count, bg_trn_count, bg_phytid_count = get_bg_gene_count("fves_phyt12_gene_fam_group")
        plant_name = "Fragaria vesca"
    elif _plant_select == 'ah':
        # query1 = "select * from kfedtschenkoi_phyt12_gene_fam_group"
        # bg_gene_count = 8027
        query1 = "select * from kfed_phyt12_gene_fam_group"
        # bg_gene_count = 9803
        # bg_trn_count = 14292
        # bg_phytid_count = 14292
        bg_gene_count, bg_trn_count, bg_phytid_count = get_bg_gene_count("kfed_phyt12_gene_fam_group")
        plant_name = "Kalanchoe fedtschenkoi"
    elif _plant_select == 'ai':
        # query1 = "select * from klaxiflora_phyt12_gene_fam_group"
        # bg_gene_count = 13990
        query1 = "select * from klax_phyt12_gene_fam_group"
        # bg_gene_count = 17208
        # bg_trn_count = 23134
        # bg_phytid_count = 23134
        bg_gene_count, bg_trn_count, bg_phytid_count = get_bg_gene_count("klax_phyt12_gene_fam_group")
        plant_name = "Kalanchoe laxiflora"
    elif _plant_select == 'aj':
        # query1 = "select * from lusitatissimum_phyt12_gene_fam_group"
        # bg_gene_count = 11634
        query1 = "select * from lusi_phyt12_gene_fam_group"
        # bg_gene_count = 14268
        # bg_trn_count = 14274
        # bg_phytid_count = 14274
        bg_gene_count, bg_trn_count, bg_phytid_count = get_bg_gene_count("lusi_phyt12_gene_fam_group")
        plant_name = "Linum usitatissimum"
    elif _plant_select == 'ak':
        # query1 = "select * from mdomestica_phyt12_gene_fam_group"
        # bg_gene_count = 16572
        query1 = "select * from mdom_phyt12_gene_fam_group"
        # bg_gene_count = 20053
        # bg_trn_count = 20054
        # bg_phytid_count = 20054
        bg_gene_count, bg_trn_count, bg_phytid_count = get_bg_gene_count("mdom_phyt12_gene_fam_group")
        plant_name = "Malus domestica"
    elif _plant_select == 'al':
        # query1 = "select * from mesculenta_phyt12_gene_fam_group"
        # bg_gene_count = 8559
        query1 = "select * from mesc_phyt12_gene_fam_group"
        # bg_gene_count = 10651
        # bg_trn_count = 13205
        # bg_phytid_count = 13205
        bg_gene_count, bg_trn_count, bg_phytid_count = get_bg_gene_count("mesc_phyt12_gene_fam_group")
        plant_name = "Manihot esculenta"
    elif _plant_select == 'an':
        # query1 = "select * from mtruncatula_phyt12_gene_fam_group"
        # bg_gene_count = 10171
        query1 = "select * from mtru_phyt12_gene_fam_group"
        # bg_gene_count = 13148
        # bg_trn_count = 16854
        # bg_phytid_count = 16854
        bg_gene_count, bg_trn_count, bg_phytid_count = get_bg_gene_count("mtru_phyt12_gene_fam_group")
        plant_name = "Medicago truncatula"
    elif _plant_select == 'ap':
        # query1 = "select * from mguttatus_phyt12_gene_fam_group"
        # bg_gene_count = 7420
        query1 = "select * from mgut_phyt12_gene_fam_group"
        # bg_gene_count = 9601
        # bg_trn_count = 11273
        # bg_phytid_count = 11273
        bg_gene_count, bg_trn_count, bg_phytid_count = get_bg_gene_count("mgut_phyt12_gene_fam_group")
        plant_name = "Mimulus guttatus"
    elif _plant_select == 'av':
        # query1 = "select * from pvulgaris_phyt12_gene_fam_group"
        # bg_gene_count = 7695
        query1 = "select * from pvul_phyt12_gene_fam_group"
        # bg_gene_count = 9544
        # bg_trn_count = 12638
        # bg_phytid_count = 12638
        bg_gene_count, bg_trn_count, bg_phytid_count = get_bg_gene_count("pvul_phyt12_gene_fam_group")
        plant_name = "Mimulus guttatus"
    elif _plant_select == 'ax':
        # query1 = "select * from pdeltoides_phyt12_gene_fam_group"
        # bg_gene_count = 11141
        query1 = "select * from pdel_phyt12_gene_fam_group"
        # bg_gene_count = 13858
        # bg_trn_count = 17973
        # bg_phytid_count = 17973
        bg_gene_count, bg_trn_count, bg_phytid_count = get_bg_gene_count("pdel_phyt12_gene_fam_group")
        plant_name = "Populus deltoides"
    elif _plant_select == 'ay':
        # query1 = "select * from ptrichocarpa_phyt12_gene_fam_group"
        # bg_gene_count = 10551
        query1 = "select * from ptri_phyt12_gene_fam_group"
        # bg_gene_count = 13120
        # bg_trn_count = 19761
        #bg_phytid_count = 19761
        bg_gene_count, bg_trn_count, bg_phytid_count = get_bg_gene_count("ptri_phyt12_gene_fam_group")
        plant_name = "Populus trichocarpa"
    elif _plant_select == 'az':
        # query1 = "select * from ppersica_phyt12_gene_fam_group"
        # bg_gene_count = 6804
        query1 = "select * from pper_phyt12_gene_fam_group"
        # bg_gene_count = 8523
        # bg_trn_count = 14720
        # bg_phytid_count = 14720
        bg_gene_count, bg_trn_count, bg_phytid_count = get_bg_gene_count("pper_phyt12_gene_fam_group")
        plant_name = "Prunus persica"
    elif _plant_select == 'ba':
        # query1 = "select * from rcommunis_phyt12_gene_fam_group"
        # bg_gene_count = 6396
        query1 = "select * from rcom_phyt12_gene_fam_group"
        # bg_gene_count = 7988
        # bg_trn_count = 7988
        # bg_phytid_count = 7988
        bg_gene_count, bg_trn_count, bg_phytid_count = get_bg_gene_count("rcom_phyt12_gene_fam_group")
        plant_name = "Ricinus communis"
    elif _plant_select == 'bb':
        # query1 = "select * from spurpurea_phyt12_gene_fam_group"
        # bg_gene_count = 9995
        query1 = "select * from spur_phyt12_gene_fam_group"
        # bg_gene_count = 12454
        # bg_trn_count = 20283
        # bg_phytid_count = 20283
        bg_gene_count, bg_trn_count, bg_phytid_count = get_bg_gene_count("spur_phyt12_gene_fam_group")
        plant_name = "Salix purpurea"
    elif _plant_select == 'bg':
        # query1 = "select * from slycopersicum_phyt12_gene_fam_group"
        # bg_gene_count = 8112
        query1 = "select * from slyc_phyt12_gene_fam_group"
        # bg_gene_count = 9979
        # bg_trn_count = 9979
        # bg_phytid_count = 9979
        bg_gene_count, bg_trn_count, bg_phytid_count = get_bg_gene_count("slyc_phyt12_gene_fam_group")
        plant_name = "Solanum lycopersicum"
    elif _plant_select == 'bf':
        # query1 = "select * from stuberosum_phyt12_gene_fam_group"
        # bg_gene_count = 8217
        query1 = "select * from stub_phyt12_gene_fam_group"
        # bg_gene_count = 10059
        # bg_trn_count = 15392
        # bg_phytid_count = 15392
        bg_gene_count, bg_trn_count, bg_phytid_count = get_bg_gene_count("stub_phyt12_gene_fam_group")
        plant_name = "Solanum tuberosum"
    elif _plant_select == 'bk':
        # query1 = "select * from tcacao_phyt12_gene_fam_group"
        # bg_gene_count = 6764
        query1 = "select * from tcac_phyt12_gene_fam_group"
        # bg_gene_count = 8349
        # bg_trn_count = 13101
        # bg_phytid_count = 13101
        bg_gene_count, bg_trn_count, bg_phytid_count = get_bg_gene_count("tcac_phyt12_gene_fam_group")
        plant_name = "Theobroma cacao"
    elif _plant_select == 'bl':
        # query1 = "select * from tpratense_phyt12_gene_fam_group"
        # bg_gene_count = 10477
        query1 = "select * from tpra_phyt12_gene_fam_group"
        # bg_gene_count = 12378
        # bg_trn_count = 12729
        # bg_phytid_count = 12729
        bg_gene_count, bg_trn_count, bg_phytid_count = get_bg_gene_count("tpra_phyt12_gene_fam_group")
        plant_name = "Trifolium pratense"
    elif _plant_select == 'bn':
        # query1 = "select * from vvinifera_phyt12_gene_fam_group"
        # bg_gene_count = 6728
        query1 = "select * from vvin_phyt12_gene_fam_group"
        # bg_gene_count = 8124
        # bg_trn_count = 8124
        # bg_phytid_count = 8124
        bg_gene_count, bg_trn_count, bg_phytid_count = get_bg_gene_count("vvin_phyt12_gene_fam_group")
        plant_name = "Vitis vinifera"
    # elif _plant_select == 'bq':
    #    query1 = "select * from acoe_phyt12_gene_fam_group"
    #    bg_gene_count = 6668
    elif _plant_select == 'bs':
        query1 = "select * from caer_phyt12_gene_fam_group"
        # bg_gene_count = 6957
        # bg_trn_count = 15238
        # bg_phytid_count = 15238
        bg_gene_count, bg_trn_count, bg_phytid_count = get_bg_gene_count("caer_phyt12_gene_fam_group")
        plant_name = "Vitis vinifera"
    elif _plant_select == 'bu':
        # query1 = "select * from lsat_phyt12_gene_fam_group"
        # bg_gene_count = 9394
        query1 = "select * from lsat_phyt12_gene_fam_group"
        # bg_gene_count = 11692
        # bg_trn_count = 19805
        # bg_phytid_count = 19805
        bg_gene_count, bg_trn_count, bg_phytid_count = get_bg_gene_count("lsat_phyt12_gene_fam_group")
        plant_name = "Lactuca sativa"
    elif _plant_select == 'bv':
        query1 = "select * from oeur_phyt12_gene_fam_group"
        # bg_gene_count = 10726
        # bg_gene_count = 13008
        # bg_trn_count = 13008
        # bg_phytid_count = 13008
        bg_gene_count, bg_trn_count, bg_phytid_count = get_bg_gene_count("oeur_phyt12_gene_fam_group")
        plant_name = "Olea europaea"
    elif _plant_select == 'by':
        # query1 = "select * from vung_phyt12_gene_fam_group"
        # bg_gene_count = 7873
        query1 = "select * from vung_phyt12_gene_fam_group"
        # bg_gene_count = 10041
        # bg_trn_count = 14110
        # bg_phytid_count = 14110
        bg_gene_count, bg_trn_count, bg_phytid_count = get_bg_gene_count("vung_phyt12_gene_fam_group")
        plant_name = "Vigna unguiculata"

    # moncots
    elif _plant_select == 'd':
        query1 = "select * from acom_phyt12_gene_fam_group"
        # bg_gene_count = 6006
        # bg_trn_count = 6006
        # bg_phytid_count = 6006
        bg_gene_count, bg_trn_count, bg_phytid_count = get_bg_gene_count("acom_phyt12_gene_fam_group")
        plant_name = "Ananas comosus"
    elif _plant_select == 'aq':
        # query1 = "select * from macuminata_phyt12_gene_fam_group"
        # bg_gene_count = 10133
        query1 = "select * from macu_phyt12_gene_fam_group"
        # bg_gene_count = 12244
        # bg_trn_count = 12244
        # bg_phytid_count = 12244
        bg_gene_count, bg_trn_count, bg_phytid_count = get_bg_gene_count("macu_phyt12_gene_fam_group")
        plant_name = "Musa acuminata"
    # based on rice
    elif _plant_select == 'j':
        # query1 = "select * from bdist_phyt12_gene_fam_group"
        # bg_gene_count = 7184
        query1 = "select * from bdis_phyt12_gene_fam_group"
        # bg_gene_count = 9439
        # bg_trn_count = 14889
        # bg_phytid_count = 14889
        bg_gene_count, bg_trn_count, bg_phytid_count = get_bg_gene_count("bdis_phyt12_gene_fam_group")
        plant_name = "Brachypodium distachyon"
    # based on rice
    elif _plant_select == 'as':
        # query1 = "select * from rice_phyt12_gene_fam_group"
        # bg_gene_count = 8034
        query1 = "select * from osat_phyt12_gene_fam_group"
        # bg_gene_count = 10340
        # bg_trn_count = 13542
        # bg_phytid_count = 13542
        bg_gene_count, bg_trn_count, bg_phytid_count = get_bg_gene_count("osat_phyt12_gene_fam_group")
        plant_name = "Oryza sativa"
    # based on rice
    elif _plant_select == 'bh':
        # query1 = "select * from sbioc_phyt12_gene_fam_group"
        # bg_gene_count = 7775
        query1 = "select * from sbio_phyt12_gene_fam_group"
        # bg_gene_count = 9966
        # bg_trn_count = 13834
        # bg_phytid_count = 13834
        bg_gene_count, bg_trn_count, bg_phytid_count = get_bg_gene_count("sbio_phyt12_gene_fam_group")
        plant_name = "Sorghum bicolor"
        # select sum(loc_len) from graimon_phyt11_gene_fam_group;
    # select sum(loc_len) from graimon_phyt11_gene_fam_group;
    elif _plant_select == 'ag':
        # query1 = "select * from hvulgare_phyt12_gene_fam_group"
        # bg_gene_count = 9563
        query1 = "select * from hvul_phyt12_gene_fam_group"
        # bg_gene_count = 9779
        # bg_trn_count = 70267
        # bg_phytid_count = 70267
        bg_gene_count, bg_trn_count, bg_phytid_count = get_bg_gene_count("hvul_phyt12_gene_fam_group")
        plant_name = "Hordeum vulgare"
    elif _plant_select == 'au':
        # query1 = "select * from pvirgatum_phyt12_gene_fam_group"
        # bg_gene_count = 17829
        query1 = "select * from pvir_phyt12_gene_fam_group"
        # bg_gene_count = 20936
        # bg_trn_count = 29063
        # bg_phytid_count = 29063
        bg_gene_count, bg_trn_count, bg_phytid_count = get_bg_gene_count("pvir_phyt12_gene_fam_group")
        plant_name = "Panicum virgatum"
    elif _plant_select == 'bp':
        # query1 = "select * from zmays_phyt12_gene_fam_group"
        # bg_gene_count = 11200
        query1 = "select * from zmay_phyt12_gene_fam_group"
        # bg_gene_count = 13324
        # bg_trn_count = 20754
        # bg_phytid_count = 20754
        bg_gene_count, bg_trn_count, bg_phytid_count = get_bg_gene_count("zmay_phyt12_gene_fam_group")
        plant_name = "Zea mays"
    elif _plant_select == 'b':
        # query1 = "select * from atrichopoda_phyt12_gene_fam_group"
        # bg_gene_count = 5365
        query1 = "select * from atri_phyt12_gene_fam_group"
        # bg_gene_count = 6437
        # bg_trn_count = 6437
        # bg_phytid_count = 6437
        bg_gene_count, bg_trn_count, bg_phytid_count = get_bg_gene_count("atri_phyt12_gene_fam_group")
        plant_name = "Amborella trichopoda"
    elif _plant_select == 'i':
        # query1 = "select * from bstricta_phyt12_gene_fam_group"
        # bg_gene_count = 6881
        query1 = "select * from bstr_phyt12_gene_fam_group"
        # bg_gene_count = 8817
        # bg_trn_count = 9529
        # bg_phytid_count = 9529
        bg_gene_count, bg_trn_count, bg_phytid_count = get_bg_gene_count("bstr_phyt12_gene_fam_group")
        plant_name = "Boechera stricta"
    elif _plant_select == 'k':
        # query1 = "select * from bstacei_phyt12_gene_fam_group"
        # bg_gene_count = 6952
        query1 = "select * from bsta_phyt12_gene_fam_group"
        # bg_gene_count = 9006
        # bg_trn_count = 11018
        # bg_phytid_count = 11018
        bg_gene_count, bg_trn_count, bg_phytid_count = get_bg_gene_count("bsta_phyt12_gene_fam_group")
        plant_name = "Brachypodium stacei"
    elif _plant_select == 'ar':
        # query1 = "select * from othomaeum_phyt12_gene_fam_group"
        # bg_gene_count = 6126
        query1 = "select * from otho_phyt12_gene_fam_group"
        # bg_gene_count = 7445
        # bg_trn_count = 7445
        # bg_phytid_count = 7445
        bg_gene_count, bg_trn_count, bg_phytid_count = get_bg_gene_count("otho_phyt12_gene_fam_group")
        plant_name = "Oropetium thomaeum"
    elif _plant_select == 'aw':
        # query1 = "select * from ppatens_phyt12_gene_fam_group"
        # bg_gene_count = 5956
        query1 = "select * from ppat_phyt12_gene_fam_group"
        # bg_gene_count = 7034
        # bg_trn_count = 22281
        # bg_phytid_count = 22281
        bg_gene_count, bg_trn_count, bg_phytid_count = get_bg_gene_count("ppat_phyt12_gene_fam_group")
        plant_name = "Physcomitrella patens"
    elif _plant_select == 'bc':
        # query1 = "select * from smoellendorffii_phyt12_gene_fam_group"
        # bg_gene_count = 4772
        query1 = "select * from smoe_phyt12_gene_fam_group"
        # bg_gene_count = 5836
        # bg_trn_count = 5836
        # bg_phytid_count = 5836
        bg_gene_count, bg_trn_count, bg_phytid_count = get_bg_gene_count("smoe_phyt12_gene_fam_group")
        plant_name = "Selaginella moellendorffii"
    elif _plant_select == 'bd':
        # query1 = "select * from sitalica_phyt12_gene_fam_group"
        # bg_gene_count = 8209
        query1 = "select * from sita_phyt12_gene_fam_group"
        # bg_gene_count = 10539
        # bg_trn_count = 13193
        # bg_phytid_count = 13193
        bg_gene_count, bg_trn_count, bg_phytid_count = get_bg_gene_count("sita_phyt12_gene_fam_group")
        plant_name = "Setaria italica"
    elif _plant_select == 'be':
        # query1 = "select * from sviridis_phyt12_gene_fam_group"
        # bg_gene_count = 7996
        query1 = "select * from svir_phyt12_gene_fam_group"
        # bg_gene_count = 9922
        # bg_trn_count = 13868
        # bg_phytid_count = 13868
        bg_gene_count, bg_trn_count, bg_phytid_count = get_bg_gene_count("svir_phyt12_gene_fam_group")
        plant_name = "Setaria viridis"
    elif _plant_select == 'bj':
        # query1 = "select * from spolyrhiza_phyt12_gene_fam_group"
        # bg_gene_count = 4861
        query1 = "select * from spol_phyt12_gene_fam_group"
        # bg_gene_count = 5995
        # bg_trn_count = 5995
        # bg_phytid_count = 5995
        bg_gene_count, bg_trn_count, bg_phytid_count = get_bg_gene_count("spol_phyt12_gene_fam_group")
        plant_name = "Spirodela polyrhiza"
    elif _plant_select == 'bm':
        # query1 = "select * from taestivum_phyt12_gene_fam_group"
        # bg_gene_count = 25487
        query1 = "select * from taes_phyt12_gene_fam_group"
        # bg_gene_count = 26390
        # bg_trn_count = 74262
        # bg_phytid_count = 74262
        bg_gene_count, bg_trn_count, bg_phytid_count = get_bg_gene_count("taes_phyt12_gene_fam_group")
        plant_name = "Triticum aestivum"
    # based on rice
    elif _plant_select == 'br':
        query1 = "select * from bsyl_phyt12_gene_fam_group"
        # bg_gene_count = 7516
        # bg_gene_count = 9932
        # bg_trn_count = 13883
        # bg_phytid_count = 13883
        bg_gene_count, bg_trn_count, bg_phytid_count = get_bg_gene_count("bsyl_phyt12_gene_fam_group")
        plant_name = "Brachypodium sylvaticum"
    elif _plant_select == 'bw':
        query1 = "select * from phal_phyt12_gene_fam_group"
        # bg_gene_count = 7414
        # bg_gene_count = 9506
        # bg_trn_count = 12832
        # bg_phytid_count = 12832
        bg_gene_count, bg_trn_count, bg_phytid_count = get_bg_gene_count("phal_phyt12_gene_fam_group")
        plant_name = "Panicum hallii"
    elif _plant_select == 'bz':
        # query1 = "select * from zmar_phyt12_gene_fam_group"
        # bg_gene_count = 5060
        query1 = "select * from zmar_phyt12_gene_fam_group"
        # bg_gene_count = 6359
        # bg_trn_count = 6359
        # bg_phytid_count = 6359
        bg_gene_count, bg_trn_count, bg_phytid_count = get_bg_gene_count("zmar_phyt12_gene_fam_group")
        plant_name = "Zostera marina"
    elif _plant_select == 'v':
        # query1 = "select * from zmar_phyt12_gene_fam_group"
        # bg_gene_count = 5060
        query1 = "select * from csat_phyt12_gene_fam_group"
        # bg_gene_count = 10062
        # bg_trn_count = 10062
        # bg_phytid_count = 0
        bg_gene_count, bg_trn_count, bg_phytid_count = get_bg_gene_count("csat_phyt12_gene_fam_group")
        plant_name = "Cucumis sativus"

    # others
    elif _plant_select == 's':
        # query1 = "select * from czofingiensis_phyt12_gene_fam_group"
        # bg_gene_count = 2068
        query1 = "select * from czof_phyt12_gene_fam_group"
        # bg_gene_count = 2774
        # bg_trn_count = 2774
        # bg_phytid_count = 2774
        bg_gene_count, bg_trn_count, bg_phytid_count = get_bg_gene_count("czof_phyt12_gene_fam_group")
        plant_name = "Chromochloris zofingiensis"
    elif _plant_select == 'x':
        # query1 = "select * from dsalina_phyt12_gene_fam_group"
        # bg_gene_count = 1613
        query1 = "select * from dsal_phyt12_gene_fam_group"
        # bg_gene_count = 2118
        # bg_trn_count = 2443
        # bg_phytid_count = 2443
        bg_gene_count, bg_trn_count, bg_phytid_count = get_bg_gene_count("dsal_phyt12_gene_fam_group")
        plant_name = "Dunaliella salina"
    elif _plant_select == 'am':
        # query1 = "select * from mpolymorpha_phyt12_gene_fam_group"
        # bg_gene_count = 3410
        query1 = "select * from mpol_phyt12_gene_fam_group"
        # bg_gene_count = 4377
        # bg_trn_count = 5732
        # bg_phytid_count = 5732
        bg_gene_count, bg_trn_count, bg_phytid_count = get_bg_gene_count("mpol_phyt12_gene_fam_group")
        plant_name = "Marchantia polymorpha"
    elif _plant_select == 'ao':
        # query1 = "select * from mpusilla_phyt12_gene_fam_group"
        # bg_gene_count = 1307
        query1 = "select * from mpus_phyt12_gene_fam_group"
        # bg_gene_count = 1797
        # bg_trn_count = 1797
        # bg_phytid_count = 1797
        bg_gene_count, bg_trn_count, bg_phytid_count = get_bg_gene_count("mpus_phyt12_gene_fam_group")
        plant_name = "Micromonas pusilla"
    elif _plant_select == 'at':
        # query1 = "select * from olucimarinus_phyt12_gene_fam_group"
        # bg_gene_count = 1127
        query1 = "select * from oluc_phyt12_gene_fam_group"
        # bg_gene_count = 1518
        # bg_trn_count = 1518
        # bg_phytid_count = 1518
        bg_gene_count, bg_trn_count, bg_phytid_count = get_bg_gene_count("oluc_phyt12_gene_fam_group")
        plant_name = "Ostreococcus lucimarinus"
    elif _plant_select == 'bx':
        query1 = "select * from pumb_phyt12_gene_fam_group"
        # bg_gene_count = 1104
        # bg_gene_count = 1462
        # bg_trn_count = 1500
        # bg_phytid_count = 1500
        bg_gene_count, bg_trn_count, bg_phytid_count = get_bg_gene_count("pumb_phyt12_gene_fam_group")
        plant_name = "Porphyra umbilicalis"
    # moss
    elif _plant_select == 'bi':
        # query1 = "select * from sfallax_phyt12_gene_fam_group"
        # bg_gene_count = 5694
        query1 = "select * from sfal_phyt12_gene_fam_group"
        # bg_gene_count = 7298
        # bg_trn_count = 8884
        # bg_phytid_count = 8884
        bg_gene_count, bg_trn_count, bg_phytid_count = get_bg_gene_count("sfal_phyt12_gene_fam_group")
        plant_name = "Sphagnum fallax"
    elif _plant_select == 'bo':
        # query1 = "select * from vcarteri_phyt12_gene_fam_group"
        # bg_gene_count = 1783
        query1 = "select * from vcar_phyt12_gene_fam_group"
        # bg_gene_count = 2369
        # bg_trn_count = 2717
        # bg_phytid_count = 2717
        bg_gene_count, bg_trn_count, bg_phytid_count = get_bg_gene_count("vcar_phyt12_gene_fam_group")
        plant_name = "Volvox carteri"
    elif _plant_select == 'bt':
        query1 = "select * from csub_phyt12_gene_fam_group"
        # bg_gene_count = 1567
        # bg_gene_count = 2054
        # bg_trn_count = 2054
        # bg_phytid_count = 2054
        bg_gene_count, bg_trn_count, bg_phytid_count = get_bg_gene_count("csub_phyt12_gene_fam_group")
        plant_name = "Coccomyxa subellipsoidea"
    elif _plant_select == 'ca':
        query1 = "select * from aoff_phyt12_gene_fam_group"
        # bg_gene_count = 7348
        # bg_trn_count = 7348
        # bg_phytid_count = 7348
        bg_gene_count, bg_trn_count, bg_phytid_count = get_bg_gene_count("aoff_phyt12_gene_fam_group")
        plant_name = "Asparagus officinalis"
    elif _plant_select == 'cb':
        query1 = "select * from bhyb_phyt12_gene_fam_group"
        # bg_gene_count = 16899
        # bg_trn_count = 19932
        # bg_phytid_count = 19932
        bg_gene_count, bg_trn_count, bg_phytid_count = get_bg_gene_count("bhyb_phyt12_gene_fam_group")
        plant_name = "Brachypodium hybridum"
    elif _plant_select == 'cc':
        query1 = "select * from crei_phyt12_gene_fam_group"
        # bg_gene_count = 2882
        # bg_trn_count = 3201
        # bg_phytid_count = 3201
        bg_gene_count, bg_trn_count, bg_phytid_count = get_bg_gene_count("crei_phyt12_gene_fam_group")
        plant_name = "Chlamydomonas reinhardtii"
    elif _plant_select == 'cd':
        query1 = "select * from cari_phyt12_gene_fam_group"
        # bg_gene_count = 8434
        # bg_trn_count = 8434
        # bg_phytid_count = 8434
        bg_gene_count, bg_trn_count, bg_phytid_count = get_bg_gene_count("cari_phyt12_gene_fam_group")
        plant_name = "Cicer arietinum"
    elif _plant_select == 'ce':
        query1 = "select * from mspr_phyt12_gene_fam_group"
        # bg_gene_count = 1837
        # bg_trn_count = 1837
        # bg_phytid_count = 1837
        bg_gene_count, bg_trn_count, bg_phytid_count = get_bg_gene_count("mspr_phyt12_gene_fam_group")
        plant_name = "Zostera marina"
    elif _plant_select == 'cf':
        query1 = "select * from msin_phyt12_gene_fam_group"
        # bg_gene_count = 18598
        # bg_trn_count = 24665
        # bg_phytid_count = 24665
        bg_gene_count, bg_trn_count, bg_phytid_count = get_bg_gene_count("msin_phyt12_gene_fam_group")
        plant_name = "Miscanthus sinensis"
        # select sum(loc_len), sum(trn_len), sum(phyt_id_len) from msin_phyt12_gene_fam_group;

    # elif _plant_select == 'x':
    #    query1 = "select * from alyrata_phyt12_gene_fam_group"
    #    bg_gene_count = 7382
    # elif _plant_select == 'ad':
    #    query1 = "select * from brapa_phyt12_gene_fam_group"
    #    bg_gene_count = 11020
    # elif _plant_select == 'f':
    #    query1 = "select * from tomato_phyt12_gene_fam_group"
    #    bg_gene_count = 8112

    cur.execute(query1)
    fam_rec = cur.fetchall()
    # gene family and their corresponding ids
    gene_fam_dict = dict()
    short_fam_dict = dict()
    gene_fam_count_dict = dict()
    get_user_id_count_for_gene_fam = dict()
    uniq_id_count_dict = dict()
    # get gene ids from user list for each gene family
    get_gene_ids_from_user_dict = dict()
    # go terms
    uniq_p = dict()
    uniq_f = dict()
    uniq_c = dict()
    # get the annotation count. number of genes from user input present in genfam database
    anot_count = 0

    flag_dict = dict()
    for item in fam_rec:
        # # phytozome locus
        if _id_type == "a":
            # item[0] is a gene family
            # item[1] is locus array_agg
            gene_fam_dict[item[0]] = [x.upper() for x in item[1]]
            short_fam_dict[item[0]] = item[7]
            get_gene_ids_from_user_dict[item[0]] = []
            # item[2] loc count
            # gene fam count for background genome
            gene_fam_count_dict[item[0]] = item[2]
            # gene fam count for user ids
            get_user_id_count_for_gene_fam[item[0]] = 0
            flag_dict[item[0]] = 0
            uniq_p[item[0]] = item[8]
            uniq_f[item[0]] = item[9]
            uniq_c[item[0]] = item[10]
            bg_gene_count = bg_gene_count
        # phytozome transcript
        elif _id_type == "b":
            # item[3] transcript id
            gene_fam_dict[item[0]] = [x.upper() for x in item[3]]
            short_fam_dict[item[0]] = item[7]
            get_gene_ids_from_user_dict[item[0]] = []
            # item[4] transcript count
            gene_fam_count_dict[item[0]] = item[4]
            get_user_id_count_for_gene_fam[item[0]] = 0
            flag_dict[item[0]] = 0
            uniq_p[item[0]] = item[8]
            uniq_f[item[0]] = item[9]
            uniq_c[item[0]] = item[10]
            # assign background count based on trn count
            bg_gene_count = bg_trn_count
        # phytozome pacId
        elif _id_type == "c":
            # item[5] PAC id
            gene_fam_dict[item[0]] = [x.upper() for x in item[5]]
            short_fam_dict[item[0]] = item[7]
            get_gene_ids_from_user_dict[item[0]] = []
            # item[6] pac count
            gene_fam_count_dict[item[0]] = item[6]
            get_user_id_count_for_gene_fam[item[0]] = 0
            flag_dict[item[0]] = 0
            uniq_p[item[0]] = item[8]
            uniq_f[item[0]] = item[9]
            uniq_c[item[0]] = item[10]
            bg_gene_count = bg_phytid_count

    # get unique id count for given file
    # id_file input user file obtained from form
    for _id in id_file:
        _id = _id.strip().upper()
        # remove the duplicate ids and keep unique
        uniq_id_count_dict[_id] = 0

    for key1 in gene_fam_dict.keys():
        for key2 in uniq_id_count_dict.keys():
            if key2.strip().upper() in gene_fam_dict[key1]:
                # if the user input id present in gene_fam_dict increment count
                get_gene_ids_from_user_dict[key1].append(key2.strip().upper())
                get_user_id_count_for_gene_fam[key1] += 1
                anot_count += 1

    enrichment_result, fdr, mapped_query_ids = enrichment_analysis(len(uniq_id_count_dict), get_user_id_count_for_gene_fam,
                                                 gene_fam_count_dict, bg_gene_count, gene_fam_dict,
                                                 _stat_sign_test, _multi_test_corr, uniq_p, uniq_f, uniq_c,
                                                 get_gene_ids_from_user_dict, short_fam_dict)


    # replace all fdr values which are greater than 1 to 1
    fdr[fdr > 1] = 1
    _fam_out_enrich_file = os.path.join(MEDIA_ROOT, "output_files", user, "fam_enrich_out.txt")
    _fam_out_all_file = os.path.join(MEDIA_ROOT, "output_files", user, "fam_all_out.txt")
    _fam_enrich_out = open(_fam_out_enrich_file, "w")
    _fam_all_out = open(_fam_out_all_file, "w")
    # Query IDs = total query ids from user (k)
    # Annotated query IDs = query ids annotated to particular gene family
    # background ids = particular gene family ids present in whole genome backround
    # background total = total genome ids

    for x in range(0, len(fdr)):
        enrichment_result[x].insert(3, anot_count)
        enrichment_result[x].insert(9, fdr[x])

    if _stat_sign_test == "d":
        # chi-s test only
        _fam_enrich_out.write("Gene Family"+"\t"+"Short Name"+"\t"+"Query total"+"\t"+"Annotated query total"+"\t"+
                          "Annotated query per family"+"\t"+"Background annotated total"+"\t"+
                          "Background annotated per family"+"\t"+"Odds ratio"+"\t"+"P-value"+"\t"+"FDR"+"\t"+
                          "GO biological process"+"\t"+"GO molecular function"+"\t"+"GO cellular component"+"\t"+
                          "Gene IDs"+"\t"+"Cells with expected frequency <5 (%)"+"\n")
        _fam_all_out.write("Gene Family"+"\t"+"Short Name"+"\t"+"Query total"+"\t"+"Annotated query total"+"\t"+"Annotated query per family"
                       +"\t"+"Background annotated total"+"\t"+"Background annotated per family"+"\t"+"Odds ratio"+
                       "\t"+"P-value"+"\t"+"FDR"+"\t"+"GO biological process"+"\t"+"GO molecular function"+"\t"+
                       "GO cellular component"+"\t"+"Gene IDs"+"\t"+"Cells with expected frequency <5 (%)"+"\n")
    else:
        _fam_enrich_out.write("Gene Family"+"\t"+"Short Name"+"\t"+"Query total"+"\t"+"Annotated query total"+"\t"+
                            "Annotated query per family"+"\t"+"Background annotated total"+"\t"+
                            "Background annotated per family" +"\t"+"Odds ratio"+"\t"+"P-value"+"\t"+"FDR"+"\t"+
                            "GO biological process" +"\t"+"GO molecular function"+"\t"+"GO cellular component"+"\t"+
                            "Gene IDs"+"\n")
        _fam_all_out.write("Gene Family"+"\t"+"Short Name"+"\t"+"Query total"+ "\t"+"Annotated query total"+"\t"+"Annotated query per family"
            +"\t"+"Background annotated total"+"\t"+"Background annotated per family"+"\t"+"Odds ratio"+
            "\t"+"P-value"+"\t"+"FDR"+"\t"+"GO biological process"+"\t"+"GO molecular function"+"\t" +
            "GO cellular component"+"\t"+"Gene IDs"+"\n")

    for x in range(0, len(enrichment_result)):
        _fam_all_out.write('\t'.join(str(v) for v in enrichment_result[x]) + "\n")
        # check pvalue less than 0.05
        if float(enrichment_result[x][8]) <= 0.05:
            _fam_enrich_out.write('\t'.join(str(v) for v in enrichment_result[x])  + "\n")
    _fam_all_out.close()
    _fam_enrich_out.close()

    # if the number input ids are less than 5, the process will stop
    # check if list is empty; this is in case of wrong input by user for gene ids or not as per phytozome format or
    # wrong species selected
    # if len(uniq_id_count_dict) >= 5 and enrichment_result and _plant_select!="z":
    if mapped_query_ids >= 5 and enrichment_result and _plant_select != "z":
        return "okay"
    # check if list is empty; this is in case of wrong input by user for gene ids or not as per phytozome format or
    # wrong species selected




