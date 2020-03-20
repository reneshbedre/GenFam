#!/usr/bin/python

# from anot.settings import BASE_DIR
import psycopg2
import scipy.stats as stats
from numpy import array, empty
import os

#BASE_DIR = imp.load_source('modulename',
#                                  '/Users/renesh/Renesh_Docs/Research/django/anot/anot/settings.py')


def correct_pvalues_for_multiple_testing(_pvalues, correction_type):
    pvalues = array(_pvalues)
    n = pvalues.shape[0]
    new_pvalues = empty(n)
    if correction_type == "Bonferroni":
        new_pvalues = n * pvalues
    elif correction_type == "Bonferroni-Holm":
        values = [ (pvalue, i) for i, pvalue in enumerate(pvalues) ]
        # sort by pvlaue
        values.sort()
        for rank, vals in enumerate(values):
            pvalue, i = vals
            new_pvalues[i] = (n-rank) * pvalue
    elif correction_type == "Benjamini-Hochberg":
        values = [ (pvalue, i) for i, pvalue in enumerate(pvalues) ]
        values.sort()
        values.reverse()
        new_values = []
        for i, vals in enumerate(values):
            rank = n - i
            pvalue, index = vals
            new_values.append((n/rank) * pvalue)
        for i in range(0, int(n)-1):
            if new_values[i] < new_values[i+1]:
                new_values[i+1] = new_values[i]
        for i, vals in enumerate(values):
            pvalue, index = vals
            new_pvalues[index] = new_values[i]
    return new_pvalues


def enrichment_analysis(_get_id_count, _get_go_id_count, _go_to_gene_count_dict, _bg_gene_count, _go_anot_dict,
                        _stat_sign_test, _multi_test_corr):
    # read https://cgrlucb.wikispaces.com/Functional+Enrichment+Analysis
    pvalues = []
    enrichment_result = []
    for key in _get_go_id_count.keys():
        # first row and first column (k)
        gene_in_category = _get_go_id_count[key]
        # first row and second column (m-k)  _get_id_count (m)
        gene_not_in_category_in_sample = _get_id_count-gene_in_category
        # second row and first column (n-k)
        gene_not_in_catgory_and_in_genome = len(_go_to_gene_count_dict[key])-gene_in_category
        # second row and second column gene_in_category+gene_not_in_catgory_and_in_genome (n)
        # N-m-n+k
        bg_in_genome = \
            _bg_gene_count-_get_id_count-(gene_in_category+gene_not_in_catgory_and_in_genome)+gene_in_category

        if _stat_sign_test == "a":
            # print key, gene_in_category, gene_not_in_category_in_sample, gene_not_in_catgory_and_in_genome, bg_in_genome
            # fisher exact test
            oddsratio, pvalue = stats.fisher_exact([[gene_in_category, gene_not_in_category_in_sample],
                                                [gene_not_in_catgory_and_in_genome, bg_in_genome]])

            enrichment_result.append([key, '\t'.join(_go_anot_dict[key]), gene_in_category,
                                      gene_not_in_category_in_sample,gene_not_in_catgory_and_in_genome, bg_in_genome,
                                      oddsratio, pvalue])
            pvalues.append(pvalue)

    # FDR Bonferroni
    if _multi_test_corr == "a":
        fdr = correct_pvalues_for_multiple_testing(pvalues, "Bonferroni")
    return enrichment_result, fdr


def go_enrich(_gene_ids_file, _plant_select, _stat_sign_test, _multi_test_corr):
    id_file = open(_gene_ids_file, "rU")
    go_dict = dict()
    go_anot_dict = dict()
    gene_go_dict = dict()
    get_go_id_count = dict()
    go_to_gene_count_dict = dict()

    # create the folder if it doesn't exist.
    try:
        os.mkdir(os.path.join(BASE_DIR, "output_files"))
    except:
        pass

    # connect to psql db
    try:
        con = psycopg2.connect(database='plant_db', user='renesh', password='927122', host="127.0.0.1", port="5432")
    except ValueError:
        print("Not able to connect database")

    cur = con.cursor()
    query1 = ""
    query2 = ""

    if _plant_select == 'a':
        query1 = "select * from rice_goslim"
        query2 = "select count(distinct gene)from rice_goslim"
    elif _plant_select == 'b':
        query1 = "select * from ath_goslim"
        query2 = "select count(distinct gene)from ath_goslim"

    cur.execute(query1)
    go_rec = cur.fetchall()
    cur.execute(query2)
    bg_gene_count = cur.fetchone()

    # get the count of each goslim term in particular plant species given by -d option
    for item in go_rec:
        go_id = item[2]
        go_class = item[3]
        go_anot = item[4]
        if go_id not in go_dict:
            go_dict[go_id] = 1
            go_anot_dict[go_id] = [go_class, go_anot]
        else:
            go_dict[go_id] += 1

    # assign GO terms to each gene id in particular plant species given by -d option
    for item in go_rec:
        go_id = item[2]
        gene_id = item[0]
        if gene_id not in gene_go_dict:
                gene_go_dict[gene_id] =[go_id]
        else:
            if go_id not in gene_go_dict[gene_id]:
                gene_go_dict[gene_id].append(go_id)

    # get the gene count for each goslim id
    for item in go_rec:
        go_id = item[2]
        gene_id = item[0]
        if go_id not in go_to_gene_count_dict:
            go_to_gene_count_dict[go_id] = [gene_id]
        else:
            if gene_id not in go_to_gene_count_dict[go_id]:
                go_to_gene_count_dict[go_id].append(gene_id)

    # assign the count of each go term to each gene/transcript in a given sample given by -i option (DEGs)
    get_id_count = 0
    for item in id_file:
        item = item.strip()
        # for total gene id count in input
        get_id_count += 1
        if item in gene_go_dict:
            # if gene is found increment the go id count for each gene
            for go_id in gene_go_dict[item]:
                if go_id not in get_go_id_count:
                    get_go_id_count[go_id] = 1
                else:
                    get_go_id_count[go_id] += 1

    enrichment_result, fdr = enrichment_analysis(get_id_count, get_go_id_count, go_to_gene_count_dict, bg_gene_count[0],
                                             go_anot_dict, _stat_sign_test, _multi_test_corr)


    # replace all fdr values which are greater than 1 to 1
    fdr[fdr > 1] = 1

    # print len(go_dict), sum(go_dict.values()), len(gene_go_dict), len(get_go_id_count), len(go_to_gene_count_dict)
    # print enrichment_result, fdr
    _go_out_file = os.path.join(BASE_DIR, "output_files", "go_enrich_out.txt")
    _go_out = open(_go_out_file, "w")
    for x in range(0, len(enrichment_result)):
        _go_out.write('\t'.join(str(v) for v in enrichment_result[x])+"\t"+str(fdr[x])+"\n")
    _go_out.close()
