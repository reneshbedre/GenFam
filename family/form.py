from django import forms
from django.utils.safestring import mark_safe


__author__ = 'renesh'


class SeqForm(forms.Form):
    # no textfield for forms
    gene_ids = forms.CharField(widget=forms.Textarea(attrs={'placeholder': 'Paste Sequence IDs...'}),
                               # label='Sequence IDs:',
                               label=mark_safe(
                                   "<b>Sequence IDs</b> "
                                   "<abbr title='Each line should contain one gene ID as per "
                                   "Phytozome v12.0 database format'><span "
                                   "class='glyphicon "
                                   "glyphicon-question-sign' aria-hidden='true'></span></abbr>"),
                               # help_text=mark_safe("<abbr title='Each line should contain one gene and provide the "
                               #                    "gene names as per Phytozome gene ID format'><span class='glyphicon "
                               #                    "glyphicon-question-sign' aria-hidden='true'></span></abbr>"),
                               required=True)
    # CATEGORY_CHOICES = [
    #    ('a', u'Reverse Complement'),
    #    ('b', u'Sequence Length'),
    #]
    # seq_radio = forms.ChoiceField(choices=CATEGORY_CHOICES, widget=forms.RadioSelect(), label='', required=False)
    # seq_file = forms.FileField(label='Upload file', required=False)
    PLANT_CHOICES = [
        ('z', '-- Select --'),
        ('a', u'Amaranthus hypochondriacus'), #1
        ('b', u'Amborella trichopoda'), #2
        #('c', u'Anacardium occidentale'),
        ('d', u'Ananas comosus'), #3
        #('h', u'Aquilegia coerulea'),
        #('f', u'Arabidopsis halleri'),
        ('g', u'Arabidopsis lyrata'), #4
        ('e', u'Arabidopsis thaliana'), #5
        ('ca', u'Asparagus officinalis'), #6
        ('i', u'Boechera stricta'), #7
        ('j', u'Brachypodium distachyon'), #8
        #('cb', u'Brachypodium hybridum'),
        #('k', u'Brachypodium stacei'),
        #('br', u'Brachypodium sylvaticum'),
        ('l', u'Brassica oleracea'),  #9
        #('m', u'Brassica rapa'),
        ('n', u'Capsella grandiflora'), #10
        ('o', u'Capsella rubella'), #11
        ('p', u'Carica papaya'), #12
        ('q', u'Chenopodium quinoa'), #13
        ('cc', u'Chlamydomonas reinhardtii'), #14
        ('s', u'Chromochloris zofingiensis'), #15
        ('cd', u'Cicer arietinum'), #16
        ('t', u'Citrus clementina'), #17
        ('u', u'Citrus sinensis'), #18
        ('bt', u'Coccomyxa subellipsoidea'), #19
        ('v', u'Cucumis sativus'), #20
        ('w', u'Daucus carota'), #21
        ('x', u'Dunaliella salina'), #22
        ('y', u'Eucalyptus grandis'), #23
        ('ab', u'Eutrema salsugineum'), #24
        ('ac', u'Fragaria vesca'), #25
        ('af', u'Glycine max Wm82.a2.v1'), #26
        #('ad', u'Gossypium hirsutum'),
        ('ae', u'Gossypium raimondii'), #27
        ('ag', u'Hordeum vulgare'), #28
        ('ah', u'Kalanchoe fedtschenkoi'), #29
        #('ai', u'Kalanchoe laxiflora'),
        ('bu', u'Lactuca sativa'), #30
        ('aj', u'Linum usitatissimum'), #31
        ('ak', u'Malus domestica'), #32
        ('al', u'Manihot esculenta'), #33
        ('am', u'Marchantia polymorpha'), #34
        ('an', u'Medicago truncatula'), #35
        ('ao', u'Micromonas pusilla'), #36
        ('ap', u'Mimulus guttatus'),  #37
        #('cf', u'Miscanthus sinensis'),
        ('aq', u'Musa acuminata'), #38
        ('bv', u'Olea europaea'),  #39
        ('ar', u'Oropetium thomaeum'),  #40
        ('as', u'Oryza sativa'),  #41
        ('at', u'Ostreococcus lucimarinus'), #42
        #('bw', u'Panicum hallii'),
       #('au', u'Panicum virgatum'),
        #('av', u'Phaseolus vulgaris'),
        ('aw', u'Physcomitrella patens'), #43
        #('ax', u'Populus deltoides'),
        ('ay', u'Populus trichocarpa'), #44
        ('bx', u'Porphyra umbilicalis'), #45
        ('az', u'Prunus persica'), #46
        ('ba', u'Ricinus communis'), #47
        #('bb', u'Salix purpurea'),
        ('bc', u'Selaginella moellendorffii'), #48
        ('bd', u'Setaria italica'), #49
        #('be', u'Setaria viridis v1.1'),
        ('bg', u'Solanum lycopersicum (ITAG2.4)'), #50
        ('bf', u'Solanum tuberosum'), #51
        ('bh', u'Sorghum bicolor'), #52
        #('bi', u'Sphagnum fallax'),
        ('bj', u'Spirodela polyrhiza'), #53
        ('bk', u'Theobroma cacao'), #54
        ('bl', u'Trifolium pratense'), #55
        ('bm', u'Triticum aestivum'), #56
        #('by', u'Vigna unguiculata'),
        ('bn', u'Vitis vinifera'), #57
        ('bo', u'Volvox carteri'), #58
        ('bp', u'Zea mays'), #59
        ('bz', u'Zostera marina'), #60
    ]
    plant_select = forms.ChoiceField(choices=PLANT_CHOICES,
                                     # label='Select Plant Species:',
                                     label=mark_safe(
                                         "<b>Select Plant Species</b> "
                                         "<abbr title='Select plant species corresponding to the gene IDs'><span class='glyphicon "
                                         "glyphicon-question-sign' aria-hidden='true'></span></abbr>"),
                                     # help_text=mark_safe(
                                     #    "<abbr title='Select the plant species for which you are studying "
                                     #    "differentially expressed genes'><span class='glyphicon "
                                     #    "glyphicon-question-sign' aria-hidden='true'></span></abbr>"),
                                     required=True)
    ID_CHOICES = [
        ('z', '-- Select --'),
        ('a', u'Phytozome locus'),
        ('b', u'Phytozome transcript'),
        ('c', u'Phytozome pacId'),
    ]
    ID_select = forms.ChoiceField(widget=forms.Select(attrs={'id': 'choic', 'class': 'reg'}),
                                  choices=ID_CHOICES,
                                  # label='Select ID Type:',
                                  label=mark_safe(
                                      "<b>Select ID Type</b> "
                                      "<abbr title='Select gene ID type for respective plant species."
                                      "&#13;"
                                      " For example allowed IDs for Rice (Oryza sativa) are, &#13;"
                                      "Phytozome locus: "
                                      "LOC_Os01g06882 &#13;"
                                      " Phytozome transcript: LOC_Os01g06882.1 &#13;"
                                      "Phytozome "
                                      " pacID: 24120792'><span class='glyphicon "
                                      "glyphicon-question-sign' aria-hidden='true'></span></abbr>"),
                                  required=True)
    STAT_TEST = [
        ('z', '-- Select --'),
        ('a', u'Fisher\'s exact'),
        ('b', u'Hypergeometric'),
        ('c', u'Binomial'),
        ('d', u'Chi-Squared'),
        # ('e', u'Z'),
        # ('f', u'Kolmogorov-Smirnov'),
        # ('g', u'Permutation'),
    ]
    stat_sign_test = forms.ChoiceField(choices=STAT_TEST,
                                       # label='Select Statistical Significance Test:',
                                       label=mark_safe(
                                           "<b>Select Statistical Significance Test</b> "
                                           "<abbr title='We recommend using Fisher exact test, chi-square test and "
                                           "hypergeometric distribution for smaller datasets (< 1000 gene IDs), and binomial "
                                           "distribution for larger datasets'><span class='glyphicon "
                                           "glyphicon-question-sign' aria-hidden='true'></span></abbr>"),
                                       required=True,
                                       initial='a')
    MULTI_TEST = [
        ('z', '-- Select --'),
        ('a', u'Bonferroni'),
        ('b', u'Bonferroni-Holm'),
        ('c', u'Benjamini-Hochberg'),
    ]
    multi_test_corr = forms.ChoiceField(choices=MULTI_TEST,
                                        # label='Select Multiple Testing Corrections:',
                                        label=mark_safe(
                                            "<b>Multiple testing correction method</b> "
                                            "<abbr title='We recommend multiple testing corrections to reduce false "
                                            "positives'><span class='glyphicon "
                                            "glyphicon-question-sign' aria-hidden='true'></span></abbr>"),
                                        required=True,
                                        initial='c')
