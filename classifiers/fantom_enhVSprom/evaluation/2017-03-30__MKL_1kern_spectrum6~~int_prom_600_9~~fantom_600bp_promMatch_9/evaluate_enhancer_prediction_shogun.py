#!/usr/bin/env python
"""
 evaluate_enhancers_shogun.py - Copyright Tony Capra 2012

 Change log: 
  - 03/11/12 started 

 -----------------------------------------------------------------------------

"""

usage = """
USAGE:

    - description coming soon.


OPTIONS:
    -d [comma separated string]
     name. blah blah blah

    -h
     help. Print this message.

"""


import os
import sys
import getopt
import random
import datetime
import commands


#import croc

import numpy as np
#np.set_printoptions(threshold='nan')

import matplotlib
matplotlib.use("Agg")
import matplotlib.pylab as plt

from scipy.interpolate import interp1d
try:
    import scipy.stats as stats
except ValueError:
    import scipy.stats as stats

from shogun.Features import *
from shogun.Classifier import *
from shogun.Kernel import *
from shogun.Evaluation import *
from shogun.ModelSelection import *

from utils_enhancer_prediction_cfeat import *


################################################################################

def main():

    # parse command line options and update variables
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hac:f:n:o:p:s:t:")
    except getopt.GetoptError:
        sys.exit(usage)

    if len(args) == 3:
        CLASSIFIER_CONFIG = args[0].split(',')
        POS_BEDS = args[1].split(',')
        NEG_BEDS = args[2].split(',')
    else:
        sys.exit(usage)
    

    INTERSECT_STAT = 'binary' # -t  binary or coverage

    now_str = datetime.datetime.now().strftime('%Y-%m-%d_%H.%M.%S')
    OUT_NAME = '%s__%s~~%s~~%s' % (now_str[:10], 
                                   CLASSIFIER_CONFIG[0].split('/')[-1].replace('.classifier_config', ''),
                                   POS_BEDS[0].split('/')[-1].replace('.bed', ''), 
                                   NEG_BEDS[0].split('/')[-1].replace('.bed', ''))  # -n
    OUT_DIR = None  # -o

    AVERAGE_KERNELS = False  # -a

    CV_ROUNDS = 10  # -c

    for opt, arg in opts:
        if opt == "-h":
            sys.exit('Use me properly.')

        elif opt == '-a':
            AVERAGE_KERNELS = True

        elif opt == '-c':
            CV_ROUNDS = int(arg)

        elif opt == '-n':
            OUT_NAME = arg

        elif opt == '-o':
            OUT_DIR = arg

        elif opt == '-t':
            INTERSECT_STAT = arg


    if OUT_NAME[:10] != now_str[:10]:
        OUT_NAME = now_str[:10] + "_" + OUT_NAME

    if OUT_DIR == None:
        OUT_DIR = 'results/svm_evaluation/%s/' % OUT_NAME.replace(' ', '_')

    if not os.path.exists(OUT_DIR): os.makedirs(OUT_DIR)
    os.system('cp %s %s' % (os.path.abspath(__file__), OUT_DIR))
    os.system('cp bin/utils_enhancer_prediction_cfeat.py %s' % (OUT_DIR))
    for config_file in CLASSIFIER_CONFIG:
        os.system('cp %s %s' % (config_file, OUT_DIR))

    PREAMBLE = '# %s %s\n#\n' % (' '.join(sys.argv), datetime.datetime.now())
    PREAMBLE += '#\n'
    PREAMBLE += '#'

    print PREAMBLE
    # Finished processing command line options

    #####################
    # PROCESS INPUT
    #####################

    print '\nProcessing input...'
    regions, raw_labels = make_shogun_labels(POS_BEDS, NEG_BEDS)
    print '\tNum Pos: %d\n\tNum Neg: %d\n' % (raw_labels.count(1), 
                                              raw_labels.count(-1))

    # combine +/- examples into a single bed file for convenience
    rand_num = np.random.randint(1000000)
    ALL_BED = 'all_regions_TEMP_%d.bed' % rand_num
    of = open(ALL_BED, 'w')
    for region in regions: 
        of.write('%s\t%d\t%d\n' % (region[0], region[1], region[2]))
    of.close()

    # This could be disasterous if the order in ALL_BED doesn't match
    # that in the regions variable. So we just make the combined bed
    # from the regions array to ensure the same order.
    #os.system('cat %s %s | cut -f 1-3 >%d.bed' % (' '.join(POS_BEDS), ' '.join(NEG_BEDS), rand_num))
    #os.system('sortBed -i %d.bed > %s' % (rand_num, ALL_BED))
    #os.remove('%d.bed' % rand_num)

    # READ CONFIG FILE
    print 'Initializing classifiers:'
    classifiers = parse_config_file(CLASSIFIER_CONFIG)

    # BUILD FEATURES
    print "\nLoading features..."
    name2feats, name2ppfeats, names2seqs = build_features(classifiers, ALL_BED, regions=regions)
    os.remove(ALL_BED)

    # COMPUTE KERNEL MATRICES
    print "\nComputing kernel matrices..."
    kernel_matrix_sets = []
    for classifier in classifiers:
        kernel_mats = []
        for kidx, kernel in enumerate(classifier.kernels):
            if classifier.seqs[kidx]: 
                #raw_features = prepare_sequence_data(classifier.seqs[kidx], regions)
                raw_features = names2seqs[str(classifier.seqs[kidx])]
            else:
                raw_features = []
                if classifier.bin_features[kidx]:
                    raw_features = [name2feats[x] for x in classifier.bin_features[kidx]]
                if classifier.con_features[kidx]:
                    raw_features += [name2feats[x] for x in classifier.con_features[kidx]]
                if classifier.pp_features[kidx]:
                    raw_features += [name2ppfeats[x] for x in classifier.pp_features[kidx]]

            km = compute_kernel_matrix(kernel, raw_features, 
                                       norm=classifier.kern_norms[kidx],
                                       order=classifier.ks[kidx],
                                       rev_comp=classifier.rev_comps[kidx])
            kernel_mats.append(km)
#            print '\n***', classifier.kern_names[kidx]
#            print 'raw_features:', raw_features, len(raw_features), len(raw_features[0])
#            print 'kernel_matrix:', km, km.shape

        kernel_matrix_sets.append(kernel_mats)

    classifier_names = [x.name for x in classifiers]
    mkl_norms = [x.mkl_norm for x in classifiers]

    c_mults=None
    #c_mults = [x.c for x in classifiers]


#    print '\nWriting full classifiers...'
#    write_full_classifiers(kernel_sets, classifier_names, raw_labels, 
#                           c_mults=c_mults, outdir=OUT_DIR)
   
    print '\nPerforming Cross Validation...'
    cv_results, subk_weights = cross_validation_shogun(kernel_matrix_sets, raw_labels, 
                                                       cv_rounds=CV_ROUNDS,
                                                       c_mults=c_mults,
                                                       mkl_norms=mkl_norms,
                                                       classifier_names=classifier_names,
                                                       region_locs=['%s\t%d\t%d' % tuple(x) for x in regions],
                                                       write_cv_scores=(OUT_DIR + OUT_NAME + "_cv_test_scores", PREAMBLE), 
                                                       average_kernels_no_mkl=AVERAGE_KERNELS)

    # create more useful Results objects
    classifier_results = [Results(cv_result) for cv_result in cv_results]

    
    # print and write classifier performance summary
    classifier_labels_roc = []
    classifier_labels_pr = []
    out_fn = '%s%s_classifier_performance.txt' % (OUT_DIR, OUT_NAME)
    of = open(out_fn, 'w')
    of.write(PREAMBLE + '\n')
    for c_idx, classifier in enumerate(classifier_names):
        print '\n', classifier
        of.write('\n' + classifier + '\n')
        print classifier_results[c_idx]
        of.write('%s\n' % classifier_results[c_idx])
        classifier_labels_roc.append(classifier_names[c_idx] + ' (%.2f)' %\
            classifier_results[c_idx].apply_to_stats(np.average)['roc_auc'])
        classifier_labels_pr.append(classifier_names[c_idx] + ' (%.2f)' %\
            classifier_results[c_idx].apply_to_stats(np.average)['pr_auc'])
    of.close()

    print '\nSubkernel Weights:'
    out_fn = '%s%s_subkernel_weights.txt' % (OUT_DIR, OUT_NAME)
    of = open(out_fn, 'w')
    of.write(PREAMBLE + '\n')
    skw_summary = ''
    for i, skws in enumerate(subk_weights):
        if len(skws[0]) == 1: 
            skw_summary = [1]
            of.write('%s_RAW\t1\n' % (classifier_names[i]))
            of.write('%s\t1\n' % (classifier_names[i]))
        else:
            skw_summary = ['%.3f' % x for x in sum(skws) / np.sum(skws)]
            of.write('%s_RAW\t%s\n' % (classifier_names[i], ' '.join(['%.3f' % x for x in sum(skws)])))
            of.write('%s\t%s\n' % (classifier_names[i], ' '.join(skw_summary)))

        print classifier_names[i], skw_summary
        for j, kernel_name in enumerate(classifiers[i].kern_names):
            print '\t', skw_summary[j], kernel_name
            of.write('\t%s\t%s\n' % (skw_summary[j], kernel_name))

    of.close()
 

    sample_points = [x / 100. for x in range(100)]

    roc_curves = []; roc_curves_ranges = []
    min_roc_curves = []; max_roc_curves = []
    for result in classifier_results:
        roc_curve, roc_curve_ranges = result.average_roc_curves(sample_points)
        roc_curves.append(roc_curve)
        roc_curves_ranges.append(roc_curve_ranges)
        #min_roc_curves.append([roc_curve[1][i] - stats.sem(x) for i, x in enumerate(roc_curve_ranges)])
        #max_roc_curves.append([roc_curve[1][i] + stats.sem(x) for i, x in enumerate(roc_curve_ranges)])
        min_roc_curves.append([sorted(x)[1] for x in roc_curve_ranges])
        max_roc_curves.append([sorted(x)[-2] for x in roc_curve_ranges])
        #min_roc_curves.append([min(x) for x in roc_curve_ranges])
        #max_roc_curves.append([max(x) for x in roc_curve_ranges])

    #roc_curves = [x.average_roc_curves(sample_points) for x in classifier_results] 


    # write ROC curves
    for i, curve in enumerate(roc_curves):
        out_fn = '%s%s__%s.roc_curve' % (OUT_DIR, OUT_NAME, classifier_names[i])
        of = open(out_fn, 'w')
        of.write(classifier_labels_roc[i] + '\n')
        of.write(' '.join([str(x) for x in curve[0]]) + '\n')
        of.write(' '.join([str(x) for x in curve[1]]) + '\n')
        of.write(' '.join([str(x) for x in min_roc_curves[i]]) + '\n')
        of.write(' '.join([str(x) for x in max_roc_curves[i]]) + '\n')
        of.close()

    # make PDF
    out_fn = '%s%s_roc.pdf' % (OUT_DIR, OUT_NAME.replace(' ', '_'))
    print '\nWriting %s...' % out_fn
    plot_curves(roc_curves, labels=classifier_labels_roc, filename=out_fn, 
                title=OUT_NAME.replace('~~', '\n'), plot_type='ROC', cis_low=min_roc_curves, cis_high=max_roc_curves)


    pr_curves = []
    min_pr_curves = []; max_pr_curves = []
    for result in classifier_results:
        pr_curve, min_pr_curve, max_pr_curve = result.average_pr_curves(sample_points)
        pr_curves.append(pr_curve)
        min_pr_curves.append(min_pr_curve)
        max_pr_curves.append(max_pr_curve)

    #pr_curves = [x.average_pr_curves(sample_points) for x in classifier_results] 

#    for c_idx, classifier in enumerate(classifier_names):
#        class_results = classifier_results[c_idx]
#        pcs = [(x[1], x[0]) for x in class_results.pr_curves]
#        for pidx, pc in enumerate(pcs):
#            plot_curves([pc], filename='poop_%s_%d.pdf' % (classifier, pidx), title=classifier+str(pidx), plot_type='PR')

    # write PR curves
    for i, curve in enumerate(pr_curves):
        out_fn = '%s%s__%s.pr_curve' % (OUT_DIR, OUT_NAME, classifier_names[i])
        of = open(out_fn, 'w')
        of.write(classifier_labels_pr[i] + '\n')
        of.write(' '.join([str(x) for x in curve[0]]) + '\n')
        of.write(' '.join([str(x) for x in curve[1]]) + '\n')
        of.write(' '.join([str(x) for x in min_pr_curves[i]]) + '\n')
        of.write(' '.join([str(x) for x in max_pr_curves[i]]) + '\n')
        of.close()

    # make PDF
    out_fn = '%s%s_pr.pdf' % (OUT_DIR, OUT_NAME.replace(' ', '_'))
    print '\nWriting %s...' % out_fn
    plot_curves(pr_curves, labels=classifier_labels_pr, filename=out_fn, 
                title=OUT_NAME.replace('~~', '\n'), plot_type='PR', cis_low=min_pr_curves, cis_high=max_pr_curves)

    #croc_curves = [zip(*croc.Curve.average(x.croc_curves)._coord) for x in classifier_results]


    # End main()


if __name__ == "__main__": 
    main()

