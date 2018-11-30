#!/usr/bin/env python
"""
 predict_enhancers_shogun.py - Copyright Tony Capra 2012

 Change log: 
  - 03/11/12 started 

 -----------------------------------------------------------------------------

"""

usage = """
USAGE:
predict_enhancers_shogun.py [options] blah blah

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
import itertools

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pylab as plt

from scipy.interpolate import interp1d

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
        opts, args = getopt.getopt(sys.argv[1:], "hc:f:n:o:p:s:t:w")
    except getopt.GetoptError:
        sys.exit(usage)

    if len(args) > 1:
        CLASSIFIER_CONFIG = args[0].split(',')
        POS_BEDS = args[1].split(',')
        NEG_BEDS = args[2].split(',')
        PREDICT_BED = args[3]
    else:
        sys.exit(usage)

    CLASSIFIER_NAME = None  # -c  # select a classifier from config file

    INTERSECT_STAT = 'binary' # -t  binary or coverage

    now_str = datetime.datetime.now().strftime('%Y-%m-%d_%H.%M.%S')
    OUT_NAME = '%s__%s__%s' % (now_str[:10], 
                               CLASSIFIER_CONFIG[0].split('/')[-1].replace('.classifier_config', ''), 
                               PREDICT_BED.split('/')[-1].replace('.bed', ''))  # -n
    OUT_DIR = None  # -o
    
    SEQ_FILES = None  # -s

    ONLY_COMPUTE_WEIGHTS = False  # -w

    for opt, arg in opts:
        if opt == "-h":
            sys.exit('Use me properly.')

        elif opt == '-c':
            CLASSIFIER_NAME = arg

        elif opt == '-n':
            OUT_NAME = arg

        elif opt == '-o':
            OUT_DIR = arg

        elif opt == '-t':
            INTERSECT_STAT = arg

        elif opt == '-w':
            ONLY_COMPUTE_WEIGHTS = True


    if OUT_NAME[:10] != now_str[:10]:
        OUT_NAME = now_str[:10] + "_" + OUT_NAME

    if OUT_DIR == None:
        OUT_DIR = 'results/predictions/%s/' % OUT_NAME.replace(' ', '_')

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
    train_regions, raw_labels = make_shogun_labels(POS_BEDS, NEG_BEDS)
    train_labels = BinaryLabels(np.array(raw_labels, dtype=np.double))


    print '\tNum Pos: %d\n\tNum Neg: %d\n' % (raw_labels.count(1), 
                                              raw_labels.count(-1))

    # combine +/- examples into a single bed file for convenience
    rand_num = np.random.randint(1000000)
    ALL_BED = 'all_regions_TEMP_%d.bed' % rand_num
    of = open(ALL_BED, 'w')
    for region in train_regions: 
        of.write('%s\t%d\t%d\n' % (region[0], region[1], region[2]))
    of.close()

    # This could be disasterous if the order in ALL_BED doesn't match
    # that in the regions variable. So we just make the combined bed
    # from the regions array to ensure the same order.
    #os.system('cat %s %s | cut -f 1-3 >%d.bed' % (' '.join(POS_BEDS), ' '.join(NEG_BEDS), rand_num))
    #os.system('sortBed -i %d.bed > %s' % (rand_num, ALL_BED))
    #os.remove('%d.bed' % rand_num)
    ##os.system('cat %s %s >%s' % (POS_BED, NEG_BED, ALL_BED))

    # READ CONFIG FILE
    print 'Initializing classifiers:'
    classifiers = parse_config_file(CLASSIFIER_CONFIG)
    classifier_names = [x.name for x in classifiers]

    classifier = None
    if CLASSIFIER_NAME == None:
        classifier = classifiers[0]
    elif CLASSIFIER_NAME in classifier_names:
        classifier = classifiers[classifier_names.index(CLASSIFIER_NAME)]
    else:
        sys.exit('ERROR! %s not defined in %s.' % (CLASSIFIER_NAME, CLASSIFIER_CONFIG))

    print '\tUsing: ', classifier.name

    # BUILD FEATURES
    print "\nLoading training features..."
    train_features, ind_train_features = build_shogun_features_for_classifier(classifier, ALL_BED)
    # When doing MKL, train_features is a combined feature object. For
    # some reason when I iterate through the feature objects,
    # get_feature_matrix() isn't found. This means I can't compute
    # within kernel feature weights when doing MKL. As a hack, I'm
    # also returning each individual feature, so that
    # get_feature_matrix() and get_features() (for StringFeatures) can
    # be called.   

    os.remove(ALL_BED)

    # TRAIN MODELS
    print "\nInitializing training kernels..."
    train_kernel = None
    if len(classifier.kernels) > 1: 
        train_kernel = CombinedKernel()
        for kernel in classifier.kernels: train_kernel.append_kernel(kernel)
    else:
        train_kernel = classifier.kernels[0]

    train_kernel.init(train_features, train_features)

    print '\nTraining models...'
    svm = None
    if len(classifier.kernels) > 1:
        svm = MKLClassification()
        svm.set_mkl_norm(classifier.mkl_norm)
        # not implemented in Classifier class:
        #svm.set_interleaved_optimization_enabled(classifier.inter_opt)
    else:
        svm = LibSVM()

    Cmult = None
    denom = float(len(raw_labels)) / Cmult if Cmult else float(len(raw_labels))
    svm.set_C(raw_labels.count(1)/denom, raw_labels.count(-1)/denom)

    svm.set_kernel(train_kernel)
    svm.set_labels(train_labels)
    svm.train()


    # COMPUTE WEIGHTS FOR KERNELS AND FEATURES
    # add subk weights to preamble
    kern2weight = {}
    if len(classifier.kernels) > 1: 
        print '\tSubkernel Weights:'
        subk_weights = train_kernel.get_subkernel_weights()
        w_sum = np.sum(subk_weights)
        k_names = classifier.kern_names
        for i, kernel in enumerate(classifier.kernels):
            s = '\t%s: %.2f (%.2f)' % (k_names[i], subk_weights[i] / w_sum, subk_weights[i])
            print s
            kern2weight[k_names[i]] = (subk_weights[i] / w_sum, subk_weights[i])
    else:
        kern2weight[classifier.kern_names[0]] = (1, 1)

    # calculate the feature weights with each kernel a la Guyon et al. 2002 
    tfs = [train_features] if not ind_train_features else ind_train_features
    kern2feat_weights = calc_within_kernel_weights(svm, classifier, tfs, raw_labels)

    # write subkernel and feature weights
    weight_fn = '%s/%s_kernel-feature-weights.txt' % (OUT_DIR, OUT_NAME)
    weight_of = open(weight_fn, 'w')
    weight_of.write(PREAMBLE + '\n')
    weight_of.write('# %d support vectors\n\n' % len(svm.get_support_vectors())) 

    for kern in sorted(kern2weight):
        weight_of.write('\n>%s\t%.3g\t%.3g\n' % (kern, kern2weight[kern][0], kern2weight[kern][1]))
        weight_of.write(kern2feat_weights[kern])

    weight_of.close()

    if ONLY_COMPUTE_WEIGHTS: sys.exit()

    svs = svm.get_support_vectors()
    #    print svs, len(svs)
    #    print svm.get_alphas(), len(svm.get_alphas())
    #    print svm.get_bias()

    km = svm.get_kernel().get_kernel_matrix()
    #    asvs =  list(km[:,svs])
    #    print len(asvs), len(asvs[0])


    # LOAD CLASSIFY REGIONS AND FEATURES
    predict_regions, predict_region_data = regions_from_bed(PREDICT_BED, return_data=True)
    print "\nLoading predict features..."
    predict_features, ind_pred_features = build_shogun_features_for_classifier(classifier, PREDICT_BED)

    # PREDICT ON CLASSIFY REGIONS
    print "\nClassifying input regions..."
    pred_labels = list(svm.apply(predict_features).get_labels())
    pred_vals = list(svm.apply(predict_features).get_values())

#    print predict_features.get_feature_matrix()
    print pred_vals

    # write predictions to a bed file
    pred_fn = '%s/%s_preds.bed' % (OUT_DIR, OUT_NAME)
    pred_of = open(pred_fn, 'w')
    pred_of.write(PREAMBLE + '\n')

    for i, pred_label in enumerate(pred_labels):
        pred_val = pred_vals[i]
        loc_str = '\t'.join([str(x) for x in predict_regions[i]])
        data_str = '\t'.join([str(x) for x in predict_region_data[i]])

        print '%s\t%s\t%.3f\t%.0f' % (loc_str, data_str, pred_val, pred_label)
        pred_of.write('%s\t%s\t%.3f\t%.0f\n' % (loc_str, data_str, pred_val, pred_label))

    pred_of.close()


    # End main()


if __name__ == "__main__": 
    main()

