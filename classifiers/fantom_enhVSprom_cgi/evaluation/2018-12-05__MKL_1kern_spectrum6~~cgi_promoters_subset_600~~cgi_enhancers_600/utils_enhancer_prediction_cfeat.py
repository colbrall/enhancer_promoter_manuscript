#!/usr/bin/env python
"""
 utils_enhancer_prediction.py - Copyright Tony Capra 2012

 Various utilities and definitions used in predicting enhancers with SHOGUN.

 -----------------------------------------------------------------------------

"""

usage = """ """


import os
import sys
import getopt
import random
import datetime
import commands
import itertools

#import pdb

import numpy as np

import matplotlib
#matplotlib.use("Agg")
import matplotlib.pylab as plt
#matplotlib.rcParams.update({'font.size': 16})

from scipy.interpolate import interp1d

from shogun.Features import *
from shogun.Classifier import *
from shogun.Kernel import *
from shogun.Evaluation import *
from shogun.ModelSelection import *


########################################
# CLASSES
########################################

class Classifier:
    """ """
    def __init__(self, name, kernels, bin_features, con_features, pp_features=None, seqs=None, 
                 kern_names=None, kern_norms=None, ks=None, c=None, mkl_norm=None, rev_comps=None):
        self.name = name
        self.kernels = kernels
        self.bin_features = bin_features
        self.con_features = con_features
        self.pp_features = pp_features
        self.kern_names = kern_names
        self.kern_norms = kern_norms
        self.rev_comps = rev_comps
        self.seqs = seqs
        self.ks = ks
        self.c = c
        self.mkl_norm = mkl_norm

    def __str__(self):
        return self.name

#    def make_names_for_kernels(self):
#        names = []
#        for i, kernel in enumerate(self.kernels):
#            files = []
#            for file in self.features[i]: 
#                files.append(file.split('/')[-1].replace('.bed', ''))
#            for file in self.pp_features[i]: 
#                files.append(file.split('/')[-1].replace('.bed', ''))
#            for file in self.seqs[i]: 
#                files.append(file.split('/')[-1].replace('.fasta', '').replace('.fa', ''))
#            names.append(';'.join(files))
            
        return names


class Results:
    """ """
    def __init__(self, result_list):
        self.stats = {
        'accuracy': [], 'error_rate': [], 'balanced_err': [], 'WRACC': [], 
        'F1': [], 'cross_corr': [],'recall': [], 'precision': [], 
        'specificity': [], 'pr_auc': [], 'roc_auc': [],
        }
        self.pr_curves = []
        self.roc_curves = []

        for result in result_list:
            ct_eval, pr_eval, roc_eval = result

            self.stats['accuracy'].append(ct_eval.get_accuracy())
            self.stats['error_rate'].append(ct_eval.get_error_rate())
            self.stats['balanced_err'].append(ct_eval.get_BAL())
            self.stats['WRACC'].append(ct_eval.get_WRACC())
            self.stats['F1'].append(ct_eval.get_F1())
            self.stats['cross_corr'].append(ct_eval.get_cross_correlation())
            self.stats['recall'].append(ct_eval.get_recall())
            self.stats['precision'].append(ct_eval.get_precision())
            self.stats['specificity'].append(ct_eval.get_specificity())

            self.stats['roc_auc'].append(roc_eval.get_auROC())
            #self.stats['pr_auc'].append(pr_eval.get_auPRC())
            self.stats['pr_auc'].append(self.interpolate_pr_auc(pr_eval.get_PRC()))

            self.roc_curves.append(roc_eval.get_ROC())
            self.pr_curves.append(pr_eval.get_PRC())


    def __str__(self):
        s = ''
        for stat in sorted(self.stats):
            data = self.stats[stat]
            s += '%s %.3f (%d)\n' % (stat.ljust(20), np.average(data), len(data))
            #s += '%s %s\n' % (stat.ljust(20), ' '.join(['%.3f' % x for x in data]))

        return s


    def apply_to_stats(self, function):
        d = {}
        for stat in self.stats:
            d[stat] = function(self.stats[stat])

        return d


    def interpolate_pr_auc(self, curve, max_prec=False):
        """ """
        r_pts = []; p_pts = [];
        for idx, prec in enumerate(curve[0]):
            rec = curve[1][idx]
            # remove 0,0 points, calculation starts with recall 1/num_pos
            if prec == 0 or rec == 0: continue 
            # make sure only one recall value is kept
            if rec in r_pts: continue

            r_pts.append(rec)
            # for each recall level r1, either append a random
            # precision or the max precision over all r2 >= r1
            p_pts.append(max(curve[0][idx:]) if max_prec else prec)

            if rec == 1.0: break

        # get rid of bit at beginning of recall range
        pr_auc = np.trapz(p_pts, r_pts) / (1. - r_pts[0])

        return pr_auc
               

    def average_pr_curves(self, sample_points, max_prec=False):
        """ """
        avg_curve = [0] * len(sample_points)
        min_curve = [1] * len(sample_points)
        max_curve = [0] * len(sample_points)

        old_settings = np.seterr(all='ignore')
        for curve in self.pr_curves:
            int_curve = [0] * len(sample_points)
# This is the line to look at if gaps in PR curve:            
            #f = interp1d(curve[1], curve[0], bounds_error=False, fill_value=max(curve[0]))
            f = interp1d(curve[1], curve[0], bounds_error=False, fill_value=max(curve[0]),kind="zero")
            for i, sp in enumerate(sample_points):
                int_curve[i] += f(sp)

            for i, sp in enumerate(sample_points):
                interp_val = max(int_curve[i:]) if max_prec else int_curve[i]
                avg_curve[i] += interp_val
                if interp_val < min_curve[i]: min_curve[i] = interp_val
                if interp_val > max_curve[i]: max_curve[i] = interp_val

        np.seterr(**old_settings)
        
        num_curves = float(len(self.pr_curves))
#        print avg_curve
        return (sample_points, [x / num_curves for x in avg_curve]), min_curve, max_curve


    def average_roc_curves(self, sample_points):
        """ """

        print >> sys.stderr, "WARNING: THIS IS FOR LAURA'S TEMPORARY USE ONLY!"
        plot_curves(self.roc_curves, filename="raw_cv_curves.pdf")

        avg_curve = [0] * len(sample_points)
        curve_ranges = [] # hold all values observed across CV runs at each sample_point
        for x in sample_points: curve_ranges.append([])

        ics = []
        
        old_settings = np.seterr(all='ignore')
        for curve in self.roc_curves:
            c0, c1 = prep_curve_for_interp(curve)
# This is the line to look at if gaps in ROC curve:
            #f = interp1d(c0, c1, kind="zero")
            f = interp1d(c0, c1, kind="linear")  
            #f = interp1d(curve[0], curve[1], kind="zero")
            ic = []
            for i, sp in enumerate(sample_points):
                interp_val = f(sp)
                avg_curve[i] += interp_val
                curve_ranges[i].append(interp_val)

                ic.append(interp_val)

            ics.append((sample_points, ic))
            #print "interp:", ic

        np.seterr(**old_settings) 

        num_curves = float(len(self.roc_curves))

        plot_curves(ics, filename="interp_cv_curves.pdf")
        for i, raw in enumerate(self.roc_curves):
#            print "\n", i
#            print raw[0]
#            print raw[1]
            plot_curves([raw, ics[i]], filename="compare_cv_curves%d.pdf" % i)

        return (sample_points, [x / num_curves for x in avg_curve]), curve_ranges

#    def average_roc_curves(self, sample_points):
#        """ """
#        avg_curve = [0] * len(sample_points)
#        curve_ranges = [] # hold all values observed across CV runs at each sample_point
#        for x in sample_points: curve_ranges.append([])
#
#        old_settings = np.seterr(all='ignore')
#        for curve in self.roc_curves:
#            f = interp1d(curve[0], curve[1])
#            for i, sp in enumerate(sample_points):
#                interp_val = f(sp)
#                avg_curve[i] += interp_val
#                curve_ranges[i].append(interp_val)
#
#        np.seterr(**old_settings) 
#
#        num_curves = float(len(self.roc_curves))
#       
#        return (sample_points, [x / num_curves for x in avg_curve]), curve_ranges
#
                


################################################################################
## GENERAL FUNCTIONS
################################################################################

def prep_curve_for_interp(curve):
    """ """
    #xs = [0]; ys = [0]
    xs = []; ys = []

    for idx, x in enumerate(curve[0]):
        if idx < (len(curve[0])-1) and  x != curve[0][idx+1]: 
            xs.append(x)
            ys.append(curve[1][idx])

    xs.append(1.); ys.append(1.)
    #print xs, ys
    return [xs,ys]


def run_bedops_intersect(bed1, bed2, bed1_regions=None):
    """ """
    if bed1_regions == None:
        bed1_regions = regions_from_bed(bed1)

    result_dict = {}
    for loc in bed1_regions: result_dict.setdefault(loc, [])

    cmd = 'bedops -e -1 %s %s' % (bed1, bed2)

    #print >> sys.stderr, cmd
    status, output = commands.getstatusoutput(cmd)

    for line in output.split('\n'):
        if line.startswith('#'): continue
        t = line.split('\t')
        if len(t) < 3: continue

        loc = (t[0], int(t[1]), int(t[2]))

        result_dict[loc].append(1)

    return result_dict


def run_intersectBed(bed1, bed2, consider_value=False):
    """ """

    result_dict = {}

    cmd = 'intersectBed -wao -a %s -b %s' % (bed1, bed2)

    #print >> sys.stderr, cmd
    status, output = commands.getstatusoutput(cmd)

    for line in output.split('\n'):
        t = line.split('\t')
        if line.startswith('#'): continue
        
        loc = (t[0], int(t[1]), int(t[2]))
        result_dict.setdefault(loc, [])

        overlap = int(t[-1])

        feature_value = t[-2] # just a string
        if consider_value:  # make a float
            if overlap == 0:
                feature_value = 0.0
            else:
                feature_value = float(feature_value) if '=' not in feature_value else float(feature_value.split('=')[-1])

        #result_dict.setdefault(loc, []).append(feature_value)
        if overlap > 0: result_dict[loc].append(feature_value)

    return result_dict


def run_coverageBed(bed1, bed2):
    """ """
    result_dict = {}

    cmd = 'coverageBed -a %s -b %s' % (bed1, bed2)

    #print >> sys.stderr, cmd
    status, output = commands.getstatusoutput(cmd)

    for line in output.split('\n'):
        t = line.split('\t')
        if line.startswith('#'): continue
        loc = (t[0], int(t[1]), int(t[2]))

        result_dict[loc] = float(t[-1])

    return result_dict


def feature_from_preprocessed(regions, pp_feat_file, handle_missing='error'):
    """ Assumes feature is in last column. """
    d = {}
    for line in open(pp_feat_file):
        t = line.split('\t')
        if line.startswith('#'): continue
        if len(t) < 4: continue

        reg = (t[0], int(t[1]), int(t[2]))
        d[reg] = float(t[-1])

    feat = []
    for region in regions:
        if region in d:
            feat.append(d[region])
        else:
            if handle_missing == 'zero':
                feat.append(0.0)
            elif handle_missing == 'error':
                s = "ERROR: %s doesn't have preprocessed features in %s." % \
                    (region, pp_feat_file)
                sys.exit(s)
            else:
                s = "ERROR: %s doesn't have preprocessed features in %s." % \
                    (region, pp_feat_file)
                sys.exit(s)
                
    return feat

def features_from_preprocessed(regions, pp_feat_file, handle_missing='error'):
    """ Return a list of feature names and a list of lists with each
    internal list containing the one feature types values for all the
    regions. All colums after the first three are treated as
    features. The first line of the pp_feat_file must have a header."""

    print "THIS IS NOT CURRENTLY IN USE!!!!!!!!!!!!"
    d = {}
    feature_names = []
    first_line = True
    for line in open(pp_feat_file):
        t = line[:-1].split('\t')

        if line.startswith('#'): continue
        if len(t) < 4: continue

        if first_line:
            first_line = False
            feature_names = t[3:]
            continue

        reg = (t[0], int(t[1]), int(t[2]))
        d[reg] = [float(x) for x in t[3:]]

    feats = []       
    for region in regions:
        if region in d:
            feats.append(d[region])
        else:
            if handle_missing == 'zero':
                feats.append([0.0] * len(feature_names))
            else:
                s = "ERROR: %s doesn't have preprocessed features in %s." % \
                    (region, pp_feat_file)
                sys.exit(s)

    return feature_names, zip(*feats)
       

def parse_fasta_file(fasta_file):
    """ Make sure to make sequence consistent case!  """
    nucleic_alphabet = ['A', 'T', 'C', 'G', 'a', 't', 'c', 'g', 'N', 'n']

    region_seq = []

    cur_name = None
    cur_seq = ''
    for line in open(fasta_file):
        line = line[:-1]
        if len(line) == 0: continue

        if line.startswith(';') or line.startswith('#'): continue

        if line.startswith('>'):
            if cur_seq != '':
                region_seq.append([cur_name, cur_seq.upper()])
                cur_seq = ''

            chrom, loc = line[1:].split(':')
            start, end = loc.split('-')
            cur_name = (chrom, int(start), int(end))


        elif line[0] in nucleic_alphabet:
            cur_seq += line
            
    region_seq.append([cur_name, cur_seq.upper()])

    return region_seq

def prepare_sequence_data(fasta_paths, regions):
    """ """
    seqs = []

    region_seq = []
    for fasta_path in fasta_paths:
        if os.path.isdir(fasta_path):
            for f in os.listdir(fasta_path):
                if f.endswith('.fa') or f.endswith('.fasta'):
                    region_seq += parse_fasta_file(fasta_path + f)
        else:
            region_seq += parse_fasta_file(fasta_path)

    region2seq = dict(region_seq)

    for region in regions:
        if region in region2seq:
            seqs.append(region2seq[region].upper().replace('N', ''))
        else:
            print 'WARNING: no sequence data for %s in %s.' % (region, fasta_paths)
            return None

    return seqs
       
#def prepare_sequence_data(fasta_files, regions):
#    """ """
#    seqs = []
#
#    region_seq = []
#    for fasta_file in fasta_files:
#        region_seq += parse_fasta_file(fasta_file)
#
#    region2seq = dict(region_seq)
#
#    for region in regions:
#        if region in region2seq:
#            seqs.append(region2seq[region].upper().replace('N', ''))
#        else:
#            print 'WARNING: no sequence data for %s in %s.' % (region, fasta_files)
#            return None
#
#    return seqs
#       

def plot_curves(curves, filename='test_curves.pdf', text=None, labels=[], title='', 
                plot_type='ROC', cis_low=None, cis_high=None):
    """ """
     
    fig = plt.figure()
    fig.set_size_inches(8, 8)
    ax1 = fig.add_subplot(111)

    plt.title(title)

    plt.grid(color='grey', linestyle=':')

    if cis_low and cis_high:
        for i, curve in enumerate(curves):
            plt.fill_between(curve[0], cis_low[i], cis_high[i], alpha=0.3, color='gray')

    plots = []
    for i, curve in enumerate(curves):
        p, = plt.plot(curve[0], curve[1], lw=2)
        plots.append(p)
        #if plot_type == 'ROC': print '\n', labels[i], curve

    if plot_type == 'ROC':
        leg_loc = 'lower right'
        plt.xlabel("False Positive Rate")
        plt.ylabel('True Positive Rate')
        plt.plot([0,1], color="grey", ls="--")
        plt.xlim(0,1)
        plt.ylim(0,1)

    elif plot_type == 'PR':
        #leg_loc = 'upper right'
        leg_loc = 'lower right'
        plt.xlabel("Recall")
        plt.ylabel('Precision')
        plt.xlim(0,1)
        plt.ylim(0,1)

    if labels != []: 
        leg = plt.legend(plots, labels, loc=leg_loc)
        leg.draw_frame(False)

    plt.savefig(filename)



################################################################################
## SHOGUN FUNCTIONS
################################################################################

def regions_from_bed(bed_file, return_data=False):
    """ """
    regions = []
    region_data = []
    for line in open(bed_file):
        t = line[:-1].split('\t')
        if line.startswith('#') or len(t) < 3: continue

        regions.append((t[0], int(t[1]), int(t[2])))
        region_data.append(t[3:])

    if return_data:
        return regions, region_data
    else:
        return regions


def make_shogun_labels(pos_beds, neg_beds):
    """ """
    l = []

    for pos_bed in pos_beds:
        for line in open(pos_bed):
            t = line.split('\t')
            if line.startswith('#') or len(t) < 3: continue
            #if line.startswith('#'): continue

            l.append((t[0], int(t[1]), int(t[2]), 1))

    for neg_bed in neg_beds:
        for line in open(neg_bed):
            t = line.split('\t')
            if line.startswith('#') or len(t) < 3: continue
            #if line.startswith('#'): continue

            l.append((t[0], int(t[1]), int(t[2]), -1))

    l.sort()

    regions = [x[:3] for x in l]
    labels = zip(*l)[-1]

    return regions, labels

def reverse_complement(seq):
    """ """
    rc_seq = ''
    rc_d = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    for bp in reversed(seq):
        rc_seq += rc_d[bp] if bp in rc_d else bp

    if len(rc_seq) != len(seq): sys.exit("Couldn't reverse complement:" + seq)
        
    return rc_seq

def build_seq_features(seqs, order=4, gap=0, reverse=False, use_sign=False, rev_comp=False):
    """ """
    charfeat = StringCharFeatures(DNA)
    charfeat.set_features(seqs)

    #seq_features = StringWordFeatures(charfeat.get_alphabet())
    seq_features = StringWordFeatures(DNA)
    seq_features.obtain_from_char(charfeat, order-1, order, gap, reverse)

    if rev_comp:
        rc_seqs = [reverse_complement(x) for x in seqs]
        rc_charfeat = StringCharFeatures(DNA)
        rc_charfeat.set_features(rc_seqs)

        rc_seq_features = StringWordFeatures(DNA)
        rc_seq_features.obtain_from_char(rc_charfeat, order-1, order,
                                         gap, reverse)

        raw_sfs = seq_features.get_features()
        raw_rcsfs = rc_seq_features.get_features()

        for idx, rsf in enumerate(raw_sfs):
            seq_features.set_feature_vector(np.concatenate((rsf, raw_rcsfs[idx]), axis=0), idx)


    preproc = SortWordString()
    preproc.init(seq_features)
    seq_features.add_preprocessor(preproc)
    seq_features.apply_preprocessor()
#    for s in seq_features: print s
    return seq_features


def eval_preds(preds, labels, verbose=False):
    """ """
    pred_labels = BinaryLabels(np.array(preds, dtype=np.double))
    #pred_labels_binary = Labels(np.array([1 if x > 0 else -1 for x in preds], dtype=np.double))

    if verbose:
        print 'True labels:\n', labels
        print 'Predicted labels:\n', pred_labels

    ct_eval = ContingencyTableEvaluation()
    ct_eval.evaluate(pred_labels, labels)

    if verbose:
        print 'accuracy:\t%.3f' % ct_eval.get_accuracy()
        print 'error_rate:\t%.3f' % ct_eval.get_error_rate()
        print 'balanced_err:\t%.3f' % ct_eval.get_BAL()
        print 'WRACC:\t%.3f' % ct_eval.get_WRACC()
        print 'F1:\t%.3f' % ct_eval.get_F1()
        print 'cross_corr:\t%.3f' % ct_eval.get_cross_correlation()
        print 'recall:\t%.3f' % ct_eval.get_recall()
        print 'precision:\t%.3f' % ct_eval.get_precision()
        print 'specificity:\t%.3f' % ct_eval.get_specificity()

    pr_eval = PRCEvaluation()
    pr_eval.evaluate(pred_labels, labels)
    if verbose: print 'PR AUC:\t%.3f' % pr_eval.get_auPRC()

    roc_eval = ROCEvaluation()
    roc_eval.evaluate(pred_labels, labels)
    if verbose: print 'ROC AUC:\t%.3f' % roc_eval.get_auROC()

    return ct_eval, pr_eval, roc_eval



def MKL_classifier_from_kernel_matrices(train_kms, train_labels, Cmult=None,
                                        mkl_norm=None, Cneg=None, Cpos=None):
    """ """
    comb_kernel = CombinedKernel()
    for km in train_kms:
        comb_kernel.append_kernel(CustomKernel(km))

    mkl = MKLClassification()
    if mkl_norm: mkl.set_mkl_norm(mkl_norm) #1,2,3
    #mkl.set_solver_type(ST_ELASTICNET)
    #mkl.set_elasticnet_lambda(1)

    mkl.set_interleaved_optimization_enabled(False)

    if Cneg == None or Cpos == None:
        raw_labels = list(train_labels.get_labels())
        num_labels = float(len(raw_labels))
        denom = num_labels / Cmult if Cmult else num_labels
        mkl.set_C(raw_labels.count(1)/denom, raw_labels.count(-1)/denom)
        #mkl.set_C(raw_labels.count(1) / num_labels, raw_labels.count(-1) / num_labels)
    else:
        mkl.set_C(Cneg, Cpos)

    mkl.set_kernel(comb_kernel)
    mkl.set_labels(train_labels)

    mkl.train()

    return mkl, comb_kernel


def SVM_classifier_from_kernel_matrix(kernel_matrix, train_labels, Cmult=None,
                                      Cneg=None, Cpos=None, svm_type=LibSVM):
    """ """
    kernel = CustomKernel()
    kernel.set_full_kernel_matrix_from_full(kernel_matrix)

    svm = svm_type(1.0, kernel, train_labels)	

    if Cneg == None or Cpos == None:
        raw_labels = list(train_labels.get_labels())
        num_labels = float(len(raw_labels))
        denom = num_labels / Cmult if Cmult else num_labels
        svm.set_C(raw_labels.count(1)/denom, raw_labels.count(-1)/denom)
        #svm.set_C(raw_labels.count(1)/num_labels, raw_labels.count(-1)/num_labels)
    else:
        svm.set_C(Cneg, Cpos)

    svm.train()

    return svm


def classify_from_kernel_matrix(trained_classifier, test_kms, subk_weights=None):
    """ If subk_weights is not None, then do MKL."""

    svs = trained_classifier.get_support_vectors()
    bias = trained_classifier.get_bias()
    alphas = trained_classifier.get_alphas() 

    if np.any(subk_weights):
        if len(test_kms) != len(subk_weights):
            sys.exit('ERROR: Subkernel weights do not match kernels!')

        km = np.zeros(test_kms[0].shape)
        for i, test_km in enumerate(test_kms):
            km += test_km * subk_weights[i]
    else:
        km = test_kms

    preds = bias + np.dot(km[:,svs], alphas.transpose())

    return preds


def cross_validation_shogun(kernel_mat_sets, raw_labels, cv_rounds=10,
                            mkl_norms=None, c_mults=None,
                            classifier_names=None,
                            write_cv_scores=False,
                            average_kernels_no_mkl=False,
                            region_locs=None):
    """Perform CV.  If not False, then write_cv_scores should be a 2-tuple
        in which the first element is the base file name and the second is
        the preamble.  If average_kernels_no_mkl, then MKL is not
        performed and the individual kernels are averaged. region_locs
        can hold a string with the chromosome location of each input
        region; must have the same length as raw_labels.

    """

    pos_idxs = []; neg_idxs = []
    for idx, label in enumerate(raw_labels):
        if label == 1:
            pos_idxs.append(idx)
        elif label == -1:
            neg_idxs.append(idx)

    bins = np.linspace(0, len(pos_idxs), num=cv_rounds+1)
    pos_fold_assignments = np.digitize(range(0, len(pos_idxs)), bins)
    bins = np.linspace(0, len(neg_idxs), num=cv_rounds+1)
    neg_fold_assignments = np.digitize(range(0, len(neg_idxs)), bins)

    # since the regions are sorted, randomize to avoid any chromosome bias
    random.seed(12345678)
    random.shuffle(pos_fold_assignments)
    random.shuffle(neg_fold_assignments)

    cv_score_files = []
    if write_cv_scores:
        for kms_idx, kms in enumerate(kernel_mat_sets):
            fn = write_cv_scores[0] + '_%s.txt' % (classifier_names[kms_idx] if classifier_names else str(kms_idx))
            cv_score_files.append(open(fn, 'w'))
            cv_score_files[-1].write(write_cv_scores[1] + '\n') # should be the preamble 

    results = [[] for k in kernel_mat_sets]  # k results for each kernel
    subk_weights = [[] for k in kernel_mat_sets]  
    for cv_k in range(cv_rounds):
        print '\tCV round: ', cv_k
        # get indices (relative to input lists) for train and test
        # sets for this run
        test_idxs = []; train_idxs = []
        for idx, fold_assign in enumerate(pos_fold_assignments): 
            if fold_assign == cv_k + 1:
                test_idxs.append(pos_idxs[idx])
            else:
                train_idxs.append(pos_idxs[idx])

        for idx, fold_assign in enumerate(neg_fold_assignments): 
            if fold_assign == cv_k + 1:
                test_idxs.append(neg_idxs[idx])
            else:
                train_idxs.append(neg_idxs[idx])

        # create shogun label and feature objects for train and test
        test_labels = BinaryLabels(np.array(raw_labels, dtype=np.double)[test_idxs])     
        train_labels = BinaryLabels(np.array(raw_labels, dtype=np.double)[train_idxs])

        for kms_idx, kms in enumerate(kernel_mat_sets):

            c_mult = c_mults[kms_idx] if c_mults else None
            mkl_norm = mkl_norms[kms_idx] if mkl_norms else None

            k_subk_weights = [1]  # default for non-MKL
            # train and predict
            #if type(kms) != type([]):
            if len(kms) == 1:
                train_km = kms[0][:,train_idxs][train_idxs,:]
                test_km  = kms[0][test_idxs,:][:,train_idxs]

                train_svm = \
                    SVM_classifier_from_kernel_matrix(train_km, train_labels,
                                                      Cmult=c_mult)

                preds = classify_from_kernel_matrix(train_svm, test_km)

            else: # got multiple kernels, so we do MKL or average
                if average_kernels_no_mkl:
                    train_kms = [km[train_idxs,:][:,train_idxs] for km in kms]
                    test_kms  = [km[test_idxs,:][:,train_idxs] for km in kms]

                    #train_km = sum(train_kms) / float(len(train_kms))
                    #test_km = sum(test_kms) / float(len(test_kms))

                    train_km = sum(train_kms) / float(len(train_kms))
                    test_km = sum(test_kms) / float(len(test_kms))

                    train_svm = \
                        SVM_classifier_from_kernel_matrix(train_km, train_labels,
                                                          Cmult=c_mult)

                    preds = classify_from_kernel_matrix(train_svm, test_km)

                    k_subk_weights = np.array([0] * len(train_kms))  # just so that this doesn't crash

                else:
                    #train_kms = [km[:,train_idxs][train_idxs,:] for km in kms]
                    #test_kms  = [km[test_idxs,:][:,train_idxs] for km in kms]

                    train_kms = [km[train_idxs,:][:,train_idxs] for km in kms]
                    test_kms  = [km[test_idxs,:][:,train_idxs] for km in kms]

                    train_mkl, train_comb_kernel = \
                        MKL_classifier_from_kernel_matrices(train_kms, train_labels,
                                                            mkl_norm=mkl_norm, 
                                                            Cmult=c_mult)

                    k_subk_weights=train_comb_kernel.get_subkernel_weights()
                    #print '\t'.join(['%.2f' % x for x in k_subk_weights / sum(k_subk_weights)])
                    #for km in kms: print km.min(), km.max()

                    #k_subk_weights = np.array([1. / len(kms)] * len(kms))  # test equal weight combo

                    preds = classify_from_kernel_matrix(train_mkl, test_kms, 
                                                        subk_weights=k_subk_weights)

            # now evaluate predictions
            result = eval_preds(preds, test_labels, verbose=False)  #ct_eval, pr_eval, roc_eval

            results[kms_idx].append(result)
            subk_weights[kms_idx].append(k_subk_weights)

            if write_cv_scores:
                relevant_labels = list(np.array(raw_labels, dtype=np.double)[test_idxs])
                relevant_locs = list(np.array(region_locs)[test_idxs])
                for x in sorted(zip(relevant_locs, list(preds), relevant_labels)):
                    cv_score_files[kms_idx].write("%s\t%g\t%d\t%d\n" % (x[0], x[1], x[2], cv_k))
        
    if write_cv_scores: 
        for f in cv_score_files: f.close()

    return results, subk_weights



def compute_kernel_matrix(kernel, raw_features, norm=None, order=4, rev_comp=False):
    """norm must be a SHOGUN normalization class like VarianceKernelNormalizer """
    features = None

    if type(raw_features[0]) == type(''):
        features = build_seq_features(raw_features, order=order, rev_comp=rev_comp)
        #features = build_seq_features(raw_features, order=k_idx+3)
    else:
        features = RealFeatures(np.array(raw_features))
        
    if norm:
        kernel.set_normalizer(norm)

    kernel.init(features, features)

    return kernel.get_kernel_matrix()


def write_full_classifiers(kernel_mat_sets, set_names, raw_labels, 
                           c_mults=None, outdir=''):
    all_labels = BinaryLabels(np.array(raw_labels, dtype=np.double))

    for km_idx, kms in enumerate(kernel_mat_sets):

        c_mult = c_mults[km_idx] if c_mults else None

        # train and write
        subk_weights = [1]
        if type(kms) != type([]):
            svm = SVM_classifier_from_kernel_matrix(kms, all_labels, Cmult=c_mult)

#            out_fn = '%s%s.shogun_classifier' % (outdir, set_names[km_idx].replace(' ', '_'))
#            print '\twriting %s...' % out_fn
#
#            fstream = SerializableAsciiFile(out_fn, "w")
#            status = svm.save_serializable(fstream)

        else: # got multiple kernels, so we do MKL
            svm, comb_kernel = \
                MKL_classifier_from_kernel_matrices(kms, all_labels, Cmult=c_mult)

            subk_weights = comb_kernel.get_subkernel_weights()

#            out_fn = '%s%s.shogun_labels' % (outdir, set_names[km_idx].replace(' ', '_'))
#            print '\twriting %s...' % out_fn
#            fstream = SerializableAsciiFile(out_fn, "w")
#            status = all_labels.save_serializable(fstream)
#            print status
#
#            out_fn = '%s%s.shogun_kernels' % (outdir, set_names[km_idx].replace(' ', '_'))
#            print '\twriting %s...' % out_fn
#            fstream = SerializableAsciiFile(out_fn, "w")
#            status = comb_kernel.save_serializable(fstream)
#            print status

#            print mkl.get_support_vectors()
#            out_fn = '%s%s.shogun_classifier' % (outdir, set_names[km_idx].replace(' ', '_'))
#            print '\twriting %s...' % out_fn
#            fstream = SerializableAsciiFile(out_fn, "w")
#            status = mkl.save_serializable(fstream)

        out_fn = '%s%s.shogun_classifier' % (outdir, set_names[km_idx].replace(' ', '_'))
        print '\twriting %s...' % out_fn
        outfile = open(out_fn, 'w')
        s = ' '.join([str(x) for x in svm.get_support_vectors()])
        outfile.write('support_vectors\t%s\n' % s)
        s = ' '.join([str(x) for x in svm.get_alphas()])
        outfile.write('alphas\t%s\n' % s)
        outfile.write('bias\t%s\n' % svm.get_bias())
        s = ' '.join([str(x) for x in subk_weights])
        outfile.write('subk_weights\t%s\n' % s)
        outfile.close()

    return None



def parse_config_file(config_files, feature_file_suffix='.bed'):
    """ Return classifiers. """
    valid_words = ['CLASSIFIER', 'KERNEL', 'KERNEL_NAME', 'KERNEL_NORM', 'FEATURE', 
                   'BIN_FEATURE', 'CON_FEATURE', 'PP_FEATURE', 'SEQ', 'END', 'C', 'MKL_NORM', 'REV_COMP']

    classifiers = []

    for config_file in config_files:
        name = None; kernels = []; c = None
        bin_features = []; con_features = []; pp_features = []
        seqs = []; ks = [];
        kern_names = []
        kern_norms = []
        rev_comps = []
        mkl_norm = 2

        line_num = 0
        for line in open(config_file):
            line_num += 1
            line = line.strip()
            if line.startswith('#') or not line: continue

            t = line.split()
            if t[0] not in valid_words: 
                sys.exit('ERROR! %s not recognized in line %d: %s.' % (t[0], line_num, line))

            elif t[0] == 'CLASSIFIER':
                name = '_'.join(t[1:])

            elif t[0] == 'KERNEL':
                bin_features.append([])
                con_features.append([])
                pp_features.append([])
                seqs.append([])
                kern_norms.append(None)
                kern_names.append(None)
                rev_comps.append(False)
                ks.append(None)
                if t[1] == 'Linear':
                    kernels.append(LinearKernel())
                elif t[1].startswith('CommWordString_') or t[1].startswith('Spectrum_'):
                    k = int(t[1].split('_')[-1])
                    ks[-1] = k
                    kernels.append(CommWordStringKernel(10, False))
                elif t[1].startswith('Gaussian_'):
                    sigma = float(t[1].split('_')[-1])
                    kernels.append(GaussianKernel(10, sigma))
                else:
                    sys.exit('ERROR! %s is not a valid kernel.' % t[1])

            elif t[0] == 'KERNEL_NORM':
                if t[1] == 'VarianceKernelNormalizer':
                    kern_norms[-1] = VarianceKernelNormalizer()
                    kernels[-1].set_normalizer(VarianceKernelNormalizer())
                elif t[1] == 'SqrtDiagKernelNormalizer':
                    kern_norms[-1] = SqrtDiagKernelNormalizer()
                    kernels[-1].set_normalizer(SqrtDiagKernelNormalizer())
                elif t[1] == 'AvgDiagKernelNormalizer':
                    kern_norms[-1] = AvgDiagKernelNormalizer()
                    kernels[-1].set_normalizer(AvgDiagKernelNormalizer())
                else:
                    sys.exit('ERROR! %s is not a recognized kernel normalizer.' % t[1])

            elif t[0] == 'KERNEL_NAME':
                kern_names[-1] = t[1]

            elif t[0] == 'REV_COMP':
                rev_comps[-1] = bool(int(t[1]))

            elif t[0] == 'BIN_FEATURE' or t[0] == 'FEATURE':  # for BW compatibility
                bf_path = t[1]
                if os.path.isdir(bf_path):
                    for fn in os.listdir(bf_path):
                        if not fn.endswith(feature_file_suffix): continue
                        full_fn = bf_path + ('/' if bf_path[-1] != '/' else '') + fn
                        if full_fn not in bin_features[-1]: bin_features[-1].append(full_fn)
                else:
                    bin_features[-1].append(bf_path)

            elif t[0] == 'CON_FEATURE':
                cf_path = t[1]
                if os.path.isdir(cf_path):
                    for fn in os.listdir(cf_path):
                        if not fn.endswith(feature_file_suffix): continue
                        full_fn = cf_path + ('/' if cf_path[-1] != '/' else '') + fn
                        if full_fn not in con_features[-1]: con_features[-1].append(full_fn)
                else:
                    con_features[-1].append(cf_path)

            elif t[0] == 'PP_FEATURE':
                pf_path = t[1]
                if os.path.isdir(pf_path):
                    for fn in os.listdir(pf_path):
                        if not fn.endswith(feature_file_suffix): continue
                        full_fn = pf_path + ('/' if pf_path[-1] != '/' else '') + fn
                        if full_fn not in pp_features[-1]: pp_features[-1].append(full_fn)
                else:
                    pp_features[-1].append(pf_path)

            elif t[0] == 'SEQ':
                seqs[-1].append(t[1])

            elif t[0] == 'C':
                c = [float(x) for x in t[1:]]

            elif t[0] == 'MKL_NORM':
                mkl_norm = int(t[1])

            elif t[0] == 'END':
                print '\t', name
                classifiers.append(Classifier(name, kernels, bin_features, con_features,
                                              pp_features=pp_features, seqs=seqs, 
                                              kern_norms=kern_norms, ks=ks, c=c, 
                                              kern_names=kern_names, 
                                              mkl_norm=mkl_norm, rev_comps=rev_comps))

                name = None; kernels = []; c = None
                bin_features = []; con_features; pp_features = []
                seqs = []; ks = []
                kern_norms = []
                kern_names = []
                rev_comps = []
                mkl_norm = 2

    return classifiers


def build_features(classifiers, region_bed, regions=None):
    """ """
    name2feats = {}
    name2pp_feats = {}
    names2seqs = {}  # NOTE THAT SEQS MAPS MULTIPLE FILES, NOT NEC. ONE

    intersect_features = []; pp_features = []; seq_sets = []
    consider_vals = []  # for each intersect_feature tells whether to consider region values
    regions = regions_from_bed(region_bed)
    for classifier in classifiers:
        for feature_set in classifier.bin_features:
            for feature in feature_set:
                if feature not in intersect_features: 
                    intersect_features.append(feature)
                    consider_vals.append(False)
        for feature_set in classifier.con_features:
            for feature in feature_set:
                if feature not in intersect_features: 
                    intersect_features.append(feature)
                    consider_vals.append(True)
        for pp_feature_set in classifier.pp_features:
            for pp_feature in pp_feature_set:
                if pp_feature not in pp_features: pp_features.append(pp_feature)
        for seq_set in classifier.seqs:
            if seq_set not in seq_sets: seq_sets.append(seq_set)

    if regions == None:
        regions = regions_from_bed(region_bed)

    for idx, feature_file in enumerate(intersect_features):
        #consider_val = True if 'phastcons' in feature_file.lower() else False
        loc2features = run_intersectBed(region_bed, feature_file, 
                                        consider_value=consider_vals[idx])
        #loc2features = run_bedops_intersect(region_bed, feature_file, bed1_regions=regions)
        #loc2features = run_coverageBed(feature_file, region_bed)   
        feature = []
        for region in regions:
            if consider_vals[idx]:
                v = np.average(loc2features[region]) if len(loc2features[region]) > 0 else 0.0
                feature.append(v)
            else:
                feature.append(float(len(loc2features[region]) > 0))
                
#            feature.append(float(len(loc2features[region]) > 0))

            ##feature.append(loc2features[region])

        print '\t%s\t%d\t%g' % (feature_file.split('/')[-1], len(feature), sum(feature))
        #print '\t%s\t%d\t%s' % (feature_file.split('/')[-1], len(feature), feature)
        name2feats[feature_file] = feature

    for pp_feat_file in pp_features:
        pp_feature = feature_from_preprocessed(regions, pp_feat_file)

        print '\t%s\t%d\t%g' % (pp_feat_file.split('/')[-1], len(pp_feature), sum(pp_feature))
        name2pp_feats[pp_feat_file] = pp_feature

## THIS SEGMENT WAS AN ATTEMPT TO ALLOW MULTIPLE PP FEATURES IN A
## SINGLE FILE. It doesn't work because of the use of the file names
## in classifier.pp_features[kidx] to find the features indexed in
## name2ppfeats in evaluate_enhancer_prediction.       
#    for pp_feat_file in pp_features: 
#        pp_feature_names, pp_feature_data = features_from_preprocessed(regions, pp_feat_file)
#        print "\tread features from %s: %s" % (pp_feat_file.split('/')[-1],
#                                               '\t'.join(pp_feature_names))
#
#        for idx, pp_feature_name in enumerate(pp_feature_names):
#            print '\t%s\t%d\t%g' % (pp_feature_name, len(pp_feature_data[idx]),
#                                    sum(pp_feature_data[idx]))
#            name2pp_feats[pp_feature_name] = pp_feature_data[idx]

    for seq_set in seq_sets:
        if seq_set:
            names2seqs[str(seq_set)] = prepare_sequence_data(seq_set, regions)

    return name2feats, name2pp_feats, names2seqs

# FIGURE OUT HOW THESE ARE DIFFERENT AND IF THEY CAN BE COMBINED!!!
# build_features() is used in evaluate_enhancer_prediction_shoguns.py
# build_shogun_features_for_classifier() is used in predict_enhancers_shoguns.py

def build_shogun_features_for_classifier(classifier, region_bed):
    """ """
    name2feats = {}
    name2pp_feats = {}
    names2seqs = {}  # NOTE THAT SEQS MAPS MULTIPLE FILES, NOT NEC. ONE

    # extract feature filenames
    intersect_features = []; pp_features = []; seq_sets = []
    consider_vals = []  # for each intersect_feature tells whether to consider region values
    regions = regions_from_bed(region_bed)
    for feature_set in classifier.bin_features:
        for feature in feature_set:
            if feature not in intersect_features: 
                intersect_features.append(feature)
                consider_vals.append(False)
    for feature_set in classifier.con_features:
        for feature in feature_set:
            if feature not in intersect_features: 
                intersect_features.append(feature)
                consider_vals.append(True)
    for pp_feature_set in classifier.pp_features:
        for pp_feature in pp_feature_set:
            if pp_feature not in pp_features: pp_features.append(pp_feature)
    for seq_set in classifier.seqs:
        if seq_set not in seq_sets: seq_sets.append(seq_set)

    # now do the intersections
    for idx, feature_file in enumerate(intersect_features):
        #consider_val = True if 'phastcons' in feature_file.lower() else False
        loc2features = run_intersectBed(region_bed, feature_file, 
                                        consider_value=consider_vals[idx])
        feature = []
        for region in regions:
            if consider_vals[idx]:
                v = np.average(loc2features[region]) if len(loc2features[region]) > 0 else 0.0
                feature.append(v)
            else:
                feature.append(float(len(loc2features[region]) > 0))

            #if len(loc2features[region]) > 0 and 'sabrina' in region_bed:
            #    print '%s\t%d\t%d\t%s' % (region[0], region[1], region[2], feature_file)
                

        print '\t', feature_file.split('/')[-1], len(feature)
        name2feats[feature_file] = feature

    for pp_feat_file in pp_features:
        pp_feature = feature_from_preprocessed(regions, pp_feat_file)

        print '\t', pp_feat_file.split('/')[-1], len(pp_feature)
        name2pp_feats[pp_feat_file] = pp_feature

    for seq_set in seq_sets:
        if seq_set:
            names2seqs[str(seq_set)] = prepare_sequence_data(seq_set, regions)


    # now build shogun feature object
    shogun_features = None
    ind_features = []
    #if len(classifier.features) > 1:
    if len(classifier.bin_features) > 1 or len(classifier.con_features) > 1:
        shogun_features = CombinedFeatures()
        for i, kernels in enumerate(classifier.kernels):
            if classifier.seqs[i] != []:
                seq_feats = build_seq_features(names2seqs[str(classifier.seqs[i])], order=classifier.ks[i], rev_comp=classifier.rev_comps[i])
                shogun_features.append_feature_obj(seq_feats)
                ind_features.append(seq_feats)
            else:
                feats = [name2feats[x] for x in classifier.bin_features[i]]
                feats += [name2feats[x] for x in classifier.con_features[i]]
                pp_feats = [name2pp_feats[x] for x in classifier.pp_features[i]]
                feat_obj = RealFeatures(np.array(feats + pp_feats))
                shogun_features.append_feature_obj(feat_obj)
                ind_features.append(feat_obj)
    else:
        if classifier.seqs[0] != []:
            seq_feats = build_seq_features(names2seqs[str(classifier.seqs[0])], 
                                           order=classifier.ks[0], rev_comp=classifier.rev_comps[0])
            #seq_feats.save_serializable()
            #sys.exit()
            shogun_features = seq_feats
        else:
            #feats = [name2feats[x] for x in classifier.features[0]]
            feats = [name2feats[x] for x in classifier.bin_features[0]]
            feats += [name2feats[x] for x in classifier.con_features[0]]
            pp_feats = [name2pp_feats[x] for x in classifier.pp_features[0]]
            feat_obj = RealFeatures(np.array(feats + pp_feats))
            shogun_features = feat_obj
    
    return shogun_features, ind_features


def calc_within_kernel_weights(svm, classifier, feature_objs, raw_labels):
    """ """
    kern2feat_weights = {}
    alphas = svm.get_alphas()
    svs = svm.get_support_vectors()


    if len(classifier.kern_names) != len(feature_objs):
        sys.exit("ERROR: Number of feature objects (%d) does not match config file (%d)." % (len(classifer.kern_names), len(feature_objs)))

    for obj_idx, feat_obj in enumerate(feature_objs):
        weight_str = ''
        feat_name = feat_obj.get_name()
        kern_name = classifier.kern_names[obj_idx]

        feats = None
        if 'String' in feat_name:
            feats = feat_obj.get_features()

            k = classifier.ks[obj_idx]
            raw_weights = np.zeros(4 ** k)

            #for feat_idx, feat in enumerate(feats):
            for alpha_idx, sv_idx in enumerate(svs):
                alpha = alphas[alpha_idx]
                # The alphas have already been multiplied by y_i
                for kmer in feats[sv_idx]:
                    raw_weights[kmer] += alpha

            #norm_weights = raw_weights / np.sum(np.absolute(raw_weights))
            norm_weights = raw_weights

            eps = 'ACGT'
            kmers = [''.join(x) for x in list(itertools.product(eps, repeat=k))]

            weights_names = zip(norm_weights, kmers)

        else:
            feat_mat = feat_obj.get_feature_matrix()[:,svs]
            labs = np.array(raw_labels, dtype=np.double)[svs]

            # since alphas for positives are positive and negative for
            #negatives, I'm wondering if they've already been
            #multiplied by the y_i
            #al = alphas * labs
            prod = alphas * feat_mat

            raw_weights = prod.sum(axis=1).transpose()
            norm_weights = raw_weights / np.sum(np.absolute(raw_weights))            
            #feat_names = [x.split('/')[-1] for x in classifier.features[obj_idx] + classifier.pp_features[obj_idx]]
            feat_names = [x.split('/')[-1] for x in classifier.bin_features[obj_idx] + classifier.con_features[obj_idx] + classifier.pp_features[obj_idx]]

            weights_names = zip(norm_weights, feat_names)

        weights_names.sort(reverse=True)
        for x in weights_names:
            weight_str += '%g\t%s\n' % (x[0], x[1])
                
        kern2feat_weights[kern_name] = weight_str

    return kern2feat_weights


################################################################################

if __name__ == "__main__": 
    sys.exit("This is just a library!")

