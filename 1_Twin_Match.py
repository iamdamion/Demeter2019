#!/usr/bin/env python
# See help option for usage
#
# Originally written by Damion Demeter & Remy Mallett 03.10.18
#
# Twin Resting State MVPA Classifier Project - This script \
# classifies twins based on resting state network connectivity. Currently hardcoded \
# for use with the UT Redux 255 ROI set. SPHERE timecourse .txt files from the DCN \
# RS Pipeline tools are expected.' (see last revision's modified date)
from __future__ import division

last_modified = 'Last Modified by Damion Demeter, 4.11.19'

import argparse,glob,itertools,os,random,sys
import numpy as np
from scipy.stats.stats import pearsonr
from scipy import stats

from sklearn.svm import SVC
from sklearn.pipeline import Pipeline
from sklearn.feature_selection import SelectPercentile, f_classif
from sklearn.model_selection import LeaveOneGroupOut, permutation_test_score, cross_validate
from sklearn import metrics

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

here = os.path.dirname(os.path.abspath(os.path.realpath(__file__)))
############
prog_descrip = 'Twin Resting State SVM Classifier' + last_modified

def main(argv=sys.argv):
    arg_parser = argparse.ArgumentParser(description=prog_descrip,formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # Check for arguments. #
    if len(sys.argv[1:])==0:
        print '\nArguments required. Use -h option to print FULL usage.\n'

    arg_parser.add_argument('-noboot', action='store_true', required=False, default=False,
                            help=('The NO bootstrap flag will STOP bootstraping classifier (shuffling labels) for better average and SE/CI calculation.'),
                            dest='noboot'
                            )
    arg_parser.add_argument('-chunk', action='store', type=int, required=False, default='10',
                            help=('TR count for each timecourse "CHUNK". (Must be an even divisor for the CROP argument)'),
                            dest='chunk'
                            )
    arg_parser.add_argument('-crop', action='store', type=int, required=False, default='150',
                            help=('TOTAL TR count for cropping ALL timecourses for even "chunking". (Will take TRs PRIOR to this count)'),
                            dest='crop'
                            )
    arg_parser.add_argument('-fs', action='store', type=int, required=False, default=5,
                            help=('Feature Slection Percentage. Top X percent features in ANOVA method.'),
                            dest='fs'
                            )
    arg_parser.add_argument('-nw', action='store', nargs='+', type=str, required=False, default=['FULL'],
                            help=('Valid network choices: UNKN, MEM, SSM_H, SSM_M, CO, AUD, DMN, VIS, FP, SAL, SUBCORT, VA, DA. (Default = FULL)'),
                            dest='network_list'
                            )
    arg_parser.add_argument('-o', metavar='OUT_PATH', action='store', required=True,
                            help=('Path to output directory. All results/pfds will be saved here.'),
                            dest='out_path'
                            )
    arg_parser.add_argument('-plot', action='store_true', required=False, default=False,
                            help=('Create plots (saved to output directory).'),
                            dest='plot'
                            )
    arg_parser.add_argument('-t', metavar='TC_DIR', action='store', required=True,
                            help=('Path to TIMECOURSE directory. Files named "LABEL_GROUP.txt". Example: "A_1.txt","A_2.txt", etc.'),
                            dest='tc_dir'
                            )
    args = arg_parser.parse_args()
    #################################################
    ## Script Argument Verification and Assignment ##
    #################################################
    print '\n--------------------- setup info ---------------------------------'
    ### CHECK CROP AND CHUNK SIZE DIVISIBLE ###
    if args.crop % args.chunk == 0:
        crop = args.crop
        chunk = args.chunk
    else:
        print 'Crop size not evenly divisible by chunk size. Check and re-run. Exiting...'
        sys.exit()
    ### VERIFY CHOSEN NETWORKS ###
    NETWORKS = dict(
        UNKN    = (range(0,24),'grey'),
        MEM     = (range(24,29),'grey'),
        SSM_H   = (range(29,59),'cyan'),
        SSM_M   = (range(59,64),'orange'),
        CO      = (range(64,78),'purple'),
        AUD     = (range(78,91),'pink'),
        DMN     = (range(91,148),'red'),
        VIS     = (range(148,179),'blue'),
        FP      = (range(179,204),'yellow'),
        SAL     = (range(204,222),'black'),
        SUBCORT = (range(222,235),'brown'),
        VA      = (range(235,244),'teal'),
        DA      = (range(244,255),'green'),
        FULL    = (range(0,255),'aliceblue'),
    )
    if 'ALL' in args.network_list:
        network_list = NETWORKS.keys()
    else: 
        for n in args.network_list:
            if n in NETWORKS.keys():
                network_list = args.network_list
            else:
                print 'The following chosen network is not in the available network list: ' + n
                print 'Check and re-run. Exiting...'
                sys.exit()
    ### VERIFY OUTPUT PATH EXISTS ###
    if os.path.isdir(args.out_path):
        out_path = args.out_path
        print 'All output files will be saved here:'
        print out_path + '\n'
    else:
        print 'Output directory could NOT be verified and is required. Check and re-run. Exiting...'
        sys.exit()
    ### VERIFY TIMECOURSE DIRECTORY PATH EXISTS ###
    if os.path.isdir(args.tc_dir):
        tc_dir = args.tc_dir
        print 'Using timecourse .txt files from:'
        print tc_dir
        print 'NOTE: File naming will NOT be verified. Be sure all files are correctly named.'
    else:
        print 'Timecourse directory could NOT be verified and is required. Check and re-run. Exiting...'
        sys.exit()
    print '--------------------------- end ---------------------------------\n'
    #################################################
    ##          Global Variable Assignment         ##
    #################################################
    sub_glob = sorted(glob.glob(os.path.join(args.tc_dir,'*.txt')))
    ### Make fs_perc variable to help with looping ###
    fs_perc = args.fs
    #################################################
    ##                  FUNCTIONS                  ##
    #################################################
    def natural_sort(l): 
        convert = lambda text: int(text) if text.isdigit() else text.lower() 
        alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
        return sorted(l, key = alphanum_key)

    def create_arrays(nw,crop,chunk):
        ## build empty arrays to hold data, labels, and group ids ##
        # each row is a vector of corr vals within a mini chunk for one sub
        n_subs = len(sub_glob)
        n_chunks = int(crop / chunk)
        n_rows = n_subs*n_chunks
        # each col is a corr val from the lower triangle of corr matrix
        n_rois = len(NETWORKS[nw][0])
        n_cols = int(n_rois*(n_rois-1) / 2)
        data = np.empty([n_rows,n_cols])
        labels = np.empty(n_rows, dtype='S10')
        groups = np.empty(n_rows, dtype='S10')
        ## iterate over subjects ##
        i = 0 # counter for inserting vectors in final data array

        for sub in sub_glob:
            # get information about subject, for sample and group labels
            # label = family | group = sib number
            label = os.path.splitext(os.path.basename(sub))[0].split('_')[0]
            group = os.path.splitext(os.path.basename(sub))[0].split('_')[1]
            # import subject timecourse (ROIs as rows, time/tr as columns)
            sub_tc = np.loadtxt(sub)
            ## CHECK FOR ZEROS AND NaNs in array ##
            if np.isnan(np.sum(sub_tc)) == False:
                pass
            else:
                print 'FOUND NANS!! In: ' + sub
                print np.isnan(np.sum(sub_tc))
                print np.isnan(sub_tc)
                print np.isinf(sub_tc)
                sys.exit()
            # extract certain network from rows,
            # and extract a subset of total timecourse from columns
            sub_network_tc = sub_tc[NETWORKS[nw][0],:crop]
            # extract "mini" timecourses from whole timecourse
            sub_minitc = np.split(sub_network_tc, n_chunks, axis=1)

            # sub_minitc is a LIST of arrays, where each list is a timecourse.
            # for each mini-timecourse, make a vector from correlation matrix
            for idx,mtc in enumerate(sub_minitc):
                # create correlation matrix for network ROIs
                sub_minicorr = np.corrcoef(mtc)
                # get upper triangle indices, WITHOUT center 1's
                indices = np.triu_indices_from(sub_minicorr,k=1)
                sub_corr_vector = sub_minicorr[indices]
                # insert data and info into corresponding arrays
                data[i,:] = sub_corr_vector
                labels[i] = label
                groups[i] = group
                # move counter forward
                i += 1

        return data, labels, groups, indices, n_rois

    def get_feat_mask(data, labels, groups, indices):
        ## BUILD CLASSIFIER ##
        # choose classifier and feature selection methods - NOTE: Using FS PIPELINE HERE
        # DEFAULT = USE ANOVA @ TOP 5% MOST "IMPORTANT" (VARIABLE) ROIs TO DEFINE FEATURES
        clf = Pipeline([
          ('feature_selection', SelectPercentile(f_classif,percentile=fs_perc)),
          ('classification', SVC(kernel='linear'))])
        top_corrs = {}
        f_of_tops = {}

        ## GRABBING FEATURES MASK ##
        for j, (train_indx, test_indx) in enumerate(LeaveOneGroupOut().split(data,labels,groups)):
            ## FIRST GET TRAINING SET FIT ##
            X, Y = data[train_indx,:], labels[train_indx]
            x, y = data[test_indx,:], labels[test_indx]
            clf.fit(X,Y)
            # get boolean vector of which features were selected
            selected_feats = clf.named_steps['feature_selection'].get_support()
            # each feature is a correlation value between two ROIs, so get those two rois for each selected feature
            feat_indices = np.argwhere(selected_feats==True).flatten()
            vals_in_triu = np.array([ [indices[0][r],indices[1][r]] for r in feat_indices ])
            top_corrs[j] = vals_in_triu
            # get the F values
            f_of_tops[j] = clf.named_steps['feature_selection'].scores_[selected_feats]
            ## THIS IS THE FEATURE MASK, LIST OF TRUE/FALSE
            feats_mask = selected_feats

        return feats_mask

    def apply_feature_mask(dataset,feature_mask):
        feats_mask = feature_mask.tolist()
        feat_count = feats_mask.count(True)
        feat_masked_data = np.zeros(shape=(dataset.shape[0],feat_count))
        for idx,row in enumerate(dataset):
            feat_masked_data[idx] = list(itertools.compress(row, feats_mask))

        return feat_masked_data

    def run_classifier(data, labels, groups):
        ## RUN CLASSIFIER ##
        # NOTE: Using *NO FS* PIPELINE HERE
        # THIS IS BECAUSE FEATURE MASK IS ALREADY APPLIED!!
        clf = Pipeline([
          ('classification', SVC(kernel='linear'))])
        # get classifier performance, and test significance with a permutation test
        score, perms, pval = permutation_test_score(
            clf, data, labels,
            groups=groups,
            cv=LeaveOneGroupOut(),
            n_permutations=1000,
        )
        # get each fold's score
        fold_scores = cross_validate(
            clf, data, labels,
            groups=groups,
            cv=LeaveOneGroupOut(),
        )

        return score, perms, pval, fold_scores

    def run_bootstrap_classifier(data, labels, groups):
        # FIRST MAKE LIST OF GROUP LABELS TO ITERATE OVER #
        n = int(groups.shape[0]/2)
        print 'Group size double check: ' + str(n)
        if n < 10:
            iter_count = 500
        else:
            iter_count = 1000
        print 'Running ' + str(iter_count) + ' bootstrap iterations of label shuffles.'
        print 'This is NOT the same as permutations. Only shuffling the train/test label between the pairs, not between individuals/families...'
        group_iterations = []
        choices = ['12','21']

        while len(group_iterations) < iter_count:
            rand = np.random.choice(choices,n)
            iter_lst = rand.tolist()
            iter_string = ''.join(iter_lst)
            iter_str_lst = [str(i) for i in iter_string]
            if not iter_str_lst in group_iterations:
                group_iterations.append(iter_str_lst)
            else:
                pass
        ## NOW WITH NON-REPEATED BOOTSTRAP GROUP LISTS MADE, RUN CLASSIFIER AND SAVE FOLD SCORES TO CALCULATE AVERAGE AND SE ##
        scores = []
        ## PERM SCORES FOR CALCULATING P VALUES!!!
        perms_for_p = []
        frozen_count = 1
        for boot in group_iterations:
            strap = np.array(boot)
            print 'iteration: ' + str(frozen_count)
            ## RUN CLASSIFIER ##
            # NOTE: Using *NO FS* PIPELINE HERE
            # THIS IS BECAUSE FEATURE MASK IS ALREADY APPLIED!!
            clf = Pipeline([
              ('classification', SVC(kernel='linear'))])
            # get classifier performance, and test significance with a permutation test
            score, perms, pval = permutation_test_score(
                clf, data, labels,
                groups=groups,
                cv=LeaveOneGroupOut(),
                n_permutations=1000,
            )
            ## ADDING PERM SCORES LIST TO FULL LIST TO CALCULATE P VAL LATER ##
            perms_for_p.append(perms)
            # get each fold's score
            fold_scores = cross_validate(
                clf, data, labels,
                groups=strap,
                cv=LeaveOneGroupOut(),
            )
            scores.append(fold_scores['test_score'][0])
            scores.append(fold_scores['test_score'][1])
            frozen_count = frozen_count + 1

        mean_score = np.mean(scores)
        se = stats.sem(scores)
        sd = np.std(scores)
        theoretical_chance = 1/n
        conf_int = stats.norm.interval(0.95, loc=mean_score, scale=sd / np.sqrt(len(scores)))

        return mean_score, se, sd, theoretical_chance, conf_int, scores, perms_for_p
    #################################################
    ##               MAIN SCRIPT ENTRY             ##
    #################################################
    ### ITERATE OVER NETWORKS ###
    group_name = os.path.basename(os.path.abspath(args.tc_dir))
    print 'Working on: ' + group_name
    print 'Looping over selected networks:'
    for idx,nw in enumerate(sorted(network_list)):
        print nw + ' (' + str((idx+1)) + '/' + str(len(network_list)) + ')'
        ## PRE-I. ADJUST DEFAULT CROP FOR HCP (TO AVOID FORGETTING)
        if 'HCP' in group_name:
            print 'Using HCP chunk/crop values...'
            crop = 420
            chunk = 28
        else:
            crop = args.crop
            chunk = args.chunk
        print 'Crop: ' + str(crop)
        print 'Chunk: ' + str(chunk)

        ## I. CREATE ORGANIZED DATAFRAMES FOR ALL SUBJECTS FOR THIS NETWORK ##
        print 'Creating CHUNKED arrays for FEATURE SELECTION...'
        data, labels, groups, indices, n_rois = create_arrays(nw,crop,chunk)
        ## II. RUN CLASSIFIER TO GET FEATURES MASK ONLY ##
        print 'Creating features mask...'
        feats_mask = get_feat_mask(data, labels, groups, indices)
        ## III. CREATE NON-CHUNKED DATAFRAME ##
        # (for the function, feeding the crop as both crop and chunk. This is for HCP/UT flexibility)
        print 'Creating NON-CHUNKED arrays for SVM...'
        data_nochunk, labels_nochunk, groups_nochunk, indices_nochunk, n_rois_nochunk = create_arrays(nw,crop,crop)
        ## IV. APPLY FEATURES MASK TO CORRELATION VECTORS IN DATA ##
        print 'Applying features mask to non-chunked data...'
        feat_masked_data = apply_feature_mask(data_nochunk,feats_mask)
        if args.noboot == True:
            ## V. RUN SVM WITH MASKED AND FULL DATASET!!
            print 'Running classifier WITHOUT label shuffle iterations...'
            print '\n----------\nMatrix Verify: ' + str(len(feat_masked_data[0]))
            score, perms, pval, fold_scores = run_classifier(feat_masked_data, labels_nochunk, groups_nochunk)
            print 'SVM (FULL Masked) SCORE: ' + str(float(score)) + ' P = ' + str(pval)
            print 'FOLD INFO:'
            print fold_scores['test_score']

        elif args.noboot == False:
            ## VI. RUN SVM WITH MASKED AND FULL DATASET - AND WITH BOOTSTRAPPED GROUP LIST!!!!!
            print 'Running TWIN-SHUFFLE version of classifier. Default 1K iterations...'
            print '\n----------\nMatrix Verify: ' + str(len(feat_masked_data[0]))
            mean_score, se, sd, theoretical_chance, conf_int, scores, perms_for_p = run_bootstrap_classifier(feat_masked_data, labels_nochunk, groups_nochunk)
            ## VIb. CALCULATE P VAL BASED ON *ALL* TWIN-SHUFFLE CHANCE PERMUTATIONS ##
            chance_scores = np.concatenate(perms_for_p)
            score = mean_score
            C = np.sum(chance_scores>score)
            N_PERMS = chance_scores.size
            pval = (C+1) / float(N_PERMS+1)
            empirical_chance = np.mean(perms_for_p)
            ## REPORT SCORES, ETC ##
            print '\nTwin-Shuffle group MEAN prediction score: ' + str(mean_score)
            print '95% Confidence Interval: ' + str(conf_int)
            print 'Empirical CHANCE: ' + str(empirical_chance)
            print 'THEORETICAL Chance: ' + str(theoretical_chance)
            print 'p = {}'.format(pval)

if __name__ == '__main__':
    sys.exit(main())