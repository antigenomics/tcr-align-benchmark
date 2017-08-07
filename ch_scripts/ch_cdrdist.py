#code used here was obtained from https://github.com/phbradley/tcr-dist and modified

import os
import pandas as pd
import numpy as np
import argparse
import ch_antigen_base as ch_w
import copy
import ch_read_blosum as ch_rbl
from Bio import pairwise2

gap_character = '.'

amino_acids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
               'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']


def compute_cdrmatrix(imatrix, lim=4):
    cdrdistmatrix = {}
    for a in amino_acids:
        for b in amino_acids:
            bab = importmatrix['matrix'][a][b]
            if a == b:
                cdrdistmatrix[(a, b)] = 0
            else:  ## different
                if bab < 0:
                    cdrdistmatrix[(a, b)] = -lim
                else:
                    cdrdistmatrix[(a, b)] = -lim + bab
    return (cdrdistmatrix)

def getmax(matrix):
    maxmatrix = copy.deepcopy(matrix)
    np.fill_diagonal(maxmatrix.values, 0)
    return(maxmatrix.values.max())

importmatrix = ch_rbl.read_blosum('../blosum_matrices/BLOSUM65')
maxmatrix = copy.deepcopy(importmatrix['matrix'])
np.fill_diagonal(maxmatrix.values, 0)
cdrdistmatrix = compute_cdrmatrix(importmatrix, maxmatrix.values.max())


class DistanceParams:
    def __init__(self, distmatrix=cdrdistmatrix, config_string=None):
        self.gap_penalty_v_region = 4
        self.gap_penalty_cdr3_region = 12  # same as gap_penalty_v_region=4 since weight_cdr3_region=3 is not applied
        self.weight_v_region = 1
        self.weight_cdr3_region = 3
        self.distance_matrix = cdrdistmatrix
        self.align_cdr3s = True
        self.trim_cdr3s = True
        self.scale_factor = 1.0
        if config_string:
            l = config_string.split(',')
            for tag, val in [x.split(':') for x in l]:
                if tag == 'gap_penalty_cdr3_region':
                    self.gap_penalty_cdr3_region = float(val)
                elif tag == 'gap_penalty_v_region':
                    self.gap_penalty_v_region = float(val)
                elif tag == 'weight_cdr3_region':
                    self.weight_cdr3_region = float(val)
                elif tag == 'weight_v_region':
                    self.weight_v_region = float(val)
                elif tag == 'scale_factor':
                    self.scale_factor = float(val)
                elif tag == 'align_cdr3s':
                    assert val in ['True', 'False']
                    self.align_cdr3s = (val == 'True')
                elif tag == 'trim_cdr3s':
                    assert val in ['True', 'False']
                    self.trim_cdr3s = (val == 'True')
                else:
                    print('unrecognized tag:', tag)
                    assert False
            print('config_string: {} self: {}'.format(config_string, self))

    def __str__(self):
        return 'DistanceParams: gap_penalty_v_region= {} gap_penalty_cdr3_region= {} weight_v_region= {} weight_cdr3_region= {} align_cdr3s= {} trim_cdr3s= {}' \
            .format(self.gap_penalty_v_region, self.gap_penalty_cdr3_region,
                    self.weight_v_region, self.weight_cdr3_region,
                    self.align_cdr3s, self.trim_cdr3s)


default_distance_params = DistanceParams()


def blosum_character_distance(a, b, gap_penalty, params):
    if a == gap_character and b == gap_character:
        return 0
    elif a == '*' and b == '*':
        return 0
    elif a == gap_character or b == gap_character or a == '*' or b == '*':
        return gap_penalty
    else:
        # assert a in amino_acids
        # assert b in amino_acids
        # maxval = min( blosum[(a,a)], blosum[(b,b)] )
        # return maxval - blosum[(a,b)]
        return -params.distance_matrix[(a, b)]


def blosum_sequence_distance(aseq, bseq, gap_penalty, params):
    assert len(aseq) == len(bseq)
    dist = 0.0
    for a, b in zip(aseq, bseq):
        if a == ' ' and b == ' ':
            continue
        dist += blosum_character_distance(a, b, gap_penalty, params)
    return dist


def align_cdr3s( a, b, gap_character, params ):
    ntrim = 3 if params.trim_cdr3s else 0
    ctrim = 2 if params.trim_cdr3s else 0
    if len(a) < 6 or len(b) < 6:
        alignments = pairwise2.align.globalds(a, b, cdrdistmatrix, -12, -12)
        s0 = alignments[0][0].replace('-', gap_character)
        s1 = alignments[0][1].replace('-', gap_character)
        assert len(s0) == len(s1)
        dist = blosum_sequence_distance(s0, s1, 0, params)
        return {'alignment': [a[0:ntrim] + s0 + a[-ctrim:], b[0:ntrim] + s1 + a[-ctrim:]], 'distance': dist}

    s0, s1 = a[ntrim:-ctrim], b[ntrim:-ctrim]

    if len(a) == len(b):
        dist = blosum_sequence_distance(s0, s1, 0, params)
        return {'alignment': [a[0:ntrim] + s0 + a[-ctrim:], b[0:ntrim] + s1 + a[-ctrim:]], 'distance': dist}

    alignments = pairwise2.align.globalds(s0, s1, cdrdistmatrix, -12, -12)
    s0 = alignments[0][0].replace('-', gap_character)
    s1 = alignments[0][1].replace('-', gap_character)

    assert len(s0) == len(s1)
    dist = blosum_sequence_distance(s0, s1, 0, params)
    return {'alignment':[a[0:ntrim]+s0+a[-ctrim:], b[0:ntrim]+s1+a[-ctrim:]], 'distance':dist}


def sequence_distance_with_gappos( shortseq, longseq, gappos, params ):
    ntrim = 3 if params.trim_cdr3s else 0
    ctrim = 2 if params.trim_cdr3s else 0
    remainder = len(shortseq)-gappos
    dist = 0.0
    count =0
    if ntrim < gappos:
        for i in range(ntrim,gappos):
            dist += -params.distance_matrix[ (shortseq[i], longseq[i] ) ]
            count += 1
    if ctrim < remainder:
        for i in range(ctrim, remainder):
            dist += -params.distance_matrix[ (shortseq[-1-i], longseq[-1-i] ) ]
            count += 1
    return dist,count

def weighted_cdr3_distance( seq1, seq2, params ):
    shortseq,longseq = (seq1,seq2) if len(seq1)<=len(seq2) else (seq2,seq1)

    ## try different positions of the gap
    #print(len(shortseq), len(longseq))
    lenshort = len(shortseq)
    lenlong = len(longseq)
    lendiff = lenlong - lenshort
#
    #assert lenshort > 1##JCC testing
    assert lendiff>=0
    #if params.trim_cdr3s:
    #    assert lenshort > 3+2 ## something to align...

    if not params.align_cdr3s:
        ## if we are not aligning, use a fixed gap position relative to the start of the CDR3
        ## that reflects the typically longer and more variable-length contributions to
        ## the CDR3 from the J than from the V. For a normal-length
        ## CDR3 this would be after the Cys+5 position (ie, gappos = 6; align 6 rsds on N-terminal side of CDR3).
        ## Use an earlier gappos if lenshort is less than 11.
        ##

        gappos = min( 6, 3 + (lenshort-5)/2 )
        best_dist, count = sequence_distance_with_gappos( shortseq, longseq, gappos, params )

    else:
        ## the CYS and the first G of the GXG are 'aligned' in the beta sheet
        ## the alignment seems to continue through roughly CYS+4
        ## ie it's hard to see how we could have an 'insertion' within that region
        ## gappos=1 would be a insertion after CYS
        ## gappos=5 would be a insertion after CYS+4 (5 rsds before the gap)
        ## the full cdr3 ends at the position before the first G
        ## so gappos of len(shortseq)-1 would be gap right before the 'G'
        ## shifting this back by 4 would be analogous to what we do on the other strand, ie len(shortseq)-1-4

        info = align_cdr3s(seq1, seq2, gap_character, params)
        lendiff = info['alignment'][0].count('.') + info['alignment'][1].count('.')
        best_dist = int(info['distance'])

    ## Note that weight_cdr3_region is not applied to the gap penalty
    ##
    return int(params.weight_cdr3_region) * best_dist + lendiff * int(params.gap_penalty_cdr3_region)

def compute_v_region_distance( v1, v2, params ):
    dist = params.weight_v_region * blosum_sequence_distance(v1, v2, params.gap_penalty_v_region, params)
    return dist

def compute_all_v_region_distances( organism, params ):
    #cdr3s_human need to be implemented
    pass

def compute_distance(t1,t2,chains,distance_params): #    t1/2 = [ va_reps, vb_reps, l['cdr3a'], l['cdr3b'] ]

    dist=0.0
    if 'A' in chains:
        dist += compute_v_region_distance(t1[0], t2[0], distance_params) +\
                weighted_cdr3_distance( t1[2], t2[2], distance_params )
    if 'B' in chains:
        dist += compute_v_region_distance(t1[1], t2[1], distance_params) +\
                weighted_cdr3_distance( t1[3], t2[3], distance_params )
    return distance_params.scale_factor * dist

def compute_auc( l0, l1, sign_factor=1 ):
    ## l0 are the true positives, l1 are the false positives
    ## if sign_factor==1 then lower scores are better, otherwise it's the opposite
    ##
    if not l0:
        return 0.0, [0,1], [0,0]
    elif not l1:
        return 1.0, [0,0,1], [0,1,1]

    l = [ (sign_factor*x,0) for x in l0 ] + [ (sign_factor*x,-1) for x in l1 ] ## in ties, take the false positive first
    l.sort()

    xvals = []
    yvals = []

    counts = [0,0]
    totals = [len(l0),len(l1)]

    area=0.0
    width = 1.0/totals[1]
    for ( score, neg_tcr_class ) in l:
        tcr_class = -1*neg_tcr_class
        counts[ tcr_class ] += 1
        xval = float( counts[1] ) / totals[1]
        yval = float( counts[0] ) / totals[0]
        xvals.append( xval )
        yvals.append( yval )
        if tcr_class==1: area += yval * width

    return area,xvals,yvals

#9
def get_rank( val, l ): ## does not require that the list l is sorted
    num_lower = 0
    num_upper = 0

    epsilon = 1e-6

    lower_neighbor = val-10000
    upper_neighbor = val+10000

    #there we can see how many lower and upper neighbours of fixed values are here
    for x in l: #1e-6 is enough to think that values are near to be equal
        if x<val-epsilon:
            num_lower += 1
            lower_neighbor = max( lower_neighbor, x )
        elif x>val+epsilon:
            num_upper += 1
            upper_neighbor = min( upper_neighbor, x )

    total = len(l)
    num_equal = total - num_lower - num_upper #0 or more, if there are x in val +- epsilon
    assert num_equal >=0

    if num_upper == 0:
        return 100.0

    elif num_lower == 0:
        return 0.0

    else:
        assert upper_neighbor>lower_neighbor
        interp = (val-lower_neighbor)/(upper_neighbor-lower_neighbor) #
        #lowest upper_neighbot - lowerneighbor is equal to 2e-6, interp is 0.5
        #if num_equal>0:print 'num_equal:',num_equal

        interp_num_lower = num_lower + interp * ( 1 + num_equal )

        return (100.0*interp_num_lower)/total


#10
## negative nbrdist_percentile means take exactly -nbrdist_percentile topn
def sort_and_compute_nbrdist_from_distances( l, nbrdist_percentile, dont_sort=False ):
    if not dont_sort: l.sort()
    #assert l[0]<=l[-1]
    if nbrdist_percentile<0:
        n = max( 1, min(len(l), -1*nbrdist_percentile ) ) # negative nbrdist_percentile means take exactly -nbrdist_percentile topn
    else:
        n = max(1, ( nbrdist_percentile * len(l) )//100 )
    #print((l))
    return sum( l[0:n])/float(n)

#11
## negative nbrdist_percentile means take exactly -nbrdist_percentile topn
def sort_and_compute_weighted_nbrdist_from_distances( l, nbrdist_percentile, dont_sort=False ):
    if not dont_sort: l.sort()
    #assert l[0]<=l[-1]
    if nbrdist_percentile<0:
        n = max( 1, min(len(l), -1*nbrdist_percentile ) ) # negative nbrdist_percentile means take exactly -nbrdist_percentile topn
    else:
        n = max(1, ( nbrdist_percentile * len(l) )//100 )

    total_wt = 0.0
    nbrdist=0.0
    for i,val in enumerate( l[:n] ):
        wt = 1.0 - float(i)/n
        total_wt += wt
        nbrdist += wt * val #val is mostly distance

    return nbrdist / total_wt
'''
t1 = ['ENSA......SN DIRSN...ME LDKKA.', 'MNH.......DT YYD....KIL P.NNSF',
      'CAAAAGNTGKLIF', 'CASSRLGASAETLYF']
t2 = ['ENSA......SN DIRSN...ME LDKKA.', 'MNH.......DT YYD....KIL P.NNSF',
      'CAAAGGNNKLTF', 'CASSSLGSSAETLYF']

print(compute_distance(t1, t2, 'AB', default_distance_params))
'''

