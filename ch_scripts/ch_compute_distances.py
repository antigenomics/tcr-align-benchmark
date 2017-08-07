"""
This file was exported from tcr-dist package with some modifications

MIT License

Copyright (c) 2017 Philip Harlan Bradley and Jeremy Chase Crawford

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""

## this script will compute distances between all the tcrs and write out a distance matrix
## the order of rows and columns in the distance matrix will match that in the input file
## it will also compute the rank scores for each tcr with respect to all the epitopes present
##

import ch_antigen_base as ch_w
import ch_cdrdist
import pandas as pd
import numpy as np
import csv
from multiprocessing import Pool, Manager

ooo = True
clones_file = '../imgt_work/vdjdb_ab.txt'
nbrdist_percentiles = [5, 10, 25]

#internal legacy hack
intrasubject_nbrdists = 1
new_nbrdists = not intrasubject_nbrdists

importmatrix = ch_cdrdist.ch_rbl.read_blosum('../blosum_matrices/BLOSUM65')
cdrdistmatrixn = ch_cdrdist.compute_cdrmatrix(importmatrix, ch_cdrdist.getmax(importmatrix['matrix']))
distance_params = ch_cdrdist.DistanceParams(distmatrix=cdrdistmatrixn)

clones_files = [clones_file]

epitope_prefixes = ['']*len(clones_files)

assert len(epitope_prefixes) == len(clones_files)

#print 'precomputing v-region distances'
#rep_dists = tcr_distances.compute_all_v_region_distances( organism, distance_params )
#print 'done precomputing v-region distances'

distfile_prefix = clones_files[0][:-len('_imgted.txt')]

outfile = '{}_nbrdists.tsv'.format(clones_files[0][:-len('_imgted.txt')])


all_tcrs = []
all_info = []

all_infields = []

for clones_file, epitope_prefix in zip(clones_files,epitope_prefixes):
    infields = []
    with open( clones_file,'r') as inp:
        reader = csv.reader(inp, delimiter='\t')
        header = next(reader)
        for row in reader:
            if not infields:
                if header[0] == '#':
                    infields = header[1:]
                else:
                    infields = header[:]
                all_infields.append( infields )
                continue
            assert infields

            l = ch_w.parse_tsv_line( row[:], infields )

            assert 'epitope' in l
            l['epitope'] = epitope_prefix + l['epitope']

            all_info.append( l )
            va_seq = l['cdr1.alpha']+' '+l['cdr2.alpha']+' '+l['cdr2.5.alpha']
            vb_seq = l['cdr1.beta'] + ' ' + l['cdr2.beta'] + ' ' + l['cdr2.5.beta']
            tcr = [ va_seq, vb_seq, l['cdr3.alpha'], l['cdr3.beta'] ]

            all_tcrs.append( tcr )


## we will add new fields here

outfields = []
for f in all_infields[0]:
    common = True
    for infields in all_infields[1:]:
        if f not in infields:
            common = False
            break
    if common:
        outfields.append(f)


def work_with_tcr(tcrinp):
    ii = tcrinp[0]
    t1 = tcrinp[1]
    if not ii % 100: print('computing {} distances {} {}'.format(chains, ii, total_tcrs))
    epitope_distances = {}
    for e in epitopes:
        epitope_distances[e] = []

    one_info = all_info[ii]
    e1 = all_info[ii]['epitope']
    m1 = all_info[ii]['subject.id']
    clone_id1 = all_info[ii]['clone.id']
    dists = []
    for jj, t2 in enumerate(all_tcrs):
        dist = ch_cdrdist.compute_distance(t1, t2, chains, distance_params)
        dists.append(dist)

        e2 = all_info[jj]['epitope']
        m2 = all_info[jj]['subject.id']
        if ii != jj:  ## dont include self distance in the nbrdists, rank scores
            if new_nbrdists:
                if m1 != m2:
                    epitope_distances[e2].append(dist)
            else:
                epitope_distances[e2].append(dist)
                # if e2==e1:
                #   my_epitope_dists.append(dist)

                ## write out a new line in the distance file
    assert len(dists) == len(all_clone_ids)

    ## get an nbrdist
    for e, edists in iter(epitope_distances.items()):
        edists.sort()
        for nbrdist_percentile in nbrdist_percentiles:
            if edists:
                nbrdist = ch_cdrdist.sort_and_compute_nbrdist_from_distances(edists, nbrdist_percentile,
                                                                             dont_sort=True)
                wtd_nbrdist = ch_cdrdist.sort_and_compute_weighted_nbrdist_from_distances(edists,
                                                                                          nbrdist_percentile,
                                                                                          dont_sort=True)
            else:
                nbrdist = 0.0
                wtd_nbrdist = 0.0
            one_info['{}_{}_nbrdist{}'.format(e, chains, nbrdist_percentile)] = nbrdist
            one_info['{}_{}_wtd_nbrdist{}'.format(e, chains, nbrdist_percentile)] = wtd_nbrdist
            if e == e1:
                epitope_self_nbrdists[''][e][nbrdist_percentile].append(nbrdist)
                epitope_self_nbrdists['wtd_'][e][nbrdist_percentile].append(wtd_nbrdist)
    return(one_info)

for chains in ['A','B','AB']:
    epitopes = list(set([x['epitope'] for x in all_info ] ) )
    epitopes.sort()

    ## for computing rank scores
    with Manager() as manager:
        epitope_self_nbrdists = {}
        for wtd in ['','wtd_']:
            epitope_self_nbrdists[wtd] = {}
            for e in epitopes:
                epitope_self_nbrdists[wtd][e] = {}
                for p in nbrdist_percentiles:
                    epitope_self_nbrdists[wtd][e][p] = manager.list()

        all_clone_ids = [x['clone.id'] for x in all_info]

        ## make a distances file
        #out.write('\t'.join(['clone_id']+all_clone_ids)+'\n')

        total_tcrs = len(all_tcrs)
        if ooo==True:
            pool = Pool(20)
            all_infos = pool.map(work_with_tcr, enumerate(all_tcrs))
            all_info = list(all_infos)
            pool.close()
            pool.join()
        else:
            for ii,t1 in enumerate(all_tcrs):
                if not ii%100: print('computing {} distances {} {}'.format(chains,ii,total_tcrs))
                epitope_distances = {}
                for e in epitopes:
                    epitope_distances[e] = []

                e1 = all_info[ii]['epitope']
                m1 = all_info[ii]['subject.id']
                clone_id1 = all_info[ii]['clone.id']

                #dists=['r']*len(all_tcrs)
                #my_epitope_dists = []
                dists = []
                for jj,t2 in enumerate(all_tcrs):
                    dist = ch_cdrdist.compute_distance( t1,t2,chains,distance_params)
                    dists[jj] = dist

                    e2 = all_info[jj]['epitope']
                    m2 = all_info[jj]['subject.id']
                    if ii!=jj: ## dont include self distance in the nbrdists, rank scores
                        if new_nbrdists:
                            if m1 != m2:
                                epitope_distances[e2].append(dist)
                        else:
                            epitope_distances[e2].append(dist)
                    #if e2==e1:
                    #   my_epitope_dists.append(dist)

                    ## write out a new line in the distance file
                assert len(dists) == len(all_clone_ids)

                ## get an nbrdist
                for e,edists in iter(epitope_distances.items()):
                    edists.sort()
                    for nbrdist_percentile in nbrdist_percentiles:
                        if edists:
                            nbrdist = ch_cdrdist.sort_and_compute_nbrdist_from_distances( edists, nbrdist_percentile, dont_sort=True )
                            wtd_nbrdist = ch_cdrdist.sort_and_compute_weighted_nbrdist_from_distances( edists, nbrdist_percentile, dont_sort=True)
                        else:
                            nbrdist = 0.0
                            wtd_nbrdist = 0.0
                        all_info[ii]['{}_{}_nbrdist{}'.format(e,chains,nbrdist_percentile)] = nbrdist
                        all_info[ii]['{}_{}_wtd_nbrdist{}'.format(e,chains,nbrdist_percentile)] = wtd_nbrdist
                        if e == e1:
                            epitope_self_nbrdists[''][e][nbrdist_percentile].append( nbrdist )
                            epitope_self_nbrdists['wtd_'][e][nbrdist_percentile].append( wtd_nbrdist )

        print('record new fields in outfields')
        ## record new fields in outfields
        for wtd in ['','wtd_']:
            for nbrdist_percentile in nbrdist_percentiles:
                for suffix in ['_{}_{}nbrdist{}'.format(chains,wtd,nbrdist_percentile),
                               '_{}_{}nbrdist{}rank'.format(chains,wtd,nbrdist_percentile) ]:
                    for epitope in epitopes:
                        outfields.append( epitope+suffix)

        for wtd in ['','wtd_']:
            for e in epitopes:
                for p in nbrdist_percentiles:
                    epitope_self_nbrdists[wtd][e][p] = list(epitope_self_nbrdists[wtd][e][p])
                    epitope_self_nbrdists[wtd][e][p].sort()

    def multiprocrank(ii):
        one_info = all_info[ii]
        for epitope in epitopes:
            for wtd in ['','wtd_']:
                for nbrdist_percentile in nbrdist_percentiles:
                    self_nbrdists = epitope_self_nbrdists[wtd][epitope][nbrdist_percentile]
                    assert self_nbrdists[0] <= self_nbrdists[-1] #confirm sorted
                    nbrdist = one_info['{}_{}_{}nbrdist{}'.format(epitope,chains,wtd,nbrdist_percentile)]
                    rank = ch_cdrdist.get_rank( nbrdist, self_nbrdists )
                    one_info['{}_{}_{}nbrdist{}rank'.format(epitope,chains,wtd,nbrdist_percentile)] \
                        = '{:.3f}'.format(rank)
        return one_info

    print('time for rank finding')
    pool = Pool(20)
    all_infos = pool.map(multiprocrank, [i for i in range(total_tcrs)])
    pool.close()
    pool.join()
    all_info = list(all_infos)
    print('done rank finding')
    # ## now rank versions of densities
    # for e in epitopes:
    #     for sdev in distance_sdevs:
    #         epitope_self_densities[e][sdev].sort()

    # for ii in range(total_tcrs):
    #     for epitope in epitopes:
    #         for distance_sdev in distance_sdevs:
    #             self_densities = epitope_self_densities[epitope][distance_sdev]
    #             assert self_densities[0] <= self_densities[-1] #confirm sorted
    #             dens = all_info[ii]['{}_{}_dens{}'.format(epitope,chains,distance_sdev)]
    #             rank = tcr_distances.get_rank( dens, self_densities )
    #             all_info[ii]['{}_{}_rdens{}'.format(epitope,chains,distance_sdev)] = '{:.3f}'.format(rank)



with open( outfile,'w') as out:
    out.write('\t'.join(outfields)+'\n')

    for outl in all_info:
        out.write( ch_w.make_tsv_line(outl, outfields)+'\n')
