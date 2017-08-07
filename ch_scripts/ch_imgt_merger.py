import os
import pandas as pd
import numpy as np
import argparse
import ch_antigen_base as ch_w
from Bio import SeqIO
import csv

#==============================
curr_version = 1.0
parser = argparse.ArgumentParser(description="This script will get information from vdjdb_full and imgt.fa and then "
                                             "it will sequences of CDR1, CDR2 and CDR3 by imgt anchor positions.")

parser.add_argument("-i", nargs=1, type=str, default="../vdjdb-2017-06-13/vdjdb_full.txt",
                    help="vdjdb path.")
parser.add_argument("-i2", nargs=1, type=str, default="../imgt_work/imgt_all.fasta", #"../vdjdb-2017-06-13/imgt.fa",
                    help="imgt.fa path.")
parser.add_argument("-o", nargs=1, type=str, default="../imgt_work/",
                    help="Output path.")
parser.add_argument("-o2", nargs=1, type=str, default="vdjdb_imgted.txt",
                    help="Output file name.")

args = parser.parse_args()

if type(args.i) is list:
    args.i = args.i[0]
input_file = os.path.abspath(args.i)

if type(args.i2) is list:
    args.i2 = args.i2[0]
input_file2 = os.path.abspath(args.i2)

if type(args.o) is list:
    args.o = args.o[0]
output_path = os.path.abspath(args.o)

ch_w.crdir(output_path)

if type(args.o2) is list:
    args.o2 = args.o2[0]
output_file = args.o2

#http://www.imgt.org/IMGTScientificChart/Nomenclature/IMGT-FRCDRdefinition.html
anchors = {'CDR1': [26, 38], 'CDR2': [55, 65], 'CDR2.5': [80, 86]} #'CDR3': [104, 116]

colnames = ['cdr1.alpha', 'cdr2.alpha', 'cdr2.5.alpha', 'cdr3.alpha', 'v.alpha', 'j.alpha', #5
            'cdr1.beta', 'cdr2.beta', 'cdr2.5.beta', 'cdr3.beta', 'v.beta', 'd.beta', 'j.beta', #11
            'species', 'mhc.a', 'mhc.b', 'mhc.class', 'epitope.seq', 'epitope', #17
            'antigen.species', 'meta.epitope.id', 'subject.id', 'clone.id'] #20

dictimgt = {}

with open(input_file2, "rU") as inp:
    for record in SeqIO.parse(inp, "fasta"):
        info = record.description.strip().split('|')
        specinfo = info[2].replace(' ', '').lower()
        if specinfo not in dictimgt:
            dictimgt[specinfo] = {}
        if info[1] in dictimgt[specinfo] and record.seq != dictimgt[specinfo][info[1]]:
            print(record.seq, info[1], dictimgt[specinfo][info[1]])
        if info[1].startswith('TRGC') or info[1].startswith('TRAC') or info[1].startswith('TRBC') or info[1].startswith('TRDC'):
            continue
        dictimgt[specinfo][info[1]] = record.seq

def changeseq(inpseq, positions):
    outdict, outdict2 = {}, {}
    for i in positions:
        outseq = str(inpseq[positions[i][0]:positions[i][1]])
        outdict[i] = outseq
        outdict2[i] = outseq.replace('.', '')
    return outdict, outdict2

result_table = []
with open(input_file) as inp:
    reader = csv.reader(inp, delimiter='\t')
    header = next(reader)
    for row in reader:
        info = ch_w.parse_tsv_line(row, header)
        cdralpha = ['', '', '', info['cdr3.alpha'], info['v.alpha'], info['j.alpha']]
        cdrbeta = ['', '', '', info['cdr3.beta'], info['v.beta'], info['d.beta'], info['j.beta']]
        otherinfo = [info['species'], info['mhc.a'], info['mhc.b'], info['mhc.class'], info['antigen.epitope'],
                     info['antigen.gene'], info['antigen.species'], info['meta.epitope.id'], info['meta.subject.id'],
                     info['meta.clone.id']]
        info['species'] = info['species'].lower()
        if info['species'] in dictimgt:
            if info['v.alpha'] in dictimgt[info['species']]:
                outdict, outdict2 = changeseq(dictimgt[info['species']][info['v.alpha']], anchors)
                cdralpha[0], cdralpha[1], cdralpha[2] = outdict['CDR1'], outdict['CDR2'], outdict['CDR2.5']
            if info['v.beta'] in dictimgt[info['species']]:
                outdict, outdict2 = changeseq(dictimgt[info['species']][info['v.beta']], anchors)
                cdrbeta[0], cdrbeta[1], cdrbeta[2] = outdict['CDR1'], outdict['CDR2'], outdict['CDR2.5']
        result_table.extend('\t'.join(cdralpha + cdrbeta + otherinfo)+'\n')

with open(os.path.join(output_path, output_file), 'w') as out:
    out.write('\t'.join(colnames)+'\n')
    out.writelines(result_table)