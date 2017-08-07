import os
import argparse
import glob
import csv
import ch_antigen_base as ch_w


#==============================

curr_version = 1.0
parser = argparse.ArgumentParser(description="This script will do additional work with vdjdb_imgted.")

parser.add_argument("-i", nargs=1, type=str, default="../imgt_work/vdjdb_imgted.txt",
                    help="input file path.")
parser.add_argument("-o", nargs=1, type=str, default="../imgt_work/",
                    help="Output path.")
parser.add_argument("-o2", nargs=1, type=str, default="vdjdb",
                    help="Output file name.")

args = parser.parse_args()

if type(args.i) is list:
    args.i = args.i[0]
input_path = os.path.abspath(args.i)

if type(args.o) is list:
    args.o = args.o[0]
output_path = os.path.abspath(args.o)

ch_w.crdir(output_path)

if type(args.o2) is list:
    args.o2 = args.o2[0]
output_file = args.o2

#t1 = ['ENSA......SN DIRSN...ME LDKKA.', 'MNH.......DT YYD....KIL P.NNSF',
#      'CAAAAGNTGKLIF', 'CASSRLGASAETLYF']
#t2 = ['ENSA......SN DIRSN...ME LDKKA.', 'MNH.......DT YYD....KIL P.NNSF',
#      'CAAAGGNNKLTF', 'CASSSLGSSAETLYF']
tcrseqs =  [] #both alpha and beta are available
alphatcr = [] #if only alpha is available
betatcr = [] #if only beta is available
tcrdeficient = [] #none is available

with open(input_path, 'r') as inp:
    reader = csv.reader(inp, delimiter='\t')
    header = next(reader)
    for row in reader:
        info = ch_w.parse_tsv_line(row, header)
        #row[0:4] – alpha seq, [4-6]: TRAV, TRAJ; [6:10] - beta seq, [10:13] - TRBV, TRBD, TRBJ
        inforow = [row[0]+' '+row[1]+' '+row[2], #cdr1-cdr2.5
                   row[6]+' '+row[7]+' '+row[8], #cdr1-cdr2.5
                   row[3], row[9],               #CDR3, CDR3
                   row[4], row[5],               #TRAV, TRAJ
                   row[10], row[11], row[12],    #TRBV, TRBD, TRBJ
                   row[17]]                      #epitope
        if (info['v.alpha'] == '') or (info['v.beta'] == ''):
            if info['v.alpha'] != '' and info['j.alpha'] != '':
                if info['cdr1.alpha'] != '' and info['cdr2.alpha'] != '' and \
                                info['cdr2.5.alpha'] !='' and info['cdr3.alpha'] != '':
                    alphatcr.extend('\t'.join(row)+'\n')

            elif info['v.beta'] != '' and info['j.beta'] != '':
                if info['cdr1.beta'] != '' and info['cdr2.beta'] != '' and \
                                info['cdr2.5.beta'] != '' and info['cdr3.beta'] != '':
                    betatcr.extend('\t'.join(row)+'\n')

            else:
                tcrdeficient.extend('\t'.join(row)+'\n')
        else:
            if info['cdr1.alpha'] == '' or info['cdr2.alpha'] == '' or \
                            info['cdr2.5.alpha'] == '' or info['cdr3.alpha'] == '':
                print('strange row: %s %s' % ('alpha', row))
            elif info['cdr1.beta'] == '' or info['cdr2.beta'] == '' or \
                                info['cdr2.5.beta'] == '' or info['cdr3.beta'] == '':
                print('strange row: %s %s' % ('beta', row))
            else:
                tcrseqs.extend('\t'.join(row)+'\n')

with open(os.path.join(output_path, '{}_{}.txt'.format(output_file, 'ab')), 'w') as out:
    out.write('\t'.join(header)+'\n')
    out.writelines(tcrseqs)


with open(os.path.join(output_path, '{}_{}.txt'.format(output_file, 'a')), 'w') as out:
    out.write('\t'.join(header)+'\n')
    out.writelines(alphatcr)

with open(os.path.join(output_path, '{}_{}.txt'.format(output_file, 'b')), 'w') as out:
    out.write('\t'.join(header)+'\n')
    out.writelines(betatcr)

with open(os.path.join(output_path, '{}_{}.txt'.format(output_file, 'o')), 'w') as out:
    out.write('\t'.join(header)+'\n')
    out.writelines(tcrdeficient)