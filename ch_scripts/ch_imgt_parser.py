import urllib.request
import urllib.parse
import argparse
import os
#http://www.imgt.org/genedb/GENElect?query=7.3+TRAV&species=Homo+sapiens
#http://www.imgt.org/genedb/GENElect?query=7.3+Group&species=Species

#==============================
curr_version = 1.0
parser = argparse.ArgumentParser(description="This script will parse aminoacid sequences for TR genes from various species from imgt")

parser.add_argument("-o", nargs=1, type=str, default="../imgt_work/imgt_all_2.fasta",
                    help="output file path.")

args = parser.parse_args()
if type(args.o) is list:
    args.o = args.o[0]
output_file = args.o

outputinfo = []
tr_species = {'Homo+sapiens': 'Homo sapiens',
              'Mus+musculus': 'Mus musculus',
              'Bos+taurus': 'Bos taurus',
              'Canis+lupus+familiaris': 'Canis lupus familiaris',
              'Macaca+mulatta': 'Macaca mulatta',
              'Mus+minutoides': 'Mus minutoides',
              'Mus+pahari': 'Mus pahari',
              'Mus+spretus': 'Mus spretus',
              'Oncorhynchus+mykiss': 'Oncorhynchus mykiss',
              'Ovis+aries': 'Ovis aries'}

tr_group = ['TRAV', 'TRAJ', 'TRAC', 'TRBV', 'TRBD', 'TRBJ', 'TRBC']

spiti = 1
for sp in tr_species:
    print(sp, str(spiti)+'/'+str(len(tr_species)))
    spiti += 1
    trgriti = 1
    for trgr in tr_group:
        print(trgr, str(trgriti)+'/'+str(len(tr_group)))
        trgriti += 1
        inpurl = urllib.request.urlopen('http://www.imgt.org/genedb/GENElect?query=7.3+'+trgr+'&species='+sp)
        inpurlheader = inpurl.headers.get_content_charset()
        for i in inpurl:
            if str(i).startswith('b\'<b>Number'): #part of the site, where it says how many entries are on the list
                for i in inpurl:
                    if str(i).startswith('b\'<pre>'): #next line will be fasta
                        for i in inpurl:
                            if str(i).startswith('b\'\\r\\n'): #end of fasta
                                break
                            else:
                                info = i.decode(inpurlheader)
                                if info.startswith('>'):
                                    info = info.split('|')
                                    info[2] = tr_species[sp]
                                    info = '|'.join(info)
                                outputinfo.append(info)
                        break
                break
    print('\n')

with open(output_file, 'w') as out:
    out.writelines(outputinfo)