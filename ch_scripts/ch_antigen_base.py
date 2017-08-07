"""
Module for functions used in order to get substitution matrices.
------- Version 1.0
"""
import os
import glob
import pandas as pd
import csv
import copy
#import graphs_var1 as ch1_gr
import random

#import substitution_matrix as ch2_subst


def crdir(path):
    if os.path.exists(path):
        return 0
    else:
        os.mkdir(path)

def parse_tsv_line(tsvline, headinfo):
    if type(tsvline) != list:
        wtsvline = tsvline.split('\t')
    else:
        wtsvline = tsvline
    wtsvline[-1] = wtsvline[-1].strip()
    assert(len(wtsvline)==len(headinfo))
    outinfo = {}
    for tag,val in zip(headinfo, wtsvline):
        outinfo[tag] = val
    return outinfo

def make_tsv_line(ininfo, headinfo, empty_string_replacement=''):
    l = []
    for tag in headinfo:
        val = ininfo[tag]
        if type(val) is str:
            if empty_string_replacement and not val:
                l.append(empty_string_replacement)
            else:
                l.append(val)
        else:
            l.append(str(val))
    return '\t'.join(l)

def getdata(path, gene=None, orgname=None):
    """This function will work with vdjdb from where it will take dictionary with antigens as keys.
input: path, gene and orgname. path to file, orgname - name of organism, from which you will get antigens."""
    seqspecies = {}
    sequences_with_antigen = {}
    if gene == None and orgname == None:
        with open(path, 'r') as inp:
            for line in inp:
                inf = line.strip().split()
                if inf[3] in seqspecies:
                    seqspecies[inf[3]].append(inf[1])
                else:
                    seqspecies[inf[3]] = [inf[1]]
                if inf[1] in sequences_with_antigen:
                    sequences_with_antigen[inf[1]].append(inf[3])
                else:
                    sequences_with_antigen[inf[1]] = [inf[3]]
    elif gene == None:
        with open(path, 'r') as inp:
            for line in inp:
                inf = line.strip().split()
                if inf[2] == orgname:
                    if inf[3] in seqspecies:
                        seqspecies[inf[3]].append(inf[1])
                    else:
                        seqspecies[inf[3]] = [inf[1]]
                    if inf[1] in sequences_with_antigen:
                        sequences_with_antigen[inf[1]].append(inf[3])
                    else:
                        sequences_with_antigen[inf[1]] = [inf[3]]
    elif orgname == None:
        with open(path, 'r') as inp:
            for line in inp:
                inf = line.strip().split()
                if inf[0] == gene:
                    if inf[3] in seqspecies:
                        seqspecies[inf[3]].append(inf[1])
                    else:
                        seqspecies[inf[3]] = [inf[1]]
                    if inf[1] in sequences_with_antigen:
                        sequences_with_antigen[inf[1]].append(inf[3])
                    else:
                        sequences_with_antigen[inf[1]] = [inf[3]]
    else:
        with open(path, 'r') as inp:
            for line in inp:
                inf = line.strip().split()
                if inf[0] == gene and inf[2] == orgname:
                    if inf[3] in seqspecies:
                        seqspecies[inf[3]].append(inf[1])
                    else:
                        seqspecies[inf[3]] = [inf[1]]
                    if inf[1] in sequences_with_antigen:
                        sequences_with_antigen[inf[1]].append(inf[3])
                    else:
                        sequences_with_antigen[inf[1]] = [inf[3]]

    return {"seqspecies":seqspecies, "sequences_with_antigen":sequences_with_antigen}


def test_and_teach(seq_dict, ant_dict, limit, times, filename):
    """seq_dict – dictionary. key = antigen. limit – how many sequences should interact with specific antigen. times –
how size of test sample is lower than quantity of all sequences"""
    seq_dict_teach = copy.deepcopy(seq_dict)
    ant_dict_teach = copy.deepcopy(ant_dict)
    ant_dict_test = {}
    for i in seq_dict_teach:
        if len(seq_dict_teach[i]) <= limit:
            pass
        else:
            for k in range(int(len(seq_dict_teach[i])/times)):
                number = random.randint(0, len(seq_dict_teach[i])-1)
                wseq = seq_dict_teach[i][number]
                if wseq in ant_dict_teach:
                    ant_dict_test[wseq] = ant_dict_teach[wseq]
                    del ant_dict_teach[wseq]


    with open(filename+'_test', 'w') as out:
        out.write('\n'.join(sorted(list(ant_dict_test))) + '\n')

    with open(filename+'_teach', 'w') as out:
        out.write('\n'.join(sorted(list(ant_dict_teach))) + '\n')

    with open(filename + '_test_antigen', 'w') as out:
        for i in ant_dict_test:
            out.write('\n'.join(list(map(lambda x: str(i) + '\t' + str(x), ant_dict_test[i]))) + '\n')
    with open(filename + '_teach_antigen', 'w') as out:
        for i in ant_dict_teach:
            out.write('\n'.join(list(map(lambda x: str(i) + '\t' + str(x), ant_dict_teach[i]))) + '\n')


    half_1 = copy.deepcopy(ant_dict_test)
    half_2 = {}
    for k in range(int(len(half_1) / 2)):
        number = random.randint(0, len(half_1) - 1)
        wseq = list(half_1)[number]
        half_2[wseq] = half_1[wseq]
        half_1.pop(wseq)


    with open(filename+'_half_first_testing', 'w') as out:
        out.write('\n'.join(half_1) + '\n')
    with open(filename+'_half_sacred_testing', 'w') as out:
        out.write('\n'.join(half_2) + '\n')
    with open(filename + '_half_first_testing_antigen', 'w') as out:
        for i in half_1:
            out.write('\n'.join(list(map(lambda x: str(i) + '\t' + str(x), half_1[i]))) + '\n')
    with open(filename + '_half_sacred_testing_antigen', 'w') as out:
        for i in half_2:
            out.write('\n'.join(list(map(lambda x: str(i) + '\t' + str(x), half_2[i]))) + '\n')

def gethseq(inp):
    hseq = []
    with open(inp, 'r') as inp:
        inf = inp.read()
        hseq += inf.split('\n')
        return hseq


def groovy_treesearch(inp, opt, out):
    os.system("groovy TreeSearch.groovy " + inp + " " + opt + " " + out)


def get_org_pairalignment(inp, inppath, outpath, hseq):
    if os.path.exists(outpath):
        answer = input(outpath + " already exists. Do you want to overwrite it? (y/n)")
        if answer == "y":
            with open(inppath, 'r') as inp:
                with open(outpath, 'w') as out:
                    for line in inp:
                        seqs = line.replace('-', '').upper().strip().split()
                        if seqs[0] in hseq and seqs[1] in hseq:
                            out.write(line)
        else:
            raise Exception("function terminated")


def get_alignments_set(wrkfolder, alignm, outfolder):
    crdir(os.path.join(wrkfolder, outfolder))
    for i in alignm:
        algn = set()
        with open(i, 'r') as inp:
            for line in inp:
                algn.add(line.strip().replace('-','').upper())
        with open(os.path.join(wrkfolder, outfolder, i.split('/')[-1]+'_set'), 'w') as out:
            out.write('\n'.join(algn))

def combine_alignments_sets(wrkfolder, alignone, aligntwo, onetwoseq, outfolder):
    crdir(os.path.join(wrkfolder, outfolder))
    algnone = {}
    algntwo = {}
    

class Sequence:
    """New type that has following attributes: seq (sequence), antigen (all antigens, to whom this sequence belongs),
inlinks and outlinks (links to other sequences). You can use following commands: add_link(new_sequence) and
outlinked(antigen, opt)"""
    def __init__(self, seq, antigen):
        self.seq = seq
        self.antigen = antigen
        self.inlinks = {}
        self.outlinks = set()
        for i in self.antigen:
            self.inlinks[i]=[]

    def add_link(self, sequence):
        """
Add inlink or outlink to sequence (depending of antigen for new sequence). Type of added sequence is Sequence.
However, by self.inlink or self.outlink you'll get only string.
        """
        inantigen = 0
        for i in sequence.antigen:
            if i in self.antigen:
                self.inlinks[i].append(sequence)
                inantigen += 1
        if inantigen == 0:
            self.outlinks.update([sequence])

    def outlinked(self, antig):
        """
This command will get all links to sequences, that do not belong to chosen antigen (by default), but do belong to other
sequence's antigens.
        """
        outlink = 0
        inlink = 0
        not_links = set(self.inlinks[antig])
        inlink += len(not_links)
        for i in self.inlinks:
            if i == antig:
                pass
            else:
                for k in self.inlinks[i]:
                    if k in not_links:
                        inlink += 1
                        pass
                    else:
                        outlink += 1
        return outlink

class AntigenGroup:
    """New type that has been created only to concentrate sequences belonging to one antigen. Following attributes are:
antigen, sequences, boundarysequences (they belong to more than one group), outlinks, inlinks, totoutlinks, totinlinks.
You can use following command: add_seq(sequence), if you want to add new sequence to group.
"""
    def __init__(self, antigen):
        self.antigen = antigen
        self.sequences = []
        self.boundarysequences = [] #they belong to more than one group
        self.outlinks = set()
        self.inlinks = set()
        self.totoutlinks = 0
        self.totinlinks = 0

    def add_seq(self, sequence):
        """Add new sequence to Antigengroup.
"""
        if len(sequence.antigen) == 1:
            if sequence not in self.sequences:
                self.sequences.append(sequence)
                self.outlinks.update(sequence.outlinks)
                self.totoutlinks += len(sequence.outlinks)
                self.inlinks.update(sequence.inlinks[self.antigen])#not all inlinks are from one antigen
                self.totinlinks += len(sequence.inlinks[self.antigen])
        elif len(sequence.antigen) > 1:
            if sequence not in self.sequences:
                self.boundarysequences.append(sequence)
                self.outlinks.update(sequence.outlinks)
                self.totoutlinks += len(sequence.outlinks)
                out = set()
                inn = set()
                inn.update(sequence.inlinks[self.antigen])
                self.totinlinks += len(sequence.inlinks[self.antigen])
                for i in sequence.antigen:
                    if i == self.antigen:
                        pass
                    else:
                        for k in sequence.inlinks[i]:
                            if k in inn:
                                pass
                                #self.totinlinks += 1
                            else:
                                out.update([k])
                                self.totoutlinks += 1
                self.outlinks.update(out)
                self.inlinks.update(inn)


def get_seqs_antigens(sequences_with_antigen):
    sequences = {}
    antigens = {}
    for i in sequences_with_antigen:
        sequences[i] = Sequence(i, sequences_with_antigen[i])
        for k in sequences_with_antigen[i]:
            if k in antigens:
                pass
            else:
                antigens[k] = AntigenGroup(k)
    return(sequences, antigens)


def get_information_from_alignment(alignmentset, sequences, antigens, sequences_with_antigen):
    sequences_2 = copy.deepcopy(sequences)
    antigens_2 = copy.deepcopy(antigens)
    with open(alignmentset, 'r') as inp:
        for line in inp:
            seq1, seq2 = line.strip().split()[0], line.strip().split()[1]
            sequences_2[seq1].add_link(sequences_2[seq2])
            sequences_2[seq2].add_link(sequences_2[seq1])
        for i in sequences_with_antigen:
            for k in sequences_with_antigen[i]:
                antigens_2[k].add_seq(sequences_2[i])
    return sequences_2, antigens_2


def write_information_from_matrices(alignmentset, outputfile, sequences, antigens, sequences_with_antigen):
    sequences_2, antigens_2 = get_information_from_alignment(alignmentset, sequences, antigens, sequences_with_antigen)

    listantigens = []
    for k in antigens_2:
        listantigens.append([k, antigens_2[k]])
    listantigens.sort(key = lambda x: len(x[1].sequences))
    with open(outputfile, 'w') as out:
        writer = csv.writer(out, delimiter='\t')
        writer.writerow(['antigen', 'number_of_sequences', 'possible_inlinks', 'real_inlinks', 'possible_outlinks',
                        'real_outlinks', 'number_of_boseq', 'inlinks_of_boseq', 'outlinks_of_boseq',
                        'possibe_boseq_links_to_other_antigens', 'real_boseq_links_to_other_antigens'])
        for k in listantigens:
            wantigen = k[0]
            wseqnumber = len(antigens_2[wantigen].sequences)+len(antigens_2[wantigen].boundarysequences)
            wposinlinks = ((wseqnumber**2)-wseqnumber)
            wrealinlinks = (antigens_2[wantigen].totinlinks)
            wposoutlinks = wseqnumber*(len(sequences_with_antigen)-wseqnumber)
            wrealoutlinks = antigens_2[wantigen].totoutlinks
            boseq = antigens_2[wantigen].boundarysequences
            if len(boseq) > 0:
                wboseqnumber = len(antigens_2[wantigen].boundarysequences)
                boseqin = 0
                boseqout = 0
                posboseqother = 0
                realboseqother = 0
                for i in boseq:
                    realboseqother += i.outlinked(wantigen)
                    boseqout += len(i.outlinks)+realboseqother
                    for antig in i.antigen:
                        if antig == wantigen:
                            pass
                        else:
                            #print(antig)
                            posboseqother += len(antigens_2[antig].sequences)
                            for seq in antigens_2[antig].boundarysequences:  #...
                                if wantigen in seq.antigen:
                                    #print(wantigen, seq.antigen, 'pass')
                                    pass
                                else:
                                    #print(wantigen, seq.antigen, 'not pass')
                                    posboseqother += 1
                    boseqin += len(i.inlinks[wantigen])
                wboseqinlinks = boseqin
                wboseqoutlinks = boseqout
                wposboseqother = posboseqother
                wrealboseqother = realboseqother
            else:
                wboseqnumber, wboseqinlinks, wboseqoutlinks, wposboseqother, wrealboseqother = 0, 0, 0, 0, 0
            writer.writerow([wantigen, wseqnumber, wposinlinks, wrealinlinks, wposoutlinks, wrealoutlinks, wboseqnumber, wboseqinlinks, wboseqoutlinks, wposboseqother, wrealboseqother])


def get_info_from_matrices(outfolder, alignms, gene='', addpar=''):
    #alignms = glob.glob(folder+'/'+glword)
    info = {}
    for name in alignms:
        with open(name, 'r') as inp:
            reader = csv.reader(inp, delimiter='\t')
            info[name] = {'realin':0,'realout':0, 'realboseqout':0, 'realboseqother':0, 'posin':0, 'posout':0, 'realboseqin':0, 'posboseqother':0}
            for row in reader:
                for row in reader:
                    info[name]['realin'] += int(row[3])
                    info[name]['realout'] += int(row[5])
                    info[name]['realboseqout'] += int(row[8])
                    info[name]['realboseqother'] += int(row[10])
                    info[name]['posin'] += int(row[2])
                    info[name]['posout'] += int(row[4])
                    info[name]['realboseqin'] += int(row[7])
                    info[name]['posboseqother'] += int(row[9])
    sortedinfo = []
    for name in info:
        sortedinfo.append([name, info[name]])
    sortedinfo.sort(key = lambda x: x[0])
    with open(os.path.join(outfolder,'{}{}sum_of_info_from_alignments.txt'.format(gene, addpar)), 'w') as out:
        writer = csv.writer(out, delimiter='\t')
        writer.writerow(['name', 'rpin/(posin+posout)', 'rpout/(posin+posout)', '(rpin-rpout)/(posin+posout)', 'realin', 'realout', 'realboseqin', 'realboseqout',
                         'realboseqother', 'posin', 'posout', 'posboseqother'])
        for inf in sortedinfo:
            name = inf[0]
            writer.writerow([name.split('/')[-1], "{0:.5f}".format(float(info[name]['realin'])/(info[name]['posin']+info[name]['posout'])),
                             "{0:.5f}".format(float(info[name]['realout'])/(info[name]['posin']+info[name]['posout'])),
                             "{0:.5f}".format((float(info[name]['realin']) / (info[name]['posin']+info[name]['posout']))-(float(info[name]['realout']) / (info[name]['posin']+info[name]['posout']))),
                             info[name]['realin'], info[name]['realout'],info[name]['realboseqin'], info[name]['realboseqout'],
                             info[name]['realboseqother'], info[name]['posin'], info[name]['posout'], info[name]['posboseqother']])


def get_occurence_matrix(occurence, alignments, substitutions, msubstitutions, antig_seq, name, folder, return_results = False):
    work_occur = copy.deepcopy(occurence)
    mwork_occur = copy.deepcopy(occurence)

    def getinfo(gisub, gioccurence): #get local aligninfo about substitutions and occurences.
        for m in range(len(pair[0])):
            a1, a2 = pair[0][m].upper(), pair[1][m].upper()
            if a1 != '-' and a2 != '-':
                # if a1 != a2:
                gisub[a1][a2] += 1
                gioccurence[a1] += 1
                gioccurence[a2] += 1
        return(gisub, gioccurence)

    for pairset in alignments: #one pair without gaps.
        loc_occurence, loc_moccurence = copy.deepcopy(occurence), copy.deepcopy(occurence) #dict

        localsub, localmsub = pd.DataFrame(0, index=occurence, columns=occurence), \
                              pd.DataFrame(0, index=occurence, columns=occurence) #square matrices

        p1, p2 = pairset.split('|')[0], pairset.split('|')[1]
        sameantig = 0
        for antig_seq1 in antig_seq[p1]: #check, whether pair is able to bind same antigen
            for antig_seq2 in antig_seq[p2]:
                if antig_seq1 == antig_seq2:
                    sameantig = 1

        for pair in alignments[pairset]: #pair with possible gaps
            if len(pair[0]) != len(pair[1]):
                print('different lens', pair[0], pair[1])
                print('ERRROROR')
                return 0
            elif sameantig == 1:
                localsub, loc_occurence = getinfo(localsub, loc_occurence)
            else:
                localmsub, loc_moccurence = getinfo(localmsub, loc_moccurence)

        substitutions += localsub / len(alignments[pairset])
        msubstitutions += localmsub / len(alignments[pairset])

        for k in work_occur:
            work_occur[k] += loc_occurence[k]/len(alignments[pairset])
        for k in mwork_occur:
            mwork_occur[k] += loc_moccurence[k]/len(alignments[pairset])

    substitutions.to_csv(os.path.join(folder, 'substitution_matrix_' + name + '.csv'))
    msubstitutions.to_csv(os.path.join(folder, 'msubstitution_matrix_' + name + '.csv'))

    out = open(os.path.join(folder, 'aminoacid_occurence_' + name), 'w')
    for i in work_occur:
        out.write(str(i) + '\t' + str(work_occur[i]) + '\n')
    out.close()

    out = open(os.path.join(folder, 'aminoacid_moccurence_' + name), 'w')
    for i in mwork_occur:
        out.write(str(i) + '\t' + str(mwork_occur[i]) + '\n')
    out.close()

    if return_results == True:
        return substitutions, msubstitutions, work_occur, mwork_occur


def get_alignments(inppath):
    alignments = {}
    with open(os.path.join(inppath), 'r') as inp:
        for line in inp:
            algn = line.strip().split()
            oseq, tseq = algn[0].replace('-', '').upper(), algn[1].replace('-', '').upper()
            if oseq + '|' + tseq in alignments:
                if [algn[0], algn[1]] in alignments[oseq + '|' + tseq]:
                    pass
                else:
                    alignments[oseq + '|' + tseq].append([algn[0], algn[1]])
            else:
                alignments[oseq + '|' + tseq] = [[algn[0], algn[1]]]
    return alignments

def get_alignments_var2(inppath):
    #here you will not get symmetry like in get_alignments(infname)
    alignments = {}
    with open(os.path.join(inppath), 'r') as inp:
        for line in inp:
            algn = line.strip().split()
            oseq, tseq = algn[0].replace('-', '').upper(), algn[1].replace('-', '').upper()
            if oseq + '|' + tseq in alignments:
                if [algn[0], algn[1]] in alignments[oseq + '|' + tseq]:
                    pass
                else:
                    alignments[oseq + '|' + tseq].append([algn[0], algn[1]])
            elif tseq + '|' + oseq in alignments:
                if [algn[1], algn[0]] in alignments[tseq + '|' + oseq]:
                    pass
                else:
                    alignments[tseq + '|' + oseq].append([algn[1], algn[0]])
            else:
                alignments[oseq + '|' + tseq] = [[algn[0], algn[1]]]
    return alignments

def get_sequences_with_antigen(inppath):
    sequences_with_antigen = {}
    with open(inppath, 'r') as inp:
        for line in inp:
            inf = line.strip().split()
            if inf[0] in sequences_with_antigen:
                sequences_with_antigen[inf[0]].append(inf[1])
            else:
                sequences_with_antigen[inf[0]] = [inf[1]]
    return sequences_with_antigen