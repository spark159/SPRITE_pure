import glob
import os, sys, subprocess, re
from argparse import ArgumentParser, FileType
import copy
import itertools
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np

chars = "ACGT"

def neighbors(pattern, d):
    assert(d <= len(pattern))

    if d == 0:
        return [pattern]

    r2 = neighbors(pattern[1:], d-1)
    r = [c + r3 for r3 in r2 for c in chars if c != pattern[0]]

    if (d < len(pattern)):
        r2 = neighbors(pattern[1:], d)
        r += [pattern[0] + r3 for r3 in r2]

    return r

def rev_comp (seq):
    dic={'A':'T', 'T':'A', 'C':'G', 'G':'C', 'N':'N'}
    output=''
    for nt in seq:
        output+=dic[nt]
    return output[::-1]

def hamming_dist (seq1, seq2):
    assert len(seq1) == len(seq2)
    dist = 0
    for i in range(len(seq1)):
        if seq1[i] != seq2[i]:
            dist +=1
    return dist

def abs_cmp (a, b):
    if abs(a) < abs(b):
        return -1
    elif abs(a) > abs(b):
        return 1
    else:
        if a < b:
            return -1
        elif a > b:
            return 1
        else:
            return 0

def all_path (N, states='ATCG'):
    if N==1:
        return list(states)
    output=[]
    for path in all_path(N-1):
        for state in states:
            output.append(path+state)
    return output

# mutate the sequence into other within the given hamming distance
def all_mutants (seq, hamming_dist, states='ATCG'):
    for nt in seq:
        assert nt in states
        
    if hamming_dist <= 0:
        return [seq]
    if hamming_dist > len(seq):
        return None

    def all_path (states_list):
        N = len(states_list)
        if N==1:
            return list(states_list[0])
        output = []
        for path in all_path(states_list[:-1]):
            for state in states_list[-1]:
                output.append(path+state)
        return output

    alts_list = []
    for i in range(len(seq)):
        alts = list(set(list(states)) - set(seq[i]))
        alts_list.append(alts)

    output = []
    for dist in range(1, hamming_dist+1):
        for pos_tuple in itertools.combinations(range(len(seq)), dist):
            pos_list = list(pos_tuple)
            for word in all_path([alts_list[pos] for pos in pos_list]):
                new_mutant = list(seq[:])
                for i in range(len(word)):
                    pos = pos_list[i]
                    new_mutant[pos] = word[i]
                new_mutant = ''.join(new_mutant)
                output.append(new_mutant)

    return output    

def read_tags (fname):
    tagID_seq, seq_tagID = {}, {}
    First = True
    for line in open(fname):
        if First:
            First = False
            continue
        well, name, seq = line.strip().split()
        cate, strand, num = name.split('_')
        seq = seq[-17:]
        tagID = cate + '_' + num
        tagID_seq[tagID] = seq
        seq_tagID[seq] = tagID
    return tagID_seq, seq_tagID
oddtagID_seq, seq_oddtagID = read_tags('Odd_Bottom_Plate.csv')
eventagID_seq, seq_eventagID = read_tags('Even_Bottom_Plate.csv')

##Sanity check
#oddtagIDs = oddtagID_seq.keys()
#min_dist = sys.maxint
#for i in range(len(oddtagIDs)-1):
#    for j in range(i+1, len(oddtagIDs)):
#        ID1, ID2 = oddtagIDs[i], oddtagIDs[j]
#        seq1, seq2 = oddtagID_seq[ID1], oddtagID_seq[ID2]
#        dist = hamming_dist(seq1, seq2)
#        if dist < min_dist:
#            min_dist = copy.deepcopy(dist)
#print min_dist

# tag ID identification
tag_format = 'TERM|SPACER|EVEN|SPACER|ODD|SPACER|EVEN|SPACER|ODD|SPACER|EVEN|SPACER|ODD'

cate_len = {'UMI':10, 'TERM':9, 'SPACER':7, 'ODD':17, 'EVEN':17, 'DPM':21}
cate_seq_tagID = {'ODD':seq_oddtagID, 'EVEN':seq_eventagID}

laxity = 2
Hamming_dist = 2

# allow sequencing error within the given hammidng distance
for cate in cate_seq_tagID:
    seq_tagID = cate_seq_tagID[cate]
    for seq, tagID in seq_tagID.items():
        new_pool = all_mutants(seq, Hamming_dist)
        for new_seq in new_pool:
            seq_tagID.update({new_seq:tagID})


tag_format_list = tag_format.split('|')
tag_format_len = 0
for cate in tag_format_list:
    tag_format_len += cate_len[cate]


path = '/home/spark159/Projects/2021_02_06_SPRITE_test_mistake/'
#fname = '1_S1_L001_R2_001.fastq'
fname = '2_S2_L001_R2_001.fastq'

temp_fname = fname.split('.')[0] + "_temp.fastq"
invalid_sort_fname = fname.split('.')[0] + "_invalid.sort"
valid_sort_fname = fname.split('.')[0] + "_valid.sort"
sort_fname = fname.split('.')[0] + ".sort"

temp_file = open(temp_fname, 'w')
invalid_sort_file = open(invalid_sort_fname, 'w')

read_count = 0
barcoded_count = 0
line_count = 0
BClen_count = {}
for line in open(path+fname):
    line = line.strip()
    if line_count % 4 == 0:
        assert line.startswith('@')
        read_id = ':'.join(re.split(':| ', line[1:]))
        read_count +=1
        line_count +=1
        continue
    elif line_count % 4 == 1:
        read_seq = line
        line_count +=1
        continue
    elif line_count % 4 == 2:
        read_opt = line
        line_count +=1
        continue
    else:
        assert line_count % 4 == 3 
        read_qual = line
        line_count +=1

    type = 'NA'
        
    #if len(read_seq) < tag_format_len-5: # discard short reads
    #    type = 'invalid:short'
    #    print >> invalid_sort_file, "%s::%s" % (read_id, type)
    #    print >> invalid_sort_file, read_seq
    #    continue

    # identify SPRITE barcodes
    tagIDs = []
    pt = 0
    offset_list = sorted(range(-laxity, laxity+1), cmp=abs_cmp)
    for cate in tag_format_list:
        if cate in ['TERM', 'SPACER']:
            pt += cate_len[cate]
            continue
        
        found = False
        for j in offset_list:
            st, ed = pt + j, pt + j + cate_len[cate]
            if st < 0 or ed > len(read_seq):
                continue
            subseq = read_seq[st:ed]
            try:
                tagID = cate_seq_tagID[cate][subseq]
                tagIDs.append(tagID)
                found = True
                break
            except:
                pass
        if not found:
            break
        pt = copy.deepcopy(ed)

    BClen = len(tagIDs)
    if BClen not in BClen_count:
        BClen_count[BClen] = 0
    BClen_count[BClen] +=1

    # discard reads with incomplete barcoding
    if len(tagIDs) < 6: 
        type = 'invalid:%s-Barcode' % (len(tagIDs))
        print >> invalid_sort_file, "@%s::%s" % (read_id, type)
        print >> invalid_sort_file, read_seq
        continue


    ## check the extra part of reads
    #read_extra = read_seq[pt:]

    ## discard reads with too short extra part
    #if len(read_extra) < cate_len['SPACER'] + cate_len['UMI'] + cate_len['DPM'] + 10:
    #    type = 'invalid:barcodeOnly'
    #    print >> invalid_sort_file, "@%s::%s" % (read_id, type)
    #    print >> invalid_sort_file, read_seq
    #    barcoded_count +=1
    #    continue

    ## naively extract UMI and insert
    #st = cate_len['SPACER']
    #ed = st + cate_len['UMI']
    #UMI = read_extra[st:ed]

    # possibly valid reads
    spriteID = ''
    for ID in tagIDs:
        spriteID += '[' + ID + ']'


    # save the remaining part of reads to further analysis
    st = pt
    ed = min(st + cate_len['SPACER'] + cate_len['UMI'] + cate_len['DPM'] + 10, len(read_seq))

    print >> temp_file, '@' + read_id + '::' + spriteID
    #print >> temp_file, read_extra[:cate_len['SPACER'] + cate_len['UMI'] + cate_len['DPM'] + 10]
    print >> temp_file, read_seq[st:ed]
    print >> temp_file, read_opt
    #print >> temp_file, read_qual[:cate_len['SPACER'] + cate_len['UMI'] + cate_len['DPM'] + 10]
    print >> temp_file, read_qual[st:ed]
    barcoded_count +=1

invalid_sort_file.close()
temp_file.close()

# Identification of the remaining part of reads by Bowtie alignment

# write a reference file for Bowtie alignment
ref_fname = "SPRITE_test"
if not glob.glob("/home/spark159/Projects/SPRITE_pure/" + ref_fname +".ref"):
    extra_format = 'SPACER|UMI|DPM'
    widom_seq = 'CTGGAGAATCCCGGTGCCGAGGCCGCTCAATTGGTCGTAGACAGCTCTAGCACCGCTTAAACGCACGTACGCGCTGTCCCCCGCGTTTTAACCGCCAAGGGGATTACTCCCTAGTCTCCAGGCACGTGTCAGATATATACATCCTGC'
    f = open(ref_fname + ".ref", 'w')
    print >> f, '>0'
    print >> f, 'TGACTTG' + '-'*10 + 'TCTTCCGATCTTGGGTGTTTT' + widom_seq
    print >> f, '>1'
    print >> f, 'TGACTTG' + '-'*10 + 'TCTTCCGATCTTGGGTGTTTT' + rev_comp(widom_seq)
    f.close()

if not glob.glob("/home/spark159/Projects/SPRITE_pure/" + ref_fname + "*.bt2"):
    subprocess.call(["bowtie2-build", ref_fname + ".ref", ref_fname])

# parameters
UMI_st = 7
UMI_ed = UMI_st + cate_len['UMI']
Widom_st = cate_len['SPACER'] + cate_len['UMI'] + cate_len['DPM']
sted_wins = {UMI_st:['UMI'], UMI_ed:['UMI'], Widom_st:['Widom']}

mm_cutoff = 10

# Reads alignment on to the template
aligner_cmd = ["bowtie2", '-x', ref_fname, '-U', temp_fname]
aligner_cmd += ['--n-ceil', 'L,'+str(cate_len['UMI'])+','+str(0.15)]
aligner_cmd += ['--score-min', 'L,' + str(-100) + ',' +str(-100)]
aligner_cmd += ['-N', str(1), '-L', str(1), '-i',  'S,1,0']
aligner_cmd += ['--np', '0']
#align_proc = subprocess.Popen(aligner_cmd, stdout=subprocess.PIPE,stderr=subprocess.PIPE)
align_proc = subprocess.Popen(aligner_cmd, stdout=subprocess.PIPE,stderr=open("/dev/null", 'w'))

# start sort the reads
invalid_sort_file = open(invalid_sort_fname,'a')
valid_sort_file = open(valid_sort_fname,'w')

spriteID_UMIs = {}
valid_count = 0
#read_count = 0
for line in align_proc.stdout:
    if line.startswith('@'):
        continue
    #print line

    cols = line.strip().split()
    read_id, flag, ref_id, pos, mapQ, cigar_str = cols[:6]
    spriteID = read_id.split('::')[-1]
    #read_id=":".join(read_id.split(':')[3:7])
    #read_count +=1
    flag, pos = int(flag), int(pos)
    pos-=1
    #print read_id

    defects = []
    
    type = 'NA'

    # invalid: mapping failure
    if pos < 0 or flag & 0x4 != 0:
        type = 'invalid:mutant'
        print >> invalid_sort_file, "@%s::%s" % (read_id, type)
        print >> invalid_sort_file, read_seq
        continue

    read_seq, qual =cols[9:11]
    ref_id = ref_id.strip()    

    AS,NM,MD = None, None, None
    for i in range(11, len(cols)):
        col = cols[i]
        if col.startswith('AS'):
            AS = int(col[5:])
        elif col.startswith('NM'):
            NM = int(col[5:])
        elif col.startswith('MD'):
            MD = col[5:]

    # invalid: too large edit distance 
    if NM > mm_cutoff + cate_len['UMI']:
        type = 'invalid:mutant'
        print >> invalid_sort_file, "@%s::%s" % (read_id, type)
        print >> invalid_sort_file, read_seq
        continue

    # get read points in the sensitive window
    ref_pt, read_pt = pos-1, -1
    win_readpos = {}
    cigar_str=re.split('(\d+)',cigar_str)[1:]

    for i in range(len(cigar_str)/2):
        s = cigar_str[2*i+1]
        num = int(cigar_str[2*i])
        for j in range(num):
            if read_pt > 0  and ref_pt in sted_wins: # not include a hit on the first nt of reads
                for win in sted_wins[ref_pt]:
                    if win not in win_readpos:
                        win_readpos[win] = []
                    if s == 'D':
                        win_readpos[win].append(read_pt+1)
                    else:
                        win_readpos[win].append(read_pt)
            if s=='M':
                ref_pt += 1
                read_pt += 1
            elif s=='I':
                read_pt += 1
            elif s=='D':
                ref_pt += 1
    # check the last aligned base-pair
    if ref_pt in sted_wins: # include a hit on the last nt of reads
        for win in sted_wins[ref_pt]:
            if win not in win_readpos:
                win_readpos[win] = []
            win_readpos[win].append(read_pt)

    read_st, read_ed = pos, ref_pt

    # UMI is missing
    if read_st > UMI_st or read_ed < UMI_ed:
        assert 'UMI' not in win_readpos or len(win_readpos['UMI']) < 2
        defects.append('UMImissing')

    else:
        assert len(win_readpos['UMI']) == 2
        st, ed = sorted(win_readpos['UMI'])
        UMI = read_seq[st:ed]

        # wrong UMI length
        if len(UMI) != cate_len['UMI']:
            defects.append('UMIlengthError')
    
    # Widom 601 DNA is missing
    if 'Widom' not in win_readpos or read_ed  - win_readpos['Widom'][0] < 5:
        defects.append('NoWidom601')

    if len(defects) > 0:
        type = 'invalid:' + '/'.join(defects)
        print >> invalid_sort_file, "@%s::%s" % (read_id, type)
        print >> invalid_sort_file, read_seq
        continue

    else:
        type = 'UMI:' + UMI
        print >> valid_sort_file, "@%s::%s" % (read_id, type)
        print >> valid_sort_file, read_seq


    if spriteID not in spriteID_UMIs:
        spriteID_UMIs[spriteID] = set([])
    spriteID_UMIs[spriteID].add(UMI)
        
    valid_count +=1

invalid_sort_file.close()
valid_sort_file.close()
subprocess.call(["rm", temp_fname])
cat_cmd = ' '.join(["cat", valid_sort_fname, invalid_sort_fname, ">", sort_fname])
subprocess.call(cat_cmd, shell=True)
#subprocess.call(["rm", invalid_sort_fname])
#subprocess.call(["rm", valid_sort_fname])

print "total read count %d" % read_count
print "barcoded reads %.3f %%" % (float(100*barcoded_count)/read_count)
print "valid reads %.3f %%" % (float(100*valid_count)/read_count)
print

# check cluster size distribution
size_num = {}
for spriteID in spriteID_UMIs:
    size = len(spriteID_UMIs[spriteID])
    #if size > 1:
    #    print spriteID
    if size not in size_num:
        size_num[size] = 0
    size_num[size] +=1

print "Cluster size and number"
for size in sorted(size_num.keys()):
    print size, size_num[size]


# check data categories
def read_sort (fname):
    cate_info_count = {}
    for line in open(fname):
        if line.startswith('@'):
            cols = line.strip().split('::')
            cate, info = cols[-1].split(':')
            if cate not in cate_info_count:
                cate_info_count[cate] = {}
            if info not in cate_info_count[cate]:
                cate_info_count[cate][info] = 0
            cate_info_count[cate][info] +=1
    return cate_info_count

cate_info_count = read_sort(sort_fname)

# valid vs invalid pie chart
valid_count, invalid_count = 0, 0
for cate in cate_info_count:
    count = sum(cate_info_count[cate].values())
    if cate == 'UMI':
        valid_count +=count
    else:
        invalid_count +=count

fig = plt.figure()
colors = ['gold', 'yellowgreen', 'lightcoral', 'lightskyblue', 'lightpink', 'lime']
plt.pie([valid_count, invalid_count],colors=colors, labels=['valid', 'invalid'], shadow=True, startangle=90, autopct='%1.1f%%', explode=[0.05, 0])
plt.axis('equal')
#plt.show()
plt.close()

# invalid information pie chart
info_list, count_list = [], []
incplt_BC_count = 0
for info, count in cate_info_count['invalid'].items():
    if info.endswith('Barcode'):
        incplt_BC_count += cate_info_count['invalid'][info]
        continue
    info_list.append(info)
    count_list.append(count)
info_list.append('Incomplete barcoding')
count_list.append(incplt_BC_count)
assert sum(count_list) == sum(cate_info_count['invalid'].values())

fig = plt.figure()
colors = ['gold', 'yellowgreen', 'lightcoral', 'lightskyblue', 'lightpink', 'lime']
plt.pie(count_list, colors=colors, labels=info_list, shadow=True, startangle=90, autopct='%1.1f%%')
plt.axis('equal')
#plt.show()
plt.close()

# barcoding length bar chart
invalid_BClen_count = [0]*7
valid_BClen_count = [0]*7
info_BClen_count = {}
for cate in cate_info_count:
    if cate == 'UMI':
        BClen = 6
        valid_BClen_count[BClen] += sum(cate_info_count[cate].values())
        continue
    for info in cate_info_count[cate]:
        count = cate_info_count[cate][info]
        if info.endswith('Barcode'):
            BClen = int(info.split('-')[0])
            invalid_BClen_count[BClen] +=count
        else:
            BClen = 6
            #invalid_BClen_count[BClen] +=count
            if info not in info_BClen_count:
                info_BClen_count[info] = [0]*7
            info_BClen_count[info][BClen] +=count


fig = plt.figure()
plt.bar(range(len(invalid_BClen_count)), invalid_BClen_count, color='black')
total_info_list = sorted([[sum(info_BClen_count[info]), info] for info in info_BClen_count], reverse=True)
info_list = [ info for total, info in total_info_list ]
for i in range(len(info_list)):
    info = info_list[i]
    BClen_count = info_BClen_count[info]
    if i ==0:
        bottom = np.asarray([0]*len(BClen_count))
    else:
        bottom += np.asarray(info_BClen_count[info_list[i-1]])
    plt.bar(range(len(BClen_count)), BClen_count, bottom=bottom, label=info)
bottom += np.asarray(info_BClen_count[info_list[i]])
plt.bar(range(len(valid_BClen_count)), valid_BClen_count, bottom=bottom, label='Valid data', hatch='x', color='yellow')
plt.ylabel('Counts')
plt.xlabel('Barcode length')
plt.legend()
#plt.show()
plt.close()

                    
                
"""    
category_seq = {'TERM':'TATTATGGT', 'ODD':'N'*16+'T',EVEN:'N'*16+'T','UMI':'N'*10, 'DPM':'TCTTCCGATCTTGGGTGTTTT'}
cate_step_seq = {'ODDtoEVEN':'TGACTTG', 'EVENtoODD':'GACAACT', 'TERMtoEVEN':'TGACTTG', 'TERMtoODD':'GACAACT', 'ODDtoUMI':'TGACTTG'}

template = ""
format_list = Format.split('|')

for i in range(len(format_list)-1):
    cate, next_cate = format_list[i], format_list[i+1]
    cate_step = cate + 'to' + next_cate
    template += category_seq[cate]
    try:
        template += cate_step_seq[cate_step]
    except:
        pass
    template += category_seq[next_cate]


########
# TAGS #
########

# Columns:
# CATEGORY NAME SEQUENCE NUM_MISMATCHES

# Valid categories are "EVEN", "ODD", "Y", "RPM", "DPM", "LIGTAG"

EVEN	Even2Bo1	ATACTGCGGCTGACG	2
EVEN	Even2Bo5	CTAGGTGGCGGTCTG	2
EVEN	Even2Bo2	GTGACATTAAGGTTG	2
"""
