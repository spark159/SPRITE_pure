import os, sys, subprocess, re
from argparse import ArgumentParser, FileType

def rev_cmp (seq):
    dic={'A':'T', 'T':'A', 'C':'G', 'G':'C', 'N':'N'}
    output=''
    for nt in seq:
        output+=dic[nt]
    return output[::-1]

# Make a template for reads alignment
Format = 'TERM|EVEN|ODD|EVEN|ODD|EVEN|ODD|UMI|DPM'

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

# Reads alignment on to the template
aligner_cmd=["bowtie2", '-x', ref_fname, '-U', read_fname ]
#aligner_cmd += ['--n-ceil', 'L,'+str(insert_len)+','+str(0.15)]
align_proc = subprocess.Popen(aligner_cmd, stdout=subprocess.PIPE,stderr=open("/dev/null", 'w'))

# start sort the reads
seq_sort_fname=open(out_fname+".sort",'w')
#seq_sort={}

read_count=0

for line in align_proc.stdout:
    if line.startswith('@'):
        continue
    #print line
    type, cut_loc, read_seq = 'NA', 'NA', 'NA'

    cols = line.strip().split()
    read_id, flag, ref_id, pos, mapQ, cigar_str = cols[:6]
    read_id=":".join(read_id.split(':')[3:7])
    read_count +=1
    flag, pos = int(flag), int(pos)
    pos-=1

    # invalid: mapping failure
    if pos < 0:
        type = 'invalid:mutant'
        #continue
    if flag & 0x4 != 0:
        type = 'invalid:mutant'
        #continue

    # invalid: ambiguous mapping
    #mapQ = float(mapQ)
    #if mapQ < 10:
    #    type = 'invalid:multimap'

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
    if NM > mm_cutoff:
        type = 'invalid:mutant'
        #continue

    # collect invalid data 
    if type != 'NA':
        assert type.startswith('invalid')
        #assert read_id not in seq_sort
        #seq_sort[read_id]=[type, insert, cut_loc, read_seq]
        if small_sort:
            continue
        print >> seq_sort_fname, "%s\t%s\t%s\t%s\t%s" % (read_id, type, ref_id, cut_loc, read_seq)
        continue

    if check_window:
        win, win_st = ref_id.split('-')
        win_st, win_size = int(win_st), len(win)
        win_ed = win_st + win_size
        if discard_multHit: # look all windows with same size
            sted_IDs = size_sted_IDs[win_size]
        else:
            sted_IDs = {} # look only aligned window
            sted_IDs[win_st] = [ref_id]
            sted_IDs[win_ed] = [ref_id]
    else: # no screen window
        sted_IDs = {}


    # find mismatch information
    #MD = re.findall('\d+|[A-Z]|\^[A-Z]+', MD)
    #pt, mismatch = 0, {}
    #for tag in MD:
    #    if re.search('\d+', tag):
    #        pt += int(tag)
    #    elif tag.startswith('^'):
    #        pt += len(tag[1:])
    #    else:
    #        mismatch[pt] = tag
    #        pt += len(tag)

    # get read points in the sensitive window
    ref_pt, read_pt = pos-1, -1
    ID_readpos = {}
    cigar_str=re.split('(\d+)',cigar_str)[1:]

    for i in range(len(cigar_str)/2):
        s = cigar_str[2*i+1]
        num = int(cigar_str[2*i])
        for j in range(num):
            if read_pt > 0  and ref_pt in sted_IDs: # not include a hit on the first nt of reads
                for ID in sted_IDs[ref_pt]:
                    if ID not in ID_readpos:
                        ID_readpos[ID] = []
                    if s == 'D':
                        ID_readpos[ID].append(read_pt+1) # account offset for D
                    else:
                        ID_readpos[ID].append(read_pt)
            if s=='M':
                ref_pt += 1
                read_pt += 1
            elif s=='I':
                read_pt += 1
            elif s=='D':
                ref_pt += 1
    # check the last aligned base-pair
    if ref_pt in sted_IDs: # include a hit on the last nt of reads
        for ID in sted_IDs[ref_pt]:
            if ID not in ID_readpos:
                ID_readpos[ID] = []
            ID_readpos[ID].append(read_pt)

    # sort by read length and alignment position
    ref_len = ref_length[ref_id]
    end_pos = min(ref_pt, ref_len - 1)        

    if pos < len_cutoff and end_pos > ref_len - len_cutoff - 1:
        type = 'freeDNA'
    elif pos < len_cutoff and end_pos <= ref_len - len_cutoff - 1: 
        cut_loc = 'L:' + str(end_pos)
    elif pos >= len_cutoff and end_pos > ref_len - len_cutoff - 1:
        cut_loc = 'R:' + str(pos)
    else:
        type = 'invalid:frag'

    #print check_window
    #print ID_readpos.keys()
    # screen the read sequences in windows
    if check_window:
        candidates = []
        for ID in ID_readpos:
            try:
                st, ed = sorted(ID_readpos[ID])
                insert = read_seq[st:ed]
                if len(insert) != win_size:
                    continue
                key = insert + '-' + ID.split('-')[1]
                if key in ref_win and key not in candidates:
                    candidates.append(key)
                #if len(candidates) > 1:
                #    break
            except:
                pass
        #print ref_id
        #print candidates
        #print 
        if len(candidates) <=0:
            type = 'invalid:instErr:' + 'NA'
        else:                
            if ref_id not in candidates: # check bowtie alignment
                type = 'invalid:instErr:'
                for insert in candidates:
                    type += ':' + insert
            elif ref_id in candidates and len(candidates) > 1: # check uniqueness
                type = 'invalid:multHit:'
                for insert in candidates:
                    type += ':' + insert

    # otherwise, valid data
    if type == 'NA':
        type = 'valid'

    #assert read_id not in seq_sort
    #seq_sort[read_id]=[type, insert, cut_loc, read_seq]
    if small_sort:
        if type !='valid':
            continue
        else:
            read_seq = 'N'

    print >> seq_sort_fname, "%s\t%s\t%s\t%s\t%s" % (read_id, type, ref_id, cut_loc, read_seq)


########
# TAGS #
########

# Columns:
# CATEGORY NAME SEQUENCE NUM_MISMATCHES

# Valid categories are "EVEN", "ODD", "Y", "RPM", "DPM", "LIGTAG"

EVEN	Even2Bo1	ATACTGCGGCTGACG	2
EVEN	Even2Bo5	CTAGGTGGCGGTCTG	2
EVEN	Even2Bo2	GTGACATTAAGGTTG	2
