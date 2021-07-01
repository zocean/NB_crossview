#!/home/yangz6/Software/anaconda3/bin/python
# Programmer : Yang Zhang 
# Contact: zocean636@gmail.com
# Last-modified: 28 Jun 2021 05:17:06 PM

import os,sys,argparse
import gzip
import urllib
from bx.intervals.intersection import Interval, Intersecter
import logging

import falcon
from io import StringIO, BytesIO
import gunicorn.app.base
from gunicorn.six import iteritems
from falcon_cors import CORS

lis = ['http://vis.nucleome.org','https://vis.nucleome.org','http://x7.andrew.cmu.edu:8080', 'http://devd.io:8000', 'https://fiddle.jshell.net', 'http://127.0.0.1:8000', 'http://127.0.0.1:8080', '*']
cors = CORS(allow_credentials_origins_list = lis,allow_origins_list= lis)
#color_class = ['#bebada','#fb8072','#80b1d3','#fdb462','#b3de69','#fccde5','#8dd3c7','#ffffb3','#d9d9d9','#bc80bd','#ccebc5','#ffed6f', '#252525']
color_class = ['#cc0000', '#cc6600', '#cccc00', '#66cc00', '#00cc00', '#00cc66', '#00cccc', '#0066cc', '#0000cc', '#6600cc', '#cc00cc', '#cc0066', '#cc0000', '#252525']
DIS_CUTOFF = 1200

def parse_arg():
    ''' This Function Parse the Argument '''
    p=argparse.ArgumentParser( description = 'Example: %(prog)s -h', epilog='Library dependency :')
    p.add_argument('-v','--version',action='version',version='%(prog)s 0.1')
    p.add_argument('-c','--chain',type=str,dest="mapping",nargs="+",help="liftover chain mapping file list")
    p.add_argument('-m','--mode',type=str,dest="mode",nargs="+",help="mode, highlight or ortho")
    p.add_argument('-l','--label',type=str,dest="label",nargs="+",help="labels for each liftover chain file")
    p.add_argument('-p','--port',dest='port',default=8556,type=int,help="port Default: %(default)s")
    p.add_argument('-a',dest='a',default=False,type=bool,help="open to 0.0.0.0 : %(default)s") 
    if len(sys.argv) < 2:
        print(p.print_help())
        exit(1)
    return p.parse_args()

def nopen(f, mode="rb"):
    if not isinstance(f, str):
        return f
    if f.startswith("|"):
        p = Popen(f[1:], stdout=PIPE, stdin=PIPE, shell=True)
        if mode[0] == "r": return p.stdout
        return p
    return {"r": sys.stdin, "w": sys.stdout}[mode[0]] if f == "-" \
        else gzip.open(f, mode) if f.endswith((".gz", ".Z", ".z")) \
        else urllib.urlopen(f) if f.startswith(("http://", "https://","ftp://")) \
        else open(f, mode)

def reader(fname):
    for l in nopen(fname):
        yield l.decode('utf8').strip().replace("\r", "")

def update_chromID(c_temp, c_target):
    '''
    Update chromsome ID styles from 'c_target' to 'c_temp'.
    Parameters
    ----------
    c_temp : str
        Template of chromsome ID
    c_target : str
        Chromosome ID that need to be updated
    Returns
    --------
    Updated chromosome ID
    Examples
    --------
    >>> update_chromID('chrX',1)
    'chr1'
    >>> update_chromID('1','chrY')
    'Y'
    '''
    c_temp = str(c_temp)
    c_target = str(c_target)
    if c_temp.startswith('chr'):
        if c_target.startswith('chr'):
            return c_target
        else:
            return ('chr' + c_target)
    else:
        if c_target.startswith('chr'):
            return c_target.replace('chr','')
        else:
            return c_target

def intersectBed(lst1, lst2):
    '''
    Return intersection of two bed regions.
    
    Parameters
    ----------
    lst1 : list
         The 1st genomic region. List of chrom, start, end.
         Example: ['chr1',10, 100]
    
    lst2 : list
          The 2nd genomic region. List of chrom, start, end.
          Example: ['chr1',50, 120]
    
    Examples
    --------
    >>> intersectBed(['chr1',10, 100],['chr1',50, 120])
    ('chr1', 50, 100)
    >>> intersectBed(['chr1',10, 100],['chr1',20, 30])
    ('chr1', 20, 30)
    
    '''
    (chr1, st1, end1) = lst1
    (chr2, st2, end2) = lst2
    if int(st1) > int(end1) or int(st2) > int(end2):
        raise Exception ("Start cannot be larger than end")
    if chr1 != chr2:
        return None
    if int(st1) > int(end2) or int(end1) < int(st2):
        return None
    return (chr1, max(st1, st2), min(end1,end2))

def read_chain_file(chain_file, print_table = False):
    '''
    Read chain file.
    
    Parameters
    ----------
    chain_file : file
        Chain format file. Input chain_file could be either plain text, compressed file
        (".gz",".Z", ".z", ".bz", ".bz2", ".bzip2"), or a URL pointing to the chain file
        ("http://","https://", "ftp://"). If url was used, chain file must be plain text.
    
    print_table : bool, optional
         Print mappings in human readable table.
    
    Returns
    -------
    maps : dict
        Dictionary with source chrom name as key, IntervalTree object as value. An
        IntervalTree contains many intervals. An interval is a start and end position
        and a value. eg. Interval(11, 12, strand="-", value = "abc")
    
    target_chromSize : dict
        Chromosome sizes of target genome
    
    source_chromSize : dict
        Chromosome sizes of source genome
    '''
    
    logging.info("Read the chain file \"%s\" " % chain_file)
    maps={}
    target_chromSize={}
    source_chromSize={}
    if print_table:
        blocks=[]
    
    for line in reader(chain_file):
        # Example: chain 4900 chrY 58368225 + 25985403 25985638 chr5 151006098 - 43257292 43257528 1
        if not line.strip():
            continue
        line=line.strip()
        if line.startswith(('#',' ')):
            continue
        fields = line.split()

        if fields[0] == 'chain' and len(fields) in [12, 13]:
            #score = int(fields[1])            # Alignment score
            source_name = fields[2]            # E.g. chrY
            source_size = int(fields[3])  # Full length of the chromosome
            source_strand = fields[4]       # Must be +
            if source_strand != '+':
                 raise Exception("Source strand in a chain file must be +. (%s)" % line)
            source_start = int(fields[5]) # Start of source region
            #source_end = int(fields[6])       # End of source region

            target_name = fields[7]            # E.g. chr5
            target_size = int(fields[8])  # Full length of the chromosome
            target_strand = fields[9]       # + or -
            target_start = int(fields[10])
            #target_end = int(fields[11])
            target_chromSize[target_name]= target_size
            source_chromSize[source_name] = source_size

            if target_strand not in ['+', '-']:
                 raise Exception("Target strand must be - or +. (%s)" % line)
            #chain_id = None if len(fields) == 12 else fields[12]
            if source_name not in maps:
                 maps[source_name] = Intersecter()
            sfrom, tfrom = source_start, target_start
         # Now read the alignment chain from the file and store it as a list (source_from, source_to) -> (target_from, target_to)

        elif fields[0] != 'chain' and len(fields) == 3:
             size, sgap, tgap = int(fields[0]), int(fields[1]), int(fields[2])
             if print_table:
                 if target_strand == '+': blocks.append((source_name,sfrom, sfrom+size, source_strand, target_name, tfrom, tfrom+size, target_strand))
                 elif  target_strand == '-': blocks.append((source_name,sfrom, sfrom+size, source_strand, target_name, target_size - (tfrom+size), target_size - tfrom, target_strand))
             if target_strand == '+':
                 maps[source_name].add_interval( Interval(sfrom, sfrom+size,(target_name,tfrom, tfrom+size,target_strand)))
             elif target_strand == '-':
                 maps[source_name].add_interval( Interval(sfrom, sfrom+size,(target_name,target_size - (tfrom+size), target_size - tfrom, target_strand)))
             sfrom += size + sgap
             tfrom += size + tgap
        elif fields[0] != 'chain' and len(fields) == 1:
            size = int(fields[0])
            if print_table:
                if target_strand == '+': blocks.append((source_name,sfrom, sfrom+size, source_strand, target_name, tfrom, tfrom+size, target_strand))
                elif  target_strand == '-': blocks.append((source_name,sfrom, sfrom+size, source_strand, target_name, target_size - (tfrom+size), target_size - tfrom, target_strand))
            if target_strand == '+':
                 maps[source_name].add_interval( Interval(sfrom, sfrom+size,(target_name,tfrom, tfrom+size,target_strand)))
            elif target_strand == '-':
                 maps[source_name].add_interval( Interval(sfrom, sfrom+size,(target_name,target_size - (tfrom+size), target_size - tfrom, target_strand)))
        else:
            raise Exception("Invalid chain format. (%s)" % line)
    #if (sfrom + size) != source_end  or (tfrom + size) != target_end:
    #      raise Exception("Alignment blocks do not match specified block sizes. (%s)" % header)
    if print_table:
        for i in blocks:
            print('\t'.join([str(n) for n in i]))
    return (maps, target_chromSize, source_chromSize)

def map_coordinates(mapping, q_chr, q_start, q_end, q_strand = '+', print_match = False):
    '''
    Map coordinates from source (i.e. original) assembly to target (i.e. new) assembly.
    
    Parameters
    ----------
    mapping : dict
         Dictionary with source chrom name as key, IntervalTree object as value.
    
    q_chr : str
         Chromosome ID of query interval
    
    q_start : int
         Start position of query interval.
    
    q_end : int
         End position of query interval.
    
    q_strand : str
         Strand of query interval.
    
    print_match : bool
         Print match table.
    '''
    
    matches = []
    complement = {'+':'-','-':'+'}
    
    if q_chr in mapping:
        targets = mapping[q_chr].find(q_start, q_end)
    elif q_chr.replace('chr','') in mapping:
        targets = mapping[q_chr.replace('chr','')].find(q_start, q_end)
    elif ('chr' + q_chr) in mapping:
        targets = mapping['chr' + q_chr].find(q_start, q_end)
    else:
        return None
    if len(targets)==0:
        return None
    elif len(targets)==1:
        s_start = targets[0].start
        s_end = targets[0].end
        t_chrom = targets[0].value[0]
        t_chrom = update_chromID(q_chr, t_chrom)
        t_start = targets[0].value[1]
        t_end = targets[0].value[2]
        t_strand = targets[0].value[3]
    
        (chr, real_start, real_end)      = intersectBed((q_chr,q_start,q_end),(q_chr,s_start,s_end))
        l_offset = abs(real_start - s_start)
        #r_offset = real_end - s_end
        size = abs(real_end - real_start)
    
        matches.append( (chr, real_start, real_end,q_strand))
        if t_strand == '+':
            i_start = t_start + l_offset
            if q_strand == '+':
                matches.append( (t_chrom, i_start, i_start + size, t_strand))
            else:
                matches.append( (t_chrom, i_start, i_start + size, complement[t_strand]))
        elif t_strand == '-':
            i_start = t_end - l_offset - size
            if q_strand == '+':
                matches.append( (t_chrom, i_start,     i_start + size, t_strand))
            else:
                matches.append( (t_chrom, i_start,     i_start + size, complement[t_strand]))
        else:
             raise Exception("Unknown strand: %s. Can only be '+' or '-'." % q_strand)
    elif len(targets) > 1:
        for t in targets:
            s_start = t.start
            s_end = t.end
            t_chrom = t.value[0]
            t_chrom = update_chromID(q_chr, t_chrom)
            t_start = t.value[1]
            t_end = t.value[2]
            t_strand = t.value[3]
    
            (chr, real_start, real_end)      = intersectBed((q_chr,q_start,q_end),(q_chr,s_start,s_end))
    
            l_offset = abs(real_start - s_start)
            #r_offset = abs(real_end - s_end)
            size = abs(real_end - real_start)
            matches.append( (chr, real_start, real_end,q_strand) )
            if t_strand == '+':
                i_start = t_start + l_offset
                if q_strand == '+':
                    matches.append( (t_chrom, i_start, i_start + size, t_strand))
                else:
                    matches.append( (t_chrom, i_start, i_start + size, complement[t_strand]))
            elif t_strand == '-':
                i_start = t_end - l_offset - size
                if q_strand == '+':
                    matches.append( (t_chrom, i_start,     i_start + size, t_strand))
                else:
                    matches.append( (t_chrom, i_start,     i_start + size, complement[t_strand]))
            else:
                raise Exception("Unknown strand: %s. Can only be '+' or '-'." % q_strand)
    
    if print_match:
        print(matches)
        # input: 'chr1',246974830,247024835
        # output: [('chr1', 246974830, 246974833, '+' ), ('chr1', 248908207, 248908210, '+' ), ('chr1', 247024833, 247024835, '+'), ('chr1', 249058210, 249058212,'+')]
        # [('chr1', 246974830, 246974833), ('chr1', 248908207, 248908210)]
    
    return matches

def crossmap_region_file(mapping, inbed, outfile=None, min_ratio = 0.95):
    '''
    Convert large genomic regions (in bed format) between assemblies.
    BED format: http://genome.ucsc.edu/FAQ/FAQformat.html#format1
    
    Parameters
    ----------
    mapping : dict
         Dictionary with source chrom name as key, IntervalTree object as value.
    
    inbed : file
         Input BED file.
    
    outfile : str, optional
         Prefix of output files.
    
    min_ratio : float, optional
         Minimum ratio of query bases that must remap
    
    '''
    
    # check if 'outfile' was set. If not set, print to screen, if set, print to file
    if outfile is not None:
        FILE_OUT = open(outfile,'w')
        UNMAP = open(outfile + '.unmap','w')
    else:
        pass
    
    for line in reader(inbed):
        if line.startswith(('#', 'track','browser')):continue
        if not line.strip():continue
        line=line.strip()
        fields=line.split()
        strand = '+'
    
        # filter out line less than 3 columns
        if len(fields)<3:
            print("Less than 3 fields. skip " + line, file=sys.stderr)
            if outfile:
                print(line + '\tInvalidBedFormat', file=UNMAP)
            continue
        try:
            int(fields[1])
        except:
            print("Start coordinate is not an integer. skip " + line, file=sys.stderr)
            if outfile:
                print(line + '\tInvalidStartPosition', file=UNMAP)
            continue
        try:
            int(fields[2])
        except:
            print("End coordinate is not an integer. skip " + line, file=sys.stderr)
            if outfile:
                print(line + '\tInvalidEndPosition', file=UNMAP)
            continue
        if int(fields[1]) > int(fields[2]):
            print("\"Start\" is larger than \"End\" coordinate is not an integer. skip " + line, file=sys.stderr)
            if outfile:
                print(line + '\tStart>End', file=UNMAP)
            continue
        # try to reset strand
        try:
            for f in fields:
                if f in ['+','-']:
                    strand = f
        except:
            pass
    
        chrom = fields[0]
        start = int(fields[1])
        end = int(fields[2])
        total_query_length = end - start     #used to calculate q_map_ratio
    
        a = map_coordinates(mapping, chrom, start, end, strand)
        # input: 'chr1',246974830,247024835
        # output: [('chr1', 246974830, 246974833, '+' ), ('chr1', 248908207, 248908210, '+' ), ('chr1', 247024833, 247024835, '+'), ('chr1', 249058210, 249058212,'+')]
        # [('chr1', 246974830, 246974833), ('chr1', 248908207, 248908210)]
    
        if (a is None) or (len(a) % 2 != 0):
            if outfile is None:
                print(line + '\tFail\tUnmap')
            else:
                print(line + '\tFail\tUnmap', file=UNMAP)
            continue
    
        #when a == 2, there is one-to-one match (i.e. 100% match)
        if len(a) == 2:
            #reset fields to target assembly
            fields[0] =  a[1][0]
            fields[1] =  a[1][1]
            fields[2] =  a[1][2]
            for i in range(0,len(fields)):     #update the strand information
                if fields[i] in ['+','-']:
                    fields[i] = a[1][3]
            if outfile is None:
                print(line + '\t->\t' + '\t'.join([str(i) for i in fields]) + "\tmap_ratio=1.0000")
            else:
                print('\t'.join([str(i) for i in fields]) + "\tmap_ratio=1.0000", file=FILE_OUT)
        #when a is an even number but bigger than 2, each segment is 100% match,
        # but the whole region is not. In this case, check *min_ratio* of the query
        if len(a) > 2 :
            a_query =  a[::2] #EVEN: [('chr1', 246974830, 246974833, '+'), ('chr1', 247024833, 247024835, '+')]
            a_query_mapped_nt = sum([i[2]-i[1] for i in a_query]) #sum([3,2])
            a_target = a[1::2] #ODDS: [('chr1', 248908207, 248908210, '+'), ('chr1', 249058210, 249058212, '+')]
            a_target_chroms = set([i[0] for i in a_target])
            a_target_chroms = set([i[0] for i in a_target])
            a_target_starts = [i[1] for i in a_target]
            a_target_ends = [i[2] for i in a_target]
            #print (a_target_ends)
            map_ratio = a_query_mapped_nt/total_query_length
    
            #map_ratio > cutoff
            if map_ratio >= min_ratio:
                if len(a_target_chroms) == 1:
                    t_chrom = a_target_chroms.pop()
                    fields[0] = t_chrom
                    fields[1] = min(a_target_starts)
                    fields[2] = max(a_target_ends)
                    if outfile is None:
                        print(line + '\t->\t' + '\t'.join([str(i) for i in fields]) + ("\tmap_ratio=%.4f" % map_ratio))
                    else:
                        print('\t'.join([str(i) for i in fields]) + ("\tmap_ratio=%.4f" % map_ratio), file=FILE_OUT)
                else:
                    if outfile is None: print(line + '\tFail\tCrossChroms')
                    else: print(line + '\tFail\tCrossChroms', file=UNMAP)
            # map_ratio > 0 but < cutoff
            elif map_ratio >0 and map_ratio < min_ratio:
                if outfile is None: print(line + '\tFail' + ("\tmap_ratio=%.4f" % map_ratio))
                else: print(line + '\tFail' + ("\tmap_ratio=%.4f" % map_ratio), file=UNMAP)

def crossmap_region_io(mapping, region_list, min_ratio = 0.95):
    '''
    Convert large genomic regions (Python list) between assemblies.
    BED format: http://genome.ucsc.edu/FAQ/FAQformat.html#format1
    
    Parameters
    ----------
    mapping : dict
         Dictionary with source chrom name as key, IntervalTree object as value.
    
    inbed : file
         Input list of genomic coordinates in chrom:start-end format
    
    min_ratio : float, optional
         Minimum ratio of query bases that must remap
    
    '''
    
    # save liftovered region into list
    out = []
    
    for region in region_list:
        # always use '+' strand
        strand = '+'
        if region == '':
            continue
        try:
            chrom,start,end = region.replace(':','-').split('-')
        except:
            continue
        start = int(start)
        end = int(end)
        total_query_length = end - start     #used to calculate q_map_ratio
        # if start is larger then end, skip it
        if start > end:
            continue
        # liftover
        a = map_coordinates(mapping, chrom, start, end, strand)
        # if no results found, skip it
        if (a is None) or (len(a) % 2 != 0):
            continue
        #when a == 2, there is one-to-one match (i.e. 100% match)
        fields = [None,None,None,None]
        if len(a) == 2:
            #reset fields to target assembly
            fields[0] =  a[1][0]
            fields[1] =  a[1][1]
            fields[2] =  a[1][2]
            for i in range(0,len(fields)):     #update the strand information
                if fields[i] in ['+','-']:
                    fields[i] = a[1][3]
            #map_ratio = 1.0
            out.append((fields[0], fields[1], fields[2]))
        #when a is an even number but bigger than 2, each segment is 100% match,
        # but the whole region is not. In this case, check *min_ratio* of the query
        if len(a) > 2 :
            a_query =  a[::2] #EVEN: [('chr1', 246974830, 246974833, '+'), ('chr1', 247024833, 247024835, '+')]
            a_query_mapped_nt = sum([i[2]-i[1] for i in a_query]) #sum([3,2])
            a_target = a[1::2] #ODDS: [('chr1', 248908207, 248908210, '+'), ('chr1', 249058210, 249058212, '+')]
            a_target_chroms = set([i[0] for i in a_target])
            a_target_chroms = set([i[0] for i in a_target])
            a_target_starts = [i[1] for i in a_target]
            a_target_ends = [i[2] for i in a_target]
            map_ratio = a_query_mapped_nt/total_query_length
            # map_ratio >= cutoff
            if map_ratio >= min_ratio:
                if len(a_target_chroms) == 1:
                    t_chrom = a_target_chroms.pop()
                    fields[0] = t_chrom
                    fields[1] = min(a_target_starts)
                    fields[2] = max(a_target_ends)
                    out.append((fields[0], fields[1], fields[2]))
                else:
                    continue
            # map_ratio > 0 but < cutoff
            elif map_ratio >0 and map_ratio < min_ratio:
                continue
    return out

def crossmap_ortho_io(mapping, region_list):
    '''
    Identify orthologous region pairs
    BED format: http://genome.ucsc.edu/FAQ/FAQformat.html#format1
    
    Parameters
    ----------
    mapping : dict
         Dictionary with source chrom name as key, IntervalTree object as value.
    
    inbed : file
         Input list of genomic coordinates in chrom:start-end format
     
    '''
    
    # save liftovered region into list
    out = []
    
    # check total length of reference region(s)
    # if exceed 1000000bp do nothing
    total_len = 0
    for region in region_list:
        try:
            chrom,start,end = region.replace(':','-').split('-')
        except:
            continue
        start = int(start)
        end = int(end)
        total_len += end - start
        if total_len > 1000000:
            return out

    for region in region_list:
        # always use '+' strand
        strand = '+'
        if region == '':
            continue
        try:
            chrom,start,end = region.replace(':','-').split('-')
        except:
            continue
        start = int(start)
        end = int(end)
        # if start is larger then end, skip it
        if start > end:
            continue
        # liftover
        a = map_coordinates(mapping, chrom, start, end, strand)
        # if no results found, skip it
        if (a is None) or (len(a) % 2 != 0):
            continue
        #when a == 2, there is one-to-one match (i.e. 100% match)
        if len(a) == 2:
            #reset fields to target assembly
            out.append([a[0], a[1], color_class[0]])
        #when a is an even number but bigger than 2, each segment is 100% match,
        # but the whole region is not. In this case, check *min_ratio* of the query
        if len(a) > 2 :
            count = 0
            block_count = 0
            last_r = [None, None, None, None]
            last_t = [None, None, None, None]
            # assume a is sorted by reference chrom then by start cooridnate
            for i in range(1, len(a), 2):
                if block_count == 14: # re-use color palette
                    count = 0
                    block_count = 0
                ref = a[i-1]
                target = a[i]
                if last_r[0] is None: # init
                    last_r = list(ref)
                    last_t = list(target)
                elif last_r[0] != ref[0] or last_t[0] != target[0]: # different chromosome => new block
                    out.append([ last_r, last_t, color_class[block_count] ])
                    last_r = list(ref)
                    last_t = list(target)
                    block_count += 1
                else: # same chromosome check difference
                    # check the distance between last block and the current block
                    dis_r = min_dis(last_r, ref)
                    dis_t = min_dis(last_t, target)
                    if abs(dis_r) < DIS_CUTOFF and abs(dis_t) < DIS_CUTOFF: # gaps between blocks are smaller then CUTOFF => merge
                        if dis_r > -0.1:
                            last_r[2] = ref[2]
                        else:
                            last_r[1] = ref[1]
                        if dis_t > 0:
                            last_t[2] = target[2]
                        elif dis_t < 0:
                            last_t[1] = target[1]
                        else:
                            last_t[1] = min(last_t[1], target[1])
                            last_t[2] = max(last_t[2], target[2])
                    else: # => new block
                        out.append([ last_r, last_t, color_class[block_count] ])
                        last_r = list(ref)
                        last_t = list(target)
                        block_count += 1
            # put the last item into 
            out.append([ last_r, last_t, color_class[block_count] ])
    return out

def min_dis(last, cur):
    """Calculate the minimal distance between the last region and the current region"""
    if cur[1] > last[2] - 0.1: # current is in the downstream of last
        return cur[1] - last[2] # positive value
    elif last[1] > cur[2] - 0.1: # current is in the upstream of last
        return cur[2] - last[1] # negative value
    else: # overlapping
        return 0

class StandaloneApplication(gunicorn.app.base.BaseApplication):
    def __init__(self, app, options=None):
        self.options = options or {}
        self.application = app
        super(StandaloneApplication, self).__init__()
    def load_config(self):
        config = dict([(key, value) for key, value in self.options.items()
                       if key in self.cfg.settings and value is not None])
        for key, value in config.items():
            self.cfg.set(key.lower(), value)
    def load(self):
        return self.application

class MappingBed(object):
    """Respond to liftover reqeust"""
    def __init__(self, mapping_table):
        self.mapping = {}
        for label in mapping_table.keys():
            chain,mode = mapping_table[label]
            if mode == 'highlight':
                self.mapping[label] = chain
    def on_get(self, req, resp):
        """Handles GET requests"""
        resp.status = falcon.HTTP_200 # This is the default status
        query_bed = req.params.get('query', None)
        label = req.params.get('label',None)
        buf = StringIO()
        if query_bed is None or label is None:
            buf.write("Need to provide query region(s) (e.g., query=chr1:1000-2000,chr2:300-500&label=hg19Tohg38")
        elif self.mapping.get(label, None) is None:
            buf.write("Can't find chain file for %s" % (label))
        else:
            result_list = crossmap_region_io(self.mapping[label], query_bed.split(','))
            for region in result_list:
                buf.write("%s\t%d\t%d\n" % (region[0], region[1], region[2]))
        resp.body = (buf.getvalue())

class MappingOrtho(object):
    """Respond to orthologous block request"""
    def __init__(self, mapping_table):
        self.mapping = {}
        for label in mapping_table.keys():
            chain, mode = mapping_table[label]
            if mode == 'ortho':
                self.mapping[label] = chain
    def on_get(self, req, resp):
        """Handles GET requests"""
        resp.status = falcon.HTTP_200 # This is the default status
        query_bed = req.params.get('query', None)
        label = req.params.get('label',None)
        buf = StringIO()
        if query_bed is None or label is None:
            buf.write("Need to provide query region(s), (e.g., query=chr1:1000-2000,chr2:300-500")
        elif self.mapping.get(label, None) is None:
            buf.write("Can't find chain file for %s" % (label))
        else:
            result_list = crossmap_ortho_io(self.mapping[label], query_bed.split(','))
            if len(result_list) > 0:
                for ref,target,color in result_list:
                    buf.write("%s\t%d\t%d\t%s\t%d\t%d\t%s\n" % (ref[0], ref[1], ref[2], target[0], target[1], target[2], color))
        resp.body = (buf.getvalue())

def generateApp(mapping_filename):
    # load liftover chain file
    # check number of chain file and labels
    assert len(args.mapping) == len(args.label) and len(args.mapping) == len(args.mode)
    mapping_table = {}
    for chain_file, label, mode in zip(args.mapping, args.label, args.mode):
        print("# load chain file: %s as %s for %s" % (chain_file, label, mode))
        mapping, target_chromSize, source_chromSize = read_chain_file(chain_file)
        mapping_table[label] = (mapping, mode)
    print("read chain file done")
    # init app
    app = falcon.API(middleware = [cors.middleware])
    # parse info
    serveLiftover = MappingBed(mapping_table) 
    serveOrtho = MappingOrtho(mapping_table)
    # add route
    app.add_route('/highlight', serveLiftover)
    app.add_route('/ortho', serveOrtho)
    return app

def run_app(args):
    app = generateApp(args.mapping)
    bindserver = "127.0.0.1"
    if args.a:
        bindserver = "0.0.0.0"
    options = {
        'bind': '%s:%d' % (bindserver, args.port),
        'workers': 1,
    }
    StandaloneApplication(app, options).run()

def main():
    global args
    args = parse_arg()
    run_app(args)
 
if __name__=="__main__":
    main()
