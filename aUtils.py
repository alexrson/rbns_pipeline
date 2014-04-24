import os
import ConfigParser
import time
import simplejson
import sys
import math
import operator
from collections import defaultdict
import random
from itertools import imap, tee, izip
import itertools
import subprocess

import numpy as np
from scipy import stats
import tabix

from table2dict import table2dict
import readAttributes
import aopen


def read_settings(settings_file):
    config = ConfigParser.ConfigParser()
    config.read(settings_file)
    settings = {}
    for section in config.sections():
        for option in config.options(section):
            settings[option] = config.get(section, option)
            settings[section] = True
    return settings


def bin_of_list(things, bns_bin, n_bins):
    return things[int(len(things) * bns_bin / n_bins):int(len(things) * (bns_bin + 1) / n_bins)]


def list_file_type(in_dir, suffix):
    
    files = os.listdir(in_dir)
    for f in files:
        if f.endswith(suffix):
            yield os.path.join(in_dir, f)
    

def monitor_percent(i, n, msg=''):
    percentage = 100 * i / n
    last_percentage = 100 * (i - 1) / n
    if percentage != last_percentage:
        print ' '.join([msg, str(percentage) + '%'])

def check_file_exists(fname):
    """
    makes sure a given file exists
    """
    if not os.path.exists(fname):
        return False
    fstats = os.stat(fname)
    if not fstats[6]:
        return False
    if not os.access(fname, os.R_OK):
        raise ValueError('Input File %s cannot be read' % fname)
    return True


# returns the sub-dictionary using the keys in keys
def extract_subdict(keys, d):
    return dict((k, d[k]) for k in keys if k in d)

def process_settings(settings_file):

    """
    reads a settings file (json)
    and returns a dict() of the info
    """
    config = ConfigParser.ConfigParser()
    config.read(settings_file)
    settings = {}
    for section in config.sections():
        for option in config.options(section):
            settings[option] = config.get(section, option)
    return settings


def launch(command, script_options=False, ppn='1', q='short', jobname='ajob', out_file=''):
    """
    launches the command on the cluster using the given script options
    """
    if not script_options:
        script_options = {'nodes': '1', 'ppn': str(ppn),
          'outf': os.path.expanduser('~/submission_log'),
          'jobname': jobname,
          'queue': q, 'workingdir': os.getcwd()}
    script_options['command'] = command
    #join_line = "#PBS -j oe"
    outtext = """#!/bin/bash
    #PBS -l nodes=%(nodes)s:ppn=%(ppn)s
    #PBS -m a
    #PBS -e /home/alexrson/outs/%(jobname)s.zerr
    #PBS -o /home/alexrson/outs/%(jobname)s.zout
    #PBS -M alexrson@mit.edu
    #PBS -N %(jobname)s
    #PBS -q %(queue)s
    #PBS -S /bin/bash
    echo $HOSTNAME
    echo Working directory is %(workingdir)s
    cd %(workingdir)s
    %(command)s
    echo "===== command finished =====" """ % script_options
    call = "qsub -"
    qsub = subprocess.Popen(call, shell=True,
      stdout=subprocess.PIPE,
      stderr=subprocess.PIPE, stdin=subprocess.PIPE)
    qsub.stdin.write(outtext)
    if out_file:
        open(out_file, 'w').write(outtext)
    output = qsub.communicate()
    if output[0].strip().endswith('.coyote.mit.edu'):
        job_id = int(output[0].split('.')[0])
        print job_id
        return job_id
    else:
        print output
        print 'Failure to launch'
        raise ValueError('Failure to launch')


def launch_many(commands, num_jobs=500, ppn='1', q='short', job_file_dir='', job_name_base='job'):
    num_jobs = min(len(commands), num_jobs)
    job_scripts = ['cd %s' % os.getcwd() for i in range(num_jobs)]
    for command_i, command in enumerate(commands):
        job_scripts[command_i % num_jobs] =\
          job_scripts[command_i % num_jobs] + '\n' + command
    job_list = []
    for i, job_script in enumerate(job_scripts):
        if job_script:
            if job_file_dir:
                #open('%s/%s.%i.sh' %
                #  (job_file_dir, job_name_base, i), 'w').write(job_script)
                out_file = '%s/%s.%i.sh' % (job_file_dir, job_name_base, i)
            else:
                out_file = ''
            job_list.append(launch(job_script,
                                   out_file=out_file,
                                   jobname='%s.%i' % (job_name_base, i),
                                   ppn=ppn,
                                   q=q))
    return job_list


def within_xpercent(a, b, max_percent=1.0):
    if not (a > 0 and b > 0):
        return False
    ratio = float(max(a, b)) / float(min(a, b))
    percent_increase = (ratio - 1.0) * 100.0
    return percent_increase < max_percent


def smooth_average(l):
    newl = [0] * len(l)
    newl[0] = l[0]
    newl[-1] = l[-1]
    for i in range(1, len(l) - 1):
        newl[i] = (l[i - 1] + l[i] + l[i + 1]) / 3
    return newl
    

def log10(f):
    f = float(f)
    if not f:
        return -5
    try:
        return math.log10(f)
    except:
        return -5

def significantly_enriched_OR(xs, ys, zthresh=2., scale='linear'):
    #assert scale in ['linear', 'log']
    #if scale =='log':
    #    xs = np.log2(xs)
    #    ys = np.log2(ys)
    #xs = stats.zscore(xs)
    #ys = stats.zscore(ys)
    #return [x > zthresh or y > zthresh for x, y in zip(xs, ys)]
    test = [xe or ye for xe,ye in zip(significantly_enriched(xs, zthresh, scale), significantly_enriched(ys, zthresh, scale))]
    test2 = np.array(significantly_enriched(xs, zthresh, scale)) + np.array(significantly_enriched(ys, zthresh, scale))
    for t1, t2 in zip(test,test2):
        assert t1==t2
    return np.array(significantly_enriched(xs, zthresh, scale)) + np.array(significantly_enriched(ys, zthresh, scale))


def significantly_enriched(xs, zthresh=2., scale='linear'):
    assert scale in ['linear', 'log']
    if scale =='log':
        xs = np.log2(xs)
    xs = stats.zscore(xs)
    return [x > zthresh for x in xs]


def significantly_unenriched(xs, ys, zthresh=2., scale='linear'):
    assert scale in ['linear', 'log']
    if scale =='log':
        xs = np.log2(xs)
        ys = np.log2(ys)
    xs = stats.zscore(xs)
    ys = stats.zscore(ys)
    return [x < -zthresh or y < -zthresh for x, y in zip(xs, ys)]


def wait_for_jobs(job_list, sleep=10):
    print len(job_list)
    while True:
        outputs = [subprocess.Popen('qstat %i' % job_id,
          shell=True, stdout=subprocess.PIPE,
          stderr=subprocess.PIPE).communicate()
          for job_id in job_list]
        completed = map(lambda output: 'Unknown Job' in output[1], outputs)
        if all(completed):
            break
        else:
            num_left = sum(map(lambda output:
              0 if 'Unknown J' in output[1] else 1, outputs))
            print num_left, 'jobs left'
            if not num_left:
                print completed
        time.sleep(sleep)
    

def fold_a_seq(seq):
    """
    Returns the folding energy for a sequence
    THis is much slower than doing it in a batch.
    """
    command = 'RNAfold'
    job = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE,
      stderr=subprocess.PIPE, stdin=subprocess.PIPE)
    job.stdin.write(seq)
    output = job.communicate()
    return float(output[0].replace('(', ' ').replace(')', ' ').split()[-1])


def get_adjacent_kmers(kmer):
    """
    returns all the k+1 mers that contain the kmer
    """
    return ['A' + kmer, 'C' + kmer, 'G' + kmer, 'T' + kmer,
      kmer + 'A', kmer + 'C', kmer + 'G', kmer + 'T']


def setToString(s, separator='\n'):
    return separator.join(s)


def remove_lambda(dd):
    pass


def iter_motif(seq, motif):
    if len(seq) - 10 < len(motif):
        return
    for start in range(0, len(seq) - len(motif) - 4, 1):
        if seq[start:start + len(motif)] == motif:
            yield start


def iter_kmers(k):
    for kmer in itertools.product('ATGC', repeat=k):
        yield ''.join(kmer)


def iterLinePairs(inFile):
    for l1, l2 in iterNlines(inFile, 2):
        yield l1, l2


def iter4Lines(inFile):
    return iterNlines(inFile, 4)


def iter_uncommented(inFile):
    for line in open(inFile):
        if not line.startswith('#'):
            yield line


def iterLines_skip_first(inFile):
    for i, line in enumerate(open(inFile)):
        if i:
            yield line

def findall(sub, string):
    """
    >>> text = "Allowed Hello Hollow"
    >>> tuple(findall('ll', text))
    (1, 10, 16)
    """
    index = 0 - len(sub)
    o = []
    if not sub in string:
        return []
    while True:
        try:
            index = string.index(sub, index + len(sub))
        except:
            break
        o.append(index)
    return o
    #except ValueError:
    #        pass


def iterLines_skip_header(inFile, num_header_lines=1):
    for i, line in enumerate(open(inFile)):
        if i >= num_header_lines:
            yield line

def overlapping_occurrences(string, sub):
    count = start = 0
    while True:
        start = string.find(sub, start) + 1
        if start > 0:
            count+=1
        else:
            return count


def iterNlines(inFile, N):
    assert N >= 1
    with aopen.open(inFile) as f:
        lines = [f.readline() for i in range(N)]
        while 1:
            yield lines
            lines = [f.readline() for i in range(N)]
            if lines[0] == '':
                break


def GC(seq):
    t = seq.upper().count('C') + seq.upper().count('G')
    return float(t) / len(seq)


def elementSum(l1, l2):
    if len(l1) != len(l2):
        print len(l1), len(l2)
        raise ValueError
    return [l1[i] + l2[i] for i in range(len(l1))]


def iterTabbed(file):
    for line in open(file):
        yield line.strip().split('\t')


def getGeneTransFromCoords(chr, ex,
  genesFile='/net/crate-04/data/burge/alexrson/finalAnalyses'
  '/long_short_exons/genelines.sorted.gff.gz'):
    start_q, end_q = min(ex), max(ex)
    chr = chrFill(chr)
    if isinstance(genesFile, str):
        gf = tabix.Tabix(genesFile)
    else:
        gf = genesFile
    regs = gf.fetch('%s:%i-%i' % (chr, start_q, end_q))
    regs = [reg for reg in regs]
    genes = set()
    trans2gene = {}
    found_trans = False
    found_CDS = False
    for reg in regs:
        regl = reg.split('\t')
        attri = readAttributes.readAttributesIntoDict(regl[8])
        start_exon, end_exon = map(int, regl[3:5])
        if regl[2] == 'exon':
            if start_exon == start_q or start_exon == end_q or \
              end_exon == start_q or end_exon == end_q:
                found_trans = attri['Parent'][0]
        elif regl[2] == 'CDS':
            if start_exon == start_q or start_exon == end_q or \
              end_exon == start_q or end_exon == end_q:
                found_CDS = attri['Parent'][0]
        elif regl[2] == 'mRNA':
            gene = attri['Parent'][0]
            trans = attri['ID'][0]
            trans2gene[trans] = gene
            genes.add(tuple([regl[1], gene]))
        elif regl[2] == 'gene':
            gene = attri['ID'][0]
            genes.add(tuple([regl[1], gene]))
    if found_CDS:
        return trans2gene[found_CDS], found_CDS
    if found_trans:
        return trans2gene[found_trans], found_trans
    if len(genes) == 1:
        return genes.pop()[1]
    for esp, gene in genes:
        if esp in ['protein_coding', 'rkb', 'liana']:
            return gene, 'unknown'
    for esp, gene in genes:
        if esp in ['ucsc.knownGene-kgXref-ensGene']:
            return gene, 'unknown'
    print genes
    return genes.pop()[1]


def getGeneFromCoords(chr, ex,
  genesFile='/net/crate-04/data/burge/alexrson/finalAnalyses'
    '/long_short_exons/genelines.sorted.gff.gz'):
    start_q, end_q = min(ex), max(ex)
    chr = chrFill(chr)
    if isinstance(genesFile, str):
        gf = tabix.Tabix(genesFile)
    else:
        gf = genesFile
    regs = gf.fetch('%s:%i-%i' % (chr, start_q, end_q))
    regs = [reg for reg in regs]
    genes = set()
    trans2gene = {}
    found_trans = False
    found_CDS = False
    for reg in regs:
        regl = reg.split('\t')
        attri = readAttributes.readAttributesIntoDict(regl[8])
        start_exon, end_exon = map(int, regl[3:5])
        if regl[2] == 'exon':
            if start_exon == start_q or start_exon == end_q or \
              end_exon == start_q or end_exon == end_q:
                found_trans = attri['Parent'][0]
        elif regl[2] == 'CDS':
            if start_exon == start_q or start_exon == end_q or \
              end_exon == start_q or end_exon == end_q:
                found_CDS = attri['Parent'][0]
        elif regl[2] == 'mRNA':
            gene = attri['Parent'][0]
            trans = attri['ID'][0]
            trans2gene[trans] = gene
            genes.add(tuple([regl[1], gene]))
        elif regl[2] == 'gene':
            gene = attri['ID'][0]
            genes.add(tuple([regl[1], gene]))
    if found_CDS:
        return trans2gene[found_CDS]
    if found_trans:
        return trans2gene[found_trans]
    if len(genes) == 1:
        return genes.pop()[1]
    for esp, gene in genes:
        if esp in ['protein_coding', 'rkb', 'liana']:
            return gene
    for esp, gene in genes:
        if esp in ['ucsc.knownGene-kgXref-ensGene']:
            return gene
    if not genes:
        print regs, chr, ex
        return None
    return genes.pop()[1]


def getGeneNameFromCoords(chr, ex,
  genesFile='/net/crate-04/data/burge/alexrson/upf1/finalAnalyses'
    '/long_short_exons/genelines.sorted.gff.gz'):
    start_q, end_q = min(ex), max(ex)
    chr = chrFill(chr)
    if isinstance(genesFile, str):
        gf = tabix.Tabix(genesFile)
    else:
        gf = genesFile
    regs = gf.fetch('%s:%i-%i' % (chr, start_q, end_q))
    regs = [reg for reg in regs]
    genes = set()
    trans2gene = {}
    found_trans = False
    found_CDS = False
    for reg in regs:
        regl = reg.split('\t')
        attri = readAttributes.readAttributesIntoDict(regl[8])
        start_exon, end_exon = map(int, regl[3:5])
        gene = attri['Name'][0]
        genes.add(tuple([regl[1], gene]))
    if found_CDS:
        return trans2gene[found_CDS]
    if found_trans:
        return trans2gene[found_trans]
    if len(genes) == 1:
        return genes.pop()[1]
    for esp, gene in genes:
        if esp in ['protein_coding', 'rkb', 'liana']:
            return gene
    for esp, gene in genes:
        if esp in ['ucsc.knownGene-kgXref-ensGene']:
            return gene
    if not genes:
        return 'n/a'
    return genes.pop()[1]

def isCoding(chr, start, end, cds_gf):
    chr = chrFill(chr)
    if isinstance(cds_gf, str):
        gf = tabix.Tabix(cds_gf)
    else:
        gf = cds_gf
    regs = gf.fetch('%s:%i-%i' % (chr, start, end))
    regs = [reg for reg in regs]
    return len(regs) >= 1


def getFirstCol(f):
    o = []
    for line in open(f):
        o.append(line.strip().split('\t')[0])
    return o


def chrFill(chr):
    return ('chr' + chr).replace('chrchr', 'chr')


def open_skip_first(f):
    first = True
    for line in open(f):
        if first:
            first = False
        else:
            yield line


# works best if r << n
# samples from [0,n)
def sampleIntsWithoutReplacement(n, r):
    chosen = set()
    assert n >= r
    while len(chosen) < r:
        i = random.randint(0, n - 1)
        chosen.add(i)
    return chosen


def getBinIndex(v, bins):
    for i in range(len(bins) - 1):
        if v >= bins[i] and v < bins[i + 1]:
            return i
    if v == bins[-1]:
        return i
    return -1


def getBinIndex_soft_upper(v, bins):
    for i in range(len(bins) - 1):
        if v > bins[i] and v <= bins[i + 1]:
            return i
    return -1


def num_above(thresh, values):
    return len([1 for v in values if v > thresh])


def monitor(num, milestone=1000000, max=0):
    if max > 0:
        percentCount = int(max / 100.0)
        if num % percentCount == 0:
            print '%i%%' % (num / percentCount)
    if num % milestone == 0 and num > 0:
        print num


def readLinesDir(dir, include=None, exclude=None, also_include=None):
    dir = os.path.expanduser(dir)
    for f in os.listdir(dir):

        if include is not None:
            if not include in f:
                continue
        if also_include is not None:
            if not also_include in f:
                continue
        if exclude is not None:
            if exclude in f:
                continue
        print include, f
        for line in aopen.open(os.path.join(dir, f)):
            yield line


def listFs_Dir(dir, include=None, exclude=None, also_include=None):
    dir = os.path.expanduser(dir)
    fs = []
    for f in os.listdir(dir):
        if include is not None:
            if not include in f:
                continue
        if also_include is not None:
            if not also_include in f:
                continue
        if exclude is not None:
            if exclude in f:
                continue
        fs.append(os.path.join(dir, f))
    return fs


def listDir(dir, include=None, exclude=None, also_include=None):
    for f in os.listdir(dir):
        if include is not None:
            if not include in f:
                continue
        if also_include is not None:
            if not also_include in f:
                continue
        if exclude is not None:
            if exclude in f:
                continue
        yield f


def exLen(ex):
    assert ex[1] >= ex[0]
    return ex[1] - ex[0]


def getMedian(numericValues):
    theValues = sorted(numericValues)
    if len(theValues) % 2 == 1:
        return theValues[(len(theValues) + 1) / 2 - 1]
    else:
        lower = theValues[len(theValues) / 2 - 1]
        upper = theValues[len(theValues) / 2]
        return (float(lower + upper)) / 2


def log2(x):
    if x == 0:
        return None
    return math.log(x, 2)


def transpose(list_of_lists):
    return zip(*list_of_lists)


def pullFPKMSfromCuffgenes(genesFile):
    t2d = table2dict(genesFile)
    out = {}
    for gene in t2d:
        if t2d[gene]['FPKM_status'] == 'OK':
            out[gene] = float(t2d[gene]['FPKM'])
            out[t2d[gene]['gene_short_name']] = float(t2d[gene]['FPKM'])
    return out


def psisFromFile(summaryFile):
    t2d = table2dict(summaryFile)
    out = {}
    for event in t2d:
        if float(t2d[event]['ci_high']) - float(t2d[event]['ci_low']) <= 0.1:
            out[event] = float(t2d[event]['miso_posterior_mean'])
        else:
            pass
    return out


def arrayNormalization(arr):
    return 1.0 * (arr - arr[0]) / (arr[-1] - arr[0] + .00000000001)


# euclidean length of one
def normalize(vec):
    length = (sum([v ** 2 for v in vec])) ** 0.5
    return [v / length for v in vec]


# adds to one
def normalizeVector(arr):
    return 1.0 * arr / sum(arr)


def normalizeList(arr):
    arr = [1.0 * e / sum(arr) for e in arr]
    return arr
# combinations not permutations


def pairs(l):
    for i, e1 in enumerate(l):
        for j, e2 in enumerate(l):
            if j > i:
                yield e1, e2


#permutation pairs
# not both the same
def permuted_pairs(l):
    for i, e1 in enumerate(l):
        for j, e2 in enumerate(l):
            if j != i:
                yield e1, e2


def withreplacement3(l):
    for i in l:
        for j in l:
            for k in l:
                yield i, j, k


def pairs2(l1, l2):
    for i, e1 in enumerate(l1):
        for j, e2 in enumerate(l2):
            yield e1, e2


def evalBool(a):
    if a.lower() == 'true':
        return True
    elif a.lower() == 'false':
        return False
    raise ValueError


def minlen(listoflists):
    return min(map(len, listoflists))


def resample(population):
    n = len(population)
    _random, _int = random.random, int  # speed hack
    result = [None] * n
    for i in xrange(n):
        j = _int(_random() * n)
        result[i] = population[j]
    return result


def bootstrap_AtoB_vs_CtoD(a_s, bs, cs, ds, n=10000):
    # calculate p value for the null hypothesis that the difference
    # between C and D (D-C) >= that of A and B (B-A)
    AtoB_diff = []
    CtoD_diff = []
    for replicate_i in range(n):
        aps, bps, cps, dps = map(resample, (a_s, bs, cs, ds))  # a primes
        AtoB_diff.append(average(aps) - average(bps))
        CtoD_diff.append(average(cps) - average(dps))
    freq_CtoD_diff_greater = 0.0
    for ab, cd in zip(AtoB_diff, CtoD_diff):
        if cd >= ab:
            freq_CtoD_diff_greater += 1
    return float(freq_CtoD_diff_greater) / n


def mode(l):
    a = defaultdict(int)
    for e in l:
        a[e] += 1
    return sorted(a.iteritems(), key=operator.itemgetter(1))[-1][0]


def lrange(l):
    return range(len(l))


def top2(l):
    return top_n(l, 2)


def top_n(l, n):
    return sorted(l)[-n:]


def listsum(l):
    o = []
    for li in l:
        o += li
    return o


def subtractScalar(l, s):
    return [e - s for e in l]


def direction(f):
    if f > 0:
        return 'up'
    elif f < 0:
        return 'down'
    else:
        return '0'


def ttest_p(a, b):
    t_stat, p = stats.ttest_ind(a, b)
    df = len(a) + len(b) - 2
    one_tailed_p_value = 1.0 - stats.t.cdf(abs(t_stat), df)
    return one_tailed_p_value


# non list version
def dictFromTableCols(tableFile, keyCol, valueCol):
    o = dict()
    for line in open(tableFile):
        k, v = line.strip(
            ).split('\t')[keyCol], line.strip().split('\t')[valueCol]
        o[k] = v
    return o


# value list version
def dictFromTableColsValueList(tableFile, keyCol, valueCol):
    o = defaultdict(list)
    for line in open(tableFile):
        try:
            k, v = line.strip().split('\t')[keyCol],\
              line.strip().split('\t')[valueCol]
        except:
            print line.strip(), keyCol, valueCol
        if v == '.':
            continue
        if not k in o:
            o[k] = []
        o[k].append(v)
    return o


def zeros(n):
    return [0.0 for i in range(n)]


def scalarMultiply(l, scalar):
    return [li * scalar for li in l]


def orderInts(a, b):
    return min(a, b), max(a, b)


def divideList(l, d):
    return [1.0 * li / d for li in l]


def listIntersect(a, b):
    return list(set(a) & set(b))


def permuteSeq(seq):
    s = list(seq)
    random.shuffle(s)
    return ''.join(s)


def getNewPermutation(seq, disallowed_seqs):
    if len(seq) < 2:
        print 'fail'
        sys.exit()
    s = list(seq)
    random.shuffle(s)
    i = 0
    disallowed_seq = ','.join(disallowed_seqs)
    while ''.join(s) in disallowed_seq:
        random.shuffle(s)
        i += 1
        if i == 5000:
            print 'failed to get new sequence from %s' % seq
            sys.exit()
    return ''.join(s)


def subseqAllowMismatch(subseq, seq):
    #return subseqAllowNMismatches(subseq,seq,1)
    if len(subseq) > len(seq):
        return False
    for i in range(len(seq) - len(subseq)):
        s = seq[i:i + len(subseq)]
        mismatches = 0
        for i in range(len(s)):
            if s[i] != subseq[i]:
                if mismatches == 1:
                    break
                mismatches = 1
        else:
            return True
    return False


def subseqAllowNMismatches(subseq, seq, n):
    if len(subseq) > len(seq):
        return False
    for i in range(len(seq) - len(subseq)):
        s = seq[i:i + len(subseq)]
        if hamming_distance(s, subseq) <= n:
            return True
    return False


def hamming_N(str1, str2):
    if not len(str1) == len(str2):
        raise(ValueError, 'lengths don\'t match')
    str1 = str1.upper()
    str2 = str2.upper()
    str1 = str1.replace('N', '#')
    return sum(imap(operator.ne, str1, str2))


# from http://code.activestate.com/recipes/499304-hamming-distance/
def hamming_distance(str1, str2):
    assert len(str1) == len(str2)
    ne = operator.ne
    return sum(imap(ne, str1, str2))


def hamming_bc(str1, str2):
    assert len(str1) == len(str2)
    ne = operator.ne
    if 'N' in str1.upper() or 'N' in str2.upper():
        return 1
    return sum(imap(ne, str1, str2))


def interleave(l1, l2):
    o = []
    for i in range(min(len(l1), len(l2))):
        o.append(l1[i])
        o.append(l2[i])
    return o


def average(l):
    if not l:
        return -1
    return 1.0 * sum(l) / len(l)


def stdev(l):
    from numpy import array
    nums = array(l)
    return nums.std()


def sqrt(f):
    return float(f) ** 0.5


def SEM(l):
    return stdev(l) / sqrt(len(l))


def binDictValues(nbins, d):
    values = d.values()
    binsize = len(d) / nbins
    values.sort()
    kv_pairs = [(key, value) for key, value in d.items()]
    kv_pairs.sort(key=lambda x: x[1])
    ks = map(operator.itemgetter(0), kv_pairs)
    bins = list(chunks(ks, binsize))
    return bins


def drange(start, stop, step):
    o = []
    r = start
    while r < stop:
        o.append(r)
        r += step
    return o


def chunks(l, n):
    """ Yield successive n-sized chunks from l.
    """
    for i in xrange(0, len(l), n):
        yield l[i:i + n]


def reverseDict(d):
    o = dict()
    for k, v in d.iteritems():
        o[v] = k
    return o


def rank_dict_by_value(d):
 return sorted(d, key=d.get)


def swap_nested_dicts(d):
    new_d = defaultdict(dict)
    for k1, sub_d in d.items():
        for k2, v in sub_d.items():
            new_d[k2][k1] = v
    return new_d

def reverseDictValueList(d):
    o = dict()
    for k, vs in d.iteritems():
        for v in vs:
            o[v] = k
    return o


def reverseDictList(d):
    o = dict()
    for k, vs in d.iteritems():
        for v in vs:
            if not v in o:
                o[v] = []
            o[v].append(k)
    return o


def reverseDictManyKey(d):
    o = {}
    for k, v in d.items():
        if not v in o:
            o[v] = []
        o[v].append(k)
    return o


def reverseDictManyValue(d):
    o = defaultdict(dict)
    for k, v in d.iteritems():
        if not v in o:
            o[v] = {}
        o[v][k] = 0
    return o


# NOT in place
def reverseList(l):
    return [r for r in reversed(l)]


def iterReverse(l):
    for i in l[::-1]:
        yield i


def transposeDictNest(d):
    o = {}
    for k1, dinner in d.iteritems():
        for k2, v in dinner.iteritems():
            if not k2 in o:
                o[k2] = {}
            o[k2][k1] = d[k1][k2]
    return o


def highest_index(l):
    max_i = 0
    max_v = l[0]
    for i, v in enumerate(l):
        if v > max_v:
            max_v = v
            max_i = i
    return max_i


def iter_fasta_gene(in_file):
    current_gene = ''
    trans2seq = {}
    trans2exons = {}
    for gene, trans, chr, strand, exons, seq in iter_fasta_trans(in_file):
        if gene != current_gene and current_gene != '':
            yield current_gene, zchr, zstrand, trans2seq, trans2exons
            trans2seq = {}
            trans2exons = {}
        assert strand in '+-'
        zchr = chr
        zstrand = strand
        current_gene = gene
        trans2seq[trans] = seq.strip()
        trans2exons[trans] = exons
    yield gene, chr, strand, trans2seq, trans2exons


def iter_fasta_trans(in_file):
    #
    for idinfo, seq in iterLinePairs(in_file):
        if 'cox_hap' in idinfo or 'hap' in idinfo:
            continue
        if 'utr_splice' in idinfo:
            continue
        try:
            gene, trans, chr, strand, exons = idinfo[1:].strip().split('_')
        except:
            print idinfo
            gene1, gene2, trans, chr, strand, exons =\
              idinfo[1:].strip().split('_')
            gene = gene1 + '_' + gene2
        exons = simplejson.loads(exons)
        yield gene, trans, chr, strand, exons, seq


def pair_list2dict_float(in_file):
    o = {}
    for line in open(in_file):
        ll = line.strip().split('\t')
        o[ll[0]] = float(ll[1])
    return o

#subsequent
def pairwise(iterable):
    """s -> (s0,s1), (s1,s2), (s2, s3), ..."""
    a, b = tee(iterable)
    next(b, None)
    return izip(a, b)


def iter_subseqs(seq, k):
    for i in range(len(seq) - k + 1):
        subseq = seq[i:i + k].upper()
        if not 'N' in subseq:
            yield subseq


def iter_RNAfold_output(energy_file):
    """
    iterates through RNAfold input and returns an energy iterator
    """
    for l1, l2 in iterLinePairs(energy_file):
        yield float(l2.split(' (')[1].replace(')', ''))


def frac_above(l, v):
    total = float(len(l))
    num_above = float(len([e for e in l if e >v]))
    return num_above / total


def generate_decoys(motif, n):
    """
    generates a set of n permutations of motif
    """
    decoys = set()
    while len(decoys) < n:
        decoy = getNewPermutation(motif, decoys)
        decoys.add(decoy)
    return decoys


def make_dir(dirname):
    """
    makes the directory.
    Doesn't throw an error if it exists.
    """
    if not os.path.exists(dirname):
        try:
            os.makedirs(dirname)
        except:
            print 'it looks like the dir was made'\
              'by another thread extremely recently'


def save_fig(fig1, path, extentions=['png', 'pdf', 'svg']):
    for ext in extentions:
        fig1.savefig(path + '.' + ext)


def isFloat(e):
    return isinstance(e, float)

if __name__ == '__main__':
    pass
