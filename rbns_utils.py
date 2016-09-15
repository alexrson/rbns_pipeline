import os
import operator
import itertools
import gzip
import numpy as np
from scipy import stats

def old_chris_formula(R, k, read_len):
    """
     Implements the formula Chris burge derived
    """
    return (4 ** k -1) * ( R * (read_len - k + 1) - (read_len - k) ) / (4 ** k + read_len - k - (R *(read_len - k +1)))


def chris_formula(R, k, read_len):
    """
     Implements the formula Chris burge derived
    """
    return (4. ** k - 1) * (read_len - k +1 ) * (R - 1) / (4. ** k -R * (read_len -k +1))


def get_adjacent_kmers(kmer):
    """
    returns all the k+1 mers that contain the kmer
    """
    return ['A' + kmer, 'C' + kmer, 'G' + kmer, 'T' + kmer,
      kmer + 'A', kmer + 'C', kmer + 'G', kmer + 'T']


def iter_RNAfold_output(energy_file):
    """
    iterates through RNAfold input and returns an energy iterator
    """
    for l1, l2 in iterLinePairs(energy_file):
        print energy_file, l1, l2
        yield float(l2.split(' (')[1].replace(')', ''))


def iterLinePairs(inFile):
    for l1, l2 in iterNlines(inFile, 2):
        yield l1, l2


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


def file_exists(fname):
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


def getBinIndex_soft_upper(v, bins):
    for i in range(len(bins) - 1):
        if v > bins[i] and v <= bins[i + 1]:
            return i
    return -1


def get_barcode(line):
    """
    extracts the barcode from the first line of a fastq quartet
    """
    return line.strip().split(':')[-1]


def get_index_from_kmer(kmer):
    """
    returns the base10 version of the base 4 DNA representation
    """
    index = 0
    base2face = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    for i, base in enumerate(kmer):
        if not base in 'ACGT':
            return -1
        power = len(kmer) - 1 - i
        index += base2face[base] * (4 ** power)
    return index


def get_kmer_from_index(kmax, index):
    """
    takes a number (essentially base 4)
    and returns the kmer it corresponds to in alphabetical order
    eg.
    AAAA = 0*1
    CA = 4*4 + 0*1
    GC = 3*4 + 1 * 1
    """
    bases = 'ACGT'
    out = ''
    for k in range(kmax - 1, -1, -1):
        face, index = divmod(index, 4 ** k)
        out += bases[face]
    return out


def yield_kmers(k):
    """
    An iterater to all kmers of length k in alphabetical order
    """
    bases = 'ACGT'
    for kmer in itertools.product(bases, repeat=k):
        yield ''.join(kmer)


def aopen(file, mode='r'):
    if file[-3:] == '.gz':
        return gzip.open(file, mode + 'b')
    else:
        return open(file, mode)


def hamming_N(str1, str2):
    if not len(str1) == len(str2):
        raise(ValueError, 'lengths don\'t match')
    str1 = str1.upper()
    str2 = str2.upper()
    str1 = str1.replace('N', '#')
    return sum(itertools.imap(operator.ne, str1, str2))


def pairwise(iterable):
    """s -> (s0,s1), (s1,s2), (s2, s3), ..."""
    a, b = itertools.tee(iterable)
    next(b, None)
    return itertools.izip(a, b)


# from http://code.activestate.com/recipes/499304-hamming-distance/
def hamming_distance(str1, str2):
    assert len(str1) == len(str2)
    ne = operator.ne
    return sum(itertools.imap(ne, str1, str2))


def iterNlines(inFile, N):
    assert N >= 1
    with aopen(inFile) as f:
        lines = [f.readline() for i in range(N)]
        while True:
            yield lines
            lines = [f.readline() for i in range(N)]
            if lines[0] == '':
                break


def save_fig(fig1, path, extentions=['png', 'pdf', 'svg']):
    for ext in extentions:
        fig1.savefig(path + '.' + ext)


def simpleaxis(sp):
    sp.spines['top'].set_visible(False)
    sp.spines['right'].set_visible(False)
    sp.get_xaxis().tick_bottom()
    sp.get_yaxis().tick_left()

def close_float_value(a, b, max_percent=1.0):
    if not (a > 0 and b > 0):
        return False
    ratio = float(max(a, b)) / float(min(a, b))
    percent_increase = (ratio - 1.0) * 100.0
    return percent_increase < max_percent


def significantly_enriched(xs, zthresh=2., scale='linear'):
    assert scale in ['linear', 'log']
    if scale =='log':
        xs = np.log2(xs)
    xs = stats.zscore(xs)
    return [x > zthresh for x in xs]


def iter4Lines(inFile):
    return iterNlines(inFile, 4)
