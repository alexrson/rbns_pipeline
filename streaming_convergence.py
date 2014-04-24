import sys
import itertools
import numpy as np
import copy
import time
import os

def main():
    try:
        k, reads_file, num_passes, out_dir = sys.argv[1:]
        k = int(k)
        assert k in range(3,12)
        assert os.path.exists(reads_file)
        num_passes = int(num_passes)
        assert num_passes in range(2,20)
        assert os.direxists(out_dir)
    except:
        print 'Input argument error!'
        print ''
        print 'Usage:\npython streaming_convergence.py k',
        print 'reads_file num_passes output_directory'
        print 'k must be in [3,11]'
        print 'num_passes must be in [2, 19]'
        print 'readsfile must exist and only contain the read sequences (not fastq or fasta).'
        print 'out_dir must exist'
    SKA_weights = stream_counts(k, reads_file, num_passes)
    of = open(os.path.join(out_dir, 'SKA_weights.txt'), 'w')
    of.write('kmer\tweight\n')
    for kmer, weight in zip(yield_kmers(k), SKA_weights):
        of.write('%s\t%g\n' % (kmer, weight))
    of.close()


def test_simulated_data():
    num_iterations = 10
    k = 6
    motif_every = int(sys.argv[1:3])
    motif_every = map(int, sys.argv[1:3])
    inFile = 'simulated.%i.fa.gz' % motif_every
    print inFile
    current_weights = np.ones(4 ** k)
    weight_history = [np.ones(4 ** k)]
    iternal_out_table = 'kmer_count_internal_convergence.%i.xls' % (motif_every)
    for iteration_i in range(num_iterations):
        print 'round', iteration_i, ' of ', num_iterations
        if iteration_i == 0:
            current_weights = stream_continual_update_with_convergence_table(k,
              current_weights, inFile, iternal_out_table)
        else:
            current_weights = stream_without_continual_update(k,
              current_weights, inFile)
        weight_history.append(current_weights)
        current_weights = copy.copy(current_weights)
    weight_history = [normalize_mean_1(ws) for ws in weight_history]
    final_out_file = 'kmer_count_convergence.%i.xls' % (motif_every)
    make_table(k, weight_history, final_out_file)


def stream_counts(k, inFile, num_iterations):
    current_weights = np.ones(4 ** k)
    for iteration_i in range(num_iterations):
        if iteration_i == 0:
            current_weights = stream_continual_update(k,
              current_weights, inFile)
        else:
            current_weights = stream_without_continual_update(k,
              current_weights, inFile)
        current_weights = copy.copy(current_weights)
    assert len(current_weights) == 4 ** k
    return current_weights


def make_table(k, weight_history, out_file):
    of = open(out_file, 'w')
    of.write('kmer\t' + '\t'.join(
      ['round_%i' % i for i in range(len(weight_history))]) + '\n')
    for kmer_i, kmer in enumerate(yield_kmers(k)):
        of.write('%s\t' % kmer)

        for col_i in range(len(weight_history)):
            assert len(weight_history[col_i]) == 4 ** k
            of.write('%g\t' % weight_history[col_i][kmer_i])
        of.write('\n')
    of.close()


def normalize_mean_1(ws):
    total = float(sum(ws))
    l = float(len(ws))
    ans = [w * l / total for w in ws]
    return ans


def normalize_sum_1(ws):
    total = float(sum(ws))
    ans = [w / total for w in ws]
    return ans


def stream_continual_update_with_convergence_table(k, weights, inFile, out_file, how_often_to_write=10000):
    internal_history = []
    import aopen
    for linei,line in enumerate(aopen.open(inFile)):
        if linei % how_often_to_write == 0:
            norm_weights = copy.copy(weights)
            norm_weights = normalize_mean_1(norm_weights)
            internal_history.append(norm_weights)
            print linei
        read_seq = line.strip()
        pk = get_kmers(read_seq, k)
        assigned_weights = assign_kmer_weights(pk, weights)
        for kmer, weight in zip(pk, assigned_weights):
            kmeri = get_index_from_kmer(kmer)
            weights[kmeri] += weight
    of = open(out_file, 'w')
    of.write('kmer\t' + '\t'.join(
      ['reads_read_%i' % (i * how_often_to_write) for i in range(len(internal_history))]) + '\n')
    for kmer_i, kmer in enumerate(yield_kmers(k)):
        of.write('%s\t' % kmer)

        for col_i in range(len(internal_history)):
            assert len(internal_history[col_i]) == 4 ** k
            of.write('%g\t' % internal_history[col_i][kmer_i])
        of.write('\n')
    of.close()
    return weights


def stream_continual_update(k, weights, inFile):
    total_lines = count_lines(inFile) * 2
    start_time = time.time()
    import aopen
    for linei,line in enumerate(aopen.open(inFile)):
        if linei % 10000 == 0 and linei:
            elapsed_time = time.time() - start_time
            print 'predicted time remaining',\
              (total_lines - linei) / linei * elapsed_time / 3600,\
              'hours'
        read_seq = line.strip()
        pk = get_kmers(read_seq, k)
        assigned_weights = assign_kmer_weights(pk, weights)
        for kmer, weight in zip(pk, assigned_weights):
            kmeri = get_index_from_kmer(kmer)
            weights[kmeri] += weight
        if linei > 100:
            pass
    return weights


def count_lines(inFile):
    t = 0
    for l in open(inFile).xreadlines():
        t += 1
    return t


def assign_kmer_weights(pk, weights):
    kmer_weight = np.array(map(
      lambda s: float(weights[get_index_from_kmer(s)]
      ), pk))
    kmer_weight /= float(sum(kmer_weight))
    return kmer_weight


def stream_without_continual_update(k, in_weights, inFile):
    new_weights = np.ones(4 ** k)
    import aopen
    for linei, line in enumerate(aopen.open(inFile)):
        read_seq = line.strip()
        pk = get_kmers(read_seq, k)
        additional_weights = assign_kmer_weights(pk, in_weights)
        assert sum(additional_weights) - 1.0 < 0.001
        for kmer, weight in zip(pk, additional_weights):
            kmeri = get_index_from_kmer(kmer)
            new_weights[kmeri] += weight
    return new_weights


def get_kmers_no_homopolymers(seq, k):
    pk = []
    for i in range(0, len(seq) - k):
        kmer = seq[i:i + k]
        if len(set(kmer)) > 1: # no homo polymers
            pk.append(kmer)
    return pk


def get_kmers(seq, k):
    pk = set()
    for i in range(0, len(seq) - k):
        kmer = seq[i:i + k]
        pk.add(kmer)
    return pk


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


def yield_kmers(k):
    """
    An iterater to all kmers of length k in alphabetical order
    """
    bases = 'ACGT'
    for kmer in itertools.product(bases, repeat=k):
        yield ''.join(kmer)

if __name__ == '__main__':
    main()
