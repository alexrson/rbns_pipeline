#!/usr/bin/env python
# -*- utf-8 -*-

import sys
import operator
import subprocess
import shutil
import os
import cPickle
import itertools
import numpy
import simplejson
import collections
import ConfigParser
import time
import scipy.stats
import matplotlib.pyplot as plt
import unittest
import math

import aopen
import aUtils
import streaming_convergence
import aColors
import table2dict
import stacked_bar_kmers
from rna import rna

class BNS_Lib:
    def __init__(self, experiment, barcode, conc=0.0, polyIC=0.0, washes=1,
      temperature='RT', is_input=False):
        """
        Constructor for Library class
        """
        self.experiment = experiment
        self.barcode = barcode
        self.conc = conc
        self.polyIC = polyIC
        self.washes = washes
        self.temperature = temperature
        self.is_input = is_input
        self.wdir = experiment.wdir
        self.rdir = experiment.rdir
        self.expname = experiment.settings['experiment_name']
        self.motif_and_decoy = {}
        self.replaceative_enrich = {}
        self.k2presence_counts = {}
        self.k2libfrac = {}
        self.k2enrichments = {}

    def compare_top_kmers(self, k, most_enriched_lib, num_top_kmers_to_comp):
        """
        compares the enrichment of the top kmers in this library to another 
        library (usually the most enriched one).

        Returns pearsonr
        """
        top_kmers = most_enriched_lib.get_top_kmers(k, num_top_kmers_to_comp)
        most_enriched_lib_enrichments =\
          [most_enriched_lib.k2enrichments[k][kmer_i] for kmer_i in top_kmers]
        this_lib_enrichments =\
          [self.k2enrichments[k][kmer_i] for kmer_i in top_kmers]
        r, p = scipy.stats.pearsonr(
          most_enriched_lib_enrichments,
          this_lib_enrichments)
        return r

    def get_top_kmers(self, k, num_top_kmers_to_compare):
        """
        returns the top kmers by enrichment in this library
        """
        top_kmer_file = os.path.join(self.rdir, 'analyses',
          'top_kmers.%i.%s.pkl' % (k, self.barcode))
        if os.path.exists(top_kmer_file):
            sorted_kmers = cPickle.load(open(top_kmer_file))
            if len(top_kmer_file) == num_top_kmers_to_compare:
                return sorted_kmers
        enrich_kmers = zip(self.k2enrichments[k], range(4 ** k))
        enrich_kmers.sort(reverse=True)
        enrich_kmers = enrich_kmers[:num_top_kmers_to_compare]
        top_enrich, sorted_kmers = zip(*enrich_kmers)
        cPickle.dump(sorted_kmers, open(top_kmer_file, 'wb'))
        return sorted_kmers

    def get_top_seq_vs_next(self, k, num_next_kmers=7):
        """
        compares the top motif of this library to the next best motifs
        returns ratio of top enrich to next n motif average
        """
        top_enrichments = sorted(self.k2enrichments[k])[-num_next_kmers - 1:]
        return top_enrichments[-1] / aUtils.average(top_enrichments[:-1])

    def get_max_enrichment(self, k_for_max_enrichment):
        """
        returns the largest enrichment kmer of size k
        """
        return max(self.k2enrichments[k_for_max_enrichment])

    def get_barcode(self):
        """ returns the library's barcode """
        return self.barcode

    def get_full_label(self):
        """
        returns a good human readable label (for say a graph)
        """
        if not 'relevant_variables' in self.experiment.settings:
            return '%s conc: %0.1f, polyIC: %0.1f, washes: %i, temp: %s' %\
              (self.barcode, self.conc, self.polyIC,
              self.washes, self.temperature)
        else:
            rv2desc = {'input_rna': 'rna conc', 'poly_ic': '[poly(IC)]',
              'temp': 'T', 'barcode': 'barcode'}
            rv2val = {'input_rna': self.input_conc, 'poly_ic': self.polyIC,
              'temp': self.temperature, 'barcode': self.barcode}
            out = '%i, ' % int(self.conc)
            for rv in self.experiment.settings['relevant_variables']:
                out += '%s: %s, ' % (rv2desc[rv], str(rv2val[rv]))
            out = out.strip(', ').strip().strip(',')
            return out.strip(', ').strip().strip(',')

    def plot_next_kmer_frac(self):
        """
        makes a graph of the specificity as k increases.
        work in progress
        """
        fig1 = plt.figure()
        sp = fig1.add_subplot(211)
        next_frac = []
        for k in self.experiment.settings['ks_to_test_naive'][:-1]:
            next_frac.append(self.calculate_next_kmer_frac_counts(k))
        sp.plot(self.experiment.settings['ks_to_test_naive'][:-1],
          next_frac, 'b-')
        sp.set_ylim([0, 1.1 * max(next_frac)])
        sp.set_title('counts')
        sp = fig1.add_subplot(212)
        next_enrich = []
        for k in self.experiment.settings['ks_to_test_naive'][:-1]:
            next_enrich.append(self.calculate_next_kmer_frac_enrich(k))
        sp.plot(self.experiment.settings['ks_to_test_naive'][:-1],
          next_enrich, 'b-')
        sp.set_ylim([0, 1.1 * max(next_enrich)])
        sp.set_title('enrichment')
        fig1.savefig(
          os.path.join(self.rdir, 'plots',
          '%s.%0.1f.kpick.pdf' % (self.barcode, self.conc)))

    def calculate_next_kmer_frac_counts(self, k):
        """
        for a given k:
        calculates the number of counts for addding a base to either side
        of the kmer and checks the counts for each of those.
        then calculates the proportion of the counts that go to the highest
        one of those and returns that fraction
        """
        assert k in self.experiment.settings['ks_to_test_naive'] and\
          k + 1 in self.experiment.settings['ks_to_test_naive']
        kmer_counts = self.k2counts[k]
        max_counts = max(kmer_counts)
        # assume only one max
        top_kmer_i = list(kmer_counts).index(max_counts)
        kmer = get_kmer_from_index(k, top_kmer_i)
        adjacent_kmers = aUtils.get_adjacent_kmers(kmer)
        ajacent_kmer_i_s = map(get_index_from_kmer, adjacent_kmers)
        adjacent_kmer_counts = [self.k2counts[k + 1][kmer_i]
          for kmer_i in ajacent_kmer_i_s]
        best_adjeacent_kmer_count = max(adjacent_kmer_counts)
        adjacent_kmer_frac = float(best_adjeacent_kmer_count)\
          / sum(adjacent_kmer_counts)
        return adjacent_kmer_frac

    def calculate_next_kmer_frac_enrich(self, k):
        """
        for a given k:
        calculates the number of counts for addding a base to either side
        of the kmer and checks the enrichments for each of those.
        then calculates the proportion of the (summed
        enrichments that go to the highest
        one of those and returns that fraction
        """
        assert k in self.experiment.settings['ks_to_test_naive'] and\
          k + 1 in self.experiment.settings['ks_to_test_naive']
        kmer_enrichments = self.k2enrichments[k]
        max_enrichment = max(kmer_enrichments)
        # assume only one max
        top_kmer_i = list(kmer_enrichments).index(max_enrichment)
        kmer = get_kmer_from_index(k, top_kmer_i)
        adjacent_kmers = aUtils.get_adjacent_kmers(kmer)
        ajacent_kmer_i_s = map(get_index_from_kmer, adjacent_kmers)
        adjacent_kmer_enrichments = [self.k2enrichments[k + 1][kmer_i]
          for kmer_i in ajacent_kmer_i_s]
        best_adjeacent_kmer_enrichment = max(adjacent_kmer_enrichments)
        adjacent_kmer_frac = float(best_adjeacent_kmer_enrichment)\
          / max_enrichment
        return adjacent_kmer_frac

    def split_reads_exist(self):
        """
        returns true if the split reads file for this library exists
        and is non empty
        does not check if it is complete
        """
        return aUtils.check_file_exists(self.get_split_reads_name())

    def slink_split_reads(self):
        """
        makes symbolic links for each of the 
        """
        if self.is_input:
            slink = os.path.join(self.rdir, 'split_reads', 
              '%s_input.reads' % (self.expname))
        else:
            slink = os.path.join(self.rdir, 'split_reads',
              '%s_%g.reads' % (self.expname, self.conc))
        try:
            os.symlink(self.get_split_reads_name(), slink)
        except:
            pass

    def get_split_reads_name(self):
        """
        returns the full path of this library's split reads file
        """
        split_reads_file = os.path.join(self.rdir, 'split_reads',
          '%s_%s.reads' % (self.expname, self.barcode))
        return split_reads_file

    def get_split_rhandle(self):
        """
        returns a read file handle to the split reads
        """
        return aopen.open(self.get_split_reads_name(), 'r')

    def get_split_whandle(self):
        """
        returns a write file handle to the split reads
        """
        return aopen.open(self.get_split_reads_name(), 'w')

    def get_naive_counts_file(self):
        """
        returns the full path of the naive counts pkl
        """
        counts_file = os.path.join(self.rdir, 'counts', 'naive',
          '%s_%s.pkl' % (self.expname, self.barcode))
        return counts_file

    def get_stream_counts_file(self, k):
        """
        returns the full path of the stream counts pkl
        """
        counts_file = os.path.join(self.rdir, 'counts', 'stream',
          '%s_%s.%i.pkl' % (self.expname, self.barcode, k))
        return counts_file

    def naive_counts_exist(self):
        """
        true if the naive counts pkl and num reads pkl exist
        """
        return aUtils.check_file_exists(self.get_naive_counts_file()) and\
          aUtils.check_file_exists(self.get_num_reads_file())

    def presence_counts_exist(self, k):
        """
        true if the stream counts pkl and num reads pkl exist
        """
        return aUtils.check_file_exists(self.get_presence_counts_file(k))

    def stream_counts_exist(self, k):
        """
        true if the stream counts pkl and num reads pkl exist
        """
        return aUtils.check_file_exists(self.get_stream_counts_file(k))

    def get_num_reads_file(self):
        """
        true if the num reads pkl exists
        """
        num_reads_file = os.path.join(self.rdir, 'counts', 'naive',
          '%s_%s.num_reads.pkl' % (self.expname, self.barcode))
        return num_reads_file

    def get_enrich_pkl(self):
        """
        returns the pkl file wich stores enrichemnts
        """
        enrich_pkl = os.path.join(self.rdir, 'tables', 
          '%s_%s.both_enrichment.pkl' % (self.expname, self.barcode))
        return enrich_pkl

    def do_naive_count(self):
        """
        NAIVE COUNTS
        if the pkls for this process exist it loads them
        other wise counts naively and saves results in pkl

        returns a dict of k -> array of counts (length is 4^k)
        """
        counts_file = self.get_naive_counts_file()
        num_reads_file = self.get_num_reads_file()
        if aUtils.check_file_exists(counts_file) and\
           aUtils.check_file_exists(num_reads_file): 
            print 'loading naive counts from %s' % counts_file
            self.k2counts = cPickle.load(open(counts_file, 'rb'))
            if not set(self.experiment.settings['ks_to_test_naive']).issubset(set(self.k2counts.keys())):
                print 'missing some k values: in naive count'
                os.remove(counts_file)
            else:
                self.num_reads = cPickle.load(open(num_reads_file, 'rb'))
                return self.k2counts
        if not (aUtils.check_file_exists(counts_file) and aUtils.check_file_exists(num_reads_file)):
            reads_file = self.get_split_reads_name()
            self.k2counts = dict()
            read_len = self.experiment.settings['read_len']
            kmer2index = {}
            for k in self.experiment.settings['ks_to_test_naive']:
                self.k2counts[k] = numpy.zeros(4 ** k, dtype=int)
                for i, kmer in enumerate(yield_kmers(k)):
                    kmer2index[kmer] = i
            print 'counting kmers in %s' % reads_file
            for i, line in enumerate(aopen.open(reads_file)):
                aUtils.monitor(i, 10000)
                if 'N' in line:
                    continue
                try:
                    for k in self.experiment.settings['ks_to_test_naive']:
                        for ki in range(0, read_len - k + 1):
                            kmer_index = kmer2index[line[ki:ki + k]]
                            self.k2counts[k][kmer_index] += 1
                except:
                    continue
            self.num_reads = i
            aUtils.make_dir(os.path.dirname(num_reads_file))
            cPickle.dump(self.num_reads, open(num_reads_file, 'wb'))
            cPickle.dump(self.k2counts, open(counts_file, 'wb'))
            print 'saving naive counts %s' % counts_file
        self.write_counts()
        # slink counts pkls
        if self.is_input:
            slink = os.path.join(self.rdir, 'counts','naive',
              '%s_naive_conc.input.pkl' % (self.expname))
        else:
            slink = os.path.join(self.rdir, 'counts','naive',
              '%s_naive_conc.%g.pkl' % (self.expname, self.conc))
        try:
            os.symlink(counts_file, slink)
            print 'symlink made', slink
        except:
            print 'making symlink failed', slink
            pass
        return self.k2counts

    def do_stream_count(self, ks='all', force_recount=False):
        """
        STREAM COUNTS
        if the pkls for this process exist it loads them
        other wise counts stream and saves results in pkl

        returns a dict of k -> array of counts (length is 4^k)
        """
        self.k2stream_weights = {}
        if ks == 'all':
            ks = self.experiment.settings['ks_to_test_streaming']
        for k in ks:
            counts_file = self.get_stream_counts_file(k)
            if aUtils.check_file_exists(counts_file) and not force_recount:
                print 'loading streaming', counts_file
                self.k2stream_weights[k] = cPickle.load(open(counts_file, 'rb'))
            else:
                reads_file = self.get_split_reads_name()
                reads_file = self.copy_to_scratch(reads_file, k)
                print 'stream counting', k, self.conc
                self.k2stream_weights[k] = streaming_convergence.stream_counts(
                  k, reads_file, 2)
                aUtils.make_dir(os.path.dirname(counts_file))
                print 'saving', counts_file
                cPickle.dump(self.k2stream_weights[k], open(counts_file, 'wb'))
            self.write_stream_counts(k)
        self.k2stream_libfracs = {}
        for k in self.k2stream_weights.keys():
            total_counts = numpy.sum(self.k2stream_weights[k])
            self.k2stream_libfracs[k] = numpy.array(
              [self.k2stream_weights[k][kmeri] / total_counts 
              for kmeri in range(4 ** k)])
            try:
                assert numpy.sum(self.k2stream_libfracs[k]) == 1.0
            except:
                pass
        return self.k2stream_weights

    def do_presence_count(self, k, force_recount=False):
        """
        Presence count
        """
        counts_file = self.get_presence_counts_file(k)
        if aUtils.check_file_exists(counts_file) and not force_recount:
            self.k2presence_counts[k] = cPickle.load(open(counts_file))
            if not set(self.experiment.settings['ks_to_test_naive']
              ).issubset(set(self.k2presence_counts.keys())):
                os.remove(counts_file)
        if not aUtils.check_file_exists(counts_file) or force_recount:
            print 'doing the presence count'
            self.k2presence_counts[k] = numpy.zeros(4 ** k)
            reads_file = self.get_split_reads_name()
            reads_file = self.copy_to_scratch(reads_file, k)
            read_len = self.experiment.settings['read_len']
            kmer2index = {}
            for i, kmer in enumerate(yield_kmers(k)):
                kmer2index[kmer] = i
            line_count = 0
            for line in aopen.open(reads_file):
                kmers_present = set()
                for ki in range(0, read_len - k + 1):
                    kmer = line[ki:ki + k]
                    if kmer in kmer2index:
                        kmers_present.add(kmer)
                for kmer in kmers_present:
                    self.k2presence_counts[k][kmer2index[kmer]] += 1.0
                line_count += 1
                aUtils.monitor(line_count)
            print 'presence counts', self.k2presence_counts[k]
            print 'line count:', line_count
            assert line_count > 0
            self.k2presence_counts[k] = self.k2presence_counts[k] / line_count
            print self.k2presence_counts[k]
            cPickle.dump(self.k2presence_counts[k], open(counts_file, 'w'))

    def get_presence_counts_file(self, k):
        """
        returns the full path of the presence counts pkl
        """
        counts_file = os.path.join(self.rdir, 'counts', 'presence',
          '%s_%s.%i.pkl' % (self.expname, self.barcode, k))
        return counts_file
    
    def copy_to_scratch(self, in_file, ind):
        out_location = os.path.join(self.wdir,
                                    'scratch_reads',
                                    os.path.basename(in_file))
        assert aUtils.check_file_exists(in_file)
        aUtils.make_dir(os.path.dirname(out_location))
        if not os.path.exists(out_location):
            print 'copying to scratch', out_location
            shutil.copy(in_file, out_location)
        return out_location

    def write_counts(self):
        """
        writes the counts in a ASCII table
        seperate file for each k
        """
        for k in self.experiment.settings['ks_to_test_naive']:
            out_file = os.path.join(self.rdir, 'counts', 'naive',
              '%s.%i.txt' % (self.barcode, k))
            of = open(out_file, 'w')
            of.write('kmer\tcounts\n')
            for kmer, count in zip(yield_kmers(k), self.k2counts[k]):
                of.write('%s\t%i\n' % (kmer, count))
            of.close()
            #make slinks
            if self.is_input:
                slink = os.path.join(self.rdir, 'counts','naive',
                  '%s.naive.conc.input.%i.txt' % (self.expname, k))
            else:
                slink = os.path.join(self.rdir, 'counts','naive',
                  '%s.naive.conc.%g.%i.txt' % (self.expname, self.conc, k))
            try:
                os.symlink(out_file, slink)
            except:
                pass


    def write_stream_counts(self, k):
        """
        writes the counts in a ASCII table
        seperate file for each k
        """
        out_file = os.path.join(self.rdir, 'counts', 'stream',
          '%s.%i.txt' % (self.barcode, k))
        of = open(out_file, 'w')
        of.write('kmer\tweight\n')
        for kmer, weight in zip(yield_kmers(k), self.k2stream_weights[k]):
            of.write('%s\t%g\n' % (kmer, weight))
        of.close()

    def count_replaceative(self, k):
        """
        Runs Albert's code and waits untils it is done if necessary
        otherwise loads from kmers.txt
        If kmers.txt is incomplete this could be a problem.
        """
        out_dir = os.path.join(self.rdir, 'counts', 'replaceative',
          str(k), self.barcode)

        aUtils.make_dir(out_dir)
        results_file = os.path.join(out_dir, 'kmers.txt')
        if aUtils.check_file_exists(results_file):
            self.replaceative_enrich[k] =\
              read_replaceative_results(results_file)
            num_kmers_count = len(self.replaceative_enrich[k])
            if num_kmers_count >= self.experiment.settings[
              'replaceative_num_top_kmers']:
                return
            print 'There aren\'t enough kmers counted for this lib. %i' %\
              num_kmers_count
        lib_seq_file_name = self.experiment.get_input_split_name()
        split_file_name = self.get_split_reads_name()
        command = 'bash DKmerFinder.sh %s %s %i %i %s 1000000 '\
          'qsub qstat %s' % \
          (split_file_name,
          lib_seq_file_name,
          k, self.experiment.settings['replaceative_num_top_kmers'],
          out_dir, self.rdir)
        print 'launching replaceative for %0.1f %i' % (self.conc, k)
        print '\t%s' % out_dir
        subprocess.call(command.split())
        self.replaceative_enrich[k] =\
          read_replaceative_results(results_file)

    def run_rna_fold(self):
        """
        Folds all the reads using RNAfold and stores the result
        will launch on cluster if requested
        """
        out_file = os.path.join(self.rdir, 'structure',
          '%s.fe' % self.barcode)
        if aUtils.check_file_exists(out_file + '.gz'):
            return 0
        print 'calculating structure for'\
          ' %s with concentration %0.1f' % (self.barcode, self.conc)
        tmp_file = os.path.join(self.wdir, '%s.reads' % self.barcode)
        tmp_file_out = os.path.join(self.wdir, '%s.fe' % self.barcode)
        tmp_file_err = os.path.join(self.wdir, '%s.err' % self.barcode)
        command = 'mkdir %s; cd %s; cp %s %s ; cat %s | RNAfold '\
          '1> %s 2> %s; cp %s %s %s ; gzip %s' %\
          (self.wdir, self.wdir, self.get_split_reads_name(),
          tmp_file, tmp_file,
          tmp_file_out, tmp_file_err,
          tmp_file_out, tmp_file_err,
          os.path.join(self.rdir, 'structure'),
          out_file)
        if not self.experiment.use_cluster:
            subprocess.Popen(command, shell=True).wait()
        else:
            script_options = {'nodes': '1', 'ppn': '1',
              'outf': self.experiment.settings['experiment_name']
              + '.submission_log',
              'jobname': self.experiment.settings['experiment_name']
              + '_' + self.barcode + '.fold',
              'queue': 'long', 'workingdir': self.rdir,
              'command': command}
            return aUtils.launch(command, script_options)

    def split_by_structure(self, run_now=False):
        """
        splits the reads for this library into seperate files
        binned by the folding free energy if they don't exist and aren't empty

        bins are defined in the settings file

        This can throw errors if there are no reads in one of the defined
        energy bins.

        will launch on cluster if requested
        """
        out_files = ['%s/structure/%s.%i.%0.1f_to_%0.1f.reads' %
          (self.rdir, self.barcode, i, conc1, conc2)
          for i, (conc1, conc2) in enumerate(self.experiment.energy_bins)]
        if all([aUtils.check_file_exists(outfile) for outfile in out_files]):
            return 0
        if not self.experiment.use_cluster or run_now:
            ofs = map(lambda x: aopen.open(x, 'w'), out_files)
            reads_file = self.get_split_reads_name()
            energy_file = os.path.join(self.rdir,
              'structure', '%s.fe.gz' % self.barcode)
            for line, energy in itertools.izip(aopen.open(reads_file),
              aUtils.iter_RNAfold_output(energy_file)):
                bin_i = aUtils.getBinIndex_soft_upper(energy,
                  self.experiment.settings['free_energy_limits'])
                ofs[bin_i].write(line)
            map(lambda of: of.close(), ofs)
        else:
            command = 'python ~alexrson/snorelax/script.py '\
              'bind_n_seq_pipeline.fold_split %s %s 1> %s 2> %s' % \
              (self.experiment.settings_file, self.barcode,
              self.rdir + '/structure/' + self.barcode + '.spl.out',
              self.rdir + '/structure/' + self.barcode + '.spl.err')
            script_options = {'nodes': '1', 'ppn': '1',
              'outf': self.experiment.settings['experiment_name']
              + '.submission_log',
              'jobname': self.experiment.settings['experiment_name']
              + '_' + self.barcode + '.splitE',
              'queue': 'long', 'workingdir': self.rdir,
              'command': command}
            return aUtils.launch(command, script_options)

    def count_motif_individual(self):
        """
        returns a dict of
        number of occurences of the motif ->
            number of reads than have that number of motifs
        """
        results_file = os.path.join(self.rdir,
          'analyses', '%s.motif_counts.pkl' % self.barcode)
        if aUtils.check_file_exists(results_file):
            self.motif_count = cPickle.load(open(results_file, 'rb'))
        else:
            self.motif_count = collections.Counter()
            for line in self.get_split_rhandle().readlines():
                count = line.count(self.experiment.known_motif)
                self.motif_count[count] += 1
            cPickle.dump(self.motif_count, self.experiment.get_rdir_fhandle(
              'analyses', '%s.motif_counts.pkl' % self.barcode))
        self.num_reads = sum(self.motif_count.values())
        return self.motif_count

    def calculate_naive_libfracs(self):
        """
        calculates the libfrac for each k and gives a dict
        for this library:

        libfrac = [kmer counts] / total_counts

        No normalization for the input library
        """
        libfrac_f = os.path.join(self.rdir, 'analyses',
          '%s.libfrac.pkl' % self.barcode)
        if aUtils.check_file_exists(libfrac_f):
            self.k2libfrac = cPickle.load(open(libfrac_f, 'rb'))
            if all([k in self.k2libfrac
              for k in self.experiment.settings['ks_to_test_naive']]):
                return
            else:
                print 'Not all naive ks have been precalculated'
        self.k2libfrac = {}
        for k in self.experiment.settings['ks_to_test_naive']:
            counts = numpy.array(self.k2counts[k], dtype=float)
            total_counts = sum(counts)
            if not total_counts:
                raise ValueError('no counts: %s' % self.barcode)
            self.k2libfrac[k] = counts / total_counts
            cPickle.dump(self.k2libfrac, open(libfrac_f, 'wb'))
        return self.k2libfrac
    

    def calculate_enrichment(self, k='all'):
        """
        calculates the dict of k-> enrichments

        enrichment = libfrac_kmer_this_library / libfrac_input_lib
        """
        enrich_file = self.get_enrich_pkl()
        if aUtils.check_file_exists(enrich_file):
            self.k2enrichments, self.k2rep_enrichments_arr = cPickle.load(open(enrich_file))
            if not set(self.experiment.settings['ks_to_test_naive']).issubset(
              set(self.k2enrichments.keys())):
                os.remove(enrich_file)
                print 'ks in enrichemnt insufficient', self.k2enrichments.keys()
        if not aUtils.check_file_exists(enrich_file):
            print 'calculating enrichemnts'
            if not self.k2libfrac:
                self.calculate_naive_libfracs()
            if not self.experiment.input_lib.k2libfrac:
                self.experiment.input_lib.calculate_naive_libfracs()
            self.k2enrichments = {}
            for k in self.experiment.settings['ks_to_test_naive']:
                self.experiment.input_lib.calculate_naive_libfracs()
                self.k2enrichments[k] =\
                  self.k2libfrac[k] / self.experiment.input_lib.k2libfrac[k]
            self.k2rep_enrichments_arr = {}
            if  self.experiment.settings['replaceative']:
                for k in self.experiment.settings['ks_to_test_replaceative']:
                    if not k in self.replaceative_enrich:
                        self.count_replaceative(k)
                    self.k2rep_enrichments_arr[k] = numpy.array(
                      [self.replaceative_enrich[k].get(kmer, 1.0)
                      for kmer in yield_kmers(k)])
            cPickle.dump((self.k2enrichments, self.k2rep_enrichments_arr),open(enrich_file, 'w'))
        # make symlinks
        if self.is_input:
            slink = os.path.join(self.rdir, 'tables', 
              '%s_input.both_enrichment.pkl' % (self.expname))
        else:
            slink = os.path.join(self.rdir, 'tables', 
              '%s_%g.both_enrichment.pkl' % (self.expname, self.conc))
        try:
            os.symlink(enrich_file, slink)
            print 'made symlink', slink
        except:
            print 'making symlink', slink, 'failed'
            pass
        if k == 'all':
            return self.k2enrichments
        else:
            return self.k2enrichments[k]

    def get_naive_enrichment_dict(self, k):
        """
        returns a dictionary of kmer -> naive enrichment
        """
        kmer2enrichment = {}
        if not self.k2enrichments:
            self.calculate_enrichment()
        for kmer_i, enrichment in enumerate(self.k2enrichments[k]):
            kmer = get_kmer_from_index(k, kmer_i)
            kmer2enrichment[kmer] = enrichment
        return kmer2enrichment

    def calcB(self, kmer):
        """
         calculates the B value 
        """
        k = len(kmer)
        kmeri = get_index_from_kmer(kmer)
        enrichment = self.k2enrichments[k][kmeri]
        return chris_formula(enrichment, 
                             k,
                             self.experiment.settings['read_len'])

    def check_libfrac_sum(self):
        """
        Verifies that the libfracs add up to one
        """

        for libfrac in self.k2libfrac.values():
            if abs(1.0 - sum(libfrac)) > 1e-7:
                print 'Lib frac warning:\nThe library fractions dont quite'\
                  'sum to 1.0. They are off by:'\
                  '%f' % abs(1.0 - sum(libfrac))

    def motif_count_frac(self, count):
        """
        returns the fraction of reads that have $count motif occurences
        """
        return 1.0 * self.motif_count[count] / self.num_reads

    def count_double_occurrence(self, decoy):
        """
        counts the number of reads with a decoy and a true motif
        """
        motif = self.experiment.known_motif
        count = 1
        for line in self.get_split_rhandle().readlines():
            if motif in line and decoy in line:
                if line.index(motif) < line.index(decoy):
                    count += 1
        self.motif_and_decoy[decoy] = count
        return count

    def get_double_motif_fname(self):
        """
        returns the filename for the double motif table
        """
        return os.path.join(self.rdir, 'analyses',
          '%s.double.motif.pkl' % self.barcode)

    def has_double_motif(self):
        """
        True iff double motif table exists
        """
        results_file = self.get_double_motif_fname()
        if aUtils.check_file_exists(results_file):
            self.motif_and_decoy = cPickle.load(open(results_file))
            return True
        return False

    def save_double_motif(self):
        """
        saves the double motif file
        """
        cPickle.dump(self.motif_and_decoy,
          open(self.get_double_motif_fname(), 'w'))


# This class contains all the data for a given Bind-N-Seq Experiment
#
class Bnse:
    def __init__(self, settings_file, use_cluster=False, run_now=True):
        if isinstance(use_cluster, str):
            use_cluster = use_cluster.lower() == 'true'
        if isinstance(run_now, str):
            run_now = run_now.lower() == 'true'
        self.settings_file = settings_file
        self.use_cluster = use_cluster
        self.run_now = run_now
        if not aUtils.check_file_exists(self.settings_file):
            raise RuntimeError(
              'settings file %s does not exist' % self.settings_file)
        self.process_settings()
        if run_now:
            self.create_rdir()
            self.run()
        elif not self.run_now:
            if use_cluster:
                self.run_on_cluster()
            elif not use_cluster:
                pass

    def run_on_cluster(self):
        """
        This will launch the Bnse class on a cluster node by piping
        a python command to call the free method, launcher into qsub.
        """
        command = 'python ~alexrson/snorelax/script.py '\
          'bind_n_seq_pipeline.Bnse %s true true 1> /net/utr/data/atf/alexr'\
          'son/bind_n_seq_pipeline/%s.out 2> /net/utr/data/atf/alexrson/b'\
          'ind_n_seq_pipeline/%s.err' % \
          (self.settings_file,
          self.settings['experiment_name'],
          self.settings['experiment_name'])
        script_options = {'nodes': '1', 'ppn': '8',
          'outf': self.settings['experiment_name'] + '.submission_log',
          'jobname': self.settings['experiment_name'],
          'queue': 'long', 'workingdir': os.getcwd(), 'command': command}
        aUtils.launch(command, script_options)

    def run(self):
        """
        Runs the requested analyses right here right now. Assumes that
        the settings file has been processed.
        """
        # copy unsplit to working dir
        print 'copy unsplit'
        self.copy_unsplit_to_wd()
        #split unsplit into the results dir
        print 'split'
        self.split_reads()
        # counts
        self.perform_kmer_counts()
        print 'make matlab inputs'
        self.create_matlab_inputs()
        # Normalize
        print 'calculate libfracs'
        self.calculate_all_libfracs()
        #
        print 'calculate enrichments'
        self.calculate_all_enrichments()
        print 'determine k'
        self.plot_k_picker()
        print 'determine concordance with most enriched lib'
        self.plot_concordance_to_most_enriched()

        print 'sort by enrichment'
        self.sort_all_kmers_by_enrichment()
        self.make_enrichment_table()
        print 'plot enrichments'
        if True:
            self.plot_enrichment_humps()
            self.plot_enrichment_humps_with_secondary()
            self.plot_relativeAffinity_humps_with_secondary()
            self.plot_enrichment_humps_with_secondary_stream()
            self.plot_enrichment_humps_with_secondary_presence()
            self.plot_enrichment_v_stream()
            self.plot_enrichment_v_presence()
            self.plot_presence_v_stream()
            print 'make enrichment table'
            self.make_enrichment_table()
            self.make_SKA_table()
            self.make_enrichment_hist()
            self.make_enrichment_hist_stacked()
            print 'verify concordances'
            self.verify_concordances()
            print 'make motif libfrac table'
            self.make_motif_libfrac_table()
            self.make_counts_tables()
            print 'plot num counts'
            self.plot_num_counts()
            self.count_motif_occurences()
            self.count_controlled_second_motif()
        print 'make table of B'
        self.make_Btable()
        # FE stuff
        print 'run rna fold'
        self.run_rna_fold()
        print 'count energy split'
        self.count_fold_split()
        print 'run plot_fold_split_hump'
        self.plot_fold_split_hump()
        self.make_fold_split_table()
        #comparisons
        print 'compare to other experiments'
        self.compare_all_other_experiments()

    def make_Btable(self):
        """
        Makes a table of the B values and R values
        """
        for k in self.settings['ks_to_test_naive']:
            t_of = self.get_rdir_fhandle('tables',
              'Btable.%i.xls' % k)
            self.make_table_header(t_of)
            for kmeri, kmer in enumerate(yield_kmers(k)):
                t_of.write(kmer + '\t')
                t_of.write('\t'.join(
                  [str(lib.calcB(kmer)) for lib in self.plibs]) + '\n')
            t_of.close() 

    def compare_all_other_experiments(self):
        """
        does all the comparisons of the experiment to the the other
        ones specified in the settings
        """
        for alt_settings_file in self.settings['experiments_to_compare']:
            if not os.path.isabs(alt_settings_file):
                alt_settings_file =\
                  os.path.join(os.getcwd(), alt_settings_file)
            assert aUtils.check_file_exists(alt_settings_file)
            self.compare_to_other_BNS(alt_settings_file)

    def compare_to_other_BNS(self, other_experiment_settings_file):
        """
        does the comparison to a single other experiment
        """
        print 'making a new BNS results instance %s' %\
          other_experiment_settings_file
        other_exp = Bnse_results(other_experiment_settings_file)
        print 'Comparing to %s' % other_exp.settings['experiment_name']
        other_exp.input_lib.calculate_enrichment()
        if True:
            # scatter of protein concentrations
            fig1 = plt.figure()
            sp = fig1.add_subplot(1, 1, 1)
            for lib in self.plibs:
                sp.plot(1.0, lib.conc, 'o', color=aColors.steel_blue)
            for lib in other_exp.plibs:
                sp.plot(2.0, lib.conc, 'o', color=aColors.bacon)
            sp.set_yscale('log')
            sp.set_xlim([0.5, 2.5])
            sp.set_xticks([1, 2])
            sp.set_xticklabels([self.settings['experiment_name'],
              other_exp.settings['experiment_name']])
            fig1.savefig(os.path.join(self.rdir,
              'comparisons',
              'protein_conc.%s.pdf'
              % other_exp.settings['experiment_name']))
            # Scatter: enrich
            k = 6
            kmers = [kmer for kmer in yield_kmers(k)]
            for lib_main, lib_other in itertools.product(self.plibs,
              other_exp.plibs):
                fig2 = plt.figure()
                sp = fig2.add_subplot(1, 1, 1)
                if not aUtils.within_xpercent(lib_main.conc, lib_other.conc, 40):
                    continue
                main_kmer2enrich = lib_main.get_naive_enrichment_dict(k)
                other_kmer2enrich = lib_other.get_naive_enrichment_dict(k)
                main_enriches = [main_kmer2enrich[kmer] for kmer in kmers]
                other_enriches = [other_kmer2enrich[kmer] for kmer in kmers]
                thresh = numpy.mean(main_enriches) + 2 * numpy.std(main_enriches)
                for x, y, kmer in zip(main_enriches, other_enriches, yield_kmers(k)):
                    sp.loglog(x,y,'.', 
                      color=aColors.protein_colors(kmer, 
                        self.settings['name_of_protein'], x > thresh))
                sig_enriched = aUtils.significantly_enriched(
                  main_enriches, zthresh=2., scale='log')
                if any(sig_enriched):
                    xs, ys = zip(*[(x, y) for x, y, e in
                      zip(main_enriches, other_enriches, sig_enriched) if e])
                    r, p = scipy.stats.pearsonr(xs, ys)
                fig2.suptitle('Scatter plot of enrichments %imers for '
                  '[protein]=%g or %g' %
                  (k, lib_main.conc, lib_other.conc))
                sp.set_xlabel('R values from ' + self.settings['experiment_name'])
                sp.set_ylabel('R values from ' + lib_other.expname)
                sp.set_aspect(1)
                simpleaxis(sp)
                fig2.savefig(os.path.join(self.rdir,
                  'comparisons',
                  'enrichment_correlation.%s.%.1f.nM.pdf'
                  % (other_exp.settings['experiment_name'],
                  lib_main.conc)))
        
        # Scatter: SKA
        if True:
            k=6
            print 'comparing SKA'
            for lib_main, lib_other in itertools.product(self.plibs,other_exp.plibs):
                fig3 = plt.figure()
                sp = fig3.add_subplot(1, 1, 1)
                if not lib_main.conc in [40,121,365,120]:
                    continue
                if not lib_other.conc in [40, 121,365,120]:
                    continue
                if not aUtils.within_xpercent(lib_main.conc, lib_other.conc, 100):
                    continue
                print 'comparing ', lib_main.conc, lib_other.conc
                main_kmer_libfracs = lib_main.k2stream_libfracs[k]
                lib_other.do_stream_count()
                other_kmer_libfracs = lib_other.k2stream_libfracs[k]
                for x, y, kmer in zip(
                  main_kmer_libfracs, other_kmer_libfracs, yield_kmers(k)):
                    sp.plot(x,y,'o',color=aColors.protein_colors(kmer, self.settings['name_of_protein'], x > 5./(4**k)))   
                sp.set_xlabel(r'SKA $F_{i}$ for ' + self.settings['experiment_name'] + ' [RBP]=%g' % lib_main.conc)
                sp.set_ylabel(r'SKA $F_{i}$ for ' + lib_other.expname + ' [RBP]=%g' % lib_other.conc)
                sp.set_aspect(1.)
                simpleaxis(sp)
                out_file = os.path.join(self.rdir, 'comparisons', 'SKA_comp_%s.%1.f.%.1f.nM.pdf' % (other_exp.settings['experiment_name'], lib_main.conc, lib_other.conc))
                print 'saving ', out_file
                fig3.savefig(out_file)
        # Scatter: B(enrich)
        k = 6
        kmers = [kmer for kmer in yield_kmers(k)]
        for lib_main, lib_other in itertools.product(self.plibs,
          other_exp.plibs):
            fig2 = plt.figure()
            sp = fig2.add_subplot(1, 1, 1)
            if not aUtils.within_xpercent(lib_main.conc, lib_other.conc, 40):
                continue
            main_kmer2enrich = lib_main.get_naive_enrichment_dict(k)
            other_kmer2enrich = lib_other.get_naive_enrichment_dict(k)
            main_enriches = [main_kmer2enrich[kmer] for kmer in kmers]
            other_enriches = [other_kmer2enrich[kmer] for kmer in kmers]
            thresh = numpy.mean(main_enriches) + 2 * numpy.std(main_enriches)
            main_read_len = self.settings['read_len']
            other_read_len = other_exp.settings['read_len']
            for x, y, kmer in zip(main_enriches, other_enriches, yield_kmers(k)):
                if chris_formula(x, k, main_read_len)< 0 or\
                  chris_formula(y, k, other_read_len) < 0:
                    continue
                sp.loglog(chris_formula(x, k, main_read_len),
                          chris_formula(y, k, other_read_len),
                          '.', 
                          color=aColors.protein_colors(
                            kmer, 
                            self.settings['name_of_protein'], 
                            x > thresh))
            sig_enriched = aUtils.significantly_enriched(
              main_enriches, zthresh=2., scale='log')
            if any(sig_enriched):
                xs, ys = zip(*[(x, y) for x, y, e in
                  zip(main_enriches, other_enriches, sig_enriched) if e])
                r, p = scipy.stats.pearsonr(xs, ys)
            fig2.suptitle('Scatter plot of B values %imers for '
              '[protein]=%g or %g' %
              (k, lib_main.conc, lib_other.conc))
            sp.set_xlabel('B values for ' + self.settings['experiment_name'])
            sp.set_ylabel('B values for ' + lib_other.expname)
            sp.set_aspect(1)
            sp.set_xlim([1,1000])
            sp.set_ylim([1,1000])
            simpleaxis(sp)
            fig2.savefig(os.path.join(self.rdir,
              'comparisons',
              'Bvalues.%s.%.1f.nM.pdf'
              % (other_exp.settings['experiment_name'],
              lib_main.conc)))

    def create_matlab_inputs(self):
        """
        makes the file that matlab needs to calculate Kds
        """
        print 'Creating MATLAB inputs'
        for k in self.settings['ks_for_matlab']:
            assert isinstance(k, int)
            kmer_list_file = os.path.join(self.rdir,
              'matlab', 'kmer_list.%i.m' % k)
            kmer_count_file = os.path.join(self.rdir,
              'matlab', 'kmer_counts.%i.m' % k)
            kmer_list_f = open(kmer_list_file, 'w')
            kmer_count_f = open(kmer_count_file, 'w')
            for kmer_i, kmer in enumerate(yield_kmers(k)):
                kmer_list_f.write(kmer + '\n')
                kmer_count_f.write('\t'.join(
                  [str(int(lib.k2counts[k][kmer_i]))
                  for lib in [self.input_lib] + self.plibs]) + '\n')
            kmer_list_f.close()
            kmer_count_f.close()
        protein_conc_file = os.path.join(self.rdir, 'matlab', 'protein_conc.m')
        protein_conc_f = open(protein_conc_file, 'w')
        protein_conc_f.write('\t'.join(
          [str(lib.conc)
          for lib in [self.input_lib] + self.plibs]) + '\n')
        protein_conc_f.close()

    def plot_num_counts(self):
        """
        creates a pdf plot with the number of reads in each multiplexed
        library. Requires the barcode_log.txt file (created during splitting)
        to exist.
        """
        fig2 = plt.figure()
        sp = fig2.add_subplot(211)
        sp.bar(range(len(self.libs)),
          [lib.num_reads for lib in self.libs])
        labels = [lib.get_full_label() for lib in self.libs]
        sp.set_xticks([i + 0.4 for i in range(len(self.libs))])
        sp.set_xticklabels(labels, rotation=90)
        fig2.savefig(self.rdir + '/plots/num_reads.pdf')

    def plot_k_picker(self):
        """
        for each protein library
        plots the kpicker plot (such as it is)
        """
        for lib in self.plibs:
            lib.plot_next_kmer_frac()
        # determine most enriched library
        k_for_max_enrichment = 6
        most_enriched_lib = self.plibs[0]
        best_enrich = 1.0
        for lib in self.plibs:
            max_enrich = lib.get_max_enrichment(k_for_max_enrichment)
            if max_enrich > best_enrich:
                most_enriched_lib = lib
                best_enrich = max_enrich
        fig1 = plt.figure()
        sp = fig1.add_subplot(111)
        top_vs_next = []
        for k in self.settings['ks_to_test_naive']:
            top_vs_next.append(most_enriched_lib.get_top_seq_vs_next(k, 2))
        sp.bar([k - 0.4 for k in self.settings['ks_to_test_naive']],
          top_vs_next)
        fig1.savefig(os.path.join(self.rdir, 'plots', 'best_v_next.pdf'))
        self.most_enriched_lib = most_enriched_lib

    def plot_concordance_to_most_enriched(self):
        """
        plots how well this library agrees with the most enriched library
        """
        most_enriched_lib = self.most_enriched_lib
        num_top_kmers_to_compare = int(self.settings.get(
          'num_top_kmers_to_compare', 50))
        k = self.settings.get('k_for_concordance_to_most_enriched', 7)
        fig1 = plt.figure()
        sp = fig1.add_subplot(211)
        for i, lib in enumerate(self.plibs):
            color = 'r' if lib is most_enriched_lib else '#19B271'
            sp.bar(i - 0.4,
              lib.compare_top_kmers(k,
              most_enriched_lib, num_top_kmers_to_compare), color=color)
        sp.set_ylabel('pearson corr of top '
          '%i %imers' % (num_top_kmers_to_compare, k))
        sp.set_xticks(range(len(self.plibs)))
        sp.set_xticklabels(
          [lib.get_full_label() for lib in self.plibs], rotation=90)
        fig1.savefig(os.path.join(self.rdir, 'plots',
          'comparison_to_top_lane.pdf'))

    def count_controlled_second_motif(self):
        """
        Counts the number of reads with two known motifs.
        Also counts the number of reads with one motif and a decoy
        Summarizes results in a table
        """
        if not all([lib.has_double_motif() for lib in self.libs]):
            self.decoys = list(aUtils.generate_decoys(self.known_motif, 25))
            for lib, decoy in itertools.product(self.libs, self.decoys):
                lib.count_double_occurrence(decoy)
            cPickle.dump(self.decoys, open(
              os.path.join(self.rdir, 'analyses', 'decoys.pkl'), 'w'))
            for lib in self.libs:
                lib.save_double_motif()
        else:
            self.decoys = cPickle.load(open(
              os.path.join(self.rdir, 'analyses', 'decoys.pkl')))
        table_f = self.get_rdir_fhandle('analyses',
          'controlled_multi_motif.count.txt')
        self.make_table_header(table_f)
        # double motif
        table_f.write('double_motif')
        for lib in self.libs:
            table_f.write('\t%f' % lib.motif_count[2])
        for i, decoy in enumerate(self.decoys):
            table_f.write('\nmotif & decoy:%s' % decoy)
            for lib in self.libs:
                table_f.write('\t%i' % lib.motif_and_decoy[decoy])
        table_f.close()

    def count_motif_occurences(self):
        """
        Counts how many times the known motif occurs in each library
        and summarizes results in table
        """
        for lib in self.libs:
            assert lib.split_reads_exist()
            lib.count_motif_individual()
        table_f = self.get_rdir_fhandle('analyses',
          'motif_occurence.count.txt')
        self.make_table_header(table_f)
        most_motifs = max([2] + [max(lib.motif_count) for lib in self.libs])
        for num_motif in range(most_motifs + 1):
            table_f.write(str(num_motif))
            for lib in self.libs:
                table_f.write('\t%i' % lib.motif_count[num_motif])
            table_f.write('\n')
        table_f.close()
        table_f = self.get_rdir_fhandle('analyses', 'motif_occurence.P.txt')
        self.make_table_header(table_f)
        for count in range(most_motifs + 1):
            table_f.write(str(count))
            for lib in self.libs:
                table_f.write('\t%f' % lib.motif_count_frac(count))
            table_f.write('\n')
        table_f.close()

    def make_table_header(self, of, include_library=False):
        """
        takes a file handle and writes a good header for it such that
        each lane is a column.
        """
        of.write('#')
        for lib in self.libs:
            of.write('\t' + lib.barcode)
        of.write('\n[%s]' % self.settings['name_of_protein'])
        for lib in self.libs:
            of.write('\t%0.3f' % lib.conc)
        of.write('\nwashes')
        for lib in self.libs:
            of.write('\t%i' % lib.washes)
        of.write('\n[poly IC]')
        for lib in self.libs:
            of.write('\t%i' % lib.polyIC)
        of.write('\nT (C)')
        for lib in self.libs:
            of.write('\t%s' % lib.temperature)
        of.write('\n')

    def make_counts_tables(self):
        """
        For each k
        Makes a table of the naive counts of each kmer.
        """
        print 'making counts tables'
        for k in self.settings['ks_to_test_naive']:
            print 'writing counts table for %i' % k
            of = self.get_rdir_fhandle('tables', 'naive_counts.%i.txt' % k)
            of_sorted = self.get_rdir_fhandle(
              'tables/naive_counts_sorted.%i.txt' % k)
            self.make_table_header(of)
            self.make_table_header(of_sorted)
            for i in range(4 ** k):
                kmer_str = get_kmer_from_index(k, i)
                of.write('\n' + kmer_str)
                for lib in self.libs:
                    of.write('\t%i' %
                     lib.k2counts[k][i])
            # sorted iteration through kemrs
            for i in self.naively_sorted_kmers[k]:
                kmer_str = get_kmer_from_index(k, i)
                of_sorted.write('\n' + kmer_str)
                for lib in self.libs:
                    of_sorted.write('\t%i' %
                      lib.k2counts[k][i])
            of.close()

    def make_motif_libfrac_table(self):
        """
        Makes a table of the library fraction for each kmer
        uses the length of the known motif for k.
        """
        print 'writing libfrac table'
        of = self.get_rdir_fhandle('tables', 'motif_libfrac.txt')
        of.write('barcode\tconcentration\t[poly_IC]\t'
          'washes\tlibfrac of known seq: %s\n' %
          self.known_motif)
        motif = self.known_motif
        motif_k = len(motif)
        motif_ind = get_index_from_kmer(motif)
        for lib in self.libs:
            if not motif_k in lib.k2libfrac:
                return
            libfrac_motif = lib.k2libfrac[motif_k][motif_ind]
            of.write('%s\t%0.3f\t%0.3f\t%i\t%f\n' %
              (lib.barcode, lib.conc, lib.polyIC, lib.washes, libfrac_motif))
        of.close()


    def plot_enrichment_humps(self):
        """
        plots the enrichment of the top kmers into a pdf
        """
        for k in self.settings['ks_to_test_naive']:
            fig1 = plt.figure()
            sp = fig1.add_subplot(221)
            sp_rep = fig1.add_subplot(222)
            for kmer_i in self.naively_sorted_kmers[k][0:
              self.settings['num_kmers_for_enrichment_humps']]:
                enrichments = [lib.k2enrichments[k][kmer_i]
                                    for lib in self.plibs]
                sp.plot(range(len(self.plibs)),
                  enrichments, 'o-', color=aColors.dodo_grey)
            if k == len(self.settings['known_motif']):
                enrichments = [lib.k2enrichments[k][
                  get_index_from_kmer(self.settings['known_motif'])]
                  for lib in self.plibs]
                sp.plot(range(len(self.plibs)), enrichments, 'o-', color=aColors.bacon)
            else:
                enrichments = []
            if k in self.replacely_sorted_kmers:
                if self.settings['replaceative']:
                    for kmer_i in self.replacely_sorted_kmers[k][0:
                      self.settings['num_kmers_for_enrichment_humps']]:
                        enrichments = [lib.k2rep_enrichments_arr[k][kmer_i]
                          for lib in self.plibs]
                        sp_rep.plot(range(len(self.plibs)),
                          enrichments, 'o-', color='b')
                if k == len(self.settings['known_motif']):
                    enrichments = [lib.k2rep_enrichments_arr[k][
                      get_index_from_kmer(self.settings['known_motif'])]
                      for lib in self.plibs]
                    sp_rep.plot(range(len(self.plibs)),
                      enrichments, 'o-', color='r')
            sp.set_xticks(range(len(self.plibs)))
            sp_rep.set_xticks(range(len(self.plibs)))
            labels = [lib.get_full_label() for lib in self.libs]
            sp.set_xticklabels(labels, rotation=45)
            sp_rep.set_xticklabels(labels, rotation=45)
            #sp.set_title('Naive')
            sp_rep.set_title('Iterative')
            #sp.set_xlabel('library')
            sp.set_ylabel('R value')
            fig1.suptitle(self.settings['experiment_name'] + (' k=%i' % k))
            print 'enrichment_humps.%i.pdf' % k
            fig1.savefig(os.path.join(self.rdir,
              'plots', 'enrichment_humps.%i.pdf' % k))
            fig1.savefig(os.path.join(self.rdir,
              'plots/enrichment_humps.%i.png' % k))

    def plot_enrichment_v_stream(self):
        ks_in_both = set(self.settings['ks_to_test_naive']) & set(self.settings['ks_to_test_streaming'])
        for k, lib in itertools.product(ks_in_both, self.plibs):
            fig1 = plt.figure(figsize=(5,7))
            sp = fig1.add_subplot(211)
            out_plot = os.path.join(self.rdir,
              'plots', 'count_compare.naive_stream.%i.%g' % (k, lib.conc))
            xs = lib.k2enrichments[k]
            ys = lib.k2stream_weights[k]
            assert len(xs) == len(ys)
            sp.loglog(xs, ys, '.')
            sp.set_xlabel('naive enrichment')
            sp.set_ylabel('stream weight')
            sp.set_title('protein: %s, conc: %g' %
              (self.settings['name_of_protein'], lib.conc))
            print out_plot
            aUtils.save_fig(fig1, out_plot)


    def plot_enrichment_v_presence(self):
        ks_in_both = set(self.settings['ks_to_test_naive']) & set(self.settings['ks_to_test_streaming'])
        for k, lib in itertools.product(ks_in_both, self.plibs):
            fig1 = plt.figure(figsize=(5,7))
            sp = fig1.add_subplot(211)
            out_plot = os.path.join(self.rdir,
              'plots', 'count_compare.naive_presence.%i.%g' % (k, lib.conc))
            xs = lib.k2enrichments[k]
            ys = lib.k2presence_counts[k]
            assert len(xs) == len(ys)
            try:
                sp.loglog(xs, ys, '.')
            except:
                print lib.barcode, k
                sys.exit()
            sp.set_xlabel('naive enrichment')
            sp.set_ylabel('Presence')
            sp.set_title('protein: %s, conc: %g' %
              (self.settings['name_of_protein'], lib.conc))
            print out_plot
            aUtils.save_fig(fig1, out_plot)

    def plot_presence_v_stream(self):
        for k, lib in itertools.product(self.settings['ks_to_test_streaming'], self.plibs):
            fig1 = plt.figure(figsize=(5,7))
            sp = fig1.add_subplot(211)
            out_plot = os.path.join(self.rdir,
              'plots', 'count_compare.presence_stream.%i.%g' % (k, lib.conc))
            xs = lib.k2presence_counts[k]
            ys = lib.k2stream_weights[k]
            assert len(xs) == len(ys)
            sp.loglog(xs, ys, '.')
            sp.set_xlabel('presence')
            sp.set_ylabel('Stream weight')
            sp.set_title('protein: %s, conc: %g' %
              (self.settings['name_of_protein'], lib.conc))
            print out_plot
            aUtils.save_fig(fig1, out_plot)


    def plot_enrichment_humps_with_secondary(self):
        """
        plots the enrichment of the top kmers into a pdf
        """
        xvals = [max(0, math.log(lib.conc + .00001)) for lib in self.plibs]
        for k in self.settings['ks_to_test_naive']:
            fig1 = plt.figure(figsize=(5,7))
            sp = fig1.add_subplot(211)
            kmers_of_interest = [kmer 
              for kmer in self.settings['motifs_of_interest'] if len(kmer) == k]
            #colors = self.get_enrich_hump_colors(k, kmers, [False] * len(kmers))
            legend_handles = []
            rand = {6:500, 7:3015, 5:129, 4:31, 8:8000, 9:32000, 10: 128000, 3:7}
            for kmer_i in range(33, 4 ** k, rand[k]):
                enrichments = [lib.k2enrichments[k][kmer_i]
                                    for lib in self.plibs]
                sp.plot(xvals,
                  enrichments, 'o-', color=aColors.background_motif)
            for kmeri, kmer in enumerate(kmers_of_interest):
                enrichments =\
                  [lib.k2enrichments[k][get_index_from_kmer(kmer)] 
                  for lib in self.plibs]
                if 'FOX' in self.settings['name_of_protein'].upper():
                    color = aColors.fox_colors(kmer, True)
                    #sp.set_ylim([0,50])
                elif 'MBNL' in self.settings['name_of_protein'].upper():
                    color = aColors.mbnl_colors(kmer, True)
                    color = aColors.colors8[kmeri]
                elif 'CUGBP' in self.settings['name_of_protein'].upper():
                    #color = aColors.cugbp_colors(kmer, True)
                    color= aColors.newcolors6[kmeri]
                    color = aColors.colors8[kmeri]
                else:
                    color = aColors.ERSP_colors(kmer, True)
                l, = sp.plot(xvals,
                        enrichments, 
                        'o-', 
                        color=color, 
                        label=kmer)
                legend_handles.append(l)
            labels = [lib.get_full_label() for lib in self.libs]
            self.sorted_labels = labels
            sp.set_xticks(xvals)
            sp.set_xticklabels(labels, rotation=45)
            sp.set_ylabel('R value')
            simpleaxis(sp)
            if 'CUGBP' in self.settings['name_of_protein'].upper() or\
              'MBNL' in self.settings['name_of_protein'].upper():
                if k == 6 or k == 7:
                    sp.set_ylim([0, 10])
                    sp.set_yticks([0,1,2,4,6,8, 10])
            fig1.suptitle(self.settings['experiment_name'] + (' k=%i' % k),
              verticalalignment='baseline')
            out_plot = os.path.join(self.rdir,
              'plots', 'enrichment_humps.kmersofinterest.%i.nolegend' % k)
            aUtils.save_fig(fig1, out_plot)
            sp.legend(legend_handles,
                      map(rna, kmers_of_interest),
                      mode='expand', 
                      loc=2,
                      borderaxespad=-2.)
            out_plot = os.path.join(self.rdir,
              'plots', 'enrichment_humps.kmersofinterest.%i.withlegend' % k)
            print 'saving', out_plot
            aUtils.save_fig(fig1, out_plot)

    def plot_relativeAffinity_humps_with_secondary(self):
        """
        plots the enrichment of the top kmers into a pdf
        """
        print 'plotting relative affinity humps'
        xvals = [max(0, math.log(lib.conc + .00001)) for lib in self.plibs]
        for k in self.settings['ks_to_test_naive']:
            fig1 = plt.figure(figsize=(5,7))
            sp = fig1.add_subplot(211)
            kmers_of_interest = [kmer 
              for kmer in self.settings['motifs_of_interest'] if len(kmer) == k]
            legend_handles = []
            for kmeri, kmer in enumerate(kmers_of_interest):
                enrichments =\
                  [lib.k2enrichments[k][get_index_from_kmer(kmer)] 
                  for lib in self.plibs]
                if 'FOX' in self.settings['name_of_protein'].upper():
                    color = aColors.fox_colors(kmer, True)
                elif 'MBNL' in self.settings['name_of_protein'].upper():
                    color = aColors.mbnl_colors(kmer, True)
                    color = aColors.colors8[kmeri]
                elif 'CUGBP' in self.settings['name_of_protein'].upper():
                    color= aColors.newcolors6[kmeri]
                    color = aColors.colors8[kmeri]
                else:
                    color = aColors.ERSP_colors(kmer, True)
                relative_affinities = [chris_formula(e, k, self.settings['read_len'])
                  for e in enrichments]
                try:
                    l, = sp.semilogy(xvals,
                            relative_affinities,
                            'o-', 
                            color=color, 
                            label=kmer)
                    legend_handles.append(l)
                except:
                    pass
            labels = [lib.get_full_label() for lib in self.libs]
            self.sorted_labels = labels
            sp.set_xticks(xvals)
            sp.set_xticklabels(labels, rotation=45)
            sp.set_ylabel('Relative Affinities')
            simpleaxis(sp)
            fig1.suptitle(self.settings['experiment_name'] + (' k=%i' % k),
              verticalalignment='baseline')
            out_plot = os.path.join(self.rdir,
              'plots', 'relative_affinity_humps.kmersofinterest.%i.nolegend' % k)
            aUtils.save_fig(fig1, out_plot)
            sp.legend(legend_handles,
                      map(rna, kmers_of_interest),
                      mode='expand', 
                      loc=2,
                      borderaxespad=-2.)
            out_plot = os.path.join(self.rdir,
              'plots', 'relative_affinity_humps.kmersofinterest.%i.withlegend' % k)
            print 'saving', out_plot
            aUtils.save_fig(fig1, out_plot)

    def plot_enrichment_humps_with_secondary_stream(self):
        """
        plots the enrichment of the top kmers into a pdf
        """
        xvals = [max(0, math.log(lib.conc + .00001)) for lib in self.plibs]
        for k in self.settings['ks_to_test_streaming']:
            fig1 = plt.figure(figsize=(5,7))
            sp = fig1.add_subplot(211)
            kmers_of_interest = [kmer 
              for kmer in self.settings['motifs_of_interest'] if len(kmer) == k]
            #colors = self.get_enrich_hump_colors(k, kmers, [False] * len(kmers))
            legend_handles = []
            rand = {6:500, 7:3015, 5:129, 4:31, 8:8000}
            for kmer_i in range(33, 4 ** k, rand[k]):
                enrichments = [lib.k2stream_weights[k][kmer_i] / 
                  numpy.sum(lib.k2stream_weights[k])
                                    for lib in self.plibs]
                sp.plot(xvals,
                  enrichments, 'o-', color=aColors.background_motif)
            for kmeri, kmer in enumerate(kmers_of_interest):
                enrichments =\
                  [lib.k2stream_weights[k][get_index_from_kmer(kmer)] / 
                  numpy.sum(lib.k2stream_weights[k]) 
                  for lib in self.plibs]
                if 'FOX' in self.settings['name_of_protein'].upper():
                    color = aColors.fox_colors(kmer, True)
                elif 'MBNL' in self.settings['name_of_protein'].upper():
                    color = aColors.mbnl_colors(kmer, True)
                    color = aColors.colors8[kmeri]
                elif 'CUGBP' in self.settings['name_of_protein'].upper():
                    color= aColors.newcolors6[kmeri]
                    color = aColors.colors8[kmeri]
                else:
                    color = aColors.ERSP_colors(kmer, True)
                l, = sp.plot(xvals,
                        enrichments, 
                        'o-', 
                        color=color, 
                        label=kmer)
                legend_handles.append(l)
            labels = [lib.get_full_label() for lib in self.libs]
            self.sorted_labels = labels
            sp.set_xticks(xvals)
            sp.set_xticklabels(labels, rotation=45)
            sp.set_ylabel('RBNS R value')
            simpleaxis(sp)
            if 'CUGBP' in self.settings['name_of_protein'].upper() or\
              'MBNL' in self.settings['name_of_protein'].upper():
                if k == 6 or k == 7:
                    pass
            fig1.suptitle(self.settings['experiment_name'] + (' k=%i' % k),
              verticalalignment='baseline')
            out_plot = os.path.join(self.rdir,
              'plots', 'stream_humps.kmersofinterest.%i.nolegend' % k)
            aUtils.save_fig(fig1, out_plot)
            sp.legend(legend_handles,
                      map(rna, kmers_of_interest),
                      mode='expand', 
                      loc=2,
                      borderaxespad=-2.)
            out_plot = os.path.join(self.rdir,
              'plots', 'stream_humps.kmersofinterest.%i.withlegend' % k)
            print 'saving', out_plot
            aUtils.save_fig(fig1, out_plot)

    def plot_enrichment_humps_with_secondary_presence(self):
        """
        plots the enrichment of the top kmers into a pdf
        """
        xvals = [max(0, math.log(lib.conc + .00001)) for lib in self.plibs]
        for k in self.settings['ks_to_test_streaming']:
            fig1 = plt.figure(figsize=(5,7))
            sp = fig1.add_subplot(211)
            kmers_of_interest = [kmer 
              for kmer in self.settings['motifs_of_interest'] if len(kmer) == k]
            legend_handles = []
            rand = {6:500, 7:3015, 5:129, 4:31, 8:8000}
            for kmer_i in range(33, 4 ** k, rand[k]):
                enrichments = [lib.k2presence_counts[k][kmer_i]
                                    for lib in self.plibs]
                sp.plot(xvals,
                  enrichments, 'o-', color=aColors.background_motif)
            for kmeri, kmer in enumerate(kmers_of_interest):
                enrichments =\
                  [lib.k2presence_counts[k][get_index_from_kmer(kmer)] 
                  for lib in self.plibs]
                if 'FOX' in self.settings['name_of_protein'].upper():
                    color = aColors.fox_colors(kmer, True)
                elif 'MBNL' in self.settings['name_of_protein'].upper():
                    color = aColors.mbnl_colors(kmer, True)
                    color = aColors.colors8[kmeri]
                elif 'CUGBP' in self.settings['name_of_protein'].upper():
                    color= aColors.newcolors6[kmeri]
                    color = aColors.colors8[kmeri]
                else:
                    color = aColors.ERSP_colors(kmer, True)
                l, = sp.plot(xvals,
                        enrichments, 
                        'o-', 
                        color=color, 
                        label=kmer)
                legend_handles.append(l)
            labels = [lib.get_full_label() for lib in self.libs]
            self.sorted_labels = labels
            sp.set_xticks(xvals)
            sp.set_xticklabels(labels, rotation=45)
            sp.set_ylabel('fraction with kmer')
            simpleaxis(sp)
            if 'CUGBP' in self.settings['name_of_protein'].upper() or\
              'MBNL' in self.settings['name_of_protein'].upper():
                if k == 6 or k == 7:
                    pass
            fig1.suptitle(self.settings['experiment_name'] + (' k=%i' % k),
              verticalalignment='baseline')
            out_plot = os.path.join(self.rdir,
              'plots', 'presence_humps.kmersofinterest.%i.nolegend' % k)
            aUtils.save_fig(fig1, out_plot)
            sp.legend(legend_handles,
                      map(rna, kmers_of_interest),
                      mode='expand', 
                      loc=2,
                      borderaxespad=-2.)
            out_plot = os.path.join(self.rdir,
              'plots', 'presence_humps.kmersofinterest.%i.withlegend' % k)
            print 'saving', out_plot
            aUtils.save_fig(fig1, out_plot)


    def verify_concordances(self):
        """
        verifies that the libraries are at least similar to each other.
        """
        print 'calculating concordances between libraries'
        lf_concordance_of = self.get_rdir_fhandle('QC/concordance_libfrac.txt')
        enrich_concordance_of =\
          self.get_rdir_fhandle('QC/concordance_enrich.txt')
        concordance_head = 'barcode 1\tbarcode 2\tconcentration1'\
          '\tconcentration 2\tPearson r\t'\
          'Pearson p\tSpearman Rho\tSpearman P\n'
        lf_concordance_of.write(concordance_head)
        enrich_concordance_of.write(concordance_head)
        for libi, libj in aUtils.pairwise(self.libs):
            pearsonr, pearsonp = self.get_libfrac_kmer_concordance(
              libi, libj, 'pearson')
            spearmanr, spearmanp = self.get_libfrac_kmer_concordance(
              libi, libj, 'spearman')
            lf_concordance_of.write(
              '%s\t%s\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\n' %
              (libi.barcode, libj.barcode, libi.conc, libj.conc,
              pearsonr, pearsonp, spearmanr, spearmanp))
            pearsonr, pearsonp = self.get_enrich_kmer_concordance(
              libi, libj, 'pearson')
            spearmanr, spearmanp = self.get_enrich_kmer_concordance(
              libi, libj, 'spearman')
            enrich_concordance_of.write(
              '%s\t%s\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\n' %
              (libi.barcode, libj.barcode, libi.conc, libj.conc,
              pearsonr, pearsonp, spearmanr, spearmanp))
            if spearmanr < 0.5:
                print 'low concordance warning', libi.barcode, libj.barcode
        lf_concordance_of.close()
        enrich_concordance_of.close()

    def get_enrich_kmer_concordance(self, libi, libj, test):
        """
        checks that two libaries are similar with regard to enrichment
        """
        k = self.settings['k_for_concordance_check']
        if libi.is_input or libj.is_input:
            return (1.0, 1.0)
        if test == 'pearson':
            #print self.experiment.proteinr
            print  self.settings['name_of_protein']
            return scipy.stats.pearsonr(
              libi.k2enrichments[k], libj.k2enrichments[k])
        elif test == 'spearman':
            return scipy.stats.spearmanr(
              libi.k2enrichments[k], libj.k2enrichments[k])
        else:
            raise ValueError('Unrecognized test %s' % test)

    def get_libfrac_kmer_concordance(self, libi, libj, test):
        """
        checks that two libaries are similar with regard to library fraction
        """
        k = self.settings['k_for_concordance_check']
        if test == 'pearson':
            return scipy.stats.pearsonr(libi.k2libfrac[k], libj.k2libfrac[k])
        elif test == 'spearman':
            return scipy.stats.spearmanr(libi.k2libfrac[k], libj.k2libfrac[k])
        else:
            raise ValueError('Unrecognized test %s' % test)

    def sort_all_kmers_by_enrichment(self):
        """
        sorts kmers by how enriched they are relative to the library
        in the specified libraries.
        """
        print 'sorting kmers by enrichment'
        self.naively_sorted_kmers = {}
        self.replacely_sorted_kmers = {}
        for k in self.settings['ks_to_test_naive']:
            self.naively_sorted_kmers[k] =\
              self.sort_kmers_by_enrichment(k, 'naive')
        if not self.settings['replaceative']:
            return
        for k in self.settings['ks_to_test_replaceative']:
            self.replacely_sorted_kmers[k] =\
              self.sort_kmers_by_enrichment(k, 'replace')

    def sort_kmers_by_enrichment(self, k, method='naive'):
        """
        sorts the kmers based on how enriched they are in the
        protein libraries.
        returns an array of the kmers (as indexes not strings) in sorted order
        """
        summed_enrichments = numpy.zeros(4 ** k)
        for barcode in self.settings['which_barcodes_to_sort_based_on']:
            lib = self.barcode2lib[barcode]
            if method == 'naive':
                print lib.k2enrichments.keys(),k
                summed_enrichments += lib.k2enrichments[k]
            elif method == 'replace':
                summed_enrichments += lib.k2rep_enrichments_arr[k]
        kmers = range(4 ** k)
        summed_enrich_kmers = zip(summed_enrichments, kmers)
        summed_enrich_kmers.sort(reverse=True)
        top_enrich, sorted_kmers = zip(*summed_enrich_kmers)
        return sorted_kmers

    def calculate_all_enrichments(self):
        """
        sorts all the kmers.
        """
        print 'calculating enrichments'
        for lib in self.plibs:
            lib.calculate_enrichment()

    def make_enrichment_table(self):
        """
        Makes a table of the enrichments for the values of k
        Also, writes a pkl with the same information
        """
        for k in self.settings['ks_to_test_naive']:
            print 'Writing Naive Enrichment Table for k=%i' % k
            #write table
            table_f = self.get_rdir_fhandle('tables/enrich_naive.%i.txt' % k)
            self.make_table_header(table_f)
            for kmer_i in self.naively_sorted_kmers[k]:
                assert isinstance(kmer_i, int)
                kmer = get_kmer_from_index(k, kmer_i)
                table_f.write(kmer)
                for lib in self.plibs:
                    table_f.write('\t%f' % lib.k2enrichments[k][kmer_i])
                table_f.write('\n')
            table_f.close()
            print 'Writing Naive Enrichment pkl for k=%i' % k
            barcode2kmer2enrich = dict()
            for lib in self.plibs:
                barcode2kmer2enrich[lib.barcode]\
                  = lib.get_naive_enrichment_dict(k)
            pkl_f = self.get_rdir_fhandle('tables/enrich_naive.%i.pkl' % k)
            cPickle.dump(barcode2kmer2enrich, pkl_f)
        print 'tables made'

        if not self.settings['replaceative']:
            return
        for k in self.settings['ks_to_test_replaceative']:
            table_f = self.get_rdir_fhandle('tables/enrich_replace.%i.txt' % k)
            self.make_table_header(table_f)
            for kmer_i in self.replacely_sorted_kmers[k]:
                table_f.write(get_kmer_from_index(k, kmer_i))
                for lib in self.plibs:
                    table_f.write('\t%f' %
                      lib.k2rep_enrichments_arr[k][kmer_i])
                table_f.write('\n')
            table_f.close()

    def make_SKA_table(self):
        for k in self.settings['ks_to_test_streaming']:
            print 'Writing SKA libfrac Table for k=%i' % k
            table_f = self.get_rdir_fhandle('tables/SKA_libfrac.%i.txt' % k)
            self.make_table_header(table_f)
            for kmer_i in self.naively_sorted_kmers[k]:
                assert isinstance(kmer_i, int)
                kmer = get_kmer_from_index(k, kmer_i)
                table_f.write(kmer)
                for lib in self.plibs:
                    table_f.write('\t%g' % lib.k2stream_libfracs[k][kmer_i])
                table_f.write('\n')
            table_f.close()

    def make_enrichment_hist_stacked(self):
        """
        make a histogram of enrichments for each barcode
        """
        protein = self.settings['name_of_protein'].upper()
        for k in self.settings['ks_to_test_naive']:
            for lib in self.plibs:
                hist_file = os.path.join(self.rdir, 'plots',
                  'stacked.hist.%i.%0.0f' % (k, lib.conc))
                print 'working on ', hist_file
                fig1 = plt.figure(figsize=(6,5))
                fig2 = plt.figure(figsize=(6,5))
                sp = fig1.add_subplot(1, 1, 1) #, frame_on=False)
                sp2 = fig2.add_subplot(1, 1, 1) #, frame_on=False)
                simpleaxis(sp)
                simpleaxis(sp2)
                bin_edges = 200
                counts, bin_edges = numpy.histogram(lib.k2enrichments[k], bin_edges)
                thresh = numpy.mean(lib.k2enrichments[k]) + 2 * numpy.std(lib.k2enrichments[k])
                sp.axvline(thresh, color='black')
                sp2.axvline(thresh, color='black')
                highest = stacked_bar_kmers.plot_stack(sp, bin_edges, lib.k2enrichments[k], k, protein)
                highest2 = stacked_bar_kmers.plot_stack(sp2, 
                                             bin_edges, 
                                             lib.k2enrichments[k], 
                                             k, 
                                             protein, 
                                             scale='log')
                assert highest == highest2
                sp.set_xlabel('RBNS R value')
                sp.set_ylabel('number of %imers' % k)
                sp2.set_xlabel('RBNS R value')
                sp2.set_ylabel('number of %imers'% k)
                if k == 5 and 'FOX' in self.settings['name_of_protein'].upper():
                    supt = 'GCATG: %0.2f, GCACG: %0.2f'\
                      % (lib.k2enrichments[k][get_index_from_kmer('GCATG')], 
                      lib.k2enrichments[k][get_index_from_kmer('GCACG')])
                    fig1.suptitle(supt)
                    fig2.suptitle(supt)
                ytick_defaults = [0, 1, 5, 10, 50, 100, 500, 1000, 2000, 5000]
                y_ts_lin = ytick_defaults[0:aUtils.getBinIndex(
                  highest, ytick_defaults) + 2]
                y_ts_log = map(lambda yt: math.log(yt + 1), y_ts_lin)
                sp2.set_ylim([y_ts_log[0], y_ts_log[-1]])
                sp2.set_yticks(y_ts_log)
                sp2.set_yticklabels(map(str, y_ts_lin))
                sp2.set_xticks(sp2.get_xticks() + [1])
                if 'CUGBP' in self.settings['name_of_protein'].upper():
                    if k == 7:
                        sp.set_xlim([0, 10])
                        sp2.set_xlim([0, 10])
                        sp.set_xticks([0, 1, 2, 4, 6,8, 10])
                        sp2.set_xticks([0, 1, 2, 4, 6,8, 10])
                if k == 6 and 'FOX' in protein:
                    supt = 'TGCATG: %0.2f, TGCACG: %0.2f'\
                      % (lib.k2enrichments[k][get_index_from_kmer('TGCATG')], 
                      lib.k2enrichments[k][get_index_from_kmer('TGCACG')])
                    sp.set_xlim([0, 25])
                    sp2.set_xlim([0, 25])
                    sp.set_xticks([0,1,5,10,15,20,25])
                    sp2.set_xticks([0,1,5,10,15,20,25])
                    fig1.suptitle(supt)
                    fig2.suptitle(supt)
                aUtils.save_fig(fig1, hist_file)
                aUtils.save_fig(fig2, hist_file + '.log')
                self.set_legend_texts(sp)
                self.set_legend_texts(sp2)
                aUtils.save_fig(fig1, hist_file+'.withlegend')
                aUtils.save_fig(fig2, hist_file + '.withlegend.log')
        print 'STACKED HISTOGRAMS MADE'

    def set_legend_texts(self, sp):
        ylim = sp.get_ylim()
        xlim = sp.get_xlim()
        a = sp.bar([1e9], [1e9], color=aColors.primary_motif)
        b = sp.bar([1e9], [1e9], color=aColors.secondary_motif)
        c = sp.bar([1e9], [1e9], color=aColors.other_significant)
        d = sp.bar([1e9], [1e9], color=aColors.background_motif)
        sp.set_ylim(ylim)
        sp.set_xlim(xlim)
        sp.legend((a,b,c,d), self.get_legend_texts())


    def get_legend_texts(self):
        if 'CUGBP' in self.settings['name_of_protein'].upper():
            return ("contains two UGU's", "contains one UGU", "other significant", "unenriched")
        if 'FOX' in self.settings['name_of_protein'].upper() or 'REG' in self.settings['name_of_protein'].upper():
            return ("contains GCAUG", "contains GCACG", "other significant", "unenriched")
        if 'MBNL' in self.settings['name_of_protein'].upper():
            return ("contains YGCU", "contains YGCC", "other significant", "unenriched")
        if 'ESRP' in self.settings['name_of_protein'].upper():
            return ("contains GGTG", "contains GGT", "other significant", "unenriched")
        if 'U1' in self.settings['name_of_protein'].upper():
            return ("a", "b", "c", "d")

        raise ValueError('not workign')

    def make_enrichment_hist(self):
        """
        make a histogram of enrichments for each barcode
        """
        for k in self.settings['ks_to_test_naive']:
            for lib in self.plibs:
                hist_file = os.path.join(self.rdir, 'plots',
                  'hist.%i.%0.1f' % (k, lib.conc))
                fig1 = plt.figure()
                fig2 = plt.figure()
                sp = fig1.add_subplot(1, 1, 1) #, frame_on=False)
                sp2 = fig2.add_subplot(1, 1, 1) #, frame_on=False)
                simpleaxis(sp)
                simpleaxis(sp2)
                bin_edges = 200
                if not 'MBNL' in self.settings['name_of_protein'].upper():
                    sp.hist(lib.k2enrichments[k], bins=bin_edges,
                     facecolor=aColors.Little_green, edgecolor='none')
                counts, bin_edges = numpy.histogram(lib.k2enrichments[k], bin_edges)
                sp.axvline(1, color=aColors.fofo_grey, lw=1.75, linestyle='-.')
                sp2.axvline(1, color=aColors.fofo_grey, lw=1.75, linestyle='-.')
                if 'CUGBP' in self.settings['name_of_protein'].upper() and k >=4:
                    fracs_UGU = fracs_passing(k,
                      lib.k2enrichments[k], 
                      bin_edges, 
                      hasUGU)
                    for x, y, color in zip(bin_edges[1:],
                                           counts, 
                                           plt.cm.jet(fracs_UGU)):
                        if not y:
                            continue
                        sp.bar(x,y,
                          width=(bin_edges[1]-bin_edges[0]), 
                          facecolor=color,
                          edgecolor='none')
                        sp2.bar(x,
                          math.log(y + 1),
                          width=(bin_edges[1]-bin_edges[0]), 
                          facecolor=color,
                          edgecolor='none')
                elif 'MBNL' in self.settings['name_of_protein'].upper() and k >=4:
                    fracs_ygcy = fracs_passing(k,
                      lib.k2enrichments[k], 
                      bin_edges, 
                      YGCY)
                    for x, y, color in zip(bin_edges[1:], counts, plt.cm.jet(fracs_ygcy)):
                        if not y:
                            continue
                        sp.bar(x,y,
                          width=(bin_edges[1]-bin_edges[0]), 
                          facecolor=color,
                          edgecolor='none')
                        sp2.bar(x,
                          math.log(y + 1),
                          width=(bin_edges[1]-bin_edges[0]), 
                          facecolor=color,
                          edgecolor='none')
                else:
                    sp2.bar(bin_edges[1:], 
                      [math.log(count + 1) for count in counts],
                      width=(bin_edges[1]-bin_edges[0]), 
                      facecolor=aColors.Little_green, 
                      edgecolor='none')

                sp.set_xlabel('%imer RBNS R value' % k)
                sp.set_ylabel('number of %imers' % k)
                sp2.set_xlabel('%imer RBNS R value' % k)
                sp2.set_ylabel('number of %imers' % k)
                if k == 5 and 'FOX' in self.settings['name_of_protein'].upper():
                    supt = 'GCATG: %0.2f, GCACG: %0.2f'\
                      % (lib.k2enrichments[k][get_index_from_kmer('GCATG')], 
                      lib.k2enrichments[k][get_index_from_kmer('GCACG')])
                    fig1.suptitle(supt)
                    fig2.suptitle(supt)
                if k == 6 and 'FOX' in self.settings['name_of_protein'].upper():
                    supt = 'TGCATG: %0.2f, TGCACG: %0.2f'\
                      % (lib.k2enrichments[k][get_index_from_kmer('TGCATG')], 
                      lib.k2enrichments[k][get_index_from_kmer('TGCACG')])
                    fig1.suptitle(supt)
                    fig2.suptitle(supt)
                    y_ts_lin = [0, 1, 5, 10, 50, 100, 500, 1000, 3000]
                    y_ts_log = map(lambda yt: math.log(yt + 1), y_ts_lin)
                    sp2.set_yticks(y_ts_log)
                    sp2.set_yticklabels(map(str, y_ts_lin))
                if 'CUGBP' in self.settings['name_of_protein'].upper():
                    if k == 7:
                        y_ts_lin = [0, 1, 5, 10, 50, 100, 500, 1000, 3000]
                        y_ts_log = map(lambda yt: math.log(yt + 1), y_ts_lin)
                        sp2.set_yticks(y_ts_log)
                        sp2.set_yticklabels(map(str, y_ts_lin))
                        sp.set_xlim([0, 10])
                        sp2.set_xlim([0, 10])
                if 'FOX' in self.settings['name_of_protein'].upper():
                    if k == 5:
                        sp2.set_xlim([0, 15])
                    if k == 6:
                        sp.set_xlim([0, 25])
                        sp.set_xticks([0,1,5,10,15,20,25])
                        sp2.set_xticks([0,1,5,10,15,20,25])
                    if k == 7:
                        sp.set_xlim([0, 30])
                        sp2.set_xlim([0, 30])

                if not 'FOX' in self.settings['name_of_protein'].upper():
                    supt = '%s, max enrich: %0.2f' % (
                      self.settings['name_of_protein'],
                      max(lib.k2enrichments[k]))
                    fig1.suptitle(supt)
                    fig2.suptitle(supt)
                aUtils.save_fig(fig1, hist_file)
                aUtils.save_fig(fig2, hist_file + '.log')
        print 'HISTOGRAMS MADE'

    def calculate_all_libfracs(self):
        """
        calculates the librfrac for each k, barcode pair and method
        """
        print 'calculating libfracs'
        for lib in self.libs:
            lib.calculate_naive_libfracs()
            lib.check_libfrac_sum()

    def perform_kmer_counts(self):
        """
        Does all the kmer counting.
        """
        # naive
        if self.settings['naive_count']:
            naive_start = time.time()
            self.count_naive()
            naive_end = time.time()
            print 'Naive count took %0.3f seconds' % (naive_end - naive_start)
        # replaceative
        if self.settings['replaceative']:
            for lib, k in itertools.product(self.plibs,
              self.settings['ks_to_test_replaceative']):
                #lib.count_replaceative(k)
                pass
        #streaming
        if self.settings['streaming_counts']:
            self.count_stream()
        print 'counting presence'
        self.count_presence()
        if self.use_cluster:
            self.wait_for_jobs_to_complete()

    def wait_for_jobs_to_complete(self, sleep=1):
        """
        Waits for the counts to be completed on the qsub.
        """
        while True:
            outputs = [subprocess.Popen('qstat %i' % job_id,
              shell=True, stdout=subprocess.PIPE,
              stderr=subprocess.PIPE).communicate()
              for job_id in self.waiting_jobs]
            completed = map(lambda output: 'Unknown Job' in output[1], outputs)
            if all(completed):
                break
            else:
                print sum(map(lambda output:
                  0 if 'Unknown J' in output[1] else 1, outputs)), 'jobs left'
            time.sleep(sleep)
        self.waiting_jobs = []

    def run_rna_fold(self):
        """
        for each barcode and readfile checks that the folding energy
        has been calculated and if not calculates it (for each read)
        This can be launched on the cluster.
        After the folding energy has been calculated
        split the reads into bins based on folding energy.
        """
        for lib in self.libs:
            if self.use_cluster:
                self.waiting_jobs.append(lib.run_rna_fold())
            else:
                lib.run_rna_fold()
        if self.use_cluster:
            self.wait_for_jobs_to_complete()
        for lib in self.libs:
            if self.use_cluster:
                self.waiting_jobs.append(lib.split_by_structure())
            else:
                lib.split_by_structure()
        if self.use_cluster:
            self.wait_for_jobs_to_complete()

    def count_fold_split(self):
        """ counts the kmers in the fold split versions
        """
        for barcode, (energy_bin_i, energy_bin) in \
          itertools.product(self.barcodes, enumerate(self.energy_bins)):
            results_file = os.path.join(self.rdir, 'structure',
              '%s.%i.%0.1f_%0.1f.ncount.pkl' %
              (barcode,
              energy_bin_i,
              energy_bin[0],
              energy_bin[1]))
            if not aUtils.check_file_exists(results_file):
                if self.use_cluster:
                    command = 'python ~alexrson/snorelax/script.py '\
                      'bind_n_seq_pipeline.fold_counter '\
                      '%s %s %i %s %s %s 1> %s 2> %s' % \
                      (self.settings_file,
                      barcode,
                      energy_bin_i,
                      '%0.1f' % energy_bin[0],
                      '%0.1f' % energy_bin[1],
                      results_file,
                      self.rdir + '/structure/' + barcode + '.fc.out',
                      self.rdir + '/structure/' + barcode + '.fc.err')
                    script_options = {'nodes': '1', 'ppn': '2',
                      'outf': self.settings['experiment_name'] + '.sub_log',
                      'jobname': self.settings['experiment_name']
                      + '_' + barcode + '.count_split',
                      'queue': 'long', 'workingdir': self.rdir,
                      'command': command}
                    job_id = aUtils.launch(command, script_options)
                    self.waiting_jobs.append(job_id)
                else:
                    self.count_kmer_fold_split(
                      barcode, energy_bin_i, energy_bin, results_file)
            else:
                print 'split by fold results exist for: ', barcode, energy_bin
        if self.use_cluster:
            self.wait_for_jobs_to_complete()
        # load job results
        self.barcode2ebini2k2counts = collections.defaultdict(dict)
        for barcode, (energy_bin_i, energy_bin) in \
          itertools.product(self.barcodes, enumerate(self.energy_bins)):
            results_file = os.path.join(self.rdir, 'structure',
              '%s.%i.%0.1f_%0.1f.ncount.pkl' %
              (barcode,
              energy_bin_i,
              energy_bin[0],
              energy_bin[1]))
            if not aUtils.check_file_exists(results_file):
                raise ValueError(results_file + ' not there')
            print 'loading', results_file
            self.barcode2ebini2k2counts[barcode][energy_bin_i] =\
              cPickle.load(open(results_file))
        print 'full fold results loaded'

    def plot_fold_split_hump(self):
        """
        this is basically the enrichment humps plotter for the
        reads split by free energy

        Only plots the known motif
        """
        fig1 = plt.figure(figsize=(7, 9))
        for ebi, energy_bin in enumerate(self.energy_bins):
            sp = fig1.add_subplot(4, 1, 1)
            sp2 = fig1.add_subplot(4, 1, 3)
            motif = self.known_motif
            k = len(motif)
            motif_i = get_index_from_kmer(motif)
            vs = []
            for lib in self.plibs:
                lib_counts = self.barcode2ebini2k2counts[
                  self.input_barc][ebi][k][motif_i]
                counts = self.barcode2ebini2k2counts[
                  lib.barcode][ebi][k][motif_i]
                input_libfrac = float(lib_counts) / sum(
                  self.barcode2ebini2k2counts[
                  self.input_barc][ebi][k])
                this_libfrac = float(counts) / sum(self.barcode2ebini2k2counts[
                  lib.barcode][ebi][k])
                vs.append(this_libfrac / input_libfrac)
                sp.set_ylim([0, 40])
                sp2.set_ylim([0, 40])
            sp.plot(range(len(self.plibs)), vs,
              color=aColors.colors_red_to_black_6[ebi],
              linewidth=2)
            sp2.plot(range(len(self.plibs)), vs,
              color=aColors.colors_red_to_black_6[ebi],
              linewidth=2)
        sp.set_xticks(range(len(self.plibs)))
        sp.set_xticklabels(self.sorted_labels, rotation=45)
        sp2.set_xticks(range(len(self.plibs)))
        sp2.set_xticklabels(self.sorted_labels, rotation=45)
        legend_labels = [('%0.1f to %0.1f' % energy_bin).replace(
          '-40.0 to ', 'less than ')
          for energy_bin in self.energy_bins]
        sp.legend(legend_labels, loc=2)
        sp.set_ylabel('enrichment')
        sp2.set_ylabel('enrichment')
        fig1.suptitle(self.settings['experiment_name'])
        fig1.savefig(os.path.join(self.rdir, 'plots', 'FE.split.pdf'))

    def get_ebin_counts(self):
        """
        returns the total number of reads in each energy bin
        """
        ebin_counts_file = os.path.join(
          self.rdir, 'analyses', 'ebin_counts.pkl')
        if aUtils.check_file_exists(ebin_counts_file):
            k_lib_ebin2total_counts = cPickle.load(open(ebin_counts_file))
            for k, ebini, lib in itertools.product(
              self.settings['ks_to_test_naive'],
              range(len(self.energy_bins)),
              self.libs):
                if not (k, lib.barcode, ebini) in k_lib_ebin2total_counts:
                    break
            else:
                return k_lib_ebin2total_counts
        k_lib_ebin2total_counts = {}
        for k, ebini, lib in itertools.product(
          self.settings['ks_to_test_naive'],
          range(len(self.energy_bins)),
          self.libs):
            counts = sum(self.barcode2ebini2k2counts[lib.barcode][ebini][k])
            k_lib_ebin2total_counts[(k, lib.barcode, ebini)] = counts
        cPickle.dump(k_lib_ebin2total_counts, open(ebin_counts_file, 'w'))
        return k_lib_ebin2total_counts

    def make_fold_split_table(self):
        """
        make the table for the structure split
        """
        k_lib_ebin2total_counts = self.get_ebin_counts()
        for k in self.settings['ks_to_test_naive']:
            results_file = os.path.join(self.rdir, 'analyses',
              'ebin_enrichments.%i.txt' % k)
            of = open(results_file, 'w')
            of.write('#barcode\t' + '\t'.join(
              [lib.barcode
              for ebin, lib in itertools.product(self.energy_bins, self.plibs)]
              ) + '\n')
            of.write('#deltaG\t' + '\t'.join(
              ['%0.1f to %0.1f' % ebin
              for ebin, lib in itertools.product(self.energy_bins, self.plibs)]
              ) + '\n')
            of.write('#conc\t' + '\t'.join(
              ['%0.1f' % lib.conc
              for ebin, lib in itertools.product(self.energy_bins, self.plibs)]
              ) + '\n')
            barcode2ebini2kmer2enrichment = {}
            for lib in self.plibs:
                barcode2ebini2kmer2enrichment[lib.barcode]\
                  = collections.defaultdict(dict)
            for kmer in yield_kmers(k):
                of.write(kmer)
                for ebini, lib in itertools.product(
                  range(len(self.energy_bins)), self.plibs):
                    print self.energy_bins
                    ebin = self.energy_bins[ebini]
                    kmer_i = get_index_from_kmer(kmer)
                    lib_counts = self.barcode2ebini2k2counts[
                      lib.barcode][ebini][k][kmer_i]
                    total_lib_counts =\
                      k_lib_ebin2total_counts[(k, lib.barcode, ebini)]
                    input_counts = self.barcode2ebini2k2counts[
                      self.input_lib.barcode][ebini][k][kmer_i]
                    total_input_counts =\
                      k_lib_ebin2total_counts[
                      (k, self.input_lib.barcode, ebini)]
                    enrichment = (1.0 * lib_counts / total_lib_counts) /\
                      (1.0 * input_counts / total_input_counts)
                    barcode2ebini2kmer2enrichment[
                      lib.barcode][ebini][kmer] = enrichment
                    of.write('\t%0.4f' % enrichment)
                of.write('\n')
            of.close()
            results_file = os.path.join(self.rdir, 'analyses',
              'energy_enrichments.%i.pkl' % k)
            cPickle.dump(barcode2ebini2kmer2enrichment,
              open(results_file, 'w'))

    def count_kmer_fold_split(self, barcode, bin_i, energy_bin, results_file):
        """
        does a naive count on the reads split by fold energy
        make this code shared with the actual naive counter
        """
        split_reads_file = os.path.join(self.rdir, 'structure',
          '%s.%i.%0.1f_to_%0.1f.reads' %
          (barcode, bin_i, energy_bin[0], energy_bin[1]))
        assert aUtils.check_file_exists(split_reads_file)
        k2counts = dict()
        kmer2index = {}
        for k in self.settings['ks_to_test_naive']:
            k2counts[k] = numpy.zeros(4 ** k, dtype=int)
            for i, kmer in enumerate(yield_kmers(k)):
                kmer2index[kmer] = i
        read_len = self.settings['read_len']
        for i, line in enumerate(aopen.open(split_reads_file)):
            aUtils.monitor(i)
            if line == '\n':
                continue
            if 'N' in line:
                continue
            try:
                for k in self.settings['ks_to_test_naive']:
                    for ki in range(0, read_len - k + 1):
                        kmer_index = kmer2index[line[ki:ki + k]]
                        k2counts[k][kmer_index] += 1
            except:
                continue
        assert sum(k2counts[6]) > 1
        cPickle.dump(k2counts, open(results_file, 'wb'))

    def count_naive(self):
        """
        Acquires the naive counts by either:
            loading from pkl
            running here
            launching on cluster for later
        """
        count_dir = os.path.join(self.rdir, 'counts', 'naive')
        if self.use_cluster:
            for lib in self.libs:
                command = 'python ~alexrson/snorelax/script.py '\
                  'bind_n_seq_pipeline.naive_counter %s %s 1> %s 2> %s' % \
                  (self.settings_file,
                  lib.barcode,
                  os.path.join(count_dir, '%s.out' % lib.barcode),
                  os.path.join(count_dir, '%s.err' % lib.barcode))
                script_options = {'nodes': '1', 'ppn': '1',
                  'outf': self.settings['experiment_name'] + '.submission_log',
                  'jobname': self.settings['experiment_name']
                  + '_' + lib.barcode + '.n',
                  'queue': 'long', 'workingdir': self.rdir,
                  'command': command}
                if not lib.naive_counts_exist():
                    self.waiting_jobs.append(
                      aUtils.launch(command, script_options))
            self.wait_for_jobs_to_complete()
        map(lambda lib: lib.do_naive_count(), self.libs)

    def count_presence(self):
        count_dir = os.path.join(self.rdir, 'counts', 'presence')
        aUtils.make_dir(count_dir)
        if self.use_cluster:
            print 'using cluster for count presence'
            for k, lib in itertools.product(self.settings['ks_to_test_streaming'], self.libs):
                command = 'python ~alexrson/snorelax/script.py '\
                  'bind_n_seq_pipeline.presence_counter %s %s %i %s 1> %s 2> %s' % \
                  (self.settings_file,
                  lib.barcode, k, self.settings['force_presence_recount'],
                  os.path.join(count_dir, '%s.out' % lib.barcode),
                  os.path.join(count_dir, '%s.err' % lib.barcode))
                if not lib.presence_counts_exist(k) or\
                  self.settings['force_presence_recount']:
                    print 'PRESENCE launching', lib.conc, k
                    self.waiting_jobs.append(
                      aUtils.launch(command, jobname='%s.%g.%i.p' % 
                      (self.settings['experiment_name'], lib.conc, k)))
            print 'wait for jobs to complete'
            self.wait_for_jobs_to_complete()
            for k in self.settings['ks_to_test_streaming']:
                map(lambda lib: lib.do_presence_count(k), self.libs)
        else:
            for k, lib in itertools.product(
              self.settings['ks_to_test_streaming'], self.libs):
                print 'counting presence', lib.conc, k
                lib.do_presence_count(k)
                

    def count_stream(self):

        """
        Acquires the stream counts by either:
            loading from pkl
            running here
            launching on cluster for later
        """
        count_dir = os.path.join(self.rdir, 'counts', 'stream')
        if self.use_cluster:
            print 'using cluster for stream'
            for k, lib in itertools.product(
              self.settings['ks_to_test_streaming'], self.libs):
                command = 'python ~alexrson/snorelax/script.py '\
                  'bind_n_seq_pipeline.stream_counter %s %s %i %s 1> %s 2> %s' % \
                  (self.settings_file,
                  lib.barcode, k, self.settings['force_stream_recount'],
                  os.path.join(count_dir, '%s.out' % lib.barcode),
                  os.path.join(count_dir, '%s.err' % lib.barcode))
                script_options = {'nodes': '1', 'ppn': '1',
                  'outf': self.settings['experiment_name'] + '.submission_log',
                  'jobname': self.settings['experiment_name']
                  + '_' + lib.barcode + '.%i.s' % k,
                  'queue': 'long', 'workingdir': self.rdir,
                  'command': command}
                if not lib.stream_counts_exist(k) or self.settings['force_stream_recount']:
                    print 'launching stream counter', k, lib.conc
                    self.waiting_jobs.append(
                      aUtils.launch(command, script_options))
                    print self.waiting_jobs
            print 'waiting'
            self.wait_for_jobs_to_complete()
        map(lambda lib: lib.do_stream_count(), self.libs)

    def get_all_barcode_handles(self):
        """
        returns a dictionary barcode-> file handles for the split reads
        """
        return dict(
          [(lib.barcode, lib.get_split_whandle()) for lib in self.libs])

    def split_reads(self):
        """
        Splits all the reads
        """
        settings = self.settings
        total_reads = 0
        bad_barcodes = collections.Counter()
        reads_per_barcode = collections.Counter()
        if all([lib.split_reads_exist() for lib in self.libs]):
            [lib.slink_split_reads() for lib in self.libs]
            return
        barcode2of = self.get_all_barcode_handles()
        for l1, l2, l3, l4 in aUtils.iter4Lines(settings['fastq']):
            aUtils.monitor(total_reads, 100000)
            barcode = get_barcode(l1)
            barcode = barcode[0:self.barcode_len]
            barcode_match = self.get_barcode_match(barcode)
            total_reads += 1
            if not barcode_match:
                bad_barcodes[barcode] += 1
                continue
            else:
                trimmed_read = self.trim_read(l2.strip())[0:int(self.settings['read_len'])]
                barcode2of[barcode_match].write(trimmed_read + '\n')
                reads_per_barcode[barcode_match] += 1
            if total_reads > settings['max_reads_to_split'] and\
              settings['max_reads_to_split'] != 0:
                break
        print 'splitting complete'
        self.write_barcode_log(reads_per_barcode,
          total_reads, bad_barcodes)
        map(lambda f: f.close(), barcode2of.values())

    def write_barcode_log(self, reads_per_barcode, total_reads, bad_barcodes):
        """
        Writes the log of which other barcodes are represented.
        """
        barcode_log_of = self.get_rdir_fhandle('split_reads/barcode_log.txt')
        barcode_log_of.write('barcode\tconc\treads assigned\n')
        for barcode, conc in itertools.izip(
          self.settings['barcodes'], self.settings['concentrations']):
            barcode_log_of.write('%s\t%s\t%i\n' %
              (barcode,
              str(conc) if barcode != self.input_barc
              else 'input library',
              reads_per_barcode[barcode]))
        barcode_log_of.write('total_reads\t%i\nbad barcode\t%i\n\n'
          % (total_reads, sum(bad_barcodes.values())))
        barcode_log_of.write('Bad Barcode Summary\n')
        for bad_barcode, count in sorted(
          bad_barcodes.items(), key=operator.itemgetter(1), reverse=True):
            barcode_log_of.write('%s\t%s\n' % (bad_barcode, count))

    def trim_read(self, full_read):
        """
        trims the read from 3' end
        """
        return full_read[0: len(full_read) - self.settings['trim_3p']]

    def process_settings(self):
        """
        reads the settings file and converts str to float or list or whatever
        stores result in self.settings as a dict()
        """
        self.waiting_jobs = []
        self.thourough_check = False
        int_keys = ['replaceative_num_top_kmers', 'read_len',
          'mismatches_allowed_in_barcode', 'trim_3p', 'max_reads_to_split',
          'k_for_concordance_check',
          'num_kmers_for_enrichment_humps']
        float_keys = []
        boolean_keys = ['naive_count', 'replaceative', 'streaming_counts', 'force_stream_recount', 'force_presence_recount']
        list_str_keys = ['barcodes', 'which_barcodes_to_sort_based_on',
          'temperature', 'relevant_variables', 'experiments_to_compare',
          'motifs_of_interest']
        list_int_keys = ['ks_to_test_replaceative', 'ks_to_test_streaming',
          'ks_to_test_naive', 'ks_for_matlab', 'washes']
        list_float_keys = ['concentrations', 'poly_ic_conc',
          'free_energy_limits', 'input_rna']
        extant_files = ['fastq']
        config = ConfigParser.ConfigParser()
        config.read(self.settings_file)
        settings = {}
        for section in config.sections():
            for option in config.options(section):
                settings[option] = config.get(section, option)
                settings[section] = True
        for k in int_keys:
            settings[k] = int(settings[k])
        for k in float_keys:
            settings[k] = float(settings[k])
        for k in boolean_keys:
            if not settings[k].lower() in ['true', 'false']:
                raise ValueError(
                  'Boolean value %s must be "true" or "false"' % k)
            settings[k] = settings[k].lower() == 'true'
        for k in list_float_keys:
            settings[k] = map(float, simplejson.loads(settings[k]))
        for k in list_int_keys:
            settings[k] = map(int, simplejson.loads(settings[k]))
        for k in list_str_keys:
            settings[k] = simplejson.loads(settings[k])
        for k in extant_files:
            assert aUtils.check_file_exists(settings[k])
        self.settings = settings
        self.wdir = settings['working_dir']
        self.rdir = settings['results_dir']
        self.check_barcode_lens()
        self.barcodes = self.settings['barcodes']
        self.barcode2washes = dict(zip(settings['barcodes'],
                          settings['washes']))
        self.barcode2conc = dict(zip(settings['barcodes'],
                          settings['concentrations']))
        self.barcode2polyIC = dict(zip(settings['barcodes'],
                          settings['poly_ic_conc']))
        self.barcode2temp = dict(zip(settings['barcodes'],
                          settings['temperature']))
        self.check_barcodes_are_separated()
        self.k2replaceative_enrich_dict = {}
        self.input_barc = self.settings['library_seq_barcode']
        for k in self.settings['ks_to_test_replaceative']:
            self.k2replaceative_enrich_dict[k] = {}
            for barcode in self.barcodes:
                if barcode == self.input_barc:
                    continue
                self.k2replaceative_enrich_dict[k][barcode] = {}
        self.energy_bins = \
          [(conc1, conc2) for conc1, conc2 in
          zip(self.settings['free_energy_limits'][0:],
          self.settings['free_energy_limits'][1:])]
        self.known_motif = self.settings['known_motif']
        self.known_motif_k = len(self.known_motif)
        self.libs = []
        self.plibs = []
        self.barcode2lib = {}
        if not 1 == len(set(map(len, [self.barcodes,
          settings['concentrations'],
          settings['poly_ic_conc'],
          settings['input_rna'],
          settings['washes'],
          settings['temperature']]))):
            print 'not all library descriptions are the same length'
            print 'barcodes: %i' % len(self.barcodes)
            print 'concentrations: %i' % len(settings['concentrations'])
            print 'polyIC: %i' % len(settings['poly_ic_conc'])
            print 'washes: %i' % len(settings['washes'])
            print 'temps: %i' % len(settings['temperature'])
            print 'input RNA: %i' % len(settings['input_rna'])
            raise ValueError('bad input')
        for barcode, conc, polyIC, washes, temp, input_rna in zip(
          self.barcodes, settings['concentrations'], settings['poly_ic_conc'],
          settings['washes'], settings['temperature'], settings['input_rna']):
            lib = BNS_Lib(self, barcode, conc, polyIC,
              washes, temp, barcode == self.input_barc)
            lib.input_rna = input_rna
            self.libs.append(lib)
            if not barcode == self.input_barc:
                self.plibs.append(lib)
            else:
                self.input_lib = lib
            self.barcode2lib[barcode] = lib
        for input_conc, lib in zip(settings['input_rna'], self.libs):
            lib.input_conc = input_conc
        labels = [lib.get_full_label() for lib in self.libs]
        self.sorted_labels = labels
        if not set(settings['ks_to_test_streaming']).issubset(
          settings['ks_to_test_naive']): 
            print self.settings_file
            raise ValueError('All ks for streaming must also be in naive')

    def create_rdir(self):
        """
        creates the results directory if it doesn't exist
        """
        for dir in ['unsplit', 'counts', 'QC', 'analyses', 'split_reads',
                    'plots', 'tables', 'structure', 'matlab', 'comparisons']:
            if not os.path.exists(os.path.join(self.rdir, dir)):
                os.makedirs(os.path.join(self.rdir, dir))

    def get_counts_file_name(self, barcode, reads_file, method):
        """
        gets the counts file name for a given reads file
        """
        reads_file_base = os.path.basename(reads_file)
        counts_file_base = reads_file_base.replace('.reads', '.pkl')
        out_dir = os.path.join(self.rdir, 'counts', method)
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
        return os.path.join(out_dir, counts_file_base)

    def check_barcode_lens(self):
        """
        verifies that all the barcodes are the same length
        """
        barcode_lens = set(map(len, self.settings['barcodes']))
        if 1 != len(barcode_lens):
            raise ValueError('all barcodes must be the same length')
        self.barcode_len = barcode_lens.pop()

    def check_barcodes_are_separated(self):
        """
        makes sure the barcodes are all totally distinguishable
        """
        for b1, b2 in itertools.combinations(self.settings['barcodes'], 2):
            hamming_dist = aUtils.hamming_distance(b1, b2)
            if hamming_dist < 2:
                raise ValueError('The barcodes supplied are not well '
                  'separated: %s-%s' % (b1, b2))

    def get_input_split_name(self):
        """
        returns the name of the split reads for the input fraction
        """
        return os.path.join(self.rdir, 'split_reads', '%s_%s.reads' %
          (self.settings['experiment_name'], self.input_barc))

    def get_barcode_match(self, barcode):
        """
        takes a barcode and returns the one it matches (hamming <= 1)
        else
        empty string
        """
        if barcode in self.barcodes:
            return barcode
        for barcode_j, conc in zip(self.settings['barcodes'],
                                   self.settings['concentrations']):
            if aUtils.hamming_N(barcode, barcode_j) <= \
              self.settings['mismatches_allowed_in_barcode']:
                return barcode_j
        return ''

    def copy_unsplit_to_wd(self):
        """
        copies and compresses the unsplit fastq to the working directory
        if needed
        """
        fastq = os.path.abspath(self.settings['fastq'])
        working_unsplit_file_name = os.path.join(
          self.rdir, 'unsplit', os.path.basename(fastq))
        if working_unsplit_file_name[-3:] != '.gz':
            working_unsplit_file_name = working_unsplit_file_name + '.gz'
        if working_unsplit_file_name == fastq:
            print 'Unsplit file is defined to be within the working dir'
            pass
        elif aUtils.check_file_exists(working_unsplit_file_name):
            print 'unsplit file in resutls'
        else:
            print 'copying unsplit to working directory'
            start_copy = time.time()
            if fastq[-3:] == '.gz':
                if not os.path.exists(os.path.join(self.rdir, 'unsplit')):
                    os.makedirs(os.path.join(self.rdir, 'unsplit'))
                shutil.copyfile(fastq, working_unsplit_file_name)
            else:
                with aopen.open(fastq, 'rb') as in_f:
                    of = self.get_rdir_fhandle(working_unsplit_file_name)
                    of.writelines(in_f)
                    of.close()
                in_f.close()
            end_copy = time.time()
            print 'Copying unsplit file took %f seconds' %\
              (end_copy - start_copy)
        self.settings['fastq'] = working_unsplit_file_name

    def get_rdir_fhandle(self, *args):
        """
        returns a filehandle to the fname in the rdir
        """
        handle_type = 'r' if args[-1] == 'r' else 'w'
        if args[-1] == 'r':
            args = args[:-1]
        out_path = os.path.join(self.rdir, *args)
        out_dir = os.path.dirname(out_path)
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
        return aopen.open(out_path, handle_type)


# Class pulls just the results from a bind n seq experiment
class Bnse_results(Bnse):
    def __init__(self, settings_file):
        self.settings_file = settings_file
        self.process_settings()


def main():
    """
    run Bnse here
    """
    args = sys.argv[1:]
    if '--use-cluster' in args:
        settings_files = args[1:]
    else:
        settings_files = args
    for settings_file in settings_files:
        if '--use-cluster' in args:
            settings_file = os.path.abspath(settings_file)
            if not aUtils.check_file_exists(settings_file):
                raise ValueError('%s doesnt exist' % settings_file)
            Bnse(settings_file, True, False)
        else:
            Bnse(settings_file)


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


def naive_counter(settings_file, barcode):
    """
    a spawned job to run only the naive count for a given library
    """
    b = Bnse(settings_file, False, False)
    print 'start naive count'
    b.barcode2lib[barcode].do_naive_count()


def stream_counter(settings_file, barcode, k, force_recount):
    """
    a spawned job to run only the naive count for a given library
    """
    b = Bnse(settings_file, False, False)
    print 'start stream count'
    b.barcode2lib[barcode].do_stream_count([int(k)], force_recount=(force_recount.lower() == 'true'))


def presence_counter(settings_file, barcode, k, force_recount):
    """
    a spawned job to run only the presence count for a given library
    """
    b = Bnse(settings_file, False, False)
    print 'start presence count'
    b.barcode2lib[barcode].do_presence_count(int(k), force_recount=(force_recount.lower() == 'true'))


def fold_counter(settings_f, barcode, bin_i, conc_lo, conc_hi, results_file):
    """
    launching free method for the counting of the folded file for a given bin
    """
    print 'starting fold counter'
    b = Bnse(settings_f, False, False)
    b.count_kmer_fold_split(barcode,
      int(bin_i), map(float, (conc_lo, conc_hi)), results_file)
    print 'done fold counter'


def fold_split(settings_file, barcode):
    b = Bnse(settings_file, False, False)
    print 'splitting on cluster by energy'
    b.barcode2lib[barcode].split_by_structure(True)
    print 'done splitting'


def get_barcode(line):
    """
    extracts the barcode from the first line of a fastq quartet
    """
    return line.split('#')[-1].split('/')[0]

def simpleaxis(sp):
    sp.spines['top'].set_visible(False)
    sp.spines['right'].set_visible(False)
    sp.get_xaxis().tick_bottom()
    sp.get_yaxis().tick_left()

def read_replaceative_results(results_file):
    """
    reads the results from the kmers.txt results_file
    """
    kmer2enrich = {}
    t2d = table2dict.table2dict(results_file)
    for kmer, d in t2d.items():
        kmer2enrich[kmer] = float(d['enrichment'])
    print 'loading replaceative results from', results_file
    return kmer2enrich

def chris_formula(R, k, read_len):
    """
     Implements the formula Chris burge derived
    """
    return (4 ** k -1) * ( R * (read_len - k + 1) - (read_len - k) ) / (4 ** k + read_len - k - (R *(read_len - k +1)))


def YGCY(seq):
    test_motif_len = 4
    assert isinstance(seq, str)
    assert seq.upper() == seq
    if len(seq) < test_motif_len:
        return 
    for starti in range(0, len(seq) - 3):
        if seq[starti] in 'CT': # Y
            if seq[starti+1:starti+3] == 'GC':
                if seq[starti+3] in 'CT':
                    return True
    return False

def fracs_passing(k, enrich_arr, bin_edges, test_func=YGCY):
    assert len(enrich_arr) == 4 ** k
    num_kmers = numpy.zeros(len(bin_edges) - 1, dtype=float)
    num_kmers_passing_test = numpy.zeros(len(bin_edges) - 1, dtype=float)
    for kmeri, (kmer, enrich) in enumerate(zip(yield_kmers(k), enrich_arr)):
        bini = aUtils.getBinIndex(enrich, bin_edges)
        if bini == -1:
            print enrich, bin_edges
            print 'bini alert'
            sys.exit()
        num_kmers[bini] += 1
        if test_func(kmer):
            num_kmers_passing_test[bini] += 1
    print num_kmers_passing_test
    print num_kmers

    return num_kmers_passing_test / num_kmers
     
def hasUGU(seq):
    assert isinstance(seq, str)
    assert seq.upper() == seq
    return 'TGT' in seq


class compare_experiments_test(unittest.TestCase):
    def setUp(self):
        self.exp_main = Bnse('settings.fox.json')
        self.exp_comparee = Bnse_results('settings.fox10.json')

    def test_bnse_scatter(self):
        self.exp_main.compare_to_other_BNS('settings.fox10_4.json')

if __name__ == '__main__':
    main()
    #unittest.main()
