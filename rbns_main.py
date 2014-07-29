#!/usr/bin/env python
import sys
import os
import argparse
import cPickle
import shutil
import time
import collections
import itertools
import subprocess
import operator
import math

import numpy as np
import scipy.stats
import matplotlib.pyplot as plt

from rna import rna
import aColors
import stacked_bar_kmers
import streaming_convergence
import rbns_utils
import rbns_settings
import rbns_cluster_utils
import rbns_lib

class Bnse:
    def __init__(self, settings, counts_on_cluster=False):
        self.settings = settings
        self.counts_on_cluster = counts_on_cluster
        self.copy_unsplit_to_wd()
        self.split_reads()
        self.do_counts()
        self.make_sure_counts_worked()
        self.initialize_libs()
        self.calculate_all_enrichments()
        self.determine_most_enriched_lib()

    def calculate_all_enrichments(self):
        """
        sorts all the kmers.
        """
        print 'calculating enrichments'
        rbns_utils.make_dir(self.rdir_path('tables'))
        for k, lib in itertools.product(self.settings.get_naiveks(), self.plibs):
            lib.calculate_enrichment(k, self.input_lib)

    def do_counts(self):
        self.waiting_jobs = []
        print 'doing counts'
        for lib_settings in self.settings.iter_lib_settings():
            for count_type in ['naive', 'stream', 'presence']:
                for k in self.settings.get_ks(count_type):
                    if self.needs_calculation(lib_settings, count_type, k):
                        if self.counts_on_cluster:
                            self.waiting_jobs.append(
                              rbns_cluster_utils.launch_counter(
                                lib_settings, count_type, k,
                                self.settings.get_property('error_dir')))
                        elif not self.counts_on_cluster:
                            self.run_count(lib_settings, count_type, k)
        if self.counts_on_cluster:
            self.wait_for_jobs_to_complete()

    def determine_most_enriched_lib(self):
        """"""
        self.k2most_enriched_lib = {}
        for k in self.settings.get_naiveks():
            most_enriched_lib = self.plibs[0]
            best_enrich = 1.0
            for lib in self.plibs:
                max_enrich = lib.get_max_enrichment(k)
                if max_enrich > best_enrich:
                    most_enriched_lib = lib
                    best_enrich = max_enrich
            self.k2most_enriched_lib[k] = most_enriched_lib

    def make_sure_counts_worked(self):
        max_reattempts = 2
        for pass_i in range(max_reattempts):
            if self.verify_counts_ran_successfully():
                return
            else:
                sys.stderr.write('Counts do not seem to all have concluded successfully')
                sys.stderr.write('Will reattempt counts')
                sys.stderr.write('This is retry number: %i' % (pass_i + 1))
                self.do_counts()
        if self.counts_on_cluster:
            self.counts_on_cluster = False
            sys.stderr.write('Counting on the cluster seems to now work well.')
            sys.stderr.write('Switching to local counting')
            self.do_counts()
        raise RuntimeError('Counts keep failing :(\nCheck the errors in: %s' %
          self.settings.get_property('error_dir'))

    def verify_counts_ran_successfully(self):
        for lib_settings in self.settings.iter_lib_settings():
            for count_type in ['naive', 'stream', 'presence']:
                for k in self.settings.get_ks(count_type):
                    if self.needs_calculation(lib_settings, count_type, k):
                        return False
        return True

    def wait_for_jobs_to_complete(self, sleep=10):
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
                jobs_left = sum(map(lambda output:
                  0 if 'Unknown J' in output[1] else 1, outputs))
                print jobs_left, 'jobs left', [
                  job for job, job_complete in zip(self.waiting_jobs, completed) if not job_complete]
            time.sleep(sleep)
        self.waiting_jobs = []

    def run_count(self, lib_settings, count_type, k):
        split_reads = lib_settings.get_split_reads()
        out_pkl = lib_settings.counts_file(count_type, k)
        rbns_utils.make_dir(os.path.dirname(out_pkl))
        print 'doing counts', count_type, k, split_reads
        if count_type == 'naive':
            count_naive(split_reads, k, out_pkl)
        elif count_type == 'stream':
            count_stream(split_reads, k, out_pkl)
        elif count_type == 'presence':
            count_presence(split_reads, k, out_pkl)
        else:
            raise ValueError('Unknown count type: %s ' % count_type)

    def initialize_libs(self):
        self.libs = []
        self.plibs = []
        for lib_settings in self.settings.iter_lib_settings():
            lib = rbns_lib.RBNS_Lib(self.settings, lib_settings)
            self.libs.append(lib)
            if lib.is_input():
                self.input_lib = lib
            else:
                self.plibs.append(lib)

    def needs_calculation(self, lib_settings, count_type, k):
        if self.settings.get_force_recount(count_type):
            return True
        return not lib_settings.counts_exist(count_type, k)

    def calculate_kds(self):
        kds_ks = self.settings.get_ks('stream')
        if kds_ks == []:
            raise ValueError('Must have stream counts to do Kd calculations')
        for k in kds_ks:
            most_enriched_lib = self.k2most_enriched_lib[k]
            relative_kds = most_enriched_lib.calculate_kds(k)
            out_file = self.rdir_path('tables', 'relativekds.%i.txt' % k)
            of = open(out_file, 'w')
            of.write('#kmer\tkd\n')
            for kmer, relative_kd in zip(rbns_utils.yield_kmers(k), relative_kds):
                of.write('%s\t%g\n' % (kmer, relative_kd))
            of.close()

    def make_tables(self):
        self.verify_concordances()
        self.spearman_concordance_table_topquartile()
        self.spearman_concordance_table_top_decile()
        #self.create_matlab_inputs()
        self.make_enrichment_table()
        if self.settings.get_ks('stream'):
            self.make_SKA_table()
            self.make_motif_libfrac_table()
        self.make_counts_tables()
        self.make_Btable()

    def make_plots(self):
        rbns_utils.make_dir(self.rdir_path('plots'))
        self.make_enrichment_hist_stacked()
        self.plot_k_picker()
        self.plot_concordance_to_most_enriched()
        self.plot_enrichment_humps()
        self.plot_enrichment_humps_with_secondary()
        self.plot_relativeAffinity_humps_with_secondary()
        self.plot_enrichment_humps_with_secondary_stream()
        self.plot_enrichment_humps_with_secondary_presence()
        self.plot_enrichment_v_stream()
        self.plot_enrichment_v_presence()
        self.plot_presence_v_stream()
        if '--all-tasks' in sys.argv or '--structure-calculations' in sys.argv:
            self.plot_fold_split_hump()

    def make_Btable(self):
        """
        Makes a table of the B values and R values
        """
        for k in self.settings.get_property('ks_to_test_naive'):
            t_of = self.get_rdir_fhandle('tables',
              'Btable.%i.xls' % k)
            self.make_table_header(t_of)
            for kmeri, kmer in enumerate(rbns_utils.yield_kmers(k)):
                t_of.write(kmer + '\t')
                t_of.write('\t'.join(
                  [str(lib.calcB(kmer)) for lib in self.plibs]) + '\n')
            t_of.close()

    def compare_all_other_experiments(self):
        """
        does all the comparisons of the experiment to the the other
        ones specified in the settings
        """
        for alt_settings_file in self.settings.get_property('experiments_to_compare'):
            if not os.path.isabs(alt_settings_file):
                alt_settings_file =\
                  os.path.join(os.getcwd(), alt_settings_file)
            assert rbns_utils.file_exists(alt_settings_file)
            alt_settings = rbns_settings.RBNS_settings(alt_settings_file)
            self.compare_to_other_BNS(alt_settings)

    def compare_to_other_BNS(self, other_experiment_settings):
        """
        does the comparison to a single other experiment
        """
        print 'making a new BNS results instance %s' %\
          other_experiment_settings
        other_exp = Bnse(other_experiment_settings)
        print 'Comparing to %s' % other_exp.settings.get_property('experiment_name')
        rbns_utils.make_dir(self.rdir_path('comparisons'))
        if True:
            # scatter of protein concentrations
            fig1 = plt.figure()
            sp = fig1.add_subplot(1, 1, 1)
            for lib in self.plibs:
                sp.plot(1.0, lib.get_conc(), 'o', color=aColors.steel_blue)
            for lib in other_exp.plibs:
                sp.plot(2.0, lib.get_conc(), 'o', color=aColors.bacon)
            sp.set_yscale('log')
            sp.set_xlim([0.5, 2.5])
            sp.set_xticks([1, 2])
            sp.set_xticklabels([self.settings.get_property('experiment_name'),
              other_exp.settings.get_property('experiment_name')])
            fig1.savefig(self.rdir_path(
              'comparisons',
              'protein_conc.%s.png'
              % other_exp.settings.get_property('experiment_name')))
            # Scatter: enrich
            print 'scatter enrich'
            k = 6
            for lib_main, lib_other in itertools.product(self.plibs,
              other_exp.plibs):
                protein_conc_closeness_thresh = 10
                if not rbns_utils.close_float_value(
                  lib_main.get_conc(),
                  lib_other.get_conc(),
                  protein_conc_closeness_thresh):
                    continue
                print lib_main.get_conc(), lib_other.get_conc()
                fig2 = plt.figure()
                sp = fig2.add_subplot(1, 1, 1)
                main_enriches = lib_main.get_enrichments(k)
                other_enriches = lib_other.get_enrichments(k)
                thresh = np.mean(main_enriches) + 2 * np.std(main_enriches)
                for x, y, kmer in zip(main_enriches, other_enriches, rbns_utils.yield_kmers(k)):
                    sp.loglog(x,y,'.',
                      color=aColors.protein_colors(kmer,
                        self.settings.get_property('protein_name'), x > thresh))
                fig2.suptitle('Scatter plot of enrichments %imers for '
                  '[protein]=%g or %g' %
                  (k, lib_main.get_conc(), lib_other.get_conc()))
                sp.set_xlabel('R values from ' + self.settings.get_property('experiment_name'))
                sp.set_ylabel('R values from ' + lib_other.get_property('experiment_name'))
                sp.set_aspect(1)
                rbns_utils.simpleaxis(sp)
                fig2.savefig(self.rdir_path(
                  'comparisons',
                  'enrichment_correlation.%s.%.1f.nM.png'
                  % (other_exp.settings.get_property('experiment_name'),
                  lib_main.get_conc())))

        # Scatter: SKA
        print 'scatter SKA'
        if True:
            k=6
            print 'comparing SKA'
            for lib_main, lib_other in itertools.product(self.plibs,other_exp.plibs):
                fig3 = plt.figure()
                sp = fig3.add_subplot(1, 1, 1)
                if not rbns_utils.close_float_value(
                  lib_main.get_conc(),
                  lib_other.get_conc(),
                  protein_conc_closeness_thresh):
                    continue
                fig3 = plt.figure()
                sp = fig3.add_subplot(1, 1, 1)
                print 'comparing ', lib_main.get_conc(), lib_other.get_conc()
                main_kmer_libfracs = lib_main.get_stream_libfracs(k)
                other_kmer_libfracs = lib_other.get_stream_libfracs(k)
                for x, y, kmer in zip(
                  main_kmer_libfracs, other_kmer_libfracs, rbns_utils.yield_kmers(k)):
                    sp.plot(x,y,'o',
                            color=aColors.protein_colors(kmer,
                            self.settings.get_property('protein_name'),
                            x > 5./(4**k)))
                sp.set_xlabel(r'SKA $F_{i}$ for '
                  + self.settings.get_property('experiment_name')
                  + ' [RBP]=%g' % lib_main.get_conc())
                sp.set_ylabel(r'SKA $F_{i}$ for '
                  + lib_other.get_property('experiment_name')
                  + ' [RBP]=%g' % lib_other.get_conc())
                sp.set_aspect(1.)
                rbns_utils.simpleaxis(sp)
                out_file = self.rdir_path('comparisons',
                  'SKA_comp_%s.%1.f.%.1f.nM.png' %
                  (other_exp.settings.get_property('experiment_name'),
                  lib_main.get_conc(),
                  lib_other.get_conc()))
                print 'saving ', out_file
                fig3.savefig(out_file)
        # Scatter: B(enrich)
        print 'scatter B'
        k = 6
        for lib_main, lib_other in itertools.product(self.plibs,
          other_exp.plibs):
            if not rbns_utils.close_float_value(lib_main.get_conc(), lib_other.get_conc(), 40):
               continue
            fig2 = plt.figure()
            sp = fig2.add_subplot(1, 1, 1)
            main_enriches = lib_main.get_enrichments(k)
            other_enriches = lib_other.get_enrichments(k)
            thresh = np.mean(main_enriches) + 2 * np.std(main_enriches)
            main_read_len = self.settings.get_property('read_len')
            other_read_len = other_exp.settings.get_property('read_len')
            for x, y, kmer in zip(main_enriches, other_enriches, rbns_utils.yield_kmers(k)):
                if rbns_utils.chris_formula(x, k, main_read_len)< 0 or\
                  rbns_utils.chris_formula(y, k, other_read_len) < 0:
                    continue
                sp.loglog(rbns_utils.chris_formula(x, k, main_read_len),
                          rbns_utils.chris_formula(y, k, other_read_len),
                          '.',
                          color=aColors.protein_colors(
                            kmer,
                            self.settings.get_property('protein_name'),
                            x > thresh))
            sig_enriched = rbns_utils.significantly_enriched(
              main_enriches, zthresh=2., scale='log')
            if any(sig_enriched):
                xs, ys = zip(*[(x, y) for x, y, e in
                  zip(main_enriches, other_enriches, sig_enriched) if e])
                r, p = scipy.stats.pearsonr(xs, ys)
            fig2.suptitle('Scatter plot of B values %imers for '
              '[protein]=%g or %g' %
              (k, lib_main.get_conc(), lib_other.get_conc()))
            sp.set_xlabel('B values for ' + self.settings.get_property('experiment_name'))
            sp.set_ylabel('B values for ' + lib_other.get_property('experiment_name'))
            sp.set_aspect(1)
            sp.set_xlim([1, 1000])
            sp.set_ylim([1, 1000])
            rbns_utils.simpleaxis(sp)
            fig2.savefig(self.rdir_path(
              'comparisons',
              'Bvalues.%s.%.1f.nM.png'
              % (other_exp.settings.get_property('experiment_name'),
              lib_main.get_conc())))

    def create_matlab_inputs(self):
        """
        makes the file that matlab needs to calculate Kds
        """
        print 'Creating MATLAB inputs'
        for k in self.settings.get_property('ks_for_matlab'):
            assert isinstance(k, int)
            kmer_list_file = self.rdir_path(
              'matlab', 'kmer_list.%i.m' % k)
            kmer_count_file = self.rdir_path(
              'matlab', 'kmer_counts.%i.m' % k)
            rbns_utils.make_dir(os.path.dirname(kmer_list_file))
            kmer_list_f = open(kmer_list_file, 'w')
            kmer_count_f = open(kmer_count_file, 'w')
            for kmer_i, kmer in enumerate(rbns_utils.yield_kmers(k)):
                kmer_list_f.write(kmer + '\n')
                kmer_count_f.write('\t'.join(
                  [str(int(lib.get_naive_counts(k)[kmer_i]))
                  for lib in [self.input_lib] + self.plibs]) + '\n')
            kmer_list_f.close()
            kmer_count_f.close()
        protein_conc_file = self.rdir_path(
          'matlab', 'protein_conc.m')
        protein_conc_f = open(protein_conc_file, 'w')
        protein_conc_f.write('\t'.join([str(lib.get_conc()) for lib in [self.input_lib] + self.plibs]) + '\n')
        protein_conc_f.close()

    def plot_k_picker(self):
        """
        for each protein library
        plots the kpicker plot (such as it is)
        """
        for lib in self.plibs:
            lib.plot_next_kmer_frac()
        # determine most enriched library
        k_for_max_enrichment = np.median(map(int, self.settings.get_naiveks()))
        if isinstance(k_for_max_enrichment, list):
            k_for_max_enrichment = k_for_max_enrichment[0]
        most_enriched_lib = self.k2most_enriched_lib[int(k_for_max_enrichment)]
        fig1 = plt.figure()
        sp = fig1.add_subplot(111)
        top_vs_next = []
        for k in self.settings.get_property('ks_to_test_naive'):
            top_vs_next.append(most_enriched_lib.get_top_seq_vs_next(k, 2))
        sp.bar([k - 0.4 for k in self.settings.get_property('ks_to_test_naive')],
          top_vs_next)
        fig1.savefig(self.rdir_path('plots', 'best_v_next.png'))

    def plot_concordance_to_most_enriched(self):
        """
        plots how well this library agrees with the most enriched library
        """
        rbns_utils.make_dir(self.rdir_path('analyses'))
        num_top_kmers_to_compare = int(self.settings.get_property(
          'num_top_kmers_to_compare', 50))
        k = self.settings.get_property('k_for_concordance_to_most_enriched', 7)
        most_enriched_lib = self.k2most_enriched_lib[k]
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
        fig1.savefig(self.rdir_path('plots',
          'comparison_to_top_lane.png'))

    def make_table_header(self, of):
        """
        takes a file handle and writes a good header for it such that
        each lane is a column.
        """
        of.write('#')
        for lib in self.libs:
            of.write('\t' + lib.get_barcode())
        of.write('\n[%s]' % self.settings.get_property('protein_name'))
        for lib in self.libs:
            of.write('\t%0.3f' % lib.get_conc())
        of.write('\nwashes')
        for lib in self.libs:
            of.write('\t%i' % lib.get_washes())
        of.write('\n[poly IC]')
        for lib in self.libs:
            of.write('\t%i' % lib.get_poly_ic())
        of.write('\nT (C)')
        for lib in self.libs:
            of.write('\t%s' % lib.get_temperature())
        of.write('\n')

    def make_counts_tables(self):
        """
        For each k
        Makes a table of the naive counts of each kmer.
        """
        print 'making counts tables'
        for k in self.settings.get_property('ks_to_test_naive'):
            print 'writing counts table for %i' % k
            of = self.get_rdir_fhandle('tables', 'naive_counts.%i.txt' % k)
            self.make_table_header(of)
            for kmeri in range(4 ** k):
                kmer_str = rbns_utils.get_kmer_from_index(k, kmeri)
                of.write('\n' + kmer_str)
                for lib in self.libs:
                    of.write('\t%i' %
                     lib.get_naive_count_by_index(k, kmeri))
            of.close()
            # sorted iteration through kmers
            of_sorted = self.get_rdir_fhandle(
              'tables/naive_counts_sorted.%i.txt' % k)
            self.make_table_header(of_sorted)
            for kmeri in self.naively_sorted_kmers[k]:
                kmer_str = rbns_utils.get_kmer_from_index(k, kmeri)
                of_sorted.write('\n' + kmer_str)
                for lib in self.libs:
                    of_sorted.write('\t%i' %
                      lib.get_naive_count_by_index(k, kmeri))
            of_sorted.close()

    def make_motif_libfrac_table(self):
        """
        Makes a table of the library fraction for each kmer
        uses the length of the known motif for k.

        uses the stream algorithm
        """
        motif = self.settings.get_property('known_motif')
        if not len(motif) in self.settings.get_naiveks():
            return
        print 'writing libfrac table'
        of = self.get_rdir_fhandle('tables', 'motif_libfrac.txt')
        of.write('barcode\tconcentration\t[poly_IC]\t'
          'washes\tlibfrac of known seq: %s\n' %
          self.settings.get_property('known_motif'))
        for lib in self.libs:
            libfrac_motif = lib.get_stream_libfrac_kmer(motif)
            of.write('%s\t%0.3f\t%0.3f\t%i\t%f\n' %
              (lib.get_barcode(),
              lib.get_conc(),
              lib.get_poly_ic(),
              lib.get_washes(),
              libfrac_motif))
        of.close()


    def plot_enrichment_humps(self):
        """
        plots the enrichment of the top kmers into a pdf
        """
        for k in self.settings.get_property('ks_to_test_naive'):
            fig1 = plt.figure()
            sp = fig1.add_subplot(111)
            for kmer_i in self.naively_sorted_kmers[k][0:
              self.settings.get_property('num_kmers_for_enrichment_humps')]:
                enrichments = [lib.get_enrichments(k)[kmer_i]
                                    for lib in self.plibs]
                sp.plot(range(len(self.plibs)),
                  enrichments, 'o-', color=aColors.dodo_grey)
            if k == len(self.settings.get_property('known_motif')):
                enrichments = [lib.get_enrichments(k)[
                  rbns_utils.get_index_from_kmer(self.settings.get_property('known_motif'))]
                  for lib in self.plibs]
                sp.plot(range(len(self.plibs)), enrichments, 'o-', color=aColors.bacon)
            else:
                enrichments = []
            sp.set_xticks(range(len(self.plibs)))
            labels = [lib.get_full_label() for lib in self.libs]
            sp.set_xticklabels(labels, rotation=45)
            sp.set_ylabel('R value')
            fig1.suptitle(self.settings.get_property('experiment_name') + (' k=%i' % k))
            rbns_utils.save_fig(fig1,self.rdir_path('plots', 'enrichment_humps.%i.png' % k))

    def plot_enrichment_v_stream(self):
        ks_in_both = set(self.settings.get_property('ks_to_test_naive')) &\
          set(self.settings.get_property('ks_to_test_stream'))
        for k, lib in itertools.product(ks_in_both, self.plibs):
            fig1 = plt.figure(figsize=(5, 7))
            sp = fig1.add_subplot(211)
            out_plot = self.rdir_path(
              'plots', 'count_compare.naive_stream.%i.%g' % (k, lib.get_conc()))
            xs = lib.get_enrichments(k)
            ys = lib.get_stream_libfracs(k)
            assert len(xs) == len(ys)
            sp.loglog(xs, ys, '.')
            sp.set_xlabel('naive enrichment')
            sp.set_ylabel('stream weight')
            sp.set_title('protein: %s, conc: %g' %
              (self.settings.get_property('protein_name'), lib.get_conc()))
            print out_plot
            rbns_utils.save_fig(fig1, out_plot)

    def plot_enrichment_v_presence(self):
        ks_in_both = set(self.settings.get_property('ks_to_test_naive')) &\
          set(self.settings.get_property('ks_to_test_stream'))
        for k, lib in itertools.product(ks_in_both, self.plibs):
            fig1 = plt.figure(figsize=(5, 7))
            sp = fig1.add_subplot(211)
            out_plot = self.rdir_path(
              'plots', 'count_compare.naive_presence.%i.%g' % (k, lib.get_conc()))
            xs = lib.get_enrichments(k)
            ys = lib.get_presence_fracs(k)
            assert len(xs) == len(ys)
            try:
                sp.loglog(xs, ys, '.')
            except:
                print lib.get_barcode(), k
                sys.exit()
            sp.set_xlabel('naive enrichment')
            sp.set_ylabel('Presence')
            sp.set_title('protein: %s, conc: %g' %
              (self.settings.get_property('protein_name'), lib.get_conc()))
            print out_plot
            rbns_utils.save_fig(fig1, out_plot)

    def plot_presence_v_stream(self):
        for k, lib in itertools.product(self.settings.get_property('ks_to_test_stream'), self.plibs):
            fig1 = plt.figure(figsize=(5, 7))
            sp = fig1.add_subplot(211)
            out_plot = self.rdir_path(
              'plots', 'count_compare.presence_stream.%i.%g' % (k, lib.get_conc()))
            xs = lib.get_presence_fracs(k)
            ys = lib.get_stream_libfracs(k)
            assert len(xs) == len(ys)
            sp.loglog(xs, ys, '.')
            sp.set_xlabel('presence')
            sp.set_ylabel('Stream weight')
            sp.set_title('protein: %s, conc: %g' %
              (self.settings.get_property('protein_name'), lib.get_conc()))
            print out_plot
            rbns_utils.save_fig(fig1, out_plot)

    def plot_enrichment_humps_with_secondary(self):
        """
        plots the enrichment of the top kmers into a pdf
        """
        xvals = [max(0, math.log(lib.get_conc() + .00001)) for lib in self.plibs]
        for k in self.settings.get_property('ks_to_test_naive'):
            fig1 = plt.figure(figsize=(7, 7))
            sp = fig1.add_subplot(211)
            kmers_of_interest = [kmer
              for kmer in self.settings.get_property('motifs_of_interest') if len(kmer) == k]
            legend_handles = []
            rand = {6:500, 7:3015, 5:129, 4:31, 8:8000, 9:32000, 10: 128000, 3:7}
            for kmer_i in range(33, 4 ** k, rand[k]):
                enrichments = [ lib.get_enrichment(k, kmer_i)
                                    for lib in self.plibs]
                sp.plot(xvals,
                  enrichments, 'o-', color=aColors.background_motif)
            for ranki, kmer in enumerate(kmers_of_interest):
                enrichments =\
                  [lib.get_enrichment_kmer(kmer)
                  for lib in self.plibs]
                if 'FOX' in self.settings.get_property('protein_name').upper():
                    color = aColors.fox_colors(kmer, True)
                elif 'MBNL' in self.settings.get_property('protein_name').upper():
                    color = aColors.mbnl_colors(kmer, True)
                    color = aColors.colors8[ranki]
                elif 'CUGBP' in self.settings.get_property('protein_name').upper():
                    color= aColors.newcolors6[ranki]
                    color = aColors.colors8[ranki]
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
            rbns_utils.simpleaxis(sp)
            if 'CUGBP' in self.settings.get_property('protein_name').upper() or\
              'MBNL' in self.settings.get_property('protein_name').upper():
                if k == 6 or k == 7:
                    sp.set_ylim([0, 10])
                    sp.set_yticks([0, 1, 2, 4, 6, 8, 10])
            fig1.suptitle(self.settings.get_property('experiment_name') + (' k=%i' % k),
              verticalalignment='baseline')
            out_plot = self.rdir_path(
              'plots', 'enrichment_humps.kmersofinterest.%i.nolegend' % k)
            rbns_utils.save_fig(fig1, out_plot)
            sp.legend(legend_handles,
                      map(rna, kmers_of_interest),
                      mode='expand',
                      loc=2,
                      borderaxespad=-2.)
            out_plot = self.rdir_path(
              'plots', 'enrichment_humps.kmersofinterest.%i.withlegend' % k)
            print 'saving', out_plot
            rbns_utils.save_fig(fig1, out_plot)

    def plot_relativeAffinity_humps_with_secondary(self):
        """
        plots the enrichment of the top kmers into a pdf
        """
        print 'plotting relative affinity humps'
        xvals = [max(0, math.log(lib.get_conc() + .00001)) for lib in self.plibs]
        for k in self.settings.get_property('ks_to_test_naive'):
            fig1 = plt.figure(figsize=(5, 7))
            sp = fig1.add_subplot(211)
            kmers_of_interest = [kmer
              for kmer in self.settings.get_property('motifs_of_interest') if len(kmer) == k]
            legend_handles = []
            for ranki, kmer in enumerate(kmers_of_interest):
                if 'FOX' in self.settings.get_property('protein_name').upper():
                    color = aColors.fox_colors(kmer, True)
                elif 'MBNL' in self.settings.get_property('protein_name').upper():
                    color = aColors.mbnl_colors(kmer, True)
                    color = aColors.colors8[ranki]
                elif 'CUGBP' in self.settings.get_property('protein_name').upper():
                    color= aColors.newcolors6[ranki]
                    color = aColors.colors8[ranki]
                else:
                    color = aColors.ERSP_colors(kmer, True)
                read_len = self.settings.get_property('read_len')
                B_values = [lib.get_B_kmer(kmer, read_len) for lib in self.plibs]
                try:
                    l, = sp.semilogy(xvals,
                            B_values,
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
            rbns_utils.simpleaxis(sp)
            fig1.suptitle(self.settings.get_property('experiment_name') + (' k=%i' % k),
              verticalalignment='baseline')
            out_plot = self.rdir_path(
              'plots', 'relative_affinity_humps.kmersofinterest.%i.nolegend' % k)
            rbns_utils.save_fig(fig1, out_plot)
            sp.legend(legend_handles,
                      map(rna, kmers_of_interest),
                      mode='expand',
                      loc=2,
                      borderaxespad=-2.)
            out_plot = self.rdir_path(
              'plots', 'relative_affinity_humps.kmersofinterest.%i.withlegend' % k)
            print 'saving', out_plot
            rbns_utils.save_fig(fig1, out_plot)

    def plot_enrichment_humps_with_secondary_stream(self):
        """
        plots the enrichment of the top kmers into a pdf
        """
        xvals = [max(0, math.log(lib.get_conc() + .00001)) for lib in self.plibs]
        for k in self.settings.get_property('ks_to_test_stream'):
            fig1 = plt.figure(figsize=(5, 7))
            sp = fig1.add_subplot(211)
            kmers_of_interest = [kmer
              for kmer in self.settings.get_property('motifs_of_interest') if len(kmer) == k]
            legend_handles = []
            rand = {6:500, 7:3015, 5:129, 4:31, 8:8000}
            for kmer_i in range(33, 4 ** k, rand[k]):
                enrichments = [lib.get_stream_libfrac(k , kmer_i)
                                    for lib in self.plibs]
                sp.plot(xvals,
                  enrichments, 'o-', color=aColors.background_motif)
            for ranki, kmer in enumerate(kmers_of_interest):
                enrichments = [lib.get_stream_libfrac(k , kmer_i)
                                    for lib in self.plibs]
                if 'FOX' in self.settings.get_property('protein_name').upper():
                    color = aColors.fox_colors(kmer, True)
                elif 'MBNL' in self.settings.get_property('protein_name').upper():
                    color = aColors.mbnl_colors(kmer, True)
                    color = aColors.colors8[ranki]
                elif 'CUGBP' in self.settings.get_property('protein_name').upper():
                    color= aColors.newcolors6[ranki]
                    color = aColors.colors8[ranki]
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
            rbns_utils.simpleaxis(sp)
            if 'CUGBP' in self.settings.get_property('protein_name').upper() or\
              'MBNL' in self.settings.get_property('protein_name').upper():
                if k == 6 or k == 7:
                    pass
            fig1.suptitle(self.settings.get_property('experiment_name') + (' k=%i' % k),
              verticalalignment='baseline')
            out_plot = self.rdir_path(
              'plots', 'stream_humps.kmersofinterest.%i.nolegend' % k)
            rbns_utils.save_fig(fig1, out_plot)
            sp.legend(legend_handles,
                      map(rna, kmers_of_interest),
                      mode='expand',
                      loc=2,
                      borderaxespad=-2.)
            out_plot = self.rdir_path(
              'plots', 'stream_humps.kmersofinterest.%i.withlegend' % k)
            print 'saving', out_plot
            rbns_utils.save_fig(fig1, out_plot)

    def plot_enrichment_humps_with_secondary_presence(self):
        """
        plots the enrichment of the top kmers into a pdf
        """
        xvals = [max(0, math.log(lib.get_conc() + .00001)) for lib in self.plibs]
        for k in self.settings.get_property('ks_to_test_stream'):
            fig1 = plt.figure(figsize=(5, 7))
            sp = fig1.add_subplot(211)
            kmers_of_interest = [kmer
              for kmer in self.settings.get_property('motifs_of_interest') if len(kmer) == k]
            legend_handles = []
            rand = {6:500, 7:3015, 5:129, 4:31, 8:8000}
            for kmer_i in range(33, 4 ** k, rand[k]):
                enrichments = [lib.get_presence_frac(k, kmer_i)
                                    for lib in self.plibs]
                sp.plot(xvals,
                  enrichments, 'o-', color=aColors.background_motif)
            for ranki, kmer in enumerate(kmers_of_interest):
                presences = [lib.get_presence_frac_kmer(kmer) for lib in self.plibs]
                if 'FOX' in self.settings.get_property('protein_name').upper():
                    color = aColors.fox_colors(kmer, True)
                elif 'MBNL' in self.settings.get_property('protein_name').upper():
                    color = aColors.mbnl_colors(kmer, True)
                    color = aColors.colors8[ranki]
                elif 'CUGBP' in self.settings.get_property('protein_name').upper():
                    color= aColors.newcolors6[ranki]
                    color = aColors.colors8[ranki]
                else:
                    color = aColors.ERSP_colors(kmer, True)
                l, = sp.plot(xvals,
                        presences,
                        'o-',
                        color=color,
                        label=kmer)
                legend_handles.append(l)
            labels = [lib.get_full_label() for lib in self.libs]
            self.sorted_labels = labels
            sp.set_xticks(xvals)
            sp.set_xticklabels(labels, rotation=45)
            sp.set_ylabel('fraction with kmer')
            rbns_utils.simpleaxis(sp)
            if 'CUGBP' in self.settings.get_property('protein_name').upper() or\
              'MBNL' in self.settings.get_property('protein_name').upper():
                if k == 6 or k == 7:
                    pass
            fig1.suptitle(self.settings.get_property('experiment_name') + (' k=%i' % k),
              verticalalignment='baseline')
            out_plot = self.rdir_path(
              'plots', 'presence_humps.kmersofinterest.%i.nolegend' % k)
            rbns_utils.save_fig(fig1, out_plot)
            sp.legend(legend_handles,
                      map(rna, kmers_of_interest),
                      mode='expand',
                      loc=2,
                      borderaxespad=-2.)
            out_plot = self.rdir_path(
              'plots', 'presence_humps.kmersofinterest.%i.withlegend' % k)
            print 'saving', out_plot
            rbns_utils.save_fig(fig1, out_plot)

    def spearman_concordance_table_top_decile(self):
        """
        """
        k = self.settings.get_property('k_for_concordance_check')
        most_enriched_lib = self.k2most_enriched_lib[k]
        num_kmers_in_top_quartile = int((4 ** k) / 10)
        upper_quartile_kmers = most_enriched_lib.get_top_kmers(
          k, num_kmers_in_top_quartile)
        new_concordance_of = self.get_rdir_fhandle(
          'QC', '%s.top_decile_enrich.txt' % self.settings.get_property('experiment_name'))
        new_concordance_of.write('barcode 1\tbarcode 2\tconcentration1'\
           '\tconcentration 2\tSpearman Rho\tSpearman P\n')
        #for libi, libj in rbns_utils.pairwise(self.plibs):
        for libi, libj in itertools.combinations(self.plibs, 2):
            libi_presences = [libi.get_presence_frac(k, kmer_i)
                                for kmer_i in upper_quartile_kmers]
            libj_presences = [libj.get_presence_frac(k, kmer_i)
                                for kmer_i in upper_quartile_kmers]
            spearmanr, spearmanp = scipy.stats.spearmanr(
              libi_presences, libj_presences)
            new_concordance_of.write('%s\t%s\t%g\t%g\t%g\t%g\n' %
              (libi.get_barcode(), libj.get_barcode(),
              libi.get_conc(), libj.get_conc(),
              spearmanr, spearmanp))
        new_concordance_of.close()

    def spearman_concordance_table_topquartile(self):
        """
        """
        k = self.settings.get_property('k_for_concordance_check')
        most_enriched_lib = self.k2most_enriched_lib[k]
        num_kmers_in_top_quartile = int((4 ** k) / 4)
        upper_quartile_kmers = most_enriched_lib.get_top_kmers(
          k, num_kmers_in_top_quartile)
        new_concordance_of = self.get_rdir_fhandle(
          'QC', '%s.top_quartile_enrich.txt' % self.settings.get_property('experiment_name'))
        new_concordance_of.write('barcode 1\tbarcode 2\tconcentration1'\
           '\tconcentration 2\nSpearman Rho\tSpearman P\n')
        #for libi, libj in rbns_utils.pairwise(self.plibs):
        for libi, libj in itertools.combinations(self.plibs, 2):
            libi_presences = [libi.get_presence_frac(k, kmer_i)
                                for kmer_i in upper_quartile_kmers]
            libj_presences = [libj.get_presence_frac(k, kmer_i)
                                for kmer_i in upper_quartile_kmers]
            spearmanr, spearmanp = scipy.stats.spearmanr(
              libi_presences, libj_presences)
            new_concordance_of.write('%s\t%s\t%g\t%g\t%g\t%g\n' %
              (libi.get_barcode(), libj.get_barcode(),
              libi.get_conc(), libj.get_conc(),
              spearmanr, spearmanp))
        new_concordance_of.close()

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
        for libi, libj in rbns_utils.pairwise(self.libs):
            pearsonr, pearsonp = self.get_libfrac_kmer_concordance(
              libi, libj, 'pearson')
            spearmanr, spearmanp = self.get_libfrac_kmer_concordance(
              libi, libj, 'spearman')
            lf_concordance_of.write(
              '%s\t%s\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\n' %
              (libi.get_barcode(), libj.get_barcode(), libi.get_conc(), libj.get_conc(),
              pearsonr, pearsonp, spearmanr, spearmanp))
            pearsonr, pearsonp = self.get_enrich_kmer_concordance(
              libi, libj, 'pearson')
            spearmanr, spearmanp = self.get_enrich_kmer_concordance(
              libi, libj, 'spearman')
            enrich_concordance_of.write(
              '%s\t%s\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\n' %
              (libi.get_barcode(), libj.get_barcode(), libi.get_conc(), libj.get_conc(),
              pearsonr, pearsonp, spearmanr, spearmanp))
            if spearmanr < 0.5:
                print 'low concordance warning', libi.get_barcode(), libj.get_barcode()
        lf_concordance_of.close()
        enrich_concordance_of.close()

    def get_enrich_kmer_concordance(self, libi, libj, test):
        """
        checks that two libaries are similar with regard to enrichment
        """
        k = self.settings.get_property('k_for_concordance_check')
        if libi.is_input or libj.is_input:
            return (1.0, 1.0)
        if test == 'pearson':
            print  self.settings.get_property('protein_name')
            return scipy.stats.pearsonr(
              libi.get_enrichments(k), libj.get_enrichments(k))
        elif test == 'spearman':
            return scipy.stats.spearmanr(
              libi.get_enrichments(k), libj.get_enrichments(k))
        else:
            raise ValueError('Unrecognized test %s' % test)

    def get_libfrac_kmer_concordance(self, libi, libj, test):
        """
        checks that two libaries are similar with regard to library fraction
        """
        k = self.settings.get_property('k_for_concordance_check')
        if test == 'pearson':
            return scipy.stats.pearsonr(
              libi.get_naive_libfracs(k), libj.get_naive_libfracs(k))
        elif test == 'spearman':
            return scipy.stats.spearmanr(
              libi.get_naive_libfracs(k), libj.get_naive_libfracs(k))
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
        for k in self.settings.get_naiveks():
            self.naively_sorted_kmers[k] =\
              self.sort_kmers_by_enrichment(k, 'naive')

    def sort_kmers_by_enrichment(self, k, method='naive'):
        """
        sorts the kmers based on how enriched they are in the
        protein libraries.
        returns an array of the kmers (as indexes not strings) in sorted order
        """
        summed_enrichments = np.zeros(4 ** k)
        barc_set =  self.settings.get_property('barcodes_for_sorting')
        for lib in self.plibs:
            if not lib.get_barcode() in barc_set:
                continue
            if method == 'naive':
                summed_enrichments += lib.get_enrichments(k)
        kmers = range(4 ** k)
        summed_enrich_kmers = zip(summed_enrichments, kmers)
        summed_enrich_kmers.sort(reverse=True)
        top_enrich, sorted_kmers = zip(*summed_enrich_kmers)
        return sorted_kmers


    def make_enrichment_table(self):
        """
        Makes a table of the enrichments for the values of k
        Also, writes a pkl with the same information
        """
        for k in self.settings.get_property('ks_to_test_naive'):
            print 'Writing Naive Enrichment Table for k=%i' % k
            #write table
            table_f = self.get_rdir_fhandle('tables/enrich_naive.%i.txt' % k)
            self.make_table_header(table_f)
            for kmer_i in self.naively_sorted_kmers[k]:
                assert isinstance(kmer_i, int)
                kmer = rbns_utils.get_kmer_from_index(k, kmer_i)
                table_f.write(kmer)
                for lib in self.plibs:
                    table_f.write('\t%f' %  lib.get_enrichment(k, kmer_i))
                table_f.write('\n')
            table_f.close()
            print 'Writing Naive Enrichment pkl for k=%i' % k
            barcode2kmer2enrich = dict()
            for lib in self.plibs:
                barcode2kmer2enrich[lib.get_barcode()]\
                  = lib.get_naive_enrichment_dict(k)
            pkl_f = self.get_rdir_fhandle('tables/enrich_naive.%i.pkl' % k)
            cPickle.dump(barcode2kmer2enrich, pkl_f)
        print 'tables made'


    def make_SKA_table(self):
        for k in self.settings.get_ks('stream'):
            print 'Writing SKA libfrac Table for k=%i' % k
            table_f = self.get_rdir_fhandle('tables', 'SKA_libfrac.%i.txt' % k)
            self.make_table_header(table_f)
            for kmer_i in self.naively_sorted_kmers[k]:
                assert isinstance(kmer_i, int)
                kmer = rbns_utils.get_kmer_from_index(k, kmer_i)
                table_f.write(kmer)
                for lib in self.plibs:
                    table_f.write('\t%g' % lib.get_stream_libfrac(k, kmer_i))
                table_f.write('\n')
            table_f.close()

    def make_enrichment_hist_stacked(self):
        """
        make a histogram of enrichments for each barcode
        """
        protein = self.settings.get_property('protein_name').upper()
        for k in self.settings.get_property('ks_to_test_naive'):
            for lib in self.plibs:
                hist_file = self.rdir_path(
                  'plots',
                  'stacked.hist.%i.%0.0f' % (k, lib.get_conc()))
                print 'working on ', hist_file
                fig1 = plt.figure(figsize=(6, 5))
                fig2 = plt.figure(figsize=(6, 5))
                sp = fig1.add_subplot(1, 1, 1)
                sp2 = fig2.add_subplot(1, 1, 1)
                rbns_utils.simpleaxis(sp)
                rbns_utils.simpleaxis(sp2)
                num_bin_edges = 200
                counts, bin_edges = np.histogram(lib.get_enrichments(k), num_bin_edges)
                thresh = np.mean(lib.get_enrichments(k)) + 2 * np.std(lib.get_enrichments(k))
                sp.axvline(thresh, color='black')
                sp2.axvline(thresh, color='black')
                highest = stacked_bar_kmers.plot_stack(sp, bin_edges, lib.get_enrichments(k), k, protein)
                highest2 = stacked_bar_kmers.plot_stack(sp2,
                                             bin_edges,
                                             lib.get_enrichments(k),
                                             k,
                                             protein,
                                             scale='log')
                assert highest == highest2
                sp.set_xlabel('RBNS R value')
                sp.set_ylabel('number of %imers' % k)
                sp2.set_xlabel('RBNS R value')
                sp2.set_ylabel('number of %imers'% k)
                if k == 5 and 'FOX' in self.settings.get_property('protein_name').upper():
                    supt = 'GCATG: %0.2f, GCACG: %0.2f'\
                      % (lib.get_enrichments(k)[rbns_utils.get_index_from_kmer('GCATG')],
                      lib.get_enrichments(k)[rbns_utils.get_index_from_kmer('GCACG')])
                    fig1.suptitle(supt)
                    fig2.suptitle(supt)
                if 'MBNL' in self.settings.get_property('protein_name').upper():
                    ytick_defaults = [0, 1, 5, 10, 50, 100, 500, 1000, 2000, 5000]
                    y_ts_lin = ytick_defaults[0:rbns_utils.getBinIndex_soft_upper(
                      highest, ytick_defaults) + 2]
                    y_ts_log = map(lambda yt: math.log(yt + 1), y_ts_lin)
                    sp2.set_ylim([y_ts_log[0], y_ts_log[-1]])
                    sp2.set_yticks(y_ts_log)
                    sp2.set_xticks(sp2.get_xticks() + [1])
                if 'CUGBP' in self.settings.get_property('protein_name').upper():
                    if k == 7:
                        sp.set_xlim([0, 10])
                        sp2.set_xlim([0, 10])
                        sp.set_xticks([0, 1, 2, 4, 6, 8, 10])
                        sp2.set_xticks([0, 1, 2, 4, 6, 8, 10])
                if k == 6 and 'FOX' in protein:
                    supt = 'TGCATG: %0.2f, TGCACG: %0.2f'\
                      % (lib.get_enrichments(k)[rbns_utils.get_index_from_kmer('TGCATG')],
                      lib.get_enrichments(k)[rbns_utils.get_index_from_kmer('TGCACG')])
                    sp.set_xlim([0, 25])
                    sp2.set_xlim([0, 25])
                    sp.set_xticks([0, 1, 5, 10, 15, 20, 25])
                    sp2.set_xticks([0, 1, 5, 10, 15, 20, 25])
                    fig1.suptitle(supt)
                    fig2.suptitle(supt)
                y_ts_lin = map(lambda ylog: math.exp(ylog) - 1, sp2.get_yticks())
                log_ylabels = [str(int(round(yl))) for yl in y_ts_lin]
                sp2.set_yticklabels(log_ylabels)
                rbns_utils.save_fig(fig1, hist_file, ['png','pdf'])
                rbns_utils.save_fig(fig2, hist_file + '.log', ['png','pdf'])
                if 'FOX' in protein or 'MBNL' in protein or 'CELF' in protein:
                    self.set_legend_texts(sp)
                    self.set_legend_texts(sp2)
                    rbns_utils.save_fig(fig1, hist_file+'.withlegend', ['png','pdf'])
                    rbns_utils.save_fig(fig2, hist_file + '.withlegend.log', ['png','pdf'])
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
        if 'CUGBP' in self.settings.get_property('protein_name').upper():
            return ("contains two UGU's", "contains one UGU", "other significant", "unenriched")
        if 'FOX' in self.settings.get_property('protein_name').upper() or\
          'REG' in self.settings.get_property('protein_name').upper():
            return ("contains GCAUG", "contains GCACG", "other significant", "unenriched")
        if 'MBNL' in self.settings.get_property('protein_name').upper():
            return ("contains UGCU", "contains UGCC", "other significant", "unenriched")
        if 'ESRP' in self.settings.get_property('protein_name').upper():
            return ("contains GGTG", "contains GGT", "other significant", "unenriched")
        if 'U1' in self.settings.get_property('protein_name').upper():
            return ("a", "b", "c", "d")
        return ("kmers")

    def analyze_structure(self):
        """
        for each barcode and readfile checks that the folding energy
        has been calculated and if not calculates it (for each read)
        This can be launched on the cluster.
        After the folding energy has been calculated
        split the reads into bins based on folding energy.
        """
        rbns_utils.make_dir(self.rdir_path('structure'))
        print 'made', self.rdir_path('structure')
        self.run_rna_fold()
        self.split_by_structure()
        self.count_fold_split()
        self.make_fold_split_table()

    def run_rna_fold(self):
        self.waiting_jobs = []
        for lib in self.libs:
            if self.counts_on_cluster:
                job_id = lib.run_rna_fold(True)
                if job_id:
                    self.waiting_jobs.append(job_id)
            else:
                lib.run_rna_fold()
        if self.counts_on_cluster:
            self.wait_for_jobs_to_complete()

    def split_by_structure(self):
        self.waiting_jobs = []
        for lib in self.libs:
            if self.counts_on_cluster:
                job_id = lib.split_by_structure(True)
                if job_id:
                    self.waiting_jobs.append(job_id)
            else:
                lib.split_by_structure(False)
        if self.counts_on_cluster:
            self.wait_for_jobs_to_complete()

    def count_fold_split(self):
        # count split reads
        for barcode, (energy_bin_i, energy_bin) in \
          itertools.product(self.settings.get_property('barcodes'),
          enumerate(self.settings.get_property('energy_bins'))):
            results_file =self.rdir_path(
              'structure',
              '%s.%i.%0.1f_%0.1f.ncount.pkl' %
              (barcode,
              energy_bin_i,
              energy_bin[0],
              energy_bin[1]))
            if not rbns_utils.file_exists(results_file):
                print 'counting', barcode, 'bin #:', energy_bin_i,
                self.count_kmer_fold_split(barcode, energy_bin_i, energy_bin, results_file)
        # load results
        self.barcode2ebini2k2counts = collections.defaultdict(dict)
        for barcode, (energy_bin_i, energy_bin) in \
          itertools.product(self.settings.get_property('barcodes'),
            enumerate(self.settings.get_property('energy_bins'))):
            results_file = self.rdir_path(
              'structure',
              '%s.%i.%0.1f_%0.1f.ncount.pkl' %
              (barcode,
              energy_bin_i,
              energy_bin[0],
              energy_bin[1]))
            if not rbns_utils.file_exists(results_file):
                raise ValueError(results_file + ' not there')
            print 'loading', results_file
            self.barcode2ebini2k2counts[barcode][energy_bin_i] =\
              cPickle.load(open(results_file))

    def make_fold_split_table(self):
        """
        make the table for the structure split
        """
        rbns_utils.make_dir(self.rdir_path('analyses'))
        print self.rdir_path('analyses')
        assert os.path.exists(self.rdir_path('analyses'))
        self.barcode2ebini2k2counts = collections.defaultdict(dict)
        for barcode, (energy_bin_i, energy_bin) in \
          itertools.product(self.settings.get_property('barcodes'),
            enumerate(self.settings.get_property('energy_bins'))):
            results_file = self.rdir_path(
              'structure',
              '%s.%i.%0.1f_%0.1f.ncount.pkl' %
              (barcode,
              energy_bin_i,
              energy_bin[0],
              energy_bin[1]))
            if not rbns_utils.file_exists(results_file):
                raise ValueError(results_file + ' not there')
            print 'loading', results_file
            self.barcode2ebini2k2counts[barcode][energy_bin_i] =\
              cPickle.load(open(results_file))
        k_lib_ebin2total_counts = self.get_ebin_counts()
        for k in self.settings.get_property('ks_to_test_naive'):
            results_file = self.rdir_path(
              'analyses',
              'ebin_enrichments.%i.txt' % k)
            of = open(results_file, 'w')
            of.write('#barcode\t' + '\t'.join(
              [lib.get_barcode()
              for ebin, lib in itertools.product(self.settings.get_property('energy_bins'), self.plibs)]
              ) + '\n')
            of.write('#deltaG\t' + '\t'.join(
              ['%0.1f to %0.1f' % ebin
              for ebin, lib in itertools.product(self.settings.get_property('energy_bins'), self.plibs)]
              ) + '\n')
            of.write('#conc\t' + '\t'.join(
              ['%0.1f' % lib.get_conc()
              for ebin, lib in itertools.product(self.settings.get_property('energy_bins'), self.plibs)]
              ) + '\n')
            barcode2ebini2kmer2enrichment = {}
            for lib in self.plibs:
                barcode2ebini2kmer2enrichment[lib.get_barcode()]\
                  = collections.defaultdict(dict)
            for kmer in rbns_utils.yield_kmers(k):
                of.write(kmer)
                for ebini, lib in itertools.product(
                  range(len(self.settings.get_property('energy_bins'))), self.plibs):
                    ebin = self.settings.get_property('energy_bins')[ebini]
                    kmer_i = rbns_utils.get_index_from_kmer(kmer)
                    lib_counts = self.barcode2ebini2k2counts[
                      lib.get_barcode()][ebini][k][kmer_i]
                    total_lib_counts =\
                      k_lib_ebin2total_counts[(k, lib.get_barcode(), ebini)]
                    input_counts = self.barcode2ebini2k2counts[
                      self.input_lib.get_barcode()][ebini][k][kmer_i]
                    total_input_counts =\
                      k_lib_ebin2total_counts[
                      (k, self.input_lib.get_barcode(), ebini)]
                    enrichment = (1.0 * lib_counts / total_lib_counts) /\
                      (1.0 * input_counts / total_input_counts)
                    barcode2ebini2kmer2enrichment[
                      lib.get_barcode()][ebini][kmer] = enrichment
                    of.write('\t%0.4f' % enrichment)
                of.write('\n')
            of.close()
            results_file = self.rdir_path(
              'analyses',
              'energy_enrichments.%i.pkl' % k)
            cPickle.dump(barcode2ebini2kmer2enrichment,
              open(results_file, 'w'))

    def get_ebin_counts(self):
        """
        returns the total number of reads in each energy bin
        """
        ebin_counts_file = self.rdir_path(
          'analyses', 'ebin_counts.pkl')
        if rbns_utils.file_exists(ebin_counts_file):
            k_lib_ebin2total_counts = cPickle.load(open(ebin_counts_file))
            for k, ebini, lib in itertools.product(
              self.settings.get_property('ks_to_test_naive'),
              range(len(self.settings.get_property('energy_bins'))),
              self.libs):
                if not (k, lib.get_barcode(), ebini) in k_lib_ebin2total_counts:
                    break
            else:
                return k_lib_ebin2total_counts
        k_lib_ebin2total_counts = {}
        for k, ebini, lib in itertools.product(
          self.settings.get_property('ks_to_test_naive'),
          range(len(self.settings.get_property('energy_bins'))),
          self.libs):
            counts = sum(self.barcode2ebini2k2counts[lib.get_barcode()][ebini][k])
            k_lib_ebin2total_counts[(k, lib.get_barcode(), ebini)] = counts
        cPickle.dump(k_lib_ebin2total_counts, open(ebin_counts_file, 'w'))
        return k_lib_ebin2total_counts

    def plot_fold_split_hump(self):
        """
        this is basically the enrichment humps plotter for the
        reads split by free energy

        Only plots the known motif
        """
        fig1 = plt.figure(figsize=(7, 9))
        for ebi, energy_bin in enumerate(self.settings.get_property('energy_bins')):
            sp = fig1.add_subplot(4, 1, 1)
            sp2 = fig1.add_subplot(4, 1, 3)
            motif = self.settings.get_property('known_motif')
            k = len(motif)
            motif_i = rbns_utils.get_index_from_kmer(motif)
            vs = []
            for lib in self.plibs:
                lib_counts = self.barcode2ebini2k2counts[
                  self.input_lib.get_barcode()][ebi][k][motif_i]
                counts = self.barcode2ebini2k2counts[
                  lib.get_barcode()][ebi][k][motif_i]
                input_libfrac = float(lib_counts) / sum(
                  self.barcode2ebini2k2counts[
                  self.input_lib.get_barcode()][ebi][k])
                this_libfrac = float(counts) / sum(self.barcode2ebini2k2counts[
                  lib.get_barcode()][ebi][k])
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
          for energy_bin in self.settings.get_property('energy_bins')]
        sp.legend(legend_labels, loc=2)
        sp.set_ylabel('enrichment')
        sp2.set_ylabel('enrichment')
        fig1.suptitle(self.settings.get_property('experiment_name'))
        fig1.savefig(self.rdir_path('plots', 'FE.split.pdf'))

    def count_kmer_fold_split(self, barcode, bin_i, energy_bin, results_file):
        """
        does a naive count on the reads split by fold energy
        make this code shared with the actual naive counter
        """
        split_reads_file = self.rdir_path(
          'structure',
          '%s.%i.%0.1f_to_%0.1f.reads' %
          (barcode, bin_i, energy_bin[0], energy_bin[1]))
        print 'checking that %s exists' % split_reads_file
        assert os.path.exists(split_reads_file)
        k2counts = dict()
        kmer2index = {}
        for k in self.settings.get_property('ks_to_test_naive'):
            k2counts[k] = np.zeros(4 ** k, dtype=int)
            for i, kmer in enumerate(rbns_utils.yield_kmers(k)):
                kmer2index[kmer] = i
        read_len = self.settings.get_property('read_len')
        for i, line in enumerate(rbns_utils.aopen(split_reads_file)):
            if line == '\n':
                continue
            if 'N' in line:
                continue
            try:
                for k in self.settings.get_naiveks():
                    for ki in range(0, read_len - k + 1):
                        kmer_index = kmer2index[line[ki:ki + k]]
                        k2counts[k][kmer_index] += 1
            except:
                continue
        for k in self.settings.get_property('ks_to_test_naive'):
            assert sum(k2counts[k]) > 1
        cPickle.dump(k2counts, open(results_file, 'wb'))

    def get_all_barcode_handles(self):
        """
        returns a dictionary barcode-> file handles for the split reads
        """
        return {
          lib_settings.barcode: open(lib_settings.get_split_reads(), 'w')
          for lib_settings in self.settings.iter_lib_settings()}

    def split_reads(self):
        """
        Splits all the reads
        """
        for lib_settings in self.settings.iter_lib_settings():
            if not lib_settings.split_reads_exist():
                break
        else:
            return
        total_reads = 0
        bad_barcodes = collections.Counter()
        reads_per_barcode = collections.Counter()
        rbns_utils.make_dir(self.rdir_path('split_reads'))
        barcode2of = self.get_all_barcode_handles()
        barcodes = self.settings.get_property('barcodes')
        print self.settings.get_fastq()
        for l1, l2, l3, l4 in rbns_utils.iter4Lines(self.settings.get_fastq()):
            barcode = rbns_utils.get_barcode(l1)
            barcode = barcode[0:self.settings.get_property('barcode_len')]
            barcode_match = self.get_barcode_match(barcode, barcodes)
            total_reads += 1
            if not barcode_match:
                bad_barcodes[barcode] += 1
                continue
            else:
                trimmed_read = self.trim_read(l2.strip())[0:int(self.settings.get_property('read_len'))]
                barcode2of[barcode_match].write(trimmed_read + '\n')
                reads_per_barcode[barcode_match] += 1
        print 'splitting complete'
        self.write_barcode_log(reads_per_barcode,
          total_reads, bad_barcodes)
        map(lambda f: f.close(), barcode2of.values())

    def write_barcode_log(self, reads_per_barcode, total_reads, bad_barcodes):
        """
        Writes the log of which other barcodes are represented.
        """
        barcode_log_of = self.get_rdir_fhandle('split_reads', 'barcode_log.txt')
        barcode_log_of.write('barcode\tconc\treads assigned\n')
        for barcode, conc in itertools.izip(
          self.settings.get_property('barcodes'), self.settings.get_property('concentrations')):
            barcode_log_of.write('%s\t%s\t%i\n' %
              (barcode,
              str(conc) if barcode != self.settings.get_property('library_seq_barcode')
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
        return full_read[0: len(full_read) - self.settings.get_property('trim_3p')]

    def get_barcode_match(self, barcode, barcodes):
        """
        takes a barcode and returns the one it matches (hamming <= 1)
        else
        empty string
        """
        if barcode in barcodes:
            return barcode
        for barcode_j in barcodes:
            if rbns_utils.hamming_N(barcode, barcode_j) <= self.settings.get_property('mismatches_allowed_in_barcode'):
                return barcode_j
        return ''

    def copy_unsplit_to_wd(self):
        """
        copies and compresses the unsplit fastq to the working directory
        if needed
        """
        fastq = self.settings.get_fastq()
        working_unsplit_file_name = self.rdir_path(
            'unsplit', os.path.basename(fastq))
        if working_unsplit_file_name[-3:] != '.gz':
            working_unsplit_file_name = working_unsplit_file_name + '.gz'
        if working_unsplit_file_name == fastq:
            print 'Unsplit file is defined to be within the working dir'
            pass
        elif rbns_utils.file_exists(working_unsplit_file_name):
            print 'unsplit file in resutls'
        else:
            print 'copying unsplit to working directory'
            start_copy = time.time()
            if fastq[-3:] == '.gz':
                if not os.path.exists(self.rdir_path('unsplit')):
                    os.makedirs(self.rdir_path('unsplit'))
                shutil.copyfile(fastq, working_unsplit_file_name)
            else:
                with rbns_utils.aopen(fastq, 'rb') as in_f:
                    of = self.get_rdir_fhandle(working_unsplit_file_name)
                    of.writelines(in_f)
                    of.close()
                in_f.close()
            end_copy = time.time()
            print 'Copying unsplit file took %f seconds' %\
                  (end_copy - start_copy)
        self.settings.set_fastq(working_unsplit_file_name)

    def rdir_path(self, *args):
        return os.path.join(self.settings.get_rdir(), *args)

    def get_rdir_fhandle(self, *args):
        """
        returns a filehandle to the fname in the rdir
        """
        out_path = self.rdir_path(*args)
        out_dir = os.path.dirname(out_path)
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
        return rbns_utils.aopen(out_path, 'w')


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("settings_file")
    parser.add_argument("--launch-onto-cluster",
                        help="launches the whole thing on the cluster",
                        action='store_true')
    parser.add_argument("--calculate-kds",
                        help="Calculates relative Kd values",
                        action='store_true')
    parser.add_argument("--counts-on-cluster",
                        help="Runs on this node, but submits"
                        "counting jobs to cluster.",
                        action='store_true')
    parser.add_argument("--make-tables",
                        help="Makes tables.",
                        action='store_true')
    parser.add_argument("--make-plots",
                        help="Makes plots.",
                        action='store_true')
    parser.add_argument("--structure-calculations",
                        help="Does the free energy of folding calculations",
                        action='store_true')
    parser.add_argument("--comparisons",
                        help="Does comparisons to other experiments",
                        action='store_true')
    parser.add_argument("--all-tasks",
                        help="Makes plots, tables, folding and comparisons",
                        action='store_true')
    args = parser.parse_args()
    if args.launch_onto_cluster and args.counts_on_cluster:
        print ('Incompatible command line arguments:\n'
               ' Choose either --counts-on-cluster XOR --launch-onto-cluster')
    return args


def launch_on_cluster(args):
    settings_file = os.path.abspath(args.settings_file)
    python_script = os.path.abspath(sys.argv[0])
    assert rbns_utils.file_exists(settings_file)
    settings = rbns_settings.RBNS_settings(settings_file)
    error_dir = settings.get_property('error_dir')
    rbns_utils.make_dir(error_dir)
    out_file = os.path.join(error_dir, 'main.out')
    err_file = os.path.join(error_dir, 'main.err')
    assert rbns_utils.file_exists(python_script)
    assert os.path.exists(error_dir)
    optional_arguments = ''
    if args.all_tasks:
        optional_arguments += ' --all-tasks'
    if args.structure_calculations:
        optional_arguments += ' --structure-calculations'
    if args.make_tables:
        optional_arguments += ' --make-tables'
    if args.make_plots:
        optional_arguments += ' --make-plots'
    if args.comparisons:
        optional_arguments += ' --comparisons'
    command = ('hostname ; python %(python_script)s'
               ' --counts-on-cluster '
               '%(optional_arguments)s '
               '%(settings_file)s '
               '1> %(out_file)s '
               '2> %(err_file)s ' % locals())
    rbns_cluster_utils.launch(
        command,
        jobname='RBNS_pipe',
        error_dir=error_dir,
        ppn='8',
        q='long')


def main():
    """
    run Bnse here
    """
    args = parse_args()
    if not rbns_utils.file_exists(args.settings_file):
        raise ValueError("Settings file %s doesn't exist" % args.settings_file)
    if args.launch_onto_cluster:
        launch_on_cluster(args)
    else:
        settings = rbns_settings.RBNS_settings(args.settings_file)
        b = Bnse(settings, args.counts_on_cluster)
        if args.make_tables or args.all_tasks or args.make_plots:
            b.sort_all_kmers_by_enrichment()
        if args.calculate_kds:
            b.calculate_kds()
        if args.make_tables or args.all_tasks:
            print 'making tables'
            b.make_tables()
        if args.structure_calculations or args.all_tasks:
            print 'starting structure calculations'
            b.analyze_structure()
        if args.make_plots or args.all_tasks:
            print 'making plots'
            b.make_plots()
        if args.comparisons or args.all_tasks:
            print 'doing comparisons'
            b.compare_all_other_experiments()


def count_naive(split_reads, k, out_pkl):
    """ Does the naive count for the given split reads and k and stores result
    """
    print 'doing naive count'
    counts = np.zeros(4 ** k, dtype=int)
    kmer2index = {}
    for i, kmer in enumerate(rbns_utils.yield_kmers(k)):
        kmer2index[kmer] = i
    with rbns_utils.aopen(split_reads) as f:
        read = f.readline()
        read_len = len(read.strip())
        while read:
            if not 'N' in read:
                for ki in range(0, read_len - k + 1):
                    kmer_index = kmer2index[read[ki:ki + k]]
                    counts[kmer_index] += 1
            read = f.readline()
    cPickle.dump(counts, open(out_pkl, 'wb'))


def count_stream(split_reads, k, out_pkl, passes=2):
    """ Does the stream count for the given split reads and k and stores result
    """
    stream_weights = streaming_convergence.stream_counts(
        k, split_reads, passes)
    cPickle.dump(stream_weights, open(out_pkl, 'wb'))


def count_presence(split_reads, k, out_pkl):
    """ Does the naive count for the given split reads and k and stores result
    """
    presence_counts = np.zeros(4 ** k, dtype=int)
    kmer2index = {}
    for i, kmer in enumerate(rbns_utils.yield_kmers(k)):
        kmer2index[kmer] = i
    with rbns_utils.aopen(split_reads) as f:
        read = f.readline()
        read_len = len(read.strip())
        i = 0
        while read:
            kmers_present = set()
            for ki in range(0, read_len - k + 1):
                kmer = read[ki:ki + k]
                if kmer in kmer2index:
                    kmers_present.add(kmer)
            for kmer in kmers_present:
                presence_counts[kmer2index[kmer]] += 1
            read = f.readline()
            i += 1
            if i % 1e6 == 0:
                print i
    cPickle.dump(presence_counts, open(out_pkl, 'wb'))


if __name__ == '__main__':
    main()
