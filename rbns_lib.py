import matplotlib.pyplot as plt
import scipy.stats
import subprocess
import os
import rbns_cluster_utils
import cPickle
import rbns_utils
import numpy as np

class RBNS_Lib:
    def __init__(self, experiment_settings, lib_settings):
        """
        Constructor for Library class
        """
        self.experiment_settings = experiment_settings
        self.lib_settings = lib_settings
        self.load_counts()
        self.get_property = self.experiment_settings.get_property
        self.get_rdir = experiment_settings.get_rdir
        self.get_wdir = experiment_settings.get_wdir

    def get_conc(self):
        return self.lib_settings.conc

    def get_washes(self):
        return self.lib_settings.washes

    def get_poly_ic(self):
        return self.lib_settings.poly_ic_conc

    def get_temperature(self):
        return self.lib_settings.temperature

    def get_rna_conc(self):
        return self.lib_settings.input_rna

    def is_input(self):
        return self.lib_settings.is_input()

    def get_barcode(self):
        """ returns the library's barcode """
        return self.lib_settings.barcode

    def get_stream_libfrac(self, k, kmer_i):
        return self.type2k2counts['stream'][k].get_libfrac(kmer_i)

    def get_stream_libfrac_kmer(self, kmer):
        k = len(kmer)
        return self.type2k2counts['stream'][k].kmer_value(kmer)

    def get_presence_frac_kmer(self, kmer):
        k = len(kmer)
        return self.type2k2counts['presence'][k].kmer_value(kmer)

    def get_presence_frac(self, k, kmer_i):
        return self.type2k2counts['presence'][k].kmeri_value(kmer_i)

    def load_counts(self):
        self.type2k2counts = {}
        for count_type in ['naive', 'stream', 'presence']:
            if not self.experiment_settings.isCountRequested(count_type):
                continue
            self.type2k2counts[count_type] = {}
            for k in self.experiment_settings.get_ks(count_type):
                self.type2k2counts[count_type][k] =\
                  RBNS_profile(self.lib_settings, k, count_type)

    def get_naive_counts(self, k):
        return self.type2k2counts['naive'][k].get_profile()

    def get_naive_libfracs(self, k):
        return self.type2k2counts['naive'][k].get_libfracs()

    def get_presence_fracs(self, k):
        return self.type2k2counts['presence'][k].get_profile()

    def get_stream_libfracs(self,k):
        return self.type2k2counts['stream'][k].get_libfracs()

    def get_naive_count_by_index(self, k, kmer_i):
        return self.type2k2counts['naive'][k].kmeri_value(kmer_i)

    def get_naive_count_kmer(self, kmer):
        k = len(kmer)
        return self.type2k2counts['naive'][k].kmer_value(kmer)

    def get_enrichments(self, k):
        return self.type2k2counts['naive'][k].get_enrichments()

    def get_enrichment(self, k, kmer_i):
        return self.type2k2counts['naive'][k].get_enrichment(kmer_i)

    def get_enrichment_kmer(self, kmer):
        k = len(kmer)
        kmer_i = rbns_utils.get_index_from_kmer(kmer)
        return self.get_enrichment(k, kmer_i)

    def calculate_enrichment(self, k, input_lib):
        enrich_pkl = os.path.join(
          self.experiment_settings.get_rdir(),
          'tables',
          '%s_%s.%i_enrichment.pkl'  %
          (self.experiment_settings.get_property('experiment_name'),
          self.lib_settings.get_barcode(), k))
        if rbns_utils.file_exists(enrich_pkl):
            self.type2k2counts['naive'][k].load_enrichments(enrich_pkl)
        else:
            input_profile = input_lib.type2k2counts['naive'][k]
            self.type2k2counts['naive'][k].calculate_enrichments(input_profile)
            self.type2k2counts['naive'][k].save_enrichments(enrich_pkl)

    def compare_top_kmers(self, k, most_enriched_lib, num_top_kmers_to_comp):
        """
        compares the enrichment of the top kmers in this library to another
        library (usually the most enriched one).

        Returns pearsonr
        """
        top_kmers = most_enriched_lib.get_top_kmers(k, num_top_kmers_to_comp)
        most_enriched_lib_enrichments =\
          [most_enriched_lib.get_enrichment(k, kmer_i) for kmer_i in top_kmers]
        this_lib_enrichments =\
          [self.get_enrichment(k, kmer_i) for kmer_i in top_kmers]
        r, p = scipy.stats.pearsonr(
          most_enriched_lib_enrichments,
          this_lib_enrichments)
        return r

    def get_top_kmers(self, k, num_top_kmers_to_compare):
        """
        returns the top kmers by enrichment in this library
        """
        top_kmer_file = os.path.join(self.get_rdir(), 'analyses',
          'top_kmers.%i.%s.%i.pkl' %
          (k, self.get_barcode(), num_top_kmers_to_compare))
        if os.path.exists(top_kmer_file):
            sorted_kmers = cPickle.load(open(top_kmer_file))
            if len(top_kmer_file) == num_top_kmers_to_compare:
                return sorted_kmers
        enrich_kmers = zip(self.get_enrichments(k), range(4 ** k))
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
        top_enrichments = sorted(self.get_enrichments(k))[-num_next_kmers - 1:]
        return top_enrichments[-1] / np.mean(top_enrichments[:-1])

    def get_max_enrichment(self, k_for_max_enrichment):
        """
        returns the largest enrichment kmer of size k
        """
        return max(self.get_enrichments(k_for_max_enrichment))

    def get_full_label(self):
        """
        returns a good human readable label (for say a graph)
        """
        out = '%i, ' % int(self.get_conc())
        rv2desc = {'input_rna': 'rna conc', 'poly_ic': '[poly(IC)]',
          'temp': 'T', 'barcode': 'barcode'}
        rv2val = {'input_rna': self.lib_settings.get_rna_conc(),
                  'poly_ic': self.get_poly_ic(),
                  'temp': self.get_temperature(),
                  'barcode': self.get_barcode()}
        for rv in self.experiment_settings.get_property('relevant_variables', []):
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
        for k in self.experiment_settings.get_property('ks_to_test_naive')[:-1]:
            next_frac.append(self.calculate_next_kmer_frac_counts(k))
        sp.plot(self.experiment_settings.get_property('ks_to_test_naive')[:-1],
          next_frac, 'b-')
        sp.set_ylim([0, 1.1 * max(next_frac)])
        sp.set_title('counts')
        sp = fig1.add_subplot(212)
        next_enrich = []
        for k in self.experiment_settings.get_property('ks_to_test_naive')[:-1]:
            next_enrich.append(self.calculate_next_kmer_frac_enrich(k))
        sp.plot(self.experiment_settings.get_property('ks_to_test_naive')[:-1],
          next_enrich, 'b-')
        sp.set_ylim([0, 1.1 * max(next_enrich)])
        sp.set_title('enrichment')
        fig1.savefig(
          os.path.join(self.experiment_settings.get_rdir(), 'plots',
          '%s.%0.1f.kpick.pdf' % (self.get_barcode(), self.get_conc())))

    def calculate_next_kmer_frac_counts(self, k):
        """
        for a given k:
        calculates the number of counts for addding a base to either side
        of the kmer and checks the counts for each of those.
        then calculates the proportion of the counts that go to the highest
        one of those and returns that fraction
        """
        assert k in self.experiment_settings.get_property('ks_to_test_naive') and\
          k + 1 in self.experiment_settings.get_property('ks_to_test_naive')
        kmer_counts = self.type2k2counts['naive'][k].get_profile()
        max_counts = max(kmer_counts)
        top_kmer_i = list(kmer_counts).index(max_counts)
        kmer = rbns_utils.get_kmer_from_index(k, top_kmer_i)
        adjacent_kmers = rbns_utils.get_adjacent_kmers(kmer)
        ajacent_kmer_i_s = map(rbns_utils.get_index_from_kmer, adjacent_kmers)
        adjacent_kmer_counts = [self.get_naive_count_by_index(k + 1, kmer_i)
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
        assert k in self.experiment_settings.get_property('ks_to_test_naive') and\
          k + 1 in self.experiment_settings.get_property('ks_to_test_naive')
        kmer_enrichments = self.get_enrichments(k)
        max_enrichment = max(kmer_enrichments)
        # assume only one max
        top_kmer_i = list(kmer_enrichments).index(max_enrichment)
        kmer = rbns_utils.get_kmer_from_index(k, top_kmer_i)
        adjacent_kmers = rbns_utils.get_adjacent_kmers(kmer)
        ajacent_kmer_i_s = map(rbns_utils.get_index_from_kmer, adjacent_kmers)
        adjacent_kmer_enrichments = [self.get_enrichment(k + 1, kmer_i)
          for kmer_i in ajacent_kmer_i_s]
        best_adjeacent_kmer_enrichment = max(adjacent_kmer_enrichments)
        adjacent_kmer_frac = float(best_adjeacent_kmer_enrichment)\
          / max_enrichment
        return adjacent_kmer_frac

    def calculate_kds(self, k):
        return
        #TODO asdf integrate kd code
        #:1
        #calculate_kds_importable.asdf


    def split_reads_exist(self):
        """
        returns true if the split reads file for this library exists
        and is non empty
        does not check if it is complete
        """
        return rbns_utils.file_exists(self.get_split_reads())

    def get_split_reads(self):
        """
        returns the full path of this library's split reads file
        """
        split_reads_file = os.path.join(
          self.get_rdir(),
          'split_reads',
          '%s_%s.reads' % (
          self.experiment_settings.get_property('protein_name'),
          self.get_barcode()))
        assert split_reads_file == self.lib_settings.get_split_reads()
        return self.lib_settings.get_split_reads()

    def get_split_whandle(self):
        """
        returns a write file handle to the split reads
        """
        return rbns_utils.aopen(self.get_split_reads(), 'w')

    def run_rna_fold(self, launch=False):
        """
        Folds all the reads using RNAfold and stores the result
        will launch on cluster if requested
        """
        out_file = os.path.join(self.get_rdir(), 'structure',
          '%s.fe' % self.get_barcode())
        if rbns_utils.file_exists(out_file + '.gz'):
            return 0
        print 'calculating structure for'\
          ' %s with concentration %0.1f' % (self.get_barcode(), self.get_conc())
        tmp_file = os.path.join(self.get_wdir(), '%s.reads' % self.get_barcode())
        tmp_file_out = os.path.join(self.get_wdir(), '%s.fe' % self.get_barcode())
        err_file = os.path.join(
          self.experiment_settings.get_property('error_dir'),
          '%s.err' % self.get_barcode())
        command = 'hostname ; mkdir -p %s; cd %s; cp %s %s ; cat %s | RNAfold '\
          '1> %s 2> %s; cp %s %s ; gzip %s' %\
          (self.get_wdir(), self.get_wdir(), self.get_split_reads(),
          tmp_file, tmp_file,
          tmp_file_out, err_file,
          tmp_file_out,
          os.path.join(self.get_rdir(), 'structure'),
          out_file)
        if not launch:
            subprocess.Popen(command, shell=True).wait() # note that this waits
        else:
            script_options = {'nodes': '1', 'ppn': '8',
              'outf': self.experiment_settings.get_property('experiment_name')
              + '.submission_log',
              'jobname': self.experiment_settings.get_property('experiment_name')
              + '_' + self.get_barcode() + '.fold',
              'queue': 'long', 'workingdir': self.get_rdir(),
              'command': command}
            return rbns_cluster_utils.launch(
              command,
              script_options,
              error_dir=self.experiment_settings.get_property('error_dir'))

    def split_by_structure(self, launch=False):
        """
        splits the reads for this library into seperate files
        binned by the folding free energy if they don't exist and aren't empty

        bins are defined in the settings file

        This can throw errors if there are no reads in one of the defined
        energy bins.

        will launch on cluster if requested
        """
        out_files = ['%s/structure/%s.%i.%0.1f_to_%0.1f.reads' %
          (self.get_rdir(), self.get_barcode(), i, dG1, dG2)
          for i, (dG1, dG2) in enumerate(self.experiment_settings.get_property('energy_bins'))]
        rbns_utils.make_dir(os.path.dirname(out_files[0]))
        if all(map(rbns_utils.file_exists, out_files)):
            return
        if launch:
            return rbns_cluster_utils.launch_energy_splitter(self.lib_settings, self.experiment_settings)
        else:
            run_energy_splitter(self.lib_settings, self.experiment_settings)

    def get_naive_enrichment_dict(self, k):
        """
        returns a dictionary of kmer -> naive enrichment
        """
        return self.type2k2counts['naive'][k].get_enrichment_dict()

    def get_B_kmer(self, kmer, read_len):
        return self.type2k2counts['naive'][len(kmer)].get_B_kmer(kmer, read_len)

    def calcB(self, kmer):
        """
         calculates the B value
        """
        read_len = self.experiment_settings.get_property('read_len')
        return self.type2k2counts['naive'][len(kmer)].get_B_kmer(kmer, read_len)


class RBNS_profile:
    def __init__(self, lib_settings, k, count_type):
        self.count_type = count_type
        self.k = k
        counts_pkl = lib_settings.counts_file(count_type, k)
        self.profile = cPickle.load(open(counts_pkl, 'rb'))
        assert len(self.profile) == 4 ** self.k
        if self.count_type == 'naive' or self.count_type == 'stream':
            self.calculate_libfracs()
            assert 0.99 < sum(self.libfrac) < 1.01

    def kmer_value(self, kmer):
        kmeri = rbns_utils.get_index_from_kmer(kmer)
        return self.kmeri_value(kmeri)

    def kmeri_value(self, kmeri):
        return self.profile[kmeri]

    def calculate_libfracs(self):
        counts = np.array(self.profile, dtype=float)
        total_counts = sum(counts)
        if not total_counts:
            raise ValueError('Naive counts are 0')
        self.libfrac = counts / total_counts

    def get_libfracs(self):
        assert self.count_type == 'naive' or self.count_type == 'stream'
        return self.libfrac

    def get_libfrac(self, kmer_i):
        assert self.count_type == 'naive' or self.count_type == 'stream'
        return self.libfrac[kmer_i]

    def get_B(self, kmer_i, read_len):
        assert kmer_i
        enrichment = self.get_enrichment(kmer_i)
        B = rbns_utils.chris_formula(enrichment, self.k, read_len)
        return B

    def get_B_kmer(self, kmer, read_len):
        assert len(kmer) == self.k
        enrichment = self.get_enrichment_kmer(kmer)
        B = rbns_utils.chris_formula(enrichment, self.k, read_len)
        return B

    def get_B_values(self, read_len):
        return [rbns_utils.chris_formula(enrich, self.k, read_len) for enrich in self.enrichments]

    def weight_dict(self):
        kmer2weight = {}
        for kmer, weight in zip(rbns_utils.yield_kmers(self.k), self.profile):
            kmer2weight[kmer] = weight
        return kmer2weight

    def get_profile(self):
        return self.profile

    def save_enrichments(self, enrich_pkl):
        cPickle.dump(self.enrichments, open(enrich_pkl, 'wb'))

    def load_enrichments(self, enrich_pkl):
        self.enrichments = cPickle.load(open(enrich_pkl, 'rb'))

    def calculate_enrichments(self, input_profile):
        self.enrichments = norm_libfracs(self, input_profile)

    def get_enrichment_kmer(self, kmer):
        kmer_i = rbns_utils.get_index_from_kmer(kmer)
        return self.get_enrichment(kmer_i)

    def get_enrichments(self):
        return self.enrichments

    def get_enrichment(self, kmer_i):
        return self.enrichments[kmer_i]

    def get_enrichment_dict(self):
        assert len(self.enrichments) == 4 ** self.k
        return {kmer: enrich
                for kmer, enrich in
                zip(rbns_utils.yield_kmers(self.k),
                self.enrichments)}


def norm_libfracs(profile, input_profile):
    assert profile.k == input_profile.k
    return profile.get_libfracs() / input_profile.get_libfracs()


def run_energy_splitter(lib_settings, experiment_settings):
    barcode = lib_settings.get_barcode()
    rdir = experiment_settings.get_rdir()
    reads_file = lib_settings.get_split_reads()
    energy_file = os.path.join(experiment_settings.get_rdir(),
      'structure', '%s.fe.gz' % lib_settings.get_barcode())
    energy_bins = str(experiment_settings.get_property('free_energy_limits'))
    error_dir = experiment_settings.get_property('error_dir')
    rbns_cluster_utils.energy_splitter_commandline(
      barcode,
      rdir,
      reads_file,
      energy_file,
      energy_bins,
      error_dir)
