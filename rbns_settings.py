import os
import ConfigParser
import simplejson
import itertools

import rbns_utils

class RBNS_settings:
    def __init__(self, settings_file):
        self.settings_file = settings_file
        self.process_settings(settings_file)
    def get_fastq(self):
        return os.path.abspath(self.settings['fastq'])
    def set_fastq(self, fastq_file):
        self.settings['fastq'] = fastq_file
    def get_ks(self, count_type):
        assert 'ks_to_test_%s' % count_type in self.settings
        settings_ks = self.settings['ks_to_test_%s' % count_type]
        if count_type == 'naive' or self.settings['stream_count'] == True:
            return settings_ks
        else:
            return []
    def get_naiveks(self):
        return self.settings['ks_to_test_naive']
    def get_force_recount(self, count_type):
        return self.settings['force_%s_recount' % count_type]
    def get_settings_file(self):
        return self.settings_file
    def get_property(self, property, default=None):
        try:
            if not property in self.settings and default != None:
                return default
            return self.settings[property]
        except:
            print self.settings
            raise  ValueError('cannot find %s' % property)
    def get_rdir(self):
        rbns_utils.make_dir(self.rdir)
        return self.rdir
    def get_wdir(self):
        rbns_utils.make_dir(self.wdir)
        return self.wdir
    def get_input_barcode(self):
        return self.settings['library_seq_barcode']

    def isCountRequested(self, count_type):
        if count_type == 'naive':
            return self.settings['naive_count']
        elif count_type == 'stream' or count_type == 'presence':
            return self.settings['stream_count']

    def iter_lib_settings(self):
        for i in range(len(self.barcodes)):
            yield RBNS_lib_settings(self,
              self.barcodes[i],
              self.settings['concentrations'][i],
              self.settings['poly_ic_conc'][i],
              self.settings['input_rna'][i],
              self.settings['washes'][i],
              self.settings['temperature'][i])

    def process_settings(self, settings_file):
        """
        reads the settings file and converts str to float or list or whatever
        stores result in self.settings as a dict()
        """
        int_keys = ['read_len',
          'mismatches_allowed_in_barcode', 'trim_3p', 'max_reads_to_split',
          'k_for_concordance_check',
          'num_kmers_for_enrichment_curves']
        float_keys = []
        boolean_keys = ['naive_count', 'stream_count', 'force_stream_recount', 'force_presence_recount', 'force_naive_recount']
        list_str_keys = ['barcodes', 'barcodes_for_sorting',
          'temperature', 'relevant_variables', 'experiments_to_compare',
          'motifs_of_interest']
        list_int_keys = ['ks_to_test_stream','ks_to_test_presence',
          'ks_to_test_naive', 'ks_for_matlab', 'washes']
        list_float_keys = ['concentrations', 'poly_ic_conc',
          'free_energy_limits', 'input_rna']
        extant_files = ['fastq']
        config = ConfigParser.ConfigParser()
        config.read(settings_file)
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
            assert rbns_utils.file_exists(settings[k])
        self.settings = settings
        self.wdir = settings['working_dir']
        self.rdir = settings['results_dir']
        self.check_barcode_lens()
        self.barcodes = self.settings['barcodes']
        self.check_barcodes_are_separated()
        self.energy_bins = \
          [(conc1, conc2) for conc1, conc2 in
          zip(self.settings['free_energy_limits'][0:],
          self.settings['free_energy_limits'][1:])]
        self.settings['energy_bins'] = self.energy_bins
        self.known_motif = self.settings['known_motif']
        self.known_motif_k = len(self.known_motif)
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
        self.protein_name = settings['protein_name']
        if settings['stream_count'] and not set(settings['ks_to_test_stream']).issubset(
          settings['ks_to_test_naive']):
            raise ValueError('All ks for stream must also be in naive')

    def check_barcode_lens(self):
        """
        verifies that all the barcodes are the same length
        """
        barcode_lens = set(map(len, self.settings['barcodes']))
        if 1 != len(barcode_lens):
            raise ValueError('all barcodes must be the same length')
        self.barcode_len = barcode_lens.pop()
        self.settings['barcode_len'] = self.barcode_len

    def check_barcodes_are_separated(self):
        """
        makes sure the barcodes are all totally distinguishable
        """
        for b1, b2 in itertools.combinations(self.settings['barcodes'], 2):
            hamming_dist = rbns_utils.hamming_distance(b1, b2)
            if hamming_dist < 2:
                raise ValueError('The barcodes supplied are not well '
                  'separated: %s-%s' % (b1, b2))



class RBNS_lib_settings:
    def __init__(self, experiment_settings, barcode, conc, poly_ic_conc, input_rna, washes, temperature):
        self.experiment_settings = experiment_settings
        self.barcode = barcode
        self.conc = conc
        self.poly_ic_conc = poly_ic_conc
        self.input_rna = input_rna
        self.washes = washes
        self.temperature = temperature

    def get_rna_conc(self):
        return self.input_rna
    def is_input(self):
        return self.barcode == self.experiment_settings.get_input_barcode()

    def get_property(self, property):
        return self.experiment_settings.get_property(property)

    def get_barcode(self):
        return self.barcode

    def get_conc(self):
        return self.conc

    def counts_file(self, count_type, k):
        counts_file = os.path.join(
          self.experiment_settings.get_rdir(),
          'counts',
          count_type,
          '%(protein)s_%(barcode)s_%(k)i.pkl' %
          {'barcode': self.barcode,
           'protein': self.experiment_settings.get_property('protein_name'),
           'k': int(k)})
        return counts_file

    def counts_exist(self, count_type, k):
        return rbns_utils.file_exists(self.counts_file(count_type, k))

    def get_split_reads(self):
        split_reads = os.path.join(
          self.experiment_settings.get_rdir(),
          'split_reads',
          '%(protein)s_%(barcode)s.reads' %
           {'protein': self.experiment_settings.get_property('protein_name'),
          'barcode': self.barcode})
        return split_reads

    def split_reads_exist(self):
        split_reads = self.get_split_reads()
        return rbns_utils.file_exists(split_reads)



