import os

import scipy.stats
import numpy as np
import matplotlib.pyplot as plt
import cPickle

import aColors
import aUtils
import table2dict
from bind_n_seq_pipeline import get_index_from_kmer, yield_kmers, simpleaxis
from rna import rna


def normalize(arr):
    """
    normalizes a numpy array such that the sum of the vector is 1
    """
    return arr / np.sum(arr)


def relatavize_kds(kmer2abs_kd):
    """
    divides the kds by the best kd and returns a new dictionary
    """
    best_kd = min(kmer2abs_kd.values())
    return {k: v / best_kd  for k, v in kmer2abs_kd.items()}


def compare_predicted_kds_to_spr(kmer2predicted_kd, plot=False):
    """
    for how well a set of kds conform to the measured auweter
    """
    relativized_predictions = relatavize_kds(kmer2predicted_kd)
    fig1, sp = plt.subplots(figsize=(6,6))
    sprs, predictions = [],[]
    for kmer, relative_kd in filter(
      lambda (kmer, kd): kd < kd_prediction_ceiling, 
      relativized_predictions.items()):
        spr_value = kmer2relative_kd_to_compare[kmer]
        xshift_factor = 1.2
        sp.loglog(relative_kd, spr_value, 'o',
          color=aColors.fox_colors(kmer, spr_value < 800))
        sp.text(relative_kd*xshift_factor, spr_value, rna(kmer))
        sprs.append(spr_value)
        predictions.append(relative_kd)
    sp.set_xlabel(r'relative $K_{d}$ Bind-n-seq')
    sp.set_ylabel(r'relative $K_{d}$ SPR')
    sp.set_xlim([.1, 1e4])
    sp.set_ylim([.1, 1e4])
    sp.set_aspect(1)
    #pearson
    r,p = scipy.stats.pearsonr(sprs, predictions)
    sp.text(.2, 5000., r'Pearson r=%.3g, p=%.2g' % (r,p))
    #sp.plot([1, 1e3],[1, 1e3], '-', marker=None)
    #sp.text(.2,3000, r'n=%i, k=%i' % (len(sprs), k))
    aUtils.make_dir(os.path.dirname(plot))
    print plot
    simpleaxis(sp)
    aUtils.save_fig(fig1, plot, ['pdf', 'png'])


def get_complex_concs(complex_conc_file='best_estimates.txt', protein_concs = [0., 1.5, 4.5, 14., 40.5, 121., 365., 1100., 3300., 9800.]):
    """
    gets the concentrations as estimated by bioanalyzer
    reads what's in the file best_estimates.txt
    """
    #pull values in pg/Ul of RNA
    t2d = table2dict.table2dict(complex_conc_file)
    complex_concs = [float(t2d['%g' % conc]['[Rcomplex]'])
                     for conc in protein_concs]
    return complex_concs


def convert_complex_units(RNA_oligo_complexed_conc,
                          conversion_factor=0.03650967):
    """
     converts RNA concentration pg/uL to nM (I hope, verify this arithmatic)
     1pg/uL * (1mol/(27390g) nano*1e9
    """
    return RNA_oligo_complexed_conc * conversion_factor


def get_libfracs_library_presence():
    """
    makes comparison based on the 'presence' metric for seq measurement
    """
    in_file = '%s/counts/presence/%s_%s.%i.pkl' % (edir, experiment_name, library_barcode, k)
    input_libfrac = np.array(cPickle.load(open(in_file)), dtype=float)
    for kmeri, ilf in enumerate(input_libfrac):
        if not ilf:
            input_libfrac[kmeri] = (40. - k + 1.) / (4. ** k)
    return input_libfrac


def get_stream_libfracs():
    """
     gets the library fractions for all the barcodes
     using the SKA algorithm
    """
    in_files = ['%s/counts/stream/%s_%s.%s.pkl' %
      (edir, experiment_name, barcode, k) for barcode in barcodes]
    return [normalize(cPickle.load(open(in_file))) for in_file in in_files]


def get_naive_libfracs():
    """
    gets the naive libfracs for all barcodes
    takes if from most enriched for size k
    
    returns list of arrays
    """
    in_files = ['%s/counts/naive/%s_%s.pkl' %
      (edir, experiment_name, barcode) for barcode in barcodes]
    lib_naive_counts = [np.array(cPickle.load(open(in_file))[k],
      dtype=float) for in_file in in_files]
    lib_naive_counts = 33 * map(normalize, lib_naive_counts) # why 33?
    return lib_naive_counts 


def calcualte_kds_with_ligand_conc_bioanalyzer(SKAs, input_presences, complex_conc):
    """
    calculates the kd based on the bioanalyzer conc, the protein conc and the complex conc.
    """
    estimated_kds = dict()
    for kmer in filter(lambda kmer: len(kmer) == k, kmer2relative_kd_to_compare):
        kmeri = get_index_from_kmer(kmer)
        Litotal = input_presences[kmeri] * RNA_oligo_concentration
        RLi = SKAs[kmeri] * complex_conc
        estimated_kds[kmer] = (Litotal - RLi) / RLi
    return estimated_kds


def plot_naive_v_stream2(naive_libfracs, stream_libfracs):
    """
    makes a plot comparing streaming libfrac and Niave libfracs
    
    labels outliers using Athma's code
    """
    for barcode, conc, stream_libfrac, naive_libfrac, kmer in \
      zip(barcodes, protein_concs, stream_libfracs, naive_libfracs, yield_kmers(k)):
        if conc != 121 and conc != 365:
            continue
        fig4 = plt.figure()
        sp = fig4.add_subplot(1, 1, 1)
        simpleaxis(sp)
        #sp.loglog(naive_libfrac, stream_libfrac, '.')
        values = np.vstack([np.log(naive_libfrac), np.log(stream_libfrac)])
        kernel = scipy.stats.gaussian_kde(values)
        sp.set_xlabel(r'Naive $F_{i}$')
        sp.set_ylabel(r'SKA $F_{i}$')
        sp.set_xlim([1e-5, 1e-1])
        sp.set_ylim([1e-7, 1.2e-1])
        sp.set_aspect(1.0)
        kerns = np.zeros(4 ** k)
        of = open('comps/kerns.%i.txt' % k, 'w')
        zkerns = kernel(values)
        enrich_thresh = np.mean(naive_libfrac) + 2 * np.std(naive_libfrac)
        for kmeri, kmer, x, y, kern in zip(range(4 ** k),
          yield_kmers(k),
          naive_libfrac, 
          stream_libfrac, zkerns):
            color = aColors.protein_colors(kmer, 'fox', x >= enrich_thresh)
            of.write('%s\t%g\n' % (kmer, zkerns[kmeri]))
            kerns[kmeri] = kern
            sp.loglog(x,y,'.',color=color)
            if y > 0.01 or 'GCATG' in kmer:
                print y, kmer
                sp.text(1.08 * x, y *.95, kmer.replace('T', 'U'), fontsize=5)
            if y > 0.1 or (zkerns[kmeri] < 0.0045 and not (
              'GAAGCA' in kmer or
              'CATGG' in kmer or
              'TGCTTG' in kmer or 
              'TTGCA' in kmer or 
              'AAAAAA' in kmer or 
              'GAACCA' in kmer)):
                sp.text(1.08 * x, y *.95, kmer.replace('T', 'U'), fontsize=5)
        simpleaxis(sp)
        aUtils.save_fig(fig4, 'comps/streamVnaive.%g' % conc,
          extentions=['png', 'pdf'])


def plot_naive_v_stream(naive_libfracs, stream_libfracs):
    """
    makes a plot comparing streaming libfrac and Niave libfracs
    
    labels outliers using Athma's code
    """
    for barcode, conc, stream_libfrac, naive_libfrac in \
      zip(barcodes, protein_concs, stream_libfracs, naive_libfracs):
        if conc != 121:
            continue
        fig4 = plt.figure()
        sp = fig4.add_subplot(1, 1, 1)
        sp.loglog(naive_libfrac, stream_libfrac, '.')
        sp.set_xlabel('naive libfracs')
        sp.set_ylabel('Streaming libfracs')
        sp.set_xlim([1e-8, 1e-1])
        sp.set_ylim([1e-8, 1])
        simpleaxis(sp)
        aUtils.save_fig(fig4, 'comps/streamVnaive.%g' % conc, extentions=['png'])
        fig4_huge = plt.figure(figsize=(40, 40))
        sph = fig4_huge.add_subplot(1, 1, 1)
        sph.loglog(naive_libfrac, stream_libfrac, '.')
        sph.set_xlabel('naive libfracs')
        sph.set_ylabel('stream libfracs')
        sph.set_xlim([1e-8, 1e-1])
        sph.set_ylim([1e-8, 1])
        for kmer, x, y in zip(yield_kmers(k), naive_libfrac, stream_libfrac):
            sph.text(x, y, kmer, fontsize=5)
        simpleaxis(sph)
        aUtils.save_fig(fig4_huge, 'comps/streamVnaive.%g_huge' % conc, extentions=['png'])


def main():
    """
    the main function. duh.
    """
    print 'getting stuff'
    RNA_oligo_complexed_concs = get_complex_concs(complex_conc_file, protein_concs)
    complex_concs = map(convert_complex_units, RNA_oligo_complexed_concs)
    input_libfrac_presence = get_libfracs_library_presence()  # its the frac of oligos with the kmer
    stream_libfracs = get_stream_libfracs()
    print 'how much complex for how much protein'
    print zip(complex_concs, protein_concs)
    for i in range(len(stream_libfracs)):
        SKAperPresence_ligand_concBioanalyzer_estimates =\
          calcualte_kds_with_ligand_conc_bioanalyzer(
          stream_libfracs[i],
          input_libfrac_presence,
          complex_concs[i])
        print experiment_name
        plot_dir_ext = '_lambertSPR' if kmer2relative_kd_to_compare\
          == kmer2relative_lambert_kd else '_auweter'
        compare_predicted_kds_to_spr(
          SKAperPresence_ligand_concBioanalyzer_estimates,
          plot='estimates_plots%s/%s/SKAperPres_bioanal/k%i.%i' % (
          plot_dir_ext, experiment_name, k, i))



if __name__ == '__main__':
    use_new_data = True
    use_old_data = False
    if use_new_data:
        protein_concs = [0., 1.5, 4.5, 14., 40.5, 121., 365., 1100., 3300., 9800.]
        kd_prediction_ceiling = 1000
        barcodes = ["ATCACG", "CGATGT", "TTAGGC", "TGACCA", "ACAGTG", "GCCAAT", "CAGATC", "ACTTGA", "GATCAG", "TAGCTT"]
        RNA_oligo_concentration = 1000
        assert len(barcodes) == len(protein_concs)
        library_barcode = 'GGCTAC'
        experiment_name = 'Fox_40'
        complex_conc_file = 'best_estimates.txt'
    elif use_old_data:
        protein_concs = [0.0, 0.5, 1.0, 2.0, 12.0, 62.0, 312.0, 1000]
        kd_prediction_ceiling = 1e9
        barcodes = ["ATCACG", "CGATGT", "TTAGGC", "TGACCA", "ACAGTG", "GCCAAT", "CAGATC", "ACTTGA"]
        RNA_oligo_concentration = 1000
        assert len(barcodes) == len(protein_concs)
        library_barcode = 'GATCAG'
        experiment_name = 'Fox_40_old'
        complex_conc_file = 'best_estimates.original.txt'

    k = 6
    edir = '/net/utr/data/atf/alexrson/bind_n_seq_pipeline/'\
      'results/%s/' % experiment_name
    #naive_libfracs = get_naive_libfracs()
    #stream_libfracs = get_stream_libfracs()
    #plot_naive_v_stream2(naive_libfracs, stream_libfracs)
    kmer2auweter_kd = { 'TGCATG': 1.6019, 
        'AGCATG': 9.264, 'TGCACG': 9.457,
        'CGCATG': 11.773, 'TACATG': 391.79,
        'TGTATG': 540.4, 'TGCATA': 3531.9,
        'TGCATGT': 1.6019, 'AGCATGT': 9.264,
        'TGCACGT': 9.457, 'CGCATGT': 11.773,
        'TACATGT': 391.79, 'TGTATGT': 540.4,
        'TGCATAT': 3531.9}
    kmer2lambert_kd = {
      'TGCATGT': 2.3,
      'GGCATGT': 10.5,
      'CGCATGT': 5.0,
      'AGCATGT': 6.6,
     # 'TGCATGG': 1.8, # more of a structural thing possibly throw out according to Nicole
      'TGCACGT': 14.8,
      'AGCACGT': 55.3,
      'TGCATTA': 640.,
      'TGCATCT': 46,
      'TGCATTT': 330,
      'TGCATG': 2.3,
      'GCATGT': 10.5,
      'CGCATG': 5.0,
      'AGCATG': 6.6,
      'TGCACG': 14.8,
      'AGCACG': 55.3,
      #'TGCATT': 640.,
      'TGCATC': 46,
      'TGCATT': 330,
      }
    kmer2relative_auweter_kd = relatavize_kds(kmer2auweter_kd)
    kmer2relative_lambert_kd = relatavize_kds(kmer2lambert_kd)
    kmer2relative_kd_to_compare = kmer2relative_lambert_kd
    #kmer2relative_kd_to_compare = kmer2relative_auweter_kd

    main()
    k = 7
    main()

    kmer2relative_kd_to_compare = kmer2relative_auweter_kd
    k = 6
    main()
    k = 7
    main()



