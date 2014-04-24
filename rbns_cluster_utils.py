import sys
import getpass
import subprocess
import time
import rbns_utils
import os
import rbns_main
import simplejson
import itertools

def launch_energy_splitter(lib_settings, experiment_settings):
    barcode = lib_settings.get_barcode()
    rdir = experiment_settings.get_rdir()
    reads_file = lib_settings.get_split_reads()
    energy_file = os.path.join(experiment_settings.get_rdir(),
      'structure', '%s.fe.gz' % lib_settings.get_barcode())
    energy_bins = str(experiment_settings.get_property('free_energy_limits')).replace(' ', '')
    cluster_python_script = os.path.abspath(__file__)
    error_dir = experiment_settings.get_property('error_dir')
    print '|'+error_dir+'|'
    out_file = os.path.join(error_dir, 'esplit.%s.out' % barcode)
    err_file = os.path.join(error_dir, 'esplit.%s.err' % barcode)
    command = ('hostname ; python %(cluster_python_script)s '
               'energy_splitter_commandline '
               '%(barcode)s '
               '%(rdir)s '
               '%(reads_file)s '
               '%(energy_file)s '
               '%(energy_bins)s '
               ' %(error_dir)s 1> %(out_file)s 2> %(err_file)s' % locals())
    return launch(command,
      jobname='splitter.'+barcode, ppn='4',
      error_dir=experiment_settings.get_property('error_dir'))


def energy_splitter_commandline(barcode, rdir, reads_file, energy_file, energy_bins, error_dir):
    energy_bins = simplejson.loads(energy_bins)
    energy_bins = map(float, energy_bins)
    energy_splitter(barcode, rdir, reads_file, energy_file, energy_bins, error_dir)


def energy_splitter(barcode, rdir, reads_file, energy_file, energy_bins, error_dir):
    energy_bin_edges = [(left_edge, right_edge) for left_edge, right_edge in rbns_utils.pairwise(energy_bins)]
    out_files = ['%s/structure/%s.%i.%0.1f_to_%0.1f.reads' %
      (rdir, barcode, i, dG1, dG2)
      for i, (dG1, dG2) in enumerate(energy_bin_edges)]
    ofs = map(lambda x: rbns_utils.aopen(x, 'w'), out_files)
    for line, energy in itertools.izip(rbns_utils.aopen(reads_file),
      rbns_utils.iter_RNAfold_output(energy_file)):
        bin_i = rbns_utils.getBinIndex_soft_upper(energy, energy_bins)
        ofs[bin_i].write(line)
    map(lambda of: of.close(), ofs)


def launch_counter(lib_settings, count_type, k, error_dir):
    split_reads = lib_settings.get_split_reads()
    out_pkl = lib_settings.counts_file(count_type, k)
    rbns_utils.make_dir(os.path.dirname(out_pkl))
    cluster_python_script = os.path.abspath(__file__)
    barcode = lib_settings.get_barcode()
    out_file = os.path.join(error_dir, 'count.%s.%s.%i.out' % (barcode, count_type, k))
    err_file = os.path.join(error_dir, 'count.%s.%s.%i.err' % (barcode, count_type, k))
    command = ('hostname ; python %(cluster_python_script)s '
               'counter '
               '%(count_type)s '
               '%(split_reads)s '
               '%(k)i '
               '%(out_pkl)s '
               '1> %(out_file)s '
               '2> %(err_file)s ' % locals())
    conc = lib_settings.get_conc()
    jobname = '%s.%s.%i.%g' % (os.path.basename(split_reads), count_type, k, conc)
    return launch(command, jobname=jobname, ppn='2', error_dir=error_dir)


def launch(command, script_options=False, ppn='1', q='long', jobname='ajob', out_file='', error_dir='~/', mem='1gb'):
    """
    launches the command on the cluster using the given script options
    """
    print 'will launch: '
    print command
    if not script_options:
        script_options = {'nodes': '1', 'ppn': str(ppn),
          'jobname': jobname,
          'queue': q, 'workingdir': os.getcwd()}
    script_options['error_dir'] = os.path.expanduser(error_dir)
    script_options['command'] = command
    script_options['username'] = getpass.getuser()
    script_options['mem'] = mem
    outtext = """#!/bin/bash
    #PBS -l nodes=%(nodes)s:ppn=%(ppn)s
    #PBS -l mem=%(mem)s
    #PBS -m a
    #PBS -e %(error_dir)s/%(jobname)s.zerr
    #PBS -o %(error_dir)s/%(jobname)s.zout
    #PBS -M %(username)s@mit.edu 
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
        print outtext
        print 'Failure to launch'
        raise ValueError('Failure to launch')


def wait_jobs(outstanding_jobs, sleep=5):
    while True:
        outputs = [subprocess.Popen('qstat %i' % job_id,
          shell=True, stdout=subprocess.PIPE,
          stderr=subprocess.PIPE).communicate()
          for job_id in outstanding_jobs]
        completed = map(lambda output: 'Unknown Job' in output[1], outputs)
        if all(completed):
            break
        else:
            print sum(map(lambda output:
              0 if 'Unknown J' in output[1] else 1, outputs)), 'jobs left'
        time.sleep(sleep)


def counter(count_type, split_reads, k, out_pkl):
    k = int(k)
    if count_type == 'naive':
        rbns_main.count_naive(split_reads, k, out_pkl)
    elif count_type == 'stream':
        rbns_main.count_stream(split_reads, k, out_pkl)
    elif count_type == 'presence':
        rbns_main.count_presence(split_reads, k, out_pkl)
    else:
        raise ValueError('Unknown count type: %s ' % count_type)
 

if __name__ == '__main__':
    fxn = sys.argv[1]
    args = ['"' + arg + '"' for arg in sys.argv[2:]]
    python_command = fxn+'('+','.join(args) + ')'
    eval(python_command)
