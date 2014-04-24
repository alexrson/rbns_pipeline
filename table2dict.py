import sys

def two_column_table2dict(inFile):
    d = {}
    for line in open(inFile):
        ll = line.split('\t')
        d[ll[0]] = ll[1].strip()
    return d


def table2idList(inFile):
    t2d = table2dict(inFile)
    return list(t2d.keys())


def table2dict(inFile):
    t2d = {}
    firstLine = True
    for line in open(inFile):
        if firstLine:
            keys = line.strip().split('\t')
            firstLine = False
            continue
        ll = line.strip().split('\t')
        if line[0] == '\t':
            continue
        id = ll[0]
        t2d[id] = {}

        for i in range(1, len(keys)):
            try:
                t2d[id][keys[i]] = ll[i]
            except:
                print 'table2dict error',line,ll,inFile
                sys.exit()
    return t2d


def readList(inFile):
    o = set()
    with open(inFile) as f:
        for line in f.readlines():
            o.add(line.strip())
    return o


def miso_comparisions_table2dict(inFile):
    t2d = table2dict(inFile)
    events = t2d.keys()
    for event in events:
        for attri in ['sample1_posterior_mean','sample1_ci_low','sample1_ci_high','sample2_posterior_mean','sample2_ci_low','sample2_ci_high','diff','bayes_factor']:
            try:
                t2d[event][attri] = float(t2d[event][attri])
            except:
                print 'failed for',attri, inFile
                print t2d[event]
                print 
                print t2d[event][attri]
                del t2d[event]
                break
    return t2d


def cuff_comparisions_table2dict(inFile):
    t2d = table2dict(inFile)
    for gene in t2d:
        for v in ['FPKM','coverage','FPKM_conf_lo','FPKM_conf_hi']:
            if t2d[gene][v] <> '-':
                t2d[gene][v] = float(t2d[gene][v])
    return t2d


def gene_to_fpkm_cuff(cuff_file):
    t2d = cuff_comparisions_table2dict(cuff_file)
    o = {}
    for gene, d in t2d.items():
        o[gene] = d['FPKM']
    return o


def loess_comparisions_table2dict(inFile):
    t2d = table2dict(inFile)
    for gene in t2d:
        for v in ['new x','new y']:
            t2d[gene][v.replace('new ','')] = float(t2d[gene][v])
    return t2d


def table2intdict(inFile):
    from collections import defaultdict
    t2d = defaultdict(int)
    firstLine = True
    for line in open(inFile):
        if firstLine:
            keys = line.strip().split('\t')
            firstLine = False
            continue
        ll = line.strip().split('\t')
        id = ll[0]
        t2d[id] = {}
        for i in range(1, len(keys)):
            t2d[id] = int(ll[i])
    return t2d


def table2dict_all_floats(inFile):
    t2d = table2dict(inFile)
    t2d_f = dict([(k, dict()) for k in t2d])
    for row_key, row_d in t2d.items():
        for col_key, cell in row_d.items():
            t2d_f[row_key][col_key] = float(t2d[row_key][col_key])
    return t2d_f
