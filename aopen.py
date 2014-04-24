import sys
import gzip
bopen = open


def open(file, mode='r'):

    if file[-3:] == '.gz':
        return gzip.open(file, mode + 'b')
    else:
        return bopen(file, mode)
