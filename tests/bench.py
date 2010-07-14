import sys
import os
sys.path.insert(0, os.path.abspath("."))
import pyfasta
from pyfasta import Fasta

import time
import random
random.seed(1234)

SEQLEN = 100000
def make_long_fasta(filename="t.fasta", nrecs=2500, seqlen=SEQLEN):
    fh = open(filename, 'w')
        #0123456789"
    s = "ACTGACTGAC"
    for i in range(nrecs):
        print >>fh, ">header%i" % i
        print >>fh, (s * (seqlen/10))

    fh.close()
    return filename

def read(f, nreads=40000, seqlen=SEQLEN):

    for k in f.keys()[:10]:
        for i in range(nreads):
            start = random.randint(0, seqlen)
            end = min(seqlen, start + random.randint(1000, 2000))
            str(f[k][start:end])
        


def main():
    fa = make_long_fasta()

    t = time.time()
    f = Fasta(fa)
    print "flatten:", time.time() - t 

    
    t = time.time()
    read(f)
    print "read:", time.time() - t 

     


if __name__ == "__main__":
    main()
