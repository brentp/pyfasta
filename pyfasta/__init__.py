import sys
from fasta import Fasta, complement, DuplicateHeaderException
from records import *
from split_fasta import split
import optparse

def main():
    help = """
    available actions:
        `extract`: extract sequences from a fasta file
        `info`: show info about the fasta file and exit.
        `split`: split a large fasta file into separate files
                 and/or into K-mers.
        `flatten`: flatten a fasta file inplace so that later
                   command-line (and programmattic) access via
                   pyfasta will use the inplace flattened version
                   rather than creating another .flat copy of the
                   sequence.

    to view the help for a particular action, use:
        pyfasta [action] --help
    e.g.:
        pyfasta extract --help
    """        
    if len(sys.argv) == 1:
        print help
        sys.exit()

    action = sys.argv[1]

    sglobals = globals()
    if not action in sglobals:
        print "%s not a valid action" % action
        print help
        sys.exit()
    
    globals()[action](sys.argv[2:])

def info(args):
    """
    >>> info(['tests/data/three_chrs.fasta'])
    <BLANKLINE>
    tests/data/three_chrs.fasta
    ===========================
    >chr3 length:3600 
    >chr2 length:80 
    >chr1 length:80 
    <BLANKLINE>
    3760 basepairs in 3 sequences
    """
    parser = optparse.OptionParser("""\
   print headers and lengths of the given fasta file in order of length. e.g.:
        pyfasta info --gc some.fasta""")

    parser.add_option("-n", "--n", type="int", dest="nseqs", 
                      help="max number of records to print. use -1 for all",
                      default=20)
    parser.add_option("--gc", dest="gc", help="show gc content",
                      action="store_true", default=False)
    options, fastas = parser.parse_args(args)
    if not (fastas):
        sys.exit(parser.print_help())
    import operator

    for fasta in fastas:
        f = Fasta(fasta)
        info = [(k, len(seq)) for k, seq in f.iteritems()]

        total_len = sum(l for k, l in info)
        nseqs = len(f)
        if options.nseqs > -1:
            info = sorted(info,  key=operator.itemgetter(1), reverse=True)
            info = info[:options.nseqs]
        else:
            info.sort()

        print "\n" + fasta
        print "=" * len(fasta)
        for k, l in info:
            gc = ""
            if options.gc:
                seq = str(f[k]).upper()
                g = seq.count('G')
                c = seq.count('C')
                gc = 100.0 * (g + c) / float(l)
                gc = "gc:%.2f%%" % gc
            print (">%s length:%i " % (k, l)) + gc

        if total_len > 1000000:
            total_len = "%.3fM" % (total_len / 1000000.)
        print
        print "%s basepairs in %i sequences" % (total_len, nseqs)

def flatten(args):
    """
    >>> flatten(['tests/data/three_chrs.fasta'])
    """
    parser = optparse.OptionParser("""flatten a fasta file *inplace* so all later access by pyfasta will use that flattend (but still viable) fasta file""")
    _, fasta = parser.parse_args(args)
    for fa in fasta:
        f = Fasta(fa, flatten_inplace=True)

def extract(args):
    """
    >>> extract(['--fasta', 'tests/data/three_chrs.fasta', 'chr2'])
    TAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAT
    """

    parser = optparse.OptionParser("""extract some sequences from a fasta file. e.g.:
               pyfasta extract --fasta some.fasta --header at2g26540 at3g45640""")
    parser.add_option("--fasta", dest="fasta", help="path to the fasta file")
    parser.add_option("--header", dest="header", help="include headers", action="store_true", default=False)
    parser.add_option("--exclude", dest="exclude", help="extract all sequences EXCEPT those listed", action="store_true", default=False)
    parser.add_option("--file", dest="file", help=\
                      "if this flag is used, the sequences to extract" \
                      " are read from the file specified in args"
                      , action="store_true", default=False)
    parser.add_option("--space", dest="space", action="store_true", help=\
                      "use the fasta identifier only up to the space as the key",
                      default=False)
    options, seqs = parser.parse_args(args)
    if not (options.fasta and len(seqs)):
        sys.exit(parser.print_help())

    key_fn = (lambda k: k.split()[0]) if options.space else None
    f = Fasta(options.fasta, key_fn=key_fn)
    if options.file:
        seqs = (x.strip() for x in open(seqs[0]))
    if options.exclude:
        seqs = sorted(frozenset(f.iterkeys()).difference(seqs))

    for seqname in seqs:
        seq = f[seqname]
        if options.header:
            print ">%s" % seqname
        print seq


if __name__ == "__main__":
    main()
