from pyfasta import Fasta
import operator
import collections
import string
import sys
import optparse
from cStringIO import StringIO


def newnames(oldname, n, kmers=None, overlap=None, header=None):
    """
    >>> newnames('some.fasta', 1)
    ['some.split.fasta']

    >>> newnames('some.fasta', 2)
    ['some.a.fasta', 'some.b.fasta']

    >>> newnames('some', 2)
    ['some.a', 'some.b']

    >>> newnames('some.fasta', 2, kmers=1000)
    ['some.a.1Kmer.fasta', 'some.b.1Kmer.fasta']

    >>> newnames('some.fasta', 2, kmers=10000, overlap=2000)
    ['some.a.10Kmer.2Koverlap.fasta', 'some.b.10Kmer.2Koverlap.fasta']

    >>> newnames('some.fasta', 1, kmers=100000, overlap=2000)
    ['some.split.100Kmer.2Koverlap.fasta']

    """
    if kmers and kmers % 1000 == 0: kmers = "%iK" % (kmers/1000)
    if overlap and overlap % 1000 == 0: overlap = "%iK" % (overlap/1000)

    p = oldname.rfind("fa")
    kstr = kmers is not None and ("%smer." % kmers) or ""
    ostr = overlap is not None and ("%soverlap." % overlap) or ""
    if p != -1:
        pattern = oldname[:p] + "%s." + kstr + ostr + oldname[p:]
    else:
        pattern = oldname + kstr + ostr + ".%s"

    
    

    if n == 1:
        names = [pattern % "split"]
    else:
        width = len(str(n))
        names = [pattern % str(i).rjust(width, '0') for i in range(n)]
    print >>sys.stderr, "creating new files:"
    print >>sys.stderr, "\n".join(names)
    return names


def print_to_fh(fh, fasta, lens, seqinfo):
    key, seqlen = seqinfo
    lens[fh.name] += seqlen
    f = fasta
    assert len(str(f[key])) == seqlen, (key, seqlen, len(str(f[key])))
    print >>fh, ">%s" % key
    print >>fh, str(f[key])


def format_kmer(seqid, start):
    """
    prints out a header with 1-based indexing.

    >>> format_kmer('chr3', 1000)
    'chr3_1001'
    """
    return "%s_%i" % (seqid, start + 1)

def split(args):
    parser = optparse.OptionParser("""\
   split a fasta file into separated files.
        pyfasta split -n 6 [-k 5000 ] some.fasta
    the output will be some.1.fasta, some.2.fasta ... some.6.fasta
    the sizes will be as even as reasonable.
   """)
    parser.add_option("--header", dest="header", metavar="FILENAME_FMT",
       help="""this overrides all other options. if specified, it will
               split the file into a separate file for each header. it
               will be a template specifying the file name for each new file.
               e.g.:    "%(fasta)s.%(seqid)s.fasta"
               where 'fasta' is the basename of the input fasta file and seqid
               is the header of each entry in the fasta file.""" ,default=None)

    parser.add_option("-n", "--n", type="int", dest="nsplits", 
                            help="number of new files to create")
    parser.add_option("-o", "--overlap", type="int", dest="overlap", 
                            help="overlap in basepairs", default=0)
    parser.add_option("-k", "--kmers", type="int", dest="kmers", default=-1,
                     help="""\
    split big files into pieces of this size in basepairs. default
    default of -1 means do not split the sequence up into k-mers, just
    split based on the headers. a reasonable value would be 10Kbp""")
    options, fasta = parser.parse_args(args)
    if not (fasta and (options.nsplits or options.header)):
        sys.exit(parser.print_help())

    if isinstance(fasta, (tuple, list)):
        assert len(fasta) == 1, fasta
        fasta = fasta[0]

    kmer = options.kmers if options.kmers != -1 else None
    overlap = options.overlap if options.overlap != 0 else None
    f = Fasta(fasta)
    if options.header:
        names = dict([(seqid, options.header % \
                      dict(fasta=f.fasta_name, seqid=seqid)) \
                                       for seqid in f.keys()])
        """
        if len(names) > 0:
            assert names[0][1] != names[1][1], ("problem with header format", options.header)
        fhs = dict([(seqid, open(fn, 'wb')) for seqid, fn in names[:200]])
        fhs.extend([(seqid, StringIO(), fn) for seqid, fn in names[200:]])
        """
        return with_header_names(f, names)
    else:
        names = newnames(fasta, options.nsplits, kmers=kmer, overlap=overlap, 
                     header=options.header)

        #fhs = [open(n, 'wb') for n in names]
    if options.kmers == -1:
        return without_kmers(f, names)
    else: 
        return with_kmers(f, names, options.kmers, options.overlap)

def with_header_names(f, names):
    """
    split the fasta into the files in fhs by headers.
    """
    for seqid, name in names.iteritems():
        fh = open(name, 'wb')
        print >>fh, ">%s" % seqid
        print >>fh, str(f[seqid])
        fh.close()

def with_kmers(f, names, k, overlap):
    """
    split the sequences in Fasta object `f` into pieces of length `k` 
    with the given `overlap` the results are written to the array of files
    `fhs`
    """
    fhs = [open(name, 'wb') for name in names]
    i = 0
    for seqid in f.keys():
        seq = f[seqid]
        for (start0, subseq) in Fasta.as_kmers(seq, k, overlap=overlap):

            fh = fhs[i % len(fhs)]
            print >>fh, ">%s" % format_kmer(seqid, start0)
            print >>fh, subseq
            i += 1

def without_kmers(f, names):
    """
    long crappy function that does not solve the bin-packing problem.
    but attempts to distribute the sequences in Fasta object `f` evenly
    among the file handles in fhs.
    """
    fhs = [open(name, 'wb') for name in names]
    name2fh = dict([(fh.name, fh) for fh in fhs])
    items = sorted([(key, len(f[key])) for key in f.keys()], 
                   key=operator.itemgetter(1))

    l1 = len(items) - 1
    l0 = 0
    lens = collections.defaultdict(int)

    n_added = 0
    while l0 < l1:
        fh = fhs[n_added % len(fhs)]
        added = False
        if n_added >= len(fhs):

            while l1 > l0:
                lmin = min(lens.itervalues())
                lmax = max(lens.itervalues())
                if float(lmin) / lmax < 0.80:
                    # it's way off, so add a large (l1)
                    name = find_name_from_len(lmin, lens)
                    fh = name2fh[name]
                    print_to_fh(fh, f, lens, items[l1])
                    l1 -= 1
                    added = True
                    n_added += 1

                elif float(lmin) / lmax < 0.94:
                    # it's just a little off so add a small (l0)
                    name = find_name_from_len(lmin, lens)
                    fh = name2fh[name]
                    print_to_fh(fh, f, lens, items[l0])
                    l0 += 1
                    added = True
                    n_added += 1
                else:
                    break

                if not added:
                    break
                # TODO: test this on glycine
                #added = False
        if added:
            continue

        print_to_fh(fh, f, lens, items[l1])
        l1 -= 1
        n_added += 1

    if l0 == l1:
        fh = fhs[l0 % len(fhs)]
        print_to_fh(fh, f, lens, items[l0])


def find_name_from_len(lmin, lens):
    """
    reverse lookup, get name from dict
    >>> lens = {'chr1': 4, 'chr2': 5, 'chr3': 4}
    >>> find_name_from_len(5, lens)
    'chr2'
    """
    for fname, l in lens.iteritems():
        if l == lmin: 
            return fname
    raise Exception('name not found')
