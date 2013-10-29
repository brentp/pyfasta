from __future__ import print_function
import string
import os.path
from collections import Mapping
import sys
import numpy as np

from records import NpyFastaRecord

# string.maketrans is bytes.maketrans in Python 3, but
# we want to deal with strings instead of bytes
try:
    # 2.x
    maketrans = string.maketrans
except AttributeError:
    # 3.x
    maketrans = str.maketrans

_complement = maketrans('ATCGatcgNnXx', 'TAGCtagcNnXx')
# Python 2: string.maketrans returns a bytes object of length 256,
#   that is used as a lookup table to translate bytes to other bytes.
# Python 3: str.maketrans returns a dict mapping Unicode code points
#   to other Unicode code points. Can't use a fully-allocated lookup
#   table since it would have to be of size `sys.maxunicode`, which
#   is equal to 1114111 on wide builds of <= 3.2 and all builds of
#   Python >= 3.3.
# In Python 2, it's safe to use a unicode object as the translation
# table; this causes str.translate to return a unicode object instead
# of a str. This is safe as long as the string that you're translating
# can be decoded as ASCII, and will fail with a UnicodeDecodeError
# otherwise.
if sys.version_info[0] < 3:
    _complement = _complement.decode('latin-1')
complement = lambda s: s.translate(_complement)

class FastaNotFound(Exception): pass

class DuplicateHeaderException(Exception):
    def __init__(self, header):
        Exception.__init__(self, 'headers must be unique: %s is duplicated' % header)

class Fasta(Mapping):
    def __init__(self, fasta_name, record_class=NpyFastaRecord,
                flatten_inplace=False, key_fn=None):
        """
            >>> from pyfasta import Fasta, FastaRecord

            >>> f = Fasta('tests/data/three_chrs.fasta',
            ...                          record_class=FastaRecord)
            >>> sorted(f.keys())
            ['chr1', 'chr2', 'chr3']

        slicing returns an object.
            >>> f['chr1']
            FastaRecord('tests/data/three_chrs.fasta.flat', 0..80)

        extract sequence with normal python syntax
            >>> print(f['chr1'][:10])
            ACTGACTGAC

        take the first basepair in each codon...
            >>> print(f['chr1'][0:10:3])
            AGTC

        """
        if not os.path.exists(fasta_name):
            raise FastaNotFound('"' + fasta_name + '"')
        self.fasta_name = fasta_name
        self.record_class = record_class
        self.key_fn = key_fn
        self.index, self.prepared = self.record_class.prepare(self,
                                              self.gen_seqs_with_headers(key_fn),
                                              flatten_inplace)

        self.chr = {}

    @classmethod
    def as_kmers(klass, seq, k, overlap=0):
        kmax = len(seq)
        assert overlap < k, ('overlap must be < kmer length')
        i = 0
        while i < kmax:
            yield i, seq[i:i + k]
            i += k - overlap

    def gen_seqs_with_headers(self, key_fn=None):
        """remove all newlines from the sequence in a fasta file
        and generate starts, stops to be used by the record class"""
        fh = open(self.fasta_name, 'r')
        # do the flattening (remove newlines)
        # check of unique-ness of headers.
        seen_headers = set()
        header = None
        seqs = None
        for line in fh:
            line = line.rstrip()
            if not line: continue
            if line[0] == ">":
                if seqs is not None:
                    if header in seen_headers:
                        raise DuplicateHeaderException(header)
                    seen_headers.add(header)
                    yield header, "".join(seqs)

                header = line[1:].strip()
                if key_fn is not None:
                    header = key_fn(header)
                seqs = []
            else:
                seqs.append(line)

        if seqs:
            if header in seen_headers:
                raise DuplicateHeaderException(header)
            yield header, "".join(seqs)
        fh.close()

    def __len__(self):
        # might not work for all backends?
        return len(self.index)

    def __iter__(self):
        return iter(self.index)

    def __getitem__(self, i):
        if i in self.chr:
            return self.chr[i]
        c = self.index[i]
        self.chr[i] = self.record_class(self.prepared, c[0], c[1])
        return self.chr[i]

    def sequence(self, f, asstring=True, auto_rc=True
            , exon_keys=None, one_based=True):
        """
        take a feature and use the start/stop or exon_keys to return
        the sequence from the assocatied fasta file:
        f: a feature
        asstring: if true, return the sequence as a string
                : if false, return as a numpy array
        auto_rc : if True and the strand of the feature == -1, return
                  the reverse complement of the sequence
        one_based: if true, query is using 1 based closed intervals, if false
                    semi-open zero based intervals

            >>> from pyfasta import Fasta
            >>> f = Fasta('tests/data/three_chrs.fasta')
            >>> print(f.sequence({'start':1, 'stop':2, 'strand':1, 'chr': 'chr1'}))
            AC

            >>> print(f.sequence({'start':1, 'stop':2, 'strand': -1, 'chr': 'chr1'}))
            GT

            >>> sorted(f.index.items())
            [('chr1', (0, 80)), ('chr2', (80, 160)), ('chr3', (160, 3760))]

        NOTE: these 2 are reverse-complement-ary because of strand
        #>>> f.sequence({'start':10, 'stop':12, 'strand': -1, 'chr': 'chr1'})
            'CAG'
            >>> print(f.sequence({'start':10, 'stop':12, 'strand': 1, 'chr': 'chr1'}))
            CTG


            >>> print(f.sequence({'start':10, 'stop':12, 'strand': -1, 'chr': 'chr3'}))
            TGC
            >>> print(f.sequence({'start':10, 'stop':12, 'strand': 1, 'chr': 'chr3'}))
            GCA

            >>> print(f['chr3'][:][-10:])
            CGCACGCTAC


        a feature can have exons:
            >>> feat = dict(start=9, stop=19, strand=1, chr='chr1'
            ...    , exons=[(9,11), (13, 15), (17, 19)])

        by default, it just returns the full sequence between start
        and stop.
            >>> print(f.sequence(feat))
            ACTGACTGACT

        but if exon_keys is set to an iterable, it will search for
        those keys and will use the first to create a sequence and
        return the concatenated result.
            >>> print(f.sequence(feat, exon_keys=('rnas', 'exons')))
            ACTACTACT

        Note that sequence is 2 characters shorter than the entire
        feature, to account for the introns at base-pairs 12 and 16.

        Also note, it looks for an item with key of 'rnas', and didn't
        fine one, so it continued on to 'exons'. If it doesn't find
        any of the exon keys, it will fall back on the start, stop of
        the feature:
            >>> print(f.sequence(feat, exon_keys=('fake', 'also_fake')))
            ACTGACTGACT
        """
        assert 'chr' in f and f['chr'] in self, (f, f['chr'], self.keys())
        fasta    = self[f['chr']]
        sequence = None
        if not exon_keys is None:
            sequence = self._seq_from_keys(f, fasta, exon_keys, one_based=one_based)

        if sequence is None:
            start = f['start'] - int(one_based)
            sequence = fasta[start: f['stop']]

        if auto_rc and f.get('strand') in (-1, '-1', '-'):
            sequence = complement(sequence)[::-1]

        if asstring: return sequence
        return np.array(sequence, dtype='c')

    def _seq_from_keys(self, f, fasta, exon_keys, base='locations', one_based=True):
        """Internal:
        f: a feature dict
        fasta: a Fasta object
        exon_keys: an iterable of keys, to look for start/stop
                   arrays to get sequence.
        base: if base ('locations') exists, look there fore the
        exon_keys, not in the base of the object:
            {'name': 'OS11G42690', 'stop': 25210251, 'locations':
            {'CDS': [(25210018, 25210251)]}, 'start': 25210018, 'chr':
            '11', 'strand': -1} set(['TRNA', 'CDS'])
        """
        fbase = f.get(base, f)
        for ek in exon_keys:
            if not ek in fbase: continue
            locs = fbase[ek]
            seq = ""
            for start, stop in locs:
                seq += fasta[start - int(one_based):stop]
            return seq
        return None
