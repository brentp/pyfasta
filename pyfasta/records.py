import cPickle
import numpy as np
import sys
import os

__all__ = ['FastaRecord', 'NpyFastaRecord', 'MemoryRecord']

MAGIC = "@flattened@"

def is_up_to_date(a, b):
    return os.path.exists(a) and os.stat(a).st_mtime >= os.stat(b).st_mtime


def ext_is_flat(ext):
    fh = open(ext)
    t = fh.read(len(MAGIC))
    fh.close()
    return MAGIC == t

class FastaRecord(object):
    __slots__ = ('fh', 'start', 'stop')
    ext = ".flat"
    idx = ".gdx"

    @classmethod
    def is_current(klass, fasta_name):
        utd = is_up_to_date(fasta_name + klass.idx, fasta_name) 
        if not utd: return False
        return is_up_to_date(fasta_name + klass.ext, fasta_name)

    def __init__(self, fh, start, stop):

        self.fh      = fh
        self.stop    = stop
        self.start   = start

    def __len__(self):
        return self.stop - self.start

    @classmethod
    def prepare(klass, fasta_obj, seqinfo_generator, flatten_inplace):
        """
        returns the __getitem__'able index. and the thing to get the seqs from.
        """
        f = fasta_obj.fasta_name
        if klass.is_current(f):
            fh = open(f + klass.idx, 'rb')
            idx = cPickle.load(fh)
            fh.close()
            if flatten_inplace or ext_is_flat(f + klass.ext): flat = klass.modify_flat(f)
            else: flat = klass.modify_flat(f + klass.ext)
            if flatten_inplace and not ext_is_flat(f + klass.ext):
                del flat
            else:
                return idx, flat

        idx = {}
        flatfh = open(f + klass.ext, 'wb')
        for i, (seqid, seq) in enumerate(seqinfo_generator):
            if flatten_inplace:
                if i == 0:
                    flatfh.write('>%s\n' % seqid)
                else:
                    flatfh.write('\n>%s\n' % seqid)
            start = flatfh.tell()
            flatfh.write(seq)
            stop = flatfh.tell() 
            idx[seqid] = (start, stop)
        flatfh.close()
            
        if flatten_inplace:
            klass.copy_inplace(flatfh.name, f)
            fh = open(f + klass.idx, 'wb')
            cPickle.dump(idx, fh, -1)
            fh.close()
            return idx, klass.modify_flat(f)

        fh = open(f + klass.idx, 'wb')
        cPickle.dump(idx, fh, -1)
        fh.close()
        return idx, klass.modify_flat(f + klass.ext)

    @classmethod
    def copy_inplace(klass, flat_name, fasta_name):
        """
        so they requested flatten_inplace, so we
        overwrite the .fasta file with the .fasta.flat
        file and save something in the .flat as a place-holder.
        """
        os.rename(flat_name, fasta_name)
        # still need the flattend file to show
        # it's current.
        flatfh = open(fasta_name + klass.ext, 'wb')
        flatfh.write(MAGIC)
        flatfh.close()


    
    @classmethod
    def modify_flat(klass, flat_file):
        return open(flat_file, 'rb')
    
    def _adjust_slice(self, islice):
        l = len(self)
        if not islice.start is None and islice.start < 0:
            istart = self.stop + islice.start
        else:
            if islice.start is None:
                istart = self.start
            else:
                istart = self.start + islice.start


        if not islice.stop is None and islice.stop < 0:
            istop = self.stop + islice.stop
        else:
            istop = islice.stop is None and self.stop or (self.start + islice.stop)

        # this will give empty string
        if istart > self.stop: return self.stop, self.stop 
        if istart < self.start: istart = self.start

        if istop  < self.start: istop = self.start
        elif istop > self.stop: istop = self.stop

        return istart, istop

    def __getitem__(self, islice):
        fh = self.fh
        fh.seek(self.start)


        if isinstance(islice, (int, long)):
            if islice < 0:
                if -islice > self.stop - self.start:
                    raise IndexError
                fh.seek(self.stop + islice)
            else:
                fh.seek(self.start + islice)
            return fh.read(1)

        # [:]
        if islice.start in (0, None) and islice.stop in (None, sys.maxint):
            if islice.step in (1, None):
                return fh.read(self.stop - self.start)
            return fh.read(self.stop - self.start)[::islice.step]
        
        istart, istop = self._adjust_slice(islice)
        if istart is None: return ""
        l = istop - istart
        if l == 0: return ""


        fh.seek(istart)

        if islice.step in (1, None):
            return fh.read(l)

        return fh.read(l)[::islice.step]


    def __str__(self):
        return self[:]

    def __repr__(self):
        return "%s('%s', %i..%i)" % (self.__class__.__name__, self.fh.name,
                                   self.start, self.stop)

    @property
    def __array_interface__(self):
        return {
            'shape': (len(self), ),
            'typestr': '|S1',
            'version': 3,
            'data': buffer(self)
        }


class NpyFastaRecord(FastaRecord):
    __slots__ = ('start', 'stop', 'mm', 'as_string')

    def __init__(self, mm, start, stop, as_string=True):
        self.mm = mm
        self.start = start
        self.stop = stop
        self.as_string = as_string

    def __repr__(self):
        return "%s(%i..%i)" % (self.__class__.__name__,
                                   self.start, self.stop)

    @classmethod
    def modify_flat(klass, flat_file):
        mm = np.memmap(flat_file, dtype="S1", mode="r")
        return mm

    def getdata(self, islice):
        if isinstance(islice, (int, long)):
            if islice >= 0:
                islice += self.start
            else:
                islice += self.stop
                if islice < 0: raise IndexError
            return self.mm[islice]

        start, stop = self._adjust_slice(islice)
        return self.mm[start:stop:islice.step]

    def __getitem__(self, islice):
        d = self.getdata(islice)
        return d.tostring() if self.as_string else d

    @property
    def __array_interface__(self):
        return {
            'shape': (len(self), ),
            'typestr': '|S1',
            'version': 3,
            'data': self[:]
        }


class MemoryRecord(FastaRecord):
    """
    dont write anything to disk, just read the whole thing
    into memory
    """
    @classmethod
    def prepare(klass, fasta_obj, seqinfo_generator, flatten_inplace=False):
        f = fasta_obj.fasta_name
        seqs = {}
        idx = {}
        for seqid, seq in seqinfo_generator:
            seqs[seqid] = (seq, None)

        return seqs, seqs

    def __init__(self, _, seq, _none):
        self.seq = seq

    def __getitem__(self, slice):
        return self.seq.__getitem__(slice)

    def __len__(self):
        return len(self.seq)



try:
    import tc
    class HDB(tc.HDB):
        def __getitem__(self, k):
            return cPickle.loads(tc.HDB.get(self, k))
        def __setitem__(self, k, v):
            tc.HDB.put(self, k, cPickle.dumps(v, -1))
        def __del__(self):
            tc.HDB.close(self)

    class TCRecord(NpyFastaRecord):
        idx = ".tct"

        @classmethod
        def prepare(klass, fasta_obj, seqinfo_generator, flatten_inplace):
            f = fasta_obj.fasta_name
            if klass.is_current(f):
                idx = HDB()
                idx.open(f + klass.idx, tc.HDBOREADER)
                if flatten_inplace or ext_is_flat(f + klass.ext): flat = klass.modify_flat(f)
                else: flat = klass.modify_flat(f + klass.ext)
                return idx, flat


            db = HDB(f + klass.idx, tc.HDBOWRITER | tc.HDBOCREAT)
            flatfh = open(f + klass.ext, 'wb')
            for i, (seqid, seq) in enumerate(seqinfo_generator):
                if flatten_inplace:
                    if i == 0:
                        flatfh.write('>%s\n' % seqid)
                    else:
                        flatfh.write('\n>%s\n' % seqid)
                start = flatfh.tell()
                flatfh.write(seq)
                stop = flatfh.tell() 
                db[seqid] = (start, stop)

            db.sync()
            flatfh.close()

            if flatten_inplace:
                klass.copy_inplace(flatfh.name, f)
                return db, klass.modify_flat(f)
            return db, klass.modify_flat(f + klass.ext)


    __all__.append('TCRecord')
except ImportError:
    pass
