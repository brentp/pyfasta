==================================================
pyfasta: pythonic access to fasta sequence files.
==================================================


:Author: Brent Pedersen (brentp)
:Email: bpederse@gmail.com
:License: MIT

.. contents ::

Implementation
==============

Requires Python >= 2.5. Stores a flattened version of the fasta file without 
spaces or headers and uses either a mmap of numpy binary format or fseek/fread so the
*sequence data is never read into memory*. Saves a pickle (.gdx) of the start, stop 
(for fseek/mmap) locations of each header in the fasta file for internal use.

Usage
=====
::
  
    >>> from pyfasta import Fasta

    >>> f = Fasta('tests/data/three_chrs.fasta')
    >>> sorted(f.keys())
    ['chr1', 'chr2', 'chr3']

    >>> f['chr1']
    NpyFastaRecord(0..80)


Slicing
-------
::

    >>> f['chr1'][:10]
    'ACTGACTGAC'

    # get the 1st basepair in every codon (it's python yo)
    >>> f['chr1'][::3]
    'AGTCAGTCAGTCAGTCAGTCAGTCAGT'

    # can query by a 'feature' dictionary (note this is one based coordinates)
    >>> f.sequence({'chr': 'chr1', 'start': 2, 'stop': 9})
    'CTGACTGA'

    # same as:
    >>> f['chr1'][1:9]
    'CTGACTGA'

    # use python, zero based coords
    >>> f.sequence({'chr': 'chr1', 'start': 2, 'stop': 9}, one_based=False)
    'TGACTGA'

    # with reverse complement (automatic for - strand)
    >>> f.sequence({'chr': 'chr1', 'start': 2, 'stop': 9, 'strand': '-'})
    'TCAGTCAG'

Key Function
------------
Sometimes your fasta will have a long header like: "AT1G51370.2 | Symbols:  | F-box family protein | chr1:19045615-19046748 FORWARD" when you only want to key off: "AT1G51370.2". In this case, specify the key_fn argument to the constructor:

::

    >>> fkey = Fasta('tests/data/key.fasta', key_fn=lambda key: key.split()[0])
    >>> sorted(fkey.keys())
    ['a', 'b', 'c']

Numpy
=====

The default is to use a memmaped numpy array as the backend. In which case it's possible to
get back an array directly...
::

    >>> f['chr1'].tostring = False
    >>> f['chr1'][:10] # doctest: +NORMALIZE_WHITESPACE
    memmap(['A', 'C', 'T', 'G', 'A', 'C', 'T', 'G', 'A', 'C'], dtype='|S1')

    >>> import numpy as np
    >>> a = np.array(f['chr2'])
    >>> a.shape[0] == len(f['chr2'])
    True

    >>> a[10:14] # doctest: +NORMALIZE_WHITESPACE
    array(['A', 'A', 'A', 'A'], dtype='|S1')

mask a sub-sequence
::

    >>> a[11:13] = np.array('N', dtype='S1')
    >>> a[10:14].tostring()
    'ANNA'


Backends (Record class)
=======================
It's also possible to specify another record class as the underlying work-horse
for slicing and reading. Currently, there's just the default: 

  * NpyFastaRecord which uses numpy memmap
  * FastaRecord, which uses using fseek/fread
  * MemoryRecord which reads everything into memory and must reparse the original
    fasta every time.
  * TCRecord which is identical to NpyFastaRecord except that it saves the index
    in a TokyoCabinet hash database, for cases when there are enough records that
    loading the entire index from a pickle into memory is unwise. (NOTE: that the
    sequence is not loaded into memory in either case).

It's possible to specify the class used with the `record_class` kwarg to the `Fasta`
constructor:
::

    >>> from pyfasta import FastaRecord # default is NpyFastaRecord
    >>> f = Fasta('tests/data/three_chrs.fasta', record_class=FastaRecord)
    >>> f['chr1']
    FastaRecord('tests/data/three_chrs.fasta.flat', 0..80)

other than the repr, it should behave exactly like the Npy record class backend

it's possible to create your own using a sub-class of FastaRecord. see the source 
in pyfasta/records.py for details.

Flattening
==========
In order to efficiently access the sequence content, pyfasta saves a separate, flattened file with all newlines and headers removed from the sequence. In the case of large fasta files, one may not wish to save 2 copies of a 5GG+ file. In that case, it's possible to flatten the file "inplace", keeping all the headers, and retaining the validity of the fasta file -- with the only change being that the new-lines are removed from each sequence. This can be specified via `flatten_inplace` = True
::
    
    >>> import os
    >>> os.unlink('tests/data/three_chrs.fasta.gdx') # cleanup non-inplace idx
    >>> f = Fasta('tests/data/three_chrs.fasta', flatten_inplace=True)
    >>> f['chr1']  # note the difference in the output from above.
    NpyFastaRecord(6..86)

    # sequence from is same as when requested from non-flat file above.
    >>> f['chr1'][1:9]
    'CTGACTGA'

    # the flattened file is kept as a place holder without the sequence data.
    >>> open('tests/data/three_chrs.fasta.flat').read()
    '@flattened@'


Command Line Interface
======================
there's also a command line interface to manipulate / view fasta files.
the `pyfasta` executable is installed via setuptools, running it will show
help text.

split a fasta file into 6 new files of relatively even size:

  $ pyfasta **split** -n 6 original.fasta

split the fasta file into one new file per header with "%(seqid)s" being filled into each filename.:
  
  $ pyfasta **split** --header "%(seqid)s.fasta" original.fasta

create 1 new fasta file with the sequence split into 10K-mers:

  $ pyfasta **split** -n 1 -k 10000 original.fasta

2 new fasta files with the sequence split into 10K-mers with 2K overlap:

  $ pyfasta **split** -n 2 -k 10000 -o 2000 original.fasta


show some info about the file (and show gc content):

  $ pyfasta **info** --gc test/data/three_chrs.fasta


**extract** sequence from the file. use the header flag to make
a new fasta file. the args are a list of sequences to extract.

  $ pyfasta **extract** --header --fasta test/data/three_chrs.fasta seqa seqb seqc

**extract** sequence from a file using a file containing the headers *not* wanted in the new file:

  $ pyfasta extract --header --fasta input.fasta --exclude --file seqids_to_exclude.txt

**extract** sequence from a fasta file with complex keys where we only want to lookup based on the part before the space.

  $ pyfasta extract --header --fasta input.with.keys.fasta --space --file seqids.txt

**flatten** a file inplace, for faster later use by pyfasta, and without creating another copy. (`Flattening`_)

  $ pyfasta flatten input.fasta 

cleanup 
=======
(though for real use these will remain for faster access)
::

    >>> os.unlink('tests/data/three_chrs.fasta.gdx')
    >>> os.unlink('tests/data/three_chrs.fasta.flat')

Testing
=======
there is currently > 99% test coverage for the 2 modules and all included 
record classes. to run the tests:
::

  $ python setup.py nosetests
