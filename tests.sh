A=`python -c "from pyfasta import Fasta; f = Fasta('tests/data/three_chrs.fasta', flatten_inplace=False); print f['chr1'][:3]"`
B=`python -c "from pyfasta import Fasta; f = Fasta('tests/data/three_chrs.fasta', flatten_inplace=True); print f['chr1'][:3]"`

if [ "$A" != "$B" ]
then
    echo "$A=$B"
fi

F=`cat tests/data/three_chrs.fasta.flat`
if [ "$F" != "@flattened@" ]
then
    echo "$F=@flattened@"
fi


rm -f tests/data/three_chrs.fasta.{flat,gdx}
C=`python -c "from pyfasta import Fasta; f = Fasta('tests/data/three_chrs.fasta', flatten_inplace=True); print f['chr1'][:3]"`
D=`python -c "from pyfasta import Fasta; f = Fasta('tests/data/three_chrs.fasta', flatten_inplace=False); print f['chr1'][:3]"`

if [ "$C" != "$D" ]
then
    echo "$C=$D"
fi

F=`cat tests/data/three_chrs.fasta.flat`
if [ "$F" != "@flattened@" ]
then
    echo "$F=@flattened@"
fi


rm -f tests/data/three_chrs.fasta.{flat,gdx}
I=`python pyfasta/__init__.py info tests/data/three_chrs.fasta`
rm -f tests/data/three_chrs.fasta.{flat,gdx}

# now make sure it works with inplace, and doesnt update.flat
C=`python -c "from pyfasta import Fasta; f = Fasta('tests/data/three_chrs.fasta', flatten_inplace=True); print f['chr1'][:3]"`
J=`python pyfasta/__init__.py info tests/data/three_chrs.fasta`
if [ "$I" != "$J" ]
then
    echo "BADDD"
fi

F=`cat tests/data/three_chrs.fasta.flat`
if [ "$F" != "@flattened@" ]
then
    echo "BADD"
fi

rm -f tests/data/three_chrs.fasta.{flat,gdx}
