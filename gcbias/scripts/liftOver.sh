#!/bin/sh

#./liftOver ~/Documentos/r/2017/lara/data/exclusion/trans/split_hg19_transMapEnsemblV4/chr1.hg19_transMapEnsemblV4.bed ~/Documentos/#r/2017/lara/data/exclusion/trans/hg19ToHg18.over.chain ~/Documentos/r/2017/lara/data/exclusion/trans/#split_hg19_transMapEnsemblV4/chr1.hg18_transMapEnsemblV4.bed unlifted.bed 

ls="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22"

cd /paula/2018/gcbias/scripts

for chr in $ls
do

sourcefolder='/paula/2018/gcbias/data/preparation/hg19/transMapEnsembl'
targetfolder='/paula/2018/gcbias/data/preparation/hg18/transMapEnsembl/lift'
unliftedfolder='/paula/2018/gcbias/data/preparation/hg18/transMapEnsembl/unlifted'
sourcefile='chr'$chr'_hg19.transMapEnsembl.bed'
targetfile='chr'$chr'_hg18.transMapEnsembl.bed'

unliftedfile='chr'$chr'_hg19.unlifted.bed' 

echo $chr
echo $sourcefile
echo $targetfile
echo $unliftedfile


~/Documentos/apps/liftOver $sourcefolder/$sourcefile ./hg19ToHg18.over.chain $targetfolder/$targetfile $unliftedfolder/$unliftedfile

done


for chr in $ls
do

sourcefolder='/paula/2018/gcbias/data/preparation/hg19/transMapRefSeq'
targetfolder='/paula/2018/gcbias/data/preparation/hg18/transMapRefSeq/lift'
unliftedfolder='/paula/2018/gcbias/data/preparation/hg18/transMapRefSeq/unlifted'
sourcefile='chr'$chr'_hg19.transMaprefSeqV4.bed'
targetfile='chr'$chr'_hg18.transMaprefSeqV4.bed'
unliftedfile='chr'$chr'.hg19_unlifted.bed' 

echo $chr
echo $sourcefile
echo $targetfile
echo $unliftedfile


~/Documentos/apps/liftOver $sourcefolder/$sourcefile ./hg19ToHg18.over.chain $targetfolder/$targetfile $unliftedfolder/$unliftedfile

done



for chr in $ls
do

sourcefolder='/paula/2018/gcbias/data/preparation/hg19/transMapRNA'
targetfolder='/paula/2018/gcbias/data/preparation/hg18/transMapRNA/lift'
unliftedfolder='/paula/2018/gcbias/data/preparation/hg18/transMapRNA/unlifted'
sourcefile='chr'$chr'_hg19.transMapRnaV4.bed'
targetfile='chr'$chr'_hg18.transMapRnaV4.bed'
unliftedfile='chr'$chr'_hg19.unlifted.bed' 

echo $chr
echo $sourcefile
echo $targetfile
echo $unliftedfile


~/Documentos/apps/liftOver $sourcefolder/$sourcefile ./hg19ToHg18.over.chain $targetfolder/$targetfile $unliftedfolder/$unliftedfile

done
