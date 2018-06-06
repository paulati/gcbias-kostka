#!/bin/bash

# 1) pseudogenes
wget -O pseudogeneHuman74 http://tables.pseudogene.org/dump.cgi?table=Human74

# 2) Ensemble gene (UCSC browser track hg19.ensGene)
wget -O hg19_ensGene.bed.gz http://genome.ucsc.edu/cgi-bin/hgTables?hgta_ctName=tb_ensGene&hgta_ctDesc=table+browser+query+on+ensGene&hgta_ctVis=full&hgta_ctUrl=&fbUpBases=200&fbQual=exon&fbExonBases=0&fbIntronBases=0&fbDownBases=200&hgta_doGetBed=get+BED

# 3) Refseq gene (UCSC browser track hg19.refGene)
wget -O hg19.refGene.bed.gz http://genome.ucsc.edu/cgi-bin/hgTables?hgta_ctName=tb_refGene&hgta_ctDesc=table+browser+query+on+refGene&hgta_ctVis=full&hgta_ctUrl=&fbUpBases=200&fbQual=exon&fbExonBases=0&fbIntronBases=0&fbDownBases=200&hgta_doGetBed=get+BED

# 4) USCS known gene (UCSC browser track hg19.knownGene)
wget -O hg19.knownGene.bed.gz https://genome.ucsc.edu/cgi-bin/hgTables?hgta_ctName=tb_knownGene&hgta_ctDesc=table+browser+query+on+knownGene&hgta_ctVis=full&hgta_ctUrl=&fbUpBases=200&fbQual=exon&fbExonBases=0&fbIntronBases=0&fbDownBases=200&hgta_doGetBed=get+BED

# 5) RNA gene (UCSC browser track hg19.rnaGene)
#bajan con el script que esta en la carpeta rna

# 6)  mRNA in GenBank (UCSC browser track hg19.mrna)
#wget -O all_mrna.txt.gz http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/all_mrna.txt.gz
wget -O all_mrna.bed.gz http://rohsdb.cmb.usc.edu/GBshape/cgi-bin/hgTables?hgta_ctName=tb_all_mrna&hgta_ctDesc=table+browser+query+on+all_mrna&hgta_ctVis=pack&hgta_ctUrl=&fbUpBases=200&fbQual=exon&fbExonBases=0&fbIntronBases=0&fbDownBases=200&hgta_doGetBed=get+BED

# 7)  Refseq gene, UCSC gene or mRNA in other species 
# (UCSC browser tracks hg18.transMapAlnRefSeq, hg18.transMapAlnUcscGenes and hg18.transMapAlnMRna)

#hg19 TransMapEnsembl
wget -O hg19.transMapEnsembl.bed.gz http://genome.ucsc.edu/cgi-bin/hgTables?hgta_ctName=tb_transMapEnsemblV4&hgta_ctDesc=table+browser+query+on+transMapEnsemblV4&hgta_ctVis=full&hgta_ctUrl=&fbUpBases=200&fbQual=exon&fbExonBases=0&fbIntronBases=0&fbDownBases=200&hgta_doGetBed=get+BED

#hg19 TransMapRNA
wget -O hg19.transMapRnaV4.bed.gz http://genome.ucsc.edu/cgi-bin/hgTables?hgta_ctName=tb_transMapRnaV4&hgta_ctDesc=table+browser+query+on+transMapRnaV4&hgta_ctVis=full&hgta_ctUrl=&fbUpBases=200&fbQual=exon&fbExonBases=0&fbIntronBases=0&fbDownBases=200&hgta_doGetBed=get+BED

#hg19 TransMap RefGene
wget -O hg19.transMaprefSeqV4.bed.gz http://genome.ucsc.edu/cgi-bin/hgTables?hgta_ctName=tb_transMapRefSeqV4&hgta_ctDesc=table+browser+query+on+transMapRefSeqV4&hgta_ctVis=full&hgta_ctUrl=&fbUpBases=200&fbQual=exon&fbExonBases=0&fbIntronBases=0&fbDownBases=200&hgta_doGetBed=get+BED




