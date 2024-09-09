<p align="center">
  <a>
    <img height="150" src="logo/BAdabouM_logo.png">
  </a>
</p>

BAdabouM is a free tool to detect Structural Varations (SVs) in Whole Genome reSequencing data.
BAdabouM is able to genotype insertions (INS), deletions (DEL), inversions (INV), copy number variations (CNVs) and inter and intra chromosomal rearrangements (CTX and ITX).
BAdabouM uses paired-ends, split-reads and aberrant depth of coverage to call SVs.

Version 2 | last update : 17/07/2018


### Instalating BAdabouM
------------------------

The easiest way to get BAdabouM is to build it from sources.

```
git clone https://github.com/cumtr/BAdabouM.git

cd BAdabouM/

make
```
In case you face some issue with the compilation, you can try to independently install (or load depending on your enviroment) samtools and then run :

```
cd src/

make
```


### Running BAdabouM
--------------------

BAdabouM only needs a sorted and indexed bam file for every sample to run. 

`./BAdabouM/src/BAdabouM Sample_1.sorted.bam > Sample_1.BAdabouM_out`


To see BAdabouM's options, type:

`./BAdabouM/src/BAdabouM -h`



### Output
----------

BAdabouM's output file consists of the following columns:

1. Chromosome 1
2. Position 1.1
3. Position 2.1
4. Chromosome 2
5. Position 1.2
6. Position 2.2
7. Type of SV
8. Size of SV

A script to convert BadabouM's output into vcf format is available in the script repertory.

### Licence
-----------

BAdabouM is distributed under CeCILL.
see [LICENCE] (https://github.com/cumtr/BAdabouM/blob/master/LICENCE) file for more informations.

### How to cite
-----------

If you used this software and enjoyed it, please cite our preprint !

Tristan Cumer, François Pompanon, Frédéric Boyer
bioRxiv 2020.04.01.018127; doi: https://doi.org/10.1101/2020.04.01.018127
