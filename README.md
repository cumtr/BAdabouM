<p align="center">
  <a>
    <img height="150" src="logo/BAdabouM_logo.png">
  </a>
</p>

BAdabouM is a free tool to detect Structural Varations (SVs) in Whole Genome reSequencing data.
BAdabouM is able to genotype insertions (INS), deletions (DEL), inversions (INV), copy number variations (CNVs) and inter and intra chromosomal rearrangements (CTX and ITX).
BAdabouM uses paired-ends, split-reads and aberrant depth of coverage to call SVs.


### Instalating BAdabouM
------------------------

The easiest way to get BAdabouM is to build it from sources.

```
git clone https://github.com/cumtr/BAdabouM.git

cd BAdabouM/

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
3. Position 1.2
4. Chromosome 2
5. Position 2.1
6. Position 2.2
7. Type of a SV
8. Size of a SV


### Licence
-----------

BAdabouM is distributed under CeCILL.
see [LICENCE] (https://github.com/cumtr/BAdabouM/blob/master/LICENCE) file for more informations.


