# ddRAD_NewWorldMonkeys
Scripts and parameters used to identify and genotype large numbers of SNP loci for taxa from across the New World monkeys  and infer phylogenetic relationships using iPYRAD (http://ipyrad.readthedocs.io/).

### Checking read quality using FASTQC (https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) 

```
fastqc filename.fastqc.zip
```

Based on FASTQC results trim your data based on:
1. Overall sequence quality
2. Per base sequence quality
3. Adapter contamination


### Adapter and quality trimming using BBDuk (https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbduk-guide/)
It is recommended to do first adapter removal and then quality filter to be able to maximize the matches between contaminants and  adapter sequence. 


**1. Adapter trimming:**
```
bbduk.sh in1=.file_R1.fq.gz in2=file_R2.fq.gz out1=R1_cl.fq out2=R2_cl.fq outm=_R1_uncl.fq outm2=R2_uncl.fq  
outs=singletons.fq stats=Stats.txt ref=adapters.fa  ktrim=r k=23 mink=11 hdist=3 hdist2=1 rcomp=t tbo tpe minlen=30 
```
- k = kmer length (at least the length of the barcode).
- tbo = trim adapters based on pair overlap detection using BBMerge (which does not require known adapter sequences).
- tpe = trim both reads to the same length (in the event that an adapter kmer was only detected in one of them).


**2. Contaminant Phix filtering**

 The PhiX genome is used in library prepartion and is used as a control during Illumina sequencing.
 It provides a quality control for cluster generation, sequencing, and alignment, and a calibration control for 
 cross-talk matrix generation, phasing, and prephasing.

```
bbduk.sh iN2=./R1_cl.fq in2=R2_cl.fq out1=R1_clphix.fq  out2=R2_clphix.fq outm=R1_unclphix.fq outm2=R2_unclphix.fq outs=singphix.fq stats=statsphix.txt ref=/home1/02202/lmv498/bin/bbmap/resources/phix174_ill.ref.fa  k=31 hdist=1 
```


### Demultiplexing using deML (https://github.com/grenaud/deML)

To be able to run deML you need three files:
1. Index file: First 5bp (assumed barcode) of your R1
2. R1 without the first 5bp
3. R2 

Use FASTX tools to trim the first 5bp (resulting in input R1) and the rest of the read (resulting in INDEX file)

```
fastx_trimmer -Q33 -f 6  -i R1_clphix.fq -o R1.fq
fastx_trimmer -Q33 -l 5  -i R1_clphix.fq -o R1_ID.fq
```

However deML might not recognize the read name of your specific output files and it might need to be modified:
- deML specific input file: @K00179:39:H7JCJBBXX:1:1101:6837:1877/1
- My file: @K00179:39:H7JCJBBXX:1:1101:6837:1877 1:N:0:CGTACG

To change names:
For R1: 
```
sed -i.bak 's/\s.*/\/1/' R1.fq *
```
For R2:
```
sed -i.bak 's/\s.*/\/2/' R2_clphix.fq
```
For ID:
```
sed -i.bak 's/\s.*//' R1_ID.fq
```

Create a barcode.txt file:
```
#Index1	Name
ATCGT	sample
```
Run deML
```
deML -i barcodes.txt -f R1.fq -r R2.fq -if1 R1_ID.fq --mm 1 -o finalreads.fq -s stats.txt 
```


### SNP Discovery using IPYRAD (http://ipyrad.readthedocs.io/)
Don't forget you can do quality trimming and demultiplex in iPYRAD if you want to.

**1. Create params file**
```ipyrad -n nwm-params```

Params file example:
```
NWM                        ## [0] [assembly_name]: Assembly name. Used to name output directories for assembly steps       
./                         ## [1] [project_dir]: Project dir (made in curdir if not present)
                           ## [2] [raw_fastq_path]: Location of raw non-demultiplexed fastq files
                           ## [3] [barcodes_path]: Location of barcodes file*
*.fq.gz                    ## [4] [sorted_fastq_path]: Location of demultiplexed/sorted fastq files
denovo                     ## [5] [assembly_method]: Assembly method (denovo, reference, denovo+reference, denovo-reference)
                           ## [6] [reference_sequence]: Location of reference sequence file
pairddrad                  ## [7] [datatype]: Datatype (see docs): rad, gbs, ddrad, etc.
CATGC, AATT                ## [8] [restriction_overhang]: Restriction overhang (cut1,) or (cut1, cut2)
7                          ## [9] [max_low_qual_bases]: Max low quality base calls (Q<20) in a read
                           ## [10] [phred_Qscore_offset]: phred Q score offset (33 is default and very standard)
6                          ## [11] [mindepth_statistical]: Min depth for statistical base calling
6                          ## [12] [mindepth_majrule]: Min depth for majority-rule base calling
10000                      ## [13] [maxdepth]: Max cluster depth within samples
0.85                       ## [14] [clust_threshold]: Clustering threshold for de novo assembly
6                          ## [15] [max_barcode_mismatch]: Max number of allowable mismatches in barcodes
6                          ## [16] [filter_adapters]: Filter for adapters/primers (1 or 2=stricter)
35                         ## [17] [filter_min_trim_len]: Min length of reads after adapter trim
2                          ## [18] [max_alleles_consens]: Max alleles per site in consensus sequences
6, 6                       ## [19] [max_Ns_consens]: Max N's (uncalled bases) in consensus (R1, R2)
8, 8                       ## [20] [max_Hs_consens]: Max Hs (heterozygotes) in consensus (R1, R2)
4                          ## [21] [min_samples_locus]: Min # samples per locus for output
35, 35                     ## [22] [max_SNPs_locus]: Max # SNPs per locus (R1, R2)
100, 100                   ## [23] [max_Indels_locus]: Max # of indels per locus (R1, R2)
1.0                        ## [24] [max_shared_Hs_locus]: Max # heterozygous sites per locus (R1, R2)
5, 0, 4, 0                 ## [25] [trim_reads]: Trim raw read edges (R1>, <R1, R2>, <R2) (see docs)
0, 0, 0, 0                 ## [26] [trim_loci]: Trim locus edges (see docs) (R1>, <R1, R2>, <R2)
*                          ## [27] [output_formats]: Output formats (see docs)
                           ## [28] [pop_assign_file]: Path to population assignment file
```

**2. Run s123**
```ipyrad -p params-file.txt -s 123 -c 24 -t 24 --ipcluster &>s123.log```

**3. Run s456**
To estimate the number of maximum Ns (uncalled bases) and Hs (heterozigous sites) in the consensus read run s45 with the following parameters:
```
100, 100                          ## [19] [max_Ns_consens]: Max N's (uncalled bases) in consensus (R1, R2)
100, 100                          ## [20] [max_Hs_consens]: Max Hs (heterozygotes) in consensus (R1, R2)*
``` 

**4. Evaluate the distribution of the number of average number of H's and N's**
``` 
cat sample_consens/*consens.gz* > merged_catcons.gz
```

Run python code to calculate the variation in Hs and Ns
```
HNs_catcons.py merged_catcons.gz 
HNs_catcons.R   
```
Look at the HNs_per_consensus.png and choose the apropriate values for the parameters based on the highest percentage of
your type of reads (i.e., merged or unmerged). Pick the 95% CI and use it as your parameter for [19] [max_Ns_consens] and [20] [max_Hs_consens].


**5. Run s456 again with the new parameters**
```
ipyrad -p params-file.txt -f  -s 45 -c 24 -t 24 --ipcluster &>s123.log
```


**6. Run s7**
To estimate the number of maximum SNPs and Indels per locus run s7 with the following parameters:
```
100, 100                          ## [22] [max_SNPs_locus]: Max # SNPs per locus (R1, R2)
100, 100                          ## [23] [max_Indels_locus]: Max # of indels per locus (R1, R2)*
```

**7. Evaluate the distribution of the distribution of the  number of SNPs and INDELS per locus**
```
SNPs_locus.py samples.loci
SNPs_locus.R
```

Look at the SNPs_per_locus.png and choose the apropriate values for the parameters. Pick the 95% CI and use it as your parameter for [22] [max_SNPs_locus] and [23] [max_Indels_locus].



## Infering phylogenies using IPYRAD

```
usage: tetrad [-h] [-v] [-f] [-s seq] [-j json] [-m method] [-q nquartets]
              [-b boots] [-l map_file] [-r resolve] [-n name] [-o outdir]
              [-t starting_tree] [-c CPUs/cores] [-x random_seed] [-d] [--MPI]
              [--ipcluster]
```






