# VarAnnoFastJ
VarAnnoFastJ is a python tool for accurate variant annotation.
## Install
### Requirements
  * Python >= 3.7
  * numpy >= 1.17.3
  * pandas >= 1.1.4
  * gffutils >= 0.10.1
   
####
    git clone https://github.com/BaoCodeLab/VarAnnoFastJ
    cd VarAnnoFastJ
    python setup.py install
    
    
## Running VarAnnoFastJ
### Usage:  VarAnnoFastJ.py [-h] {preprocess,snp,indel}

#### Command line usage                        
    -h, --help          Show the help message.
    preprocess          Preprocessing of gene annotation files.
    snp                 Annotation of SNPs.
    indel               Annotation of Indels.

### Preprocessing of gene annotation files
#### Usage:  VarAnnoFastJ.py preprocess [-h] -g GBK -o1 GTFF -o2 OUTSEQ -t TAX -r REF -p prefix

#### {arguments}
    -h, --help                          Show this help message and exit.
    -g GBK, --gbk GBK                   The gbk annotation file of the reference genome.
    -o1 GTFF, --gtff GTFF               The tab-delimited gtff-format converted from gbk-format file. 
    -o2 OUTSEQ, --outSeq OUTSEQ         The output of protein sequences extracted from gene annotation file of gbk format.
    -t TAX, --tax TAX                   The taxonomy or organism name of the reference genome.
    -r REF, --ref REF                   The reference genome file of fasta format.
    -p PREFIX, --prefix PREFIX          The prefix of dictionary output files.
   
#### The preprocess command processes the gbk annotation file of the reference genome and generates three sets of output files. The output includes the tab-delimited file (general tab file format), the protein sequences encoded by the reference genome, and a series of python dictionary files representing distinct attributes of the genome annotation. The dictionary files will be used as input for variant annotation. An example of the dictionary output files is as following:
    gtff_parse.contg.dic         The contigs in the reference genome.
    gtff_parse.gene.dic              All genes encoded by the reference genome.
    gtff_parse.mol.dic           The molecular types of each gene encoded by the reference genome (i.e., gene, CDS, rRNA, tRNA).
    gtff_parse.CDS.dic           The coding genes encoded by the reference genome.
    gtff_parse.CDSseq.dic        The sequence of coding genes. 
    gtff_parse.pseudo.dic        The pseudogenes of the reference genome.
    gtff_parse.pseudoSeq.dic     The sequences of pseudogenes.
    gtff_parse.RNA.dic           The RNAs of the reference genome (i.e., rRNA, tRNA).
    gtff_parse.rnaSeq.dic        The sequences of RNAs.
    
 ### Annotation of SNPs
 #### Usage:    VarAnnoFastJ.py snp [-h] -c CONTG -g GENE -C CDS [-m MOL] -s CDS_SEQ -t [TABLE [TABLE ...]]
 #### Example:  VarAnnoFastJ.py snp -c gtff_parse.contg.dic -g gtff_parse.gene.dic -C gtff_parse.CDS.dic -m gtff_parse.mol.dic -s gtff_parse.CDSseq.dic -t chr1.snp.table       chr2.snp.table chr3.snp.table
 
 #### {arguments}
     -h, --help            show this help message and exit
     -c CONTG, --contg CONTG      dictionary of contig file parsed from parse_gpff
     -g GENE, --gene GENE  dictionary of gene file parsed from parse_gpff
     -C CDS, --CDS CDS     dictionary of CDS file parsed from parse_gpff
     -m MOL, --mol MOL     dictionary of molecular file parsed from parse_gpff
     -s CDS_SEQ, --CDS_seq CDS_SEQ        dictionary of CDS sequence file parsed from parse_gpff
     -t [TABLE [TABLE ...]], --table [TABLE [TABLE ...]]       tab-delimited file containing SNPs

#### The snp command makes annotations for SNPs using the dictionary output files from the command preprocess. The snp file should be in vcf format or user-prepared tab-delimited format with each line for a single SNP. The first, second, and third column of the user-prepared SNP file should be contig ID/chromosome ID, SNP start position on contig/chromosome, SNP end position/chromosome, respectively. SNP end position is only needed for fragments of consecutive SNP.   
      

      
      
      
      
      
      
  
                        
   
  
  
  
  
