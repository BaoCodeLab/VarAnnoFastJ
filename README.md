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
   
#### The preprocess command processes the gbk annotation file of the reference genome and generates three sets of output files for variant annotation. The output includes the tab-delimited file (general tab file format), the protein sequences encoded by the reference genome, and a series of python dictionary files representing distinct attributes of the genome annotation. An example of the dictionary output files is as following:
      gtff_parse.CDS.dic
      gtff_parse.CDSseq.dic
      gtff_parse.contg.dic
      gtff_parse.gene.dic
      gtff_parse.mol.dic
      gtff_parse.pseudo.dic
      gtff_parse.pseudoSeq.dic
      gtff_parse.RNA.dic
     gtff_parse.rnaSeq.dic
      
      
      
      
      
      
      
      
  
                        
   
  
  
  
  
