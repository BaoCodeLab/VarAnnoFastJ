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
    gtff_parse.gene.dic          All genes encoded by the reference genome.
    gtff_parse.mol.dic           The molecular types of each gene encoded by the reference genome (i.e., gene, CDS, rRNA, tRNA).
    gtff_parse.CDS.dic           The coding genes encoded by the reference genome.
    gtff_parse.CDSseq.dic        The sequence of coding genes. 
    gtff_parse.pseudo.dic        The pseudogenes of the reference genome.
    gtff_parse.pseudoSeq.dic     The sequences of pseudogenes.
    gtff_parse.RNA.dic           The RNAs of the reference genome (i.e., rRNA, tRNA).
    gtff_parse.rnaSeq.dic        The sequences of RNAs.
    
 ### Annotation of SNPs
 #### Usage:  VarAnnoFastJ.py snp [-h] 
 
 #### {arguments}
 
      
      
      
      
      
      
  
                        
   
  
  
  
  
