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

### Preprocessing of reference genome annotation files
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
 #### Usage:    VarAnnoFastJ.py snp [-h] -c CONTG -g GENE -C CDS [-m MOL] -s CDS_SEQ -t [TABLE [TABLE ...]]
 #### Example:  VarAnnoFastJ.py snp -c gtff_parse.contg.dic -g gtff_parse.gene.dic -C gtff_parse.CDS.dic -m gtff_parse.mol.dic -s gtff_parse.CDSseq.dic -t chr1.snp.table       chr2.snp.table chr3.snp.table
 
 #### {arguments}
     -h, --help                               Show this help message and exit
     -c CONTG, --contg CONTG                  Dictionary file of contigs parsed from parse_gpff
     -g GENE, --gene GENE                     Dictionary file of genes parsed from parse_gpff
     -C CDS, --CDS CDS                        Dictionary file of CDS parsed from parse_gpff
     -m MOL, --mol MOL                        Dictionary file of molecular type parsed from parse_gpff
     -s CDS_SEQ, --CDS_seq CDS_SEQ            Dictionary file of CDS sequences parsed from parse_gpff
     -t [TABLE ...], --table [TABLE ...]      Tab-delimited files containing SNPs

#### The snp command makes annotations for SNPs using the dictionary output files from the command preprocess. The snp file should be in vcf format or user-prepared tab-delimited format with each line for a single SNP. The SNP file could contains as many columns as you want. But the first, second, and third column of the user-prepared SNP file should be contig ID/chromosome ID, SNP start position on contig/chromosome, and SNP end position/chromosome, respectively. SNP end position is only needed for fragments of consecutive SNP. Multiple SNP files can be provided and a separated SNP annotation file with suffix .anno will be generated for each SNP input file. An example of the output looks like:
       Chromosome  Start    End   Ref  Align    MutType        GeneSymbol     GeneID                         Strand   AAchange    AAposition
       FM252032.1    157    157    C     T     CDS_synon       dnaA           SSUBM407_0001                     +                     40
       FM252032.1    204    204    G     T     CDS_nonSynon    dnaA           SSUBM407_0001                     +     Ser56Ile        56
       FM252032.1   1384   1384    C     G     inter-gene-3-5  dnaA||dnaN||   SSUBM407_0001||SSUBM407_0002||    0
       FM252032.1   2571   2571    G    T,G    CDS_nonSynon    dnaN           SSUBM407_0002                     +     Ser349Ile      349 
#### The fourth and fifth column are Reference allele and Aligned allele, respectively. The sixth column "MutType" refers to the mutation type, including CDS_synon,  CDS_nonSynon, inter-gene-3-5, inter-gene-5-3, inter-gene-3-3, inter-gene-5-5, UTR, intron-splice-5, intron-splice-3, intron, and pseudo-gene|ncRNA. The number "5" and "3" indicate the genic region of the mutation relative to the upstream or downstream gene. For example, "inter-gene-3-5" means the mutation is located between the 3'-end of the upstream gene and 5'-end of the downstream gene, "intron-splice-5" means the mutation is located in the splicing sites of the 5'-end of a gene.      
      
### Annotation of Indels
#### Usage:  VarAnnoFastJ.py indel [-h] -c CONTG -g GENE -C CDS [-m] -s CDS_SEQ -t [TABLE [TABLE ...]]
#### Example:  VarAnnoFastJ.py indel -c gtff_parse.contg.dic -g gtff_parse.gene.dic -C gtff_parse.CDS.dic -m gtff_parse.mol.dic -s gtff_parse.CDSseq.dic -t chr1.indel.table chr2.indel.table chr3.indel.table


 #### {arguments}
     -h, --help                               Show this help message and exit
     -c CONTG, --contg CONTG                  Dictionary file of contigs parsed from parse_gpff
     -g GENE, --gene GENE                     Dictionary file of genes parsed from parse_gpff
     -C CDS, --CDS CDS                        Dictionary file of CDS parsed from parse_gpff
     -m MOL, --mol MOL                        Dictionary file of molecular type parsed from parse_gpff
     -s CDS_SEQ, --CDS_seq CDS_SEQ            Dictionary file of CDS sequences parsed from parse_gpff
     -t [TABLE ...], --table [TABLE ...]      Tab-delimited files containing Indels

#### The indel command makes annotations for Indels using the dictionary output files from the command preprocess. The indel file should be in vcf format or user-prepared tab-delimited format with each line for a single Indel. The Indel file could contains as many columns as you want. However, the first, second, and third column of the user-prepared Indel file should be contig ID/chromosome ID, Indel start position on contig/chromosome, and Indel end position/chromosome, respectively. Indel start and end position indicate the start and end position of the fragment of insertion/deletion. The users should make sure uniform alignment mode (left or right alignment) for indels of tandom repeat sequences in the alignment stage. Multiple Indel files can be provided and a separated Indel annotation file with suffix .anno will be generated for each Indel input file. An example of the output looks like:
     Chromosome     Start      End   Ref  Align    MutType        GeneSymbol         GeneID               Strand   ShiftPosition
       AP53        708559   708562   AAT    ---    in-frame              mac         AP53_713                -
       AP53        708569   708572   ---    TTA    in-frame-Stp          mac         AP53_713                -
       AP53        743775   743775    T      -     frame-shift          amiC         AP53_750                +          86
       AP53        755767   755768    --     CT    inter-gene-3-5    parE||parC||    AP53_756||AP53_757||    0
 ####  Similar to the annotation file of SNPs, the fourth and fifth column are Reference allele and Aligned allele, respectively. The sixth column "MutType" refers to the mutation type, including in-frame,  in-frame-Stp, frame-shift, inter-gene-3-5, inter-gene-5-3, inter-gene-3-3, inter-gene-5-5, UTR, intron-splice-5, intron-splice-3, intron, and pseudo-gene|ncRNA. The number "5" and "3" indicate the genic region of the mutation relative to the upstream and downstream gene. For example, "inter-gene-3-5" means the mutation is located between the 3'-end of the upstream gene and 5'-end of the downstream gene, "intron-splice-5" means the mutation is located in the splicing sites of the 5'-end of a gene.   
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
      
      
      
      
      
  
                        
   
  
  
  
  
