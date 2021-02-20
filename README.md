# VarAnnoFastJ
VarAnnoFastJ is a python tool for accurate variant annotation.
## Install
### Requirements
   Python >= 3.7
   numpy >= 1.17.3
   pandas >= 1.1.4
   gffutils >= 0.10.1
   
####
    git clone https://github.com/BaoCodeLab/VarAnnoFastJ
    cd VarAnnoFastJ
    python setup.py install
    
    
## Running VarAnnoFastJ
### Usage:  VarAnnoFastJ.py [-h] {preprocess,snp,indel}

#### Command line usage                        
    -h, --help          show the help message.
    preprocess             Calculate the SNP density on the genome.
    snp             Perform SNP clustering.
    indel                Calculate the significance p-value of SNP clustering.
    
