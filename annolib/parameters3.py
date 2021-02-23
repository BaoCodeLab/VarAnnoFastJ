import argparse
from annolib.snps2 import snps
from annolib.indels2 import indels
from annolib.preprocess import preprocess


pars=argparse.ArgumentParser(description="VarAnnoFastJ, a command line tool for variant annotation.")
pars_subparse=pars.add_subparsers()


def parse_arguments(args=None):
    snpPars=pars_subparse.add_parser("snp",help="annotation for SNPs")
    snpPars.add_argument('-c','--contg',required=True,help="dictionary of contig file parsed from parse_gtff")
    snpPars.add_argument('-g','--gene',required=True,help="dictionary of gene file parsed from parse_gtff")
    snpPars.add_argument('-C','--CDS',required=True,help="dictionary of CDS file parsed from parse_gtff")
    snpPars.add_argument('-m','--mol',help="dictionary of molecular file parsed from parse_gtff")
    snpPars.add_argument('-s','--CDS_seq',required=True,help="dictionary of CDS sequence file parsed from parse_gtff")

    snpPars_grp=snpPars.add_mutually_exclusive_group()
    snpPars_grp.add_argument('--table',nargs='*',help="tab-delimited files containing SNPs")
    snpPars_grp.add_argument('--vcf',nargs='*',help="vcf files containing SNPs")
#    snpPars.add_argument('-t','--table',nargs='*',required=True,help="tab-delimited file containing SNPs")
    snpPars.set_defaults(func=snps)

    indelPars=pars_subparse.add_parser("indel",help="annotation for INDELs")
    indelPars.add_argument('-c','--contg',required=True,help="dictionary of contig file parsed from parse_gtff")
    indelPars.add_argument('-g','--gene',required=True,help="dictionary of gene file parsed from parse_gtff")
    indelPars.add_argument('-C','--CDS',required=True,help="dictionary of CDS file parsed from parse_gtff")
    indelPars.add_argument('-m','--mol',help="dictionary of molecular file parsed from parse_gtff")
    indelPars.add_argument('-s','--CDS_seq',required=True,help="dictionary of CDS sequence file parsed from parse_gtff")

    indelPars_grp=indelPars.add_mutually_exclusive_group(required=True)
    indelPars_grp.add_argument('-t','--table',nargs='*',help="tab-delimited file containing INDELs")
    indelPars_grp.add_argument('-v','--vcf',nargs='*',help="vcf files containing INDELs")
#    indelPars.add_argument('-t','--table',nargs='*',required=True,help="tab-delimited file containing INDELs")
    indelPars.set_defaults(func=indels)

    annoPars=pars_subparse.add_parser("preprocess",help="preprocessing of annotation files")
    annoPars.add_argument('-g','--gbk',required=True,help="The gbk annotation file of the reference genome")
    annoPars.add_argument('-o1','--gtff',required=True,help="The tab-delimited gtff-format output")
    annoPars.add_argument('-o2','--outSeq',required=True,help="The output of protein sequences extracted from gbk annotation file")
    annoPars.add_argument('-t','--tax',required=True,help="The taxonomy or organism name of the reference genome")

    annoPars.add_argument('-r','--ref',required=True,help="The reference genome file")
    annoPars.add_argument('-p','--prefix',required=True,help="The prefix of dictionary output files")
    annoPars.set_defaults(func=preprocess)

    return pars.parse_args(args=args)



