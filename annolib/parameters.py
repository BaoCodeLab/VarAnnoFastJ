import argparse
import os
import sys
from annolib.snps import snp
from annolib.indels import indel

pars=argparse.ArgumentParser(description="VarAnnoFastJ, a command line tool for variant annotation.")
pars_subparse=pars.add_subparsers()


def parse_arguments(args=None):
    snpPars=pars_subparse.add_parser("snp",help="annotation for SNPs")
    snpPars.add_argument('-c','--contg',required=True,help="dictionary of contig file parsed from parse_gpff")
    snpPars.add_argument('-g','--gene',required=True,help="dictionary of gene file parsed from parse_gpff")
    snpPars.add_argument('-C','--CDS',required=True,help="dictionary of CDS file parsed from parse_gpff")
    snpPars.add_argument('-m','--mol',help="dictionary of molecular file parsed from parse_gpff")
    snpPars.add_argument('-s','--CDS_seq',required=True,help="dictionary of CDS sequence file parsed from parse_gpff")
    snpPars.add_argument('-t','--table',nargs='*',required=True,help="tab-delimited file containing SNPs")
    snpPars.set_defaults(func=snp)

    indelPars=pars_subparse.add_parser("indel",help="annotation for INDELs")
    indelPars.add_argument('-c','--contg',required=True,help="dictionary of contig file parsed from parse_gpff")
    indelPars.add_argument('-g','--gene',required=True,help="dictionary of gene file parsed from parse_gpff")
    indelPars.add_argument('-C','--CDS',required=True,help="dictionary of CDS file parsed from parse_gpff")
    indelPars.add_argument('-m','--mol',action="store_false",help="dictionary of molecular file parsed from parse_gpff")
    indelPars.add_argument('-s','--CDS_seq',required=True,help="dictionary of CDS sequence file parsed from parse_gpff")
    indelPars.add_argument('-t','--table',nargs='*',required=True,help="tab-delimited file containing INDELs")
    indelPars.set_defaults(func=indel)


    return pars.parse_args(args=args)

