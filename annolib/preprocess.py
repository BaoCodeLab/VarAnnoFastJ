import argparse
from annolib import gbk2seqGene2
from annolib import parseGtff


def preprocess(args):
    gbk2seqGene2.gbkConvert(args)
    parseGtff.parseGTFF(args)





