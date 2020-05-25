#!/usr/bin/python
import sys
import argparse
import gffutils
import pandas
import numpy

parser=argparse.ArgumentParser(description="transform the gff format to tab format")
parser.add_argument("-gff",required=True,help="The annotation .gff file of the reference")
parser.add_argument("-out",required=True,help="The output file")
args=parser.parse_args()

if args.gff:
    gff_file=args.gff
else:
    raise Exception("Please provide the .gff file")

geneID,st_ed,strand,geneName,locusTag,codetype=[],[],[],[],[],[]
myDB=gffutils.create_db(gff_file,dbfn="gff.db",force=True,keep_order=True,merge_strategy='merge',sort_attribute_values=False)
myDB=gffutils.FeatureDB('gff.db')
features=myDB.all_features(featuretype='gene')
print ('\nThe type of features is: ', type(features))
for ii in features:
    geneID+=ii.attributes['ID']
    st_ed+=[[ii.start,ii.end]]
    strand+=[ii.strand]
    geneName+=ii.attributes['Name']
#    locusTag+=ii.attributes['locus_tag']
    locusTag+=['-']
    codetype+=ii.attributes['gene_biotype']

geneFrame=pandas.DataFrame(st_ed,columns=['start','end'],index=geneID)
geneFrame['strand']=strand
geneFrame['length']=geneFrame.max(axis=1)-geneFrame.min(axis=1)+1
geneFrame['geneID']=geneID
geneFrame['geneSymb']=geneName
geneFrame['locus']=locusTag
geneFrame['code_type']=codetype
geneFrame['gene_type']='GENE'

dist=[0]
start=list(geneFrame['start'])
end=list(geneFrame['end'])
for i in range(1,len(start)):
    dist+=[min(start[i],end[i])-max(start[i-1],end[i-1])]

geneFrame['dist']=dist

#print(dist[0:10])

geneID2,product,protID=[],[],[]
featuresCDS=myDB.all_features(featuretype='CDS')
print ('\nThe type of features is: ', type(featuresCDS))
for jj in featuresCDS:
    geneID2+=jj.attributes['Parent']
    product+=jj.attributes['product']
    protID+=jj.attributes['protein_id']


cdsFrame=pandas.DataFrame(protID,columns=['protein_id'],index=geneID2)
cdsFrame['geneID']=geneID2
cdsFrame['product']=product
cdsFrame['cds_type']='CDS'


#totalFrame=pandas.merge(geneFrame,cdsFrame,on=["geneID","cod_type"],how="outer",sort=False)
totalFrame=pandas.merge(geneFrame,cdsFrame,on="geneID",how="outer",sort=False)

newGeneFrame=totalFrame[['start','end','strand','length','geneID','geneSymb','locus','gene_type','dist','product','protein_id']]
newCdsFrame=totalFrame[['start','end','strand','length','geneID','geneSymb','locus','cds_type','dist','product','protein_id']]
newGeneFrame.rename(columns={'gene_type':'mol_type'},inplace=True)
newCdsFrame.rename(columns={'cds_type':'mol_type'},inplace=True)

finalFrame=pandas.concat([newGeneFrame,newCdsFrame])

#print(finalFrame)

if args.out:
    finalFrame.to_csv(args.out,sep="\t",index=True,float_format='%.4f')
else:
    raise Exception("Please provide the output file name")




