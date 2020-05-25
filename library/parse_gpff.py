#!/usr/bin/python
##usage: xxx.py gbk-annotation-file output-prefix genome-ref-file
##This script was from parse_M23ND_gpff.py
##fixing the bug of lacking pseudogene/RNA location
import sys
import pickle

anno_file=open(sys.argv[1],'r')
out_contgDic=open(sys.argv[2]+'.contg.dic','w')
out_molDic=open(sys.argv[2]+'.mol.dic','w')
out_geneDic=open(sys.argv[2]+'.gene.dic','w')
out_CDSDic=open(sys.argv[2]+'.CDS.dic','w')
out_CDSseqDic=open(sys.argv[2]+'.CDSseq.dic','w')
out_pseudoDic=open(sys.argv[2]+'.pseudo.dic','w')
out_rnaDic=open(sys.argv[2]+'.RNA.dic','w')
out_pseudoSeqDic=open(sys.argv[2]+'.pseudoSeq.dic','w')
out_rnaSeqDic=open(sys.argv[2]+'.rnaSeq.dic','w')

big_sz=1000000000
gene_ref_file=open(sys.argv[3],'r')
gene_ref0,cunt='',0
contg_ref={}

while True:
	gene_ref1=gene_ref_file.read(big_sz)
	cunt+=1
	print cunt
	if gene_ref1=='':
		gene_ref=gene_ref0+'\n>'
	else:
		gene_ref=gene_ref0+gene_ref1
	st=gene_ref.find('>',0)
	while True:
                ed=gene_ref.find('\n>',st+1)
		if ed==-1:
			gene_ref0=gene_ref[st:]
			break
		the_blck=gene_ref[st:ed].split('\n')
#		print the_blck[0],the_blck[1]
#		print the_blck[0].split('|')[3]
		tit=the_blck[0].split()[0][1:]
		seq=''.join(the_blck[1:])
		contg_ref[tit]=seq
		
		st=ed+1
	if gene_ref1=='':
		break
gene_ref_file.close()


chr0,chrSt,chrEd,chrOri=1,2,3,4
contg,contgSt,contgEd,contgOri=5,6,7,8
gene,geneID,region,acc=9,10,11,13
nearRNAacc,RNAacc='',{}
contg_dic,mol_dic,gene_dic,CDS_dic,pseudo_dic,RNA_dic={},{},{},{},{},{}
while True:
	lin=anno_file.readline().rstrip()
	if lin=='':
		break
	lst=lin.split('\t')
	theRegion=lst[region]
	theContg=lst[contg].split('.')[0]
	theGene,theGeneID=lst[gene].split('.')[0],lst[geneID]
	theAcc=lst[acc].split('.')[0]

	ST,ED,DIR=int(lst[contgSt]),int(lst[contgEd]),lst[contgOri]

#	if theRegion=='RNA':
#		nearRNAacc=theAcc
#	if theAcc=='-':
#		theAcc=nearRNAacc
	if theRegion=='CDS' and theAcc=='-':
		continue
		
	if theRegion=='GENE' and contg_dic.has_key(theContg)==0:
		contg_dic[theContg] =[(ST,ED,DIR,theGeneID,theGene)]  # each contig contains all genes on it.
	elif theRegion=='GENE' and contg_dic.has_key(theContg)==1:
		contg_dic[theContg]+=[(ST,ED,DIR,theGeneID,theGene)]

        elif theRegion=='PSEUDO' and mol_dic.has_key(theGeneID)==0:
		mol_dic[theGeneID] =['PSEUDO']
                pseudo_dic[(theContg,theGeneID,theAcc)]=[(ST,ED,DIR,theGene)]  # This dictionary contains all pseudogenes and ranges
	elif theRegion=='PSEUDO' and mol_dic.has_key(theGeneID)==1:
#		mol_dic[theGeneID]+=['PSEUDO']
                pseudo_dic[(theContg,theGeneID,theAcc)]+=[(ST,ED,DIR,theGene)]
	elif theRegion.find('RNA')!=-1 and mol_dic.has_key(theGeneID)==0: 
		mol_dic[theGeneID]=['RNA']
                RNA_dic[(theContg,theGeneID,theAcc)]=[(ST,ED,DIR,theGene)]  # This dictionary contains all non-coding RNAs and ranges
#		RNAacc[theGeneID]=[theAcc]
	elif theRegion.find('RNA')!=-1 and mol_dic.has_key(theGeneID)==1:
#		mol_dic[theGeneID]+=['RNA']
                RNA_dic[(theContg,theGeneID,theAcc)]=+[(ST,ED,DIR,theGene)]
#		RNAacc[theGeneID]+=[theAcc]

	elif (theRegion=='UTR' or theRegion=='CDS') and gene_dic.has_key(theGeneID)==0: # This dictionary contains all coding genes and ranges
		gene_dic[theGeneID] =[(ST,ED,DIR,theRegion,theAcc)]
	elif (theRegion=='UTR' or theRegion=='CDS') and (ST,ED,DIR,theRegion,theAcc) not in gene_dic[theGeneID]:
		gene_dic[theGeneID]+=[(ST,ED,DIR,theRegion,theAcc)]

	if   theRegion=='CDS' and CDS_dic.has_key((theContg,theGeneID,theAcc))==0: # This dictionary contains all coding genes and coding regions
		CDS_dic[(theContg,theGeneID,theAcc)] =[(ST,ED,DIR,theGene)]
	elif theRegion=='CDS' and CDS_dic.has_key((theContg,theGeneID,theAcc))==1:
		CDS_dic[(theContg,theGeneID,theAcc)]+=[(ST,ED,DIR,theGene)]

anno_file.close()
#for kys in sorted(CDS_dic.keys()):
#	print kys,CDS_dic[kys]

##output the gene-sequence in a dictionary
#gene_refDic=open(sys.argv[3]+'.gene.dic','w')
CDSseq_dic={}
for CDS_kys in sorted(CDS_dic.keys()):
	geneSeq=''
	seqContg=CDS_kys[0]

	if contg_ref.has_key(seqContg)==0:
		print 'the contig %s does not exist!' % seqContg
		break

	CDS_pos=sorted(CDS_dic[CDS_kys])
	for echPos in CDS_pos:
		geneSeq+=contg_ref[seqContg][echPos[0]-1:echPos[1]]
	CDSseq_dic[CDS_kys]=geneSeq
	if len(geneSeq)%3!=0:
#		print CDS_pos
                print 'the protein is incomplete, recheck your annotation: ', CDS_kys

pseudoSeq_dic={}
for pseudo_kys in sorted(pseudo_dic.keys()):
        pseudoSeq=''
        pseudoContg=pseudo_kys[0]

        if contg_ref.has_key(pseudoContg)==0:
                print 'the contig %s does not exist!' % pseudoContg
                break
        pseudo_pos=sorted(pseudo_dic[pseudo_kys])
        for ech_pseudoPos in pseudo_pos:
                pseudoSeq+=contg_ref[pseudoContg][ech_pseudoPos[0]-1:ech_pseudoPos[1]]
        pseudoSeq_dic[pseudo_kys]=pseudoSeq
        if len(pseudoSeq)%3==0:
                print 'the pseudogene seems a normal coding gene,recheck your annotation: ', pseudo_kys

rnaSeq_dic={}
for RNA_kys in sorted(RNA_dic.keys()):
        rnaSeq=''
        rnaContg=RNA_kys[0]

        if contg_ref.has_key(rnaContg)==0:
                print 'the contig %s does not exist!' % rnaContg
                break
        rna_pos=sorted(RNA_dic[RNA_kys])
        for ech_rnaPos in rna_pos:
                rnaSeq+=contg_ref[rnaContg][ech_rnaPos[0]-1:ech_rnaPos[1]]
        rnaSeq_dic[RNA_kys]=rnaSeq


#pickle.dump(gene_ref_dic,out_ProtDic)
#out_ProtDic.close()

#sort is not needed here!!!
pickle.dump(contg_dic,out_contgDic)
pickle.dump(mol_dic,out_molDic)
pickle.dump(gene_dic,out_geneDic)
pickle.dump(CDS_dic,out_CDSDic)
pickle.dump(CDSseq_dic,out_CDSseqDic)
pickle.dump(pseudo_dic,out_pseudoDic)
pickle.dump(RNA_dic,out_rnaDic)
pickle.dump(pseudoSeq_dic,out_pseudoSeqDic)
pickle.dump(rnaSeq_dic,out_rnaSeqDic)


out_contgDic.close()
out_molDic.close()
out_geneDic.close()
out_CDSDic.close()
out_CDSseqDic.close()
out_pseudoDic.close()
out_rnaDic.close()
out_pseudoSeqDic.close()
out_rnaSeqDic.close()
