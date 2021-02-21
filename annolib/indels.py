## author: yjbao
import pickle


contg_dic,gene_dic,CDS_dic,mol_dic,CDS_seq_dic={},{},{},{},{}


class StrToBytes:
    def __init__(self,fileobj):
        self.fileobj = fileobj
    def read(self,size):
        return self.fileobj.read(size).encode()
    def readline(self,size=-1):
        return self.fileobj.readline(size).encode()

def create_contg_dic(contg_filefile):
    global contg_dic
    contg_file=open(contg_filefile,'r')
    contg_dic=pickle.load(StrToBytes(contg_file))

    contg_file.close()
#    return contg_dic


def create_gene_dic(gene_filefile):
    global gene_dic
    gene_file=open(gene_filefile,'r')
    gene_dic=pickle.load(StrToBytes(gene_file))
    gene_file.close()
#   return gene_dic


def create_CDS_dic(CDS_filefile):
    global CDS_dic
    CDS_file=open(CDS_filefile,'r')
    CDS_dic=pickle.load(StrToBytes(CDS_file))
    CDS_file.close()
#   return CDS_dic


def create_mol_dic(mol_filefile):
    global mol_dic
    mol_file=open(mol_filefile,'r')
    mol_dic=pickle.load(StrToBytes(mol_file))
    mol_file.close()
#   return mol_dic


## input CDS sequences for determination of CDS location
def create_CDSseq_dic(CDS_seq_filefile):
    global CDS_seq_dic
    CDS_seq_file=open(CDS_seq_filefile,'r')
    CDS_seq_dic=pickle.load(StrToBytes(CDS_seq_file))
    CDS_seq_file.close()
#   return CDS_seq_dic



cod=    {'frame-shift':-2,'in-frame':-1,'in-frame-Stp':-3,'UTR':1,'intron-splice-5':2,'intron-splice-3':3,'intron':4,'inter-gene-3-5':5,'inter-gene-5':6,'inter-gene-3':7,'inter-gene-5-3':8,'inter-gene':9,'pseudo-gene|ncRNA':10}
antiCod={-2:'frame-shift',-1:'in-frame',-3:'in-frame-Stp',1:'UTR',2:'intron-splice-5',3:'intron-splice-3',4:'intron',5:'inter-gene-3-5',6:'inter-gene-5',7:'inter-gene-3',8:'inter-gene-5-3',9:'inter-gene',10:'pseudo-gene|ncRNA'}

nt=['TTT','TTC','TTA','TTG','TCT','TCC','TCA','TCG','TAT','TAC','TAA','TAG','TGT','TGC','TGA','TGG',\
    'CTT','CTC','CTA','CTG','CCT','CCC','CCA','CCG','CAT','CAC','CAA','CAG','CGT','CGC','CGA','CGG',\
    'ATT','ATC','ATA','ATG','ACT','ACC','ACA','ACG','AAT','AAC','AAA','AAG','AGT','AGC','AGA','AGG',\
    'GTT','GTC','GTA','GTG','GCT','GCC','GCA','GCG','GAT','GAC','GAA','GAG','GGT','GGC','GGA','GGG']

aa=['Phe','Phe','Leu','Leu','Ser','Ser','Ser','Ser','Tyr','Tyr','Stp','Stp','Cys','Cys','Stp','Trp',\
    'Leu','Leu','Leu','Leu','Pro','Pro','Pro','Pro','His','His','Gln','Gln','Arg','Arg','Arg','Arg',\
    'Ile','Ile','Ile','Met','Thr','Thr','Thr','Thr','Asn','Asn','Lys','Lys','Ser','Ser','Arg','Arg',\
    'Val','Val','Val','Val','Ala','Ala','Ala','Ala','Asp','Asp','Glu','Glu','Gly','Gly','Gly','Gly']

aa_single=['F','F','L','L','S','S','S','S','Y','Y','*','*','C','C','*','W',\
    'L','L','L','L','P','P','P','P','H','H','Q','Q','R','R','R','R',\
    'I','I','I','M','T','T','T','T','N','N','K','K','S','S','R','R',\
    'V','V','V','V','A','A','A','A','D','D','E','E','G','G','G','G']

coden={}

def coden_dic():
    global coden
    for ii in range(0,len(aa)):
        coden[nt[ii]]=aa[ii]


def char_cmp(str0):
    if str0=='A':
        return 'T'
    elif str0=='T':
        return 'A'
    elif str0=='G':
        return 'C'
    elif str0=='C':
        return 'G'
    else:
        return 'N'


def find_cod(str1,str2,str3,snp_pos,snp_allel):
    snpCDSpos=0
    CDS_ky=(str1,str2,str3)
    CDS_seq=CDS_seq_dic[CDS_ky]
    CDS_lg=len(CDS_seq)

    if CDS_lg%3!=0:
        snp_anno,snp_prot='CDS_unanno',''
        return snp_anno,snp_prot


    CDS_dir=CDS_dic[CDS_ky][0][2]
    for CDS_val in sorted(CDS_dic[CDS_ky]):
        if snp_pos < CDS_val[0]:
            snp_anno,snp_prot='UTR',''
            break
        snpCDSpos+=min(snp_pos,CDS_val[1])-CDS_val[0]+1  ## the SNP position in the CDS

    trip_st,trip_pos=(snpCDSpos-1)//3*3,snpCDSpos%3 # the coden starting position in the CDS and the SNP frame in the coden
    aa_pos=(snpCDSpos-1)//3+1  ## amino acid position in the protein


    trip_ref=CDS_seq[trip_st:trip_st+3]
    trip_ref_lst=list(trip_ref)
    trip_allel_lst=list(trip_ref)
    trip_allel_lst[trip_pos-1]=snp_allel

    if CDS_dir=='-':
        aa_pos=CDS_lg/3-(snpCDSpos-1)//3
        trip_ref_lst.reverse()
        trip_allel_lst.reverse()

        trip_ref_lst=map(char_cmp,trip_ref_lst)
        trip_allel_lst=map(char_cmp,trip_allel_lst)


    ref_cod=coden[''.join(trip_ref_lst)]
    allel_cod=coden[''.join(trip_allel_lst)]
	
    if ref_cod == allel_cod:
        snp_anno,snp_prot='CDS_synon',''
    else:
        snp_anno,snp_prot='CDS_nonSynon','%s%d%s'  % (ref_cod,aa_pos,allel_cod)

    return snp_anno,snp_prot


def find_cod_indel(str1,str2,str3,snp_pos,snp_allel):

    snpCDSpos=0
    CDS_ky=(str1,str2,str3)
    CDS_seq=CDS_seq_dic[CDS_ky]
    CDS_lg=len(CDS_seq)

    if CDS_lg%3!=0:
#       snp_anno,snp_prot='CDS_unanno',''
        aa_pos=0
        return aa_pos

    CDS_dir=CDS_dic[CDS_ky][0][2]
    for CDS_val in sorted(CDS_dic[CDS_ky]):
        if snp_pos < CDS_val[0]:
            break
        snpCDSpos+=min(snp_pos,CDS_val[1])-CDS_val[0]+1  ## the SNP position in the CDS

#   trip_st=(snpCDSpos-1)//3*3  ## the coden starting position in the CDS
#   trip_pos=snpCDSpos%3
    aa_pos=(snpCDSpos-1)//3+1  ## amino acid position in the protein


    if CDS_dir == '-' and snp_allel.find('-') != -1: ## if it is a deletion in reverse direction
        aa_pos=CDS_lg/3-(snpCDSpos+len(snp_allel)-1-1)//3
    elif CDS_dir == '-' and snp_allel.find('-') == -1:  ## if it is an insertion in reverse direction. Cautious here: an insertion in reverse direction need a location shift.
        aa_pos=CDS_lg/3-(snpCDSpos-2)//3

### cautious again: indels upstream start codon or downstream stop codon may be annotated mistakenly as frame-shift
	
    return int(aa_pos)


def trans_aa(snp_allel,ori):

    snpLg=len(snp_allel)
    aaTrans=''

    if ori == '-':
        snp_lst=list(snp_allel)
        snp_lst.reverse()
        snp_lst=map(char_cmp,snp_lst)
        snp_allel=''.join(snp_lst)


    for kk in range(0,snpLg,3):
        trip_cod=snp_allel[kk:kk+3]
        if trip_cod not in coden:
            aaTrans+='X'
        else:
            aaTrans+=coden[trip_cod]

    return aaTrans


#cds_dx,gene_dx=10,11
def anno_indel(snp_table):

    snp_file=open(snp_table,'r')
    out_file=open(snp_table+'.anno','w')
#   snp_file.readline()

    while True:
        lin=snp_file.readline().rstrip()
        if lin=='':
            break

        lst=lin.split('\t')
        contig,pos=lst[0].split('.')[0],int(lst[1])
        allel_ref,allel=lst[3],lst[4].split(',')

        try:
            allel.remove(allel_ref)
            allel_var=allel[0]
        except ValueError:
            allel_var=allel[0]


        anno,gene,geneIDlst,snp_dir,aa_chg,aa_chgPos=[],[],[],[],[],[]
        ed0,geneSymb0,geneID0=0,'',''
#       allel_key = lambda xx: allel_ref if '-' in allel_var else allel_var
#        allel_var = lambda allel_ref: allel_var*len(allel_ref) if '-' in allel_var
        if '-' in allel_var:
            allel_var = allel_var*len(allel_ref)

        if contig not in contg_dic:

            anno,gene,geneIDlst,snp_dir,aa_chg,aa_chgPos=['inter-gene'],[''],[''],['0'],[''],[''] # if the dic has no the contig, the SNP is inter-genic
            funx=anno[0]
            funx_gene=gene[0]
            funx_geneID=geneIDlst[0]
            map_dir=snp_dir[0]
            funx_aa=aa_chg[0]
            funx_aaPos=aa_chgPos[0]
            out_file.write('\t'.join(lst)+'\t'+funx+'\t'+funx_gene+'\t'+funx_geneID+'\t'+map_dir+'\t'+funx_aa+'\t'+funx_aaPos+'\n')


#           out_file.write('\t'.join(lst)+'\n')
#           out_file.write('%s\t%d\t%s\t%s\n' % (contig,pos,funx,funx_gene) )
            continue

#       sort_contgDic=sorted(contg_dic[contig])+[]]
        for contg_val in sorted(contg_dic[contig]):

            flag,flag0=0,0
            st,ed,ori=contg_val[0],contg_val[1],contg_val[2]
            geneID,geneSymb=contg_val[3],contg_val[4]

            if ed0 < pos < st: # intergenic

                the_anno,the_gene,the_geneID='inter-gene','',''

                if pos-ed0 <= 500 and ori == '+':
                    the_anno+='-3'
                    the_gene+=geneSymb0+'||'
                    the_geneID+=geneID0+'||'

                if pos-ed0 <= 1000 and ori == '-':
                    the_anno+='-5'
                    the_gene+=geneSymb0+'||'
                    the_geneID+=geneID0+'||'

                if st-pos  <= 1000 and ori == '+':
                    the_anno+='-5'
                    the_gene+=geneSymb+'||'
                    the_geneID+=geneID+'||'

                if st-pos <= 500 and ori == '-':
                    the_anno+='-3'
                    the_gene+=geneSymb+'||'
                    the_geneID+=geneID+'||'

                anno+=[cod[the_anno]]
                #anno+=[the_anno]
                gene+=[the_gene]
                geneIDlst+=[the_geneID]
                snp_dir+=['0']
                aa_chg+=['']
                aa_chgPos+=['']

                break # if intergenic,stop searching

            if st <= pos <= ed and geneID not in gene_dic: # the SNP is in some psudo-gene or non-coding RNA
                anno+=[cod['pseudo-gene|ncRNA']]
                #anno+=['pseudo-gene']
                gene+=[geneSymb]
                geneIDlst+=[geneID]
                snp_dir+=[ori]
                aa_chg+=['']
                aa_chgPos+=['']
                flag,flag0=1,1
				
            elif st <= pos <= ed and geneID in gene_dic: # in some gene

                gene_ed0=st-1
                sort_geneDic=sorted(gene_dic[geneID])+[[ed,ed,ori,'UTR','-']] # include the last UTR-fragment

                for geneID_val in sort_geneDic:
                    gene_st,gene_ed,gene_dir=geneID_val[0],geneID_val[1],geneID_val[2]
                    gene_reg,gene_Acc=geneID_val[3],geneID_val[4]

                    if gene_ed0 < pos < gene_st: # intronic

                        the_anno,the_gene,the_geneID='intron',geneSymb,geneID

                        if gene_st-pos  <= 2:
                            the_anno+='-splice-3'
                        if pos-gene_ed0 <= 2:
                            the_anno+='-splice-5'

                        anno+=[cod[the_anno]]
                        #anno+=[the_anno]
                        gene+=[the_gene]
                        geneIDlst+=[the_geneID]
                        snp_dir+=[gene_dir]
                        aa_chg+=['']
                        aa_chgPos+=['']

                        #break # if intronic, stop searching
						
                    if gene_st <= pos <= gene_ed: # exonic

                        the_anno,the_gene,the_geneID=gene_reg,geneSymb,geneID
						
                        if the_anno != 'CDS':
                            anno+=[cod[the_anno]]
                            #anno+=[the_anno]
                            gene+=[the_gene]
                            geneIDlst+=[the_geneID]
                            snp_dir+=[gene_dir]
                            aa_chg+=['']
                            aa_chgPos+=['']

                        elif the_anno == 'CDS' and len(allel_var)%3 == 0 and find_cod_indel(contig,geneID,gene_Acc,pos,allel_var) == find_cod_indel(contig,geneID,gene_Acc,pos+2,allel_var) and trans_aa(allel_var,gene_dir).find('Stp')==-1:
                            anno+=[cod['in-frame']]
                            gene+=[the_gene]
                            geneIDlst+=[the_geneID]
                            snp_dir+=[gene_dir]
                            aa_chg+=['']
                            aa_chgPos+=[find_cod_indel(contig,geneID,gene_Acc,pos,allel_var)]

                        elif the_anno == 'CDS' and len(allel_var)%3 == 0 and find_cod_indel(contig,geneID,gene_Acc,pos,allel_var) == find_cod_indel(contig,geneID,gene_Acc,pos+2,allel_var) and trans_aa(allel_var,gene_dir).find('Stp')!=-1:
                            anno+=[cod['in-frame-Stp']]
                            gene+=[the_gene]
                            geneIDlst+=[the_geneID]
                            snp_dir+=[gene_dir]
                            aa_chg+=['']
                            aa_chgPos+=[find_cod_indel(contig,geneID,gene_Acc,pos,allel_var)]

                        elif the_anno == 'CDS' and len(allel_var)%3 == 0 and find_cod_indel(contig,geneID,gene_Acc,pos,allel_var) != find_cod_indel(contig,geneID,gene_Acc,pos+2,allel_var):
                            anno+=[cod['in-frame']]
                            gene+=[the_gene]
                            geneIDlst+=[the_geneID]
                            snp_dir+=[gene_dir]
                            aa_chg+=['']
                            aa_chgPos+=[find_cod_indel(contig,geneID,gene_Acc,pos,allel_var)]

                        elif the_anno == 'CDS' and len(allel_var)%3 > 0:
                            anno+=[cod['frame-shift']]
                            gene+=[the_gene]
                            geneIDlst+=[the_geneID]
                            snp_dir+=[gene_dir]
                            aa_chg+=['']
                            aa_chgPos+=[find_cod_indel(contig,geneID,gene_Acc,pos,allel_var)]
							

                    gene_ed0=gene_ed

                flag,flag0=1,1
			
            elif pos > ed: # if the location is in the last inter-gene region of that contig,stop
                flag=0
					
            if flag==0 and flag0==1:
                break
					

            ed0=ed
            geneSymb0=geneSymb
            geneID0=geneID
	

        zipped=zip(anno,gene,geneIDlst,snp_dir,aa_chg,aa_chgPos)
        zipped.sort()

        if zipped!=[]:
            zipped0=zipped[0]
            funx=antiCod[zipped0[0]]
            funx_gene=zipped0[1]
            funx_geneID=zipped0[2]
            map_dir=zipped0[3]
            funx_aa=zipped0[4]
            funx_aaPos=zipped0[5]

        else:

            funx='inter-gene'
            funx_gene=''
            funx_geneID=''
            map_dir='0'
            funx_aa=''
            funx_aaPos=''
			

#       out_file.write('%s\t%d\t%s\t%s\t%s\t%s\n' % (contig,pos,funx,funx_gene,funx_geneID,map_dir) )
        out_file.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % ('\t'.join(lst),funx,funx_gene,funx_geneID,map_dir,funx_aa,str(funx_aaPos)) )


#       if lst[cds_dx].find('CDS')==-1:
#           lst[cds_dx]=funx
#           lst[gene_dx]=funx_gene
#           lst+=[map_dir]
#       elif lst[cds_dx].find('CDS')!=-1:
#           lst+=[map_dir]
#       out_file.write('\t'.join(lst)+'\n')		


    snp_file.close()
    out_file.close()



def indels(args):

    if args.contg:
        contg_filefile = args.contg
        create_contg_dic(contg_filefile)
#       print(contg_dic.keys())
    else:
        raise Exception("Please provide the contig dictionary file !")


    if args.gene:
        gene_filefile = args.gene
        create_gene_dic(gene_filefile)
    else:
        raise Exception("Please provide the gene dictionary file !")


    if args.CDS:
        CDS_filefile = args.CDS
        create_CDS_dic(CDS_filefile)
    else:
        raise Exception("Please provide the CDS dictionary file !")


    if args.CDS_seq:
        CDS_seq_filefile = args.CDS_seq
        create_CDSseq_dic(CDS_seq_filefile)
    else:
        raise Exception("Please provide the CDS-Seq dictionary file !")


    if args.mol:
        mol_filefile = args.mol
        create_mol_dic(mol_filefile)
    else:
        pass


    if args.table:
        tables = args.table
    else:
        raise Exception("Please provide mutation table file !")



    coden_dic()
    for ech_table in tables:
        anno_indel(ech_table)

