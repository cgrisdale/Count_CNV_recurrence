#!/usr/bin/python
#Cameron Grisdale
#Dec 27, 2017

import sys
import re
from collections import Counter
import resource

#Compare output of Assign_CNVs_SingleFile.py to find gene-CNV recurrence, counts, etc.

def Read_in_cnv(File):
  '''Read in filtered CNV file, return list'''

  cnvl,cnvd,chrm,start,end,copyn,typea,tml,chrmd=[],{},'','','','','',[],{}

  with open(File, 'r') as f:

    fname=str(File).split('.')[2]

    for line in f:
      line=line.strip().split('\t')
      if int(line[5])!=0:
        chrm,start,end,copyn,typea,ngen,genes=line[0],line[1],line[2],line[3],line[4],line[5],line[6]
        tml=[fname,chrm,copyn,typea,ngen,genes]
        #Load biglist
        cnvl.append(tml)
      else:
        pass

  return cnvl

def Read_file_names(myfile):
  '''Read file names from file to be used to distinguish control from experimental data'''
  fnamelist,curfile=[],''
  with open(myfile, 'r') as f:
    for line in f:
      line=line.strip().split('.')[0]
      curfile=str(line)
      fnamelist.append(curfile)
  #print fnamelist
  return fnamelist

def Load_list(alist,gens):
  gens=gens.split(';')
  for y in gens:
    alist.append(y)
  return alist

def Get_stats(cnvs,expm,cntr):
  '''Go through CNVs to get recurrence of exp vs control'''
  #cnvs=list of lists: tml=[fname,chrm,copyn,typea,ngen,genes]
  exp,cnt=Read_file_names(expm),Read_file_names(cntr)
  expgenelist,cntgenelist,expexclusive,cntexclusive=[],[],[],[]

  for x in cnvs:
    fname,chrm,copyn,copyt,gcount,genes=x[0],x[1],x[2],x[3],x[4],x[5]

    if fname in exp:
      expgenelist=Load_list(expgenelist,genes)
    elif fname in cnt:
      cntgenelist=Load_list(cntgenelist,genes)
    else:
      sys.exit("File name not in experimental or control file name list; exiting")

  #Go through list of exp and cnt genes and check if they're in the other list or if they're exclusive to one
  for gene in expgenelist:
    if gene not in cntgenelist:
      expexclusive.append(gene)
  for gene in cntgenelist:
    if gene not in expgenelist:
      cntexclusive.append(gene)

  print "Total experimental:",len(expgenelist),"Total control:",len(cntgenelist),'\n'
  print "Unique experimental",len(Counter(expgenelist).keys()),"Unique control",len(Counter(cntgenelist).keys()),'\n'
  print "Experimental exclusive",len(expexclusive),"Control exclusive",len(cntexclusive),'\n'

  expcounter=Counter(expgenelist)
  cntcounter=Counter(cntgenelist)
  expexclcounter=Counter(expexclusive)
  cntexclcounter=Counter(cntexclusive)
  expreccount=Counter(expcounter.values())
  cntreccount=Counter(cntcounter.values())
  expexclreccount=Counter(expexclcounter.values())
  cntexclreccount=Counter(cntexclcounter.values())

  print "Number of recurrent Exp genes",expreccount
  print "Number of recurrent Cnt genes",cntreccount
  print "Number of recurrent Exp exclusive genes",expexclreccount
  print "Number of recurrent Cnt exclusive genes",cntexclreccount

  return expcounter,cntcounter,expexclcounter,cntexclcounter

def Countout(cd,cnv,outfile,outsample):
  '''Take Counter dict and match with CNV info list of lists to output gene,chrm,type,count'''
  cnvout,fout={},{}

  for k,v in cd.items(): #k=gene, v=count/recurrence
    curchr,curtype,prevchr,prevtype,cdcount='','','','',0

    for c in cnv: #list of lists
      fname,chrm,copyn,typea,ngen,genes1=c[0],c[1],c[2],c[3],c[4],c[5]
      genes=genes1.split(';') #make into list of genes instead of string of geneA;geneB;etc.
      types,fnames,temp=[],[],[]

      if ngen>0:

        if k in genes:
          curchr,curtype=chrm,typea

          if k in cnvout:
            types,fnames=cnvout[k][2],cnvout[k][3]
            types.append(typea)
            fnames.append(fname)
            temp=[cnvout[k][0],cnvout[k][1],types,fnames]
            cnvout[k]=temp
          else:
            cnvout[k]=[chrm,v,[typea],[fname]]

      #Keep names to be checked in next iteration
      #prevchr,prevtype=curchr,curtype

  print "cnvout dict:",len(cnvout)

  for r,t in cnvout.items():
    gls=t[2] #list of gains/losses
    gains,losses=Counter(gls)['gain'],Counter(gls)['loss']
    fout[r]=[t[0],t[1],gains,losses]
    #total=gains+losses
    #if total > t[1]:
    #  print r,t

  header="Gene"+"\t"+"Chrm"+"\t"+"Recurrence"+"\t"+"Gains"+"\t"+"Losses"+"\n"
  outfile.write(header)

  for l,m in fout.items():
    myline=l+'\t'+'\t'.join(str(n) for n in m)
    outfile.write(myline)
    outfile.write('\n')

  cnvc=0
  header1="Gene"+"\t"+"Chrm"+"\t"+"Sample"+"\t"+"Type"+"\n"
  outsample.write(header1)

  #output a line for each samplename for each gene, ie. 3 lines for GeneA if it occurs in 3 samples
  for w,e in cnvout.items():
    for x,y in zip(e[2],e[3]): #types,fnames
      if cnvc<10:
        print w,x,y
        cnvc+=1
      myline1=w+"\t"+x+"\t"+y #Gene type fname
      outsample.write(myline1)
      outsample.write('\n')

  return cnvout


if __name__ == "__main__":

  if len(sys.argv)>4:
    experimentals=sys.argv[1] #file names for "experimental" condition
    controls=sys.argv[2] #file names for "control" condition
    Files=sys.argv[3:] #CNV.outfiles

  else:
    sys.exit("Script works for multiple filtered CNV files only; exiting")

  CV,QL=[],[]
  #filename=str(Files)

  expfiles,cntfiles=Read_file_names(experimentals),Read_file_names(controls)

  #Read in multiple filtered CNV files, return list of lists and add to master list
  for i in range(len(Files)):
    filename=Files[i].split('.')[2]
    if filename in expfiles or filename in cntfiles: #it's one of the files I want to analyze
      cv=Read_in_cnv(Files[i])
      CV.append(cv)
    else:
      pass

  #Flatten list of lists of lists down one level to list of lists
  QL=[item for sublist in CV for item in sublist] #[[fname,chrm,copynumber,gain/loss,numgenes,geneA;geneB;geneC],[]]
  print "Number of CNVs examined: ",len(QL),'\n'


  expgenes,cntgenes,expexcl,cntexcl=Get_stats(QL,experimentals,controls)

  outf=open('CNV.merged.nonexcl.out.tsv', 'w')
  outs=open('CNV.merged.nonexcl.sample.tsv', 'w')
  gened=Countout(expgenes,QL,outf,outs)
  outf.close()
  outs.close()




'''

      #if any([n in cd for n in genes]): #are any genes present in count dict
      #  for i in genes:
      #    gcount=cd[i]


  print "Unique Exp genes",len(Counter(expgenes).keys())
  print "Unique Cnt genes",len(Counter(cntgenes).keys())
  print "Exp:",Counter(expgenes).most_common(10)
  print "Cnt:",Counter(cntgenes).most_common(10)
  expcounter=Counter(expgenes)
  cntcounter=Counter(cntgenes)
  expexclcounter=Counter(expexcl)
  cntexclcounter=Counter(cntexcl)
  print "Number of recurrent Exp genes",Counter(expcounter.values())
  print "Number of recurrent Cnt genes",Counter(cntcounter.values())
  print "Number of recurrent Exp exclusive genes",Counter(expexclcounter.values())
  print "Number of recurrent Cnt exclusive genes",Counter(cntexclcounter.values())
  print ''

  recgenes=[]
  for k,v in expexclcounter.items():
    if v>12:
      recgenes.append(k)
  print recgenes


  #cvoutstats='\t'.join(str(z) for z in CVstats)+'\n'
  #header="G_range_sum"+'\t'+"L_range_sum"+'\t'+"G_range_avg"+'\t'+"L_range_avg"+'\t'+"G_copy_avg"+'\t'+"L_copy_avg"+'\n'
  #outstats=open(filename+'.stats.tsv', 'w')
  #outstats.write(header)
  #outstats.write(cvoutstats)
  #outstats.close()


  outf=open('CNV.outfile.'+filename+'.tsv', 'w')
  y=0
  for l,m in c.items():
    #y+=1
    #print k,v
    #if y>9:
    #  sys.exit(0)
    for k,v in m.items():
      tmpx=[v[0],v[1][0],v[1][-1],v[2],v[3]]
      genec=len(v[4])
      genenm=';'.join(str(z) for z in v[4])
      tmpx.append(genec) #append gene count and string of gene names separated by ;
      tmpx.append(genenm)
      myline='\t'.join(str(y) for y in tmpx)
      outf.write(myline)
      outf.write('\n')
  outf.close()
'''






