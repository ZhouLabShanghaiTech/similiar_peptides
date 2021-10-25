from __future__ import division
import sys,getopt
import pandas as pd
import re
from typing import Sequence

def Percentage(peptide_sequence):
    charged=len(re.findall(r'[EDRK]',peptide_sequence))
    polar=len(re.findall(r'[STQNYC]',peptide_sequence))
    PG=peptide_sequence.count('P')+peptide_sequence.count('G')
    aliphatic=len(re.findall(r'[MFWVILA]',peptide_sequence))
    return charged,polar,PG,aliphatic

argv=sys.argv[3:]
try:
    opts,args=getopt.getopt(argv,"w:",["--window="])
except getopt.GetoptError:
    print("ERROR! Usage: search.py <inputfile> <outputfile> -w <read_window length>")
    sys.exit(2)
for opt,arg in opts:
    if opt in ['-w','--window']:
        window=arg
        window=int(window)

proteome={}
name_list=""
c_thres=float(input('please input the minimum number of charged amino acids\n'))
pg_thres=float(input('please input the minimum number of Pro and Gly\n'))
hy_min,hy_max=(input('please input the range of hydrophobic amino acids\n').split('-'))
hy_max=float(hy_max)
hy_min=float(hy_min)

try:
    f=open(sys.argv[1], 'r')
    list1=f.readlines()
    for line in range(0,len(list1)):
        list1[line]=list1[line].strip('\n')
        i=0
        str=''.join(list1[line])
        if str.startswith('>'):
            i+=1
            name_list=str
            aasequence=""
        else:
            aasequence=aasequence+str
            proteome[name_list]=aasequence
      
finally:
    if f:
        f.close()

protein_ID=[]
protein_seq=[]
for key in proteome.keys():
    sequence=proteome[key]
    length=len(sequence)
    n=length-window
    if n>=0:
	    i=0
	    while i<n:
             fragment=sequence[i:(i+window)]
             i+=1
             charged,polar,PG,aliphatic=Percentage(fragment)
             c=charged/window
             a=aliphatic/window
             p=PG/window
             if c>=c_thres and p>=pg_thres:
                     if a>=hy_min and a<=hy_max:
                             key_vaule=key.split('/n')
                             protein_ID.append(key_vaule)
                             protein_seq.append(fragment)

d={'Uniprot ID':protein_ID,'AA sequence':protein_seq}
df1=pd.DataFrame(d)
df1.to_csv(sys.argv[2],index=False,header=False)
