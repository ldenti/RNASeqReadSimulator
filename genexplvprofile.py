#!/usr/bin/env python3
'''
This file randomly assigns weights for each transcript, and gets the transcript statistics by a given transcript annotation file (BED File).

Usage:

genexplvprofile.py {OPTIONS} <BED-File> 

Options:

	-h		Print this message

	-e mu,sigma	Specify the mean and variance of the lognormial distribution used to assign expression levels. Default -4,4
	
	-f		Print the statistics only; do not assign expression levels.

Note:

	1. To get a good group information, the BED file is suggested to sort according to the chromosome name and start position. In Unix systems, use "sort -k 1,1 -k 2,2n in.BED > out.BED" to get a sorted version (out.BED) of the bed file (in.BED).  

	2. The weight is at the 8th column, if -f option is not specified. The expression level of each transcript (RPKM) can be calculated as column[8]*10^9/column[2]/sum(column[8]).

History:

	02/08/2012
		Initialization.
'''

import sys;
import pydoc;
import os;
import re;
import fileinput;
import random;

argvi=1;
mindist=50;
minscore=2;
mu=-4;
sigma=4;
assignexplv=True;

# Parse one line in count data
def parsebed(lines):
  fd=lines.strip().split('\t');
  if len(fd)!=12:
    return ['',-1,-1,0];
  if fd[10].endswith(','):
    fd[10]=fd[10][:-1];
  if fd[11].endswith(','):
    fd[11]=fd[11][:-1];
  seglen=[int(x) for x in fd[10].split(',')];
  segstart=[int(x) for x in fd[11].split(',')];
  #jstart=int(fd[1])+seglen[0]+1;
  #jend=int(fd[1])+segstart[1]+1;
  jstart=int(fd[1])+1; # start is 0-base; increase 1 to convert to 1-base
  jend=int(fd[2]);
  # jscore=int(fd[4]);
  #seg1=[jstart+segstart[i] for i in range(len(segstart))];
  #seg2=[jstart+segstart[i]+seglen[i]-1 for i in range(len(segstart))]; 
  # [seg1,seg2] are now 1-base inclusive
  return [fd[0],jstart,jend,fd[3],sum(seglen),fd[5],fd[9]];

allfile=[];

while argvi <(len(sys.argv)):
  if sys.argv[argvi]=="-h":
    print(pydoc.render_doc(sys.modules[__name__]),file=sys.stderr);
    sys.exit();
  elif sys.argv[argvi]=="-f":
    assignexplv=False;
  elif sys.argv[argvi]=="-e":
    ms=sys.argv[argvi+1].split(",");
    argvi=argvi+1;
    if len(ms)!=2:
      print('Error: incorrect parameter for -e.',file=sys.stderr);
      sys.exit(); 
    try:
      mu=float(ms[0]);
      sigma=float(ms[1]);
    except ValueError:
      print('Error: incorrect parameter for -e.',file=sys.stderr);
      sys.exit(); 
    print('Mean and variance for lognormal distribution: '+str(mu)+','+str(sigma),file=sys.stderr);
  else:
    allfile.append(sys.argv[argvi]);
  argvi=argvi+1;


allid={};

prevchr="";
prevrange=[0,0];
rangeid=0;

nline=0;

currentgene=[];
groupid=0;

print('#ID\tLength\tDir\tExons\tPosition\tGroupID\tNIsoformInGroup',end='');
if assignexplv==True:
  print('\tExplv');
else:
  print();

for lines in fileinput.input(allfile):
  nline=nline+1;
  pf=parsebed(lines);
  chrname=pf[0];jstart=pf[1];jend=pf[2];id=pf[3];length=pf[4];direction=pf[5];nexon=pf[6];
  if len(chrname)==0 and jstart<0:
    continue;
  if chrname!=prevchr or jstart-prevrange[1]>0:
    if len(prevchr)!=0:
      groupid=groupid+1;
      for item in currentgene:
        print(item[0]+"\t"+str(groupid)+"\t"+str(len(currentgene)),end='');
        if assignexplv==True:
          weight=random.lognormvariate(mu,sigma)*item[1];
          print("\t"+str(weight));
        else:
          print();
    prevrange[0]=jstart;
    prevrange[1]=jend;
    prevchr=chrname;
    rangeid=rangeid+1;
    currentgene=[];
  elif jstart<prevrange[0]:
    print('Warning: the range is not sorted at line '+str(nline),file=sys.stderr);
  else:
    if jend>prevrange[1]:
      prevrange[1]=jend;
  currentgene.append((id+"\t"+str(length)+"\t"+direction+"\t"+str(nexon)+"\t"+chrname+":"+str(jstart)+"-"+str(jend),length));


if len(prevchr)!=0:
  groupid=groupid+1;
  for item in currentgene:
    print(item[0]+"\t"+str(groupid)+"\t"+str(len(currentgene)),end='');
    if assignexplv==True:
      weight=random.lognormvariate(mu,sigma)*item[1];
      print("\t"+str(weight));
    else:
      print();
