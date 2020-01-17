#!/usr/bin/env python
'''
This file randomly assigns weights for each transcript, and gets the transcript statistics by a given transcript annotation file (BED File).

USAGE 

	genexplvprofile.py {OPTIONS} <BED-File|-> 

OPTIONS

	-h/--help\tPrint this message

	-e/--lognormal\tmu,sigma	Specify the mean and variance of the lognormal distribution used to assign expression levels. Default -4,4
        --geometric\tmu			Use geometric distribution with parameter mu instead of lognormal distribution to assign expression levels.
	
	-f/--statonly\tPrint the statistics only; do not assign expression levels.

NOTE

	1. To get a good group information, the BED file is suggested to sort according to the chromosome name and start position. In Unix systems, use "sort -k 1,1 -k 2,2n in.BED > out.BED" to get a sorted version (out.BED) of the bed file (in.BED).  

	2. The weight is at the 8th column, if -f option is not specified. The expression level of each transcript (RPKM) can be calculated as column[8]*10^9/column[2]/sum(column[8]).

HISTORY

	07/24/2012
	  Enable geometric distribution for expression level assignments. Require numpy package.

	02/16/2012
	  Run on Python 2.7

	02/08/2012
	  Initialization.
'''

import sys, os
import argparse
import pydoc
import fileinput
import random

import numpy

def parsebed(lines):
  # Parse one line in count data
  fd = lines.strip().split('\t')

  if len(fd) != 12:
    return ['',-1,-1,0,-1,-1,-1]

  jstart = int(fd[1])+1 # start is 0-base; increase 1 to convert to 1-base
  jend = int(fd[2])
  fd[10] = fd[10][:-1] if fd[10][-1] == ',' else fd[10]
  seglen = [int(x) for x in fd[10].split(',')]

  # fd[11] = fd[11][:-1] if fd[11][-1] == ',' else fd[11]
  # segstart = [int(x) for x in fd[11].split(',')]
  # jstart = int(fd[1])+seglen[0]+1
  # jend = int(fd[1])+segstart[1]+1
  # jscore = int(fd[4])
  # seg1 = [jstart+segstart[i] for i in range(len(segstart))]
  # seg2 = [jstart+segstart[i]+seglen[i]-1 for i in range(len(segstart))]
  # [seg1,seg2] are now 1-base inclusive

  return [fd[0],jstart,jend,fd[3],sum(seglen),fd[5],fd[9]]

def is_valid_geometric(parser, arg):
  try:
    arg = float(arg)
    if not 0.0 <= arg <= 1.0:
      raise ValueError
  except ValueError:
    parser.error("Geometric parameter prob must be a float value in range [0.0, 1.0]")
  return arg

def main():
  parser = argparse.ArgumentParser(description="Assign expression level to a set of transcripts.")

  # Optional arguments
  parser.add_argument("-f", "--statonly",
                      dest="assignexplv",
                      help="Print the statistics only; do not assign expression levels.",
                      required=False,
                      action="store_false",
                      default=True)

  parser.add_argument("-e", "--lognormal",
                      dest="lognormal",
                      help="Specify the mean mu and variance sigma of the lognormal distribution used to assign expression levels. Default -4.0 4.0",
                      required=False,
                      nargs=2,
                      metavar=('mu', 'sigma'),
                      type = float,
                      action="store",
                      default=[-4,4])

  parser.add_argument("-g", "--geometric",
                      dest="geometric",
                      help="Use geometric distribution with parameter prob (instead of lognormal) to assign expression levels. Default: None",
                      required=False,
                      metavar='prob',
                      type = lambda x: is_valid_geometric(parser, x),
                      action="store",
                      default=None)

  # Positional arguments
  parser.add_argument("bed_paths", nargs='+', help="Path to bed file")

  args = parser.parse_args()

  distype="lognormal"
  mu, sigma = args.lognormal
  p = args.geometric
  if p is not None:
    distype="geometric"
    print("Using geometric distribution with p={}".format(p), file=sys.stderr, end=' ')
  else:
    print("Using lognormal distribution with mu={} and sigma={}".format(mu, sigma), file=sys.stderr, end=' ')
  print("on {}.".format(','.join(args.bed_paths)), file=sys.stderr)

  print("#ID\tLength\tDir\tExons\tPosition\tGroupID\tNIsoformInGroup", end='')
  if args.assignexplv:
    print("\tExplv")
  else:
    print("")

  prevchr, prevstart, prevend = "", 0, 0

  rangeid = 0
  currentgene = []
  groupid = 0
  nline = 0
  for line in fileinput.input(args.bed_paths):
    nline += 1
    chrname, jstart, jend, idx, length, direction, nexon = parsebed(line)
    if len(chrname)==0 and jstart<0:
      print("Skipping line {} (malformed)".format(nline), file = sys.stderr)
      continue

    # A group is a set of overlapping isoforms from the same chromosome
    if chrname != prevchr or jstart > prevend:
      # We parse the current group
      if len(prevchr) != 0:
        groupid += 1
        for item in currentgene:
          print(item[0] + '\t' + str(groupid) + '\t' + str(len(currentgene)), end='')
          if args.assignexplv:
            if distype == "geometric":
              weight = numpy.random.geometric(mu)*item[1]
            else:
              weight = random.lognormvariate(mu,sigma)*item[1]
            print("\t" + str(weight))
          else:
            print("")
      prevstart = jstart
      prevend = jend
      prevchr = chrname
      rangeid = rangeid+1
      currentgene = []
    elif jstart < prevstart:
      print("Warning: the range is not sorted at line {}".format(nline), file=sys.stderr)
    else:
      if jend > prevend:
        prevend = jend
    currentgene.append((idx+"\t"+str(length)+"\t"+direction+"\t"+str(nexon)+"\t"+chrname+":"+str(jstart)+"-"+str(jend),length))


  # Here we parse the last group not parsed in the for
  if len(prevchr)!=0:
    groupid += 1;
    for item in currentgene:
      print(item[0]+"\t"+str(groupid)+"\t"+str(len(currentgene)),end='')
      if args.assignexplv:
        if distype=="geometric":
          weight=numpy.random.geometric(mu)*item[1]
        else:
          weight=random.lognormvariate(mu,sigma)*item[1]
        print("\t"+str(weight))
      else:
        print("")

if __name__ == "__main__":
  main()
