# MAKE A BIB FILE FROM a LATEX CODE
# IN THE FUTURE IMPLEMENT AUTOMATIC BIBTEX DOWNLOAD FROM ADS
# FOR NOW THE CODE WILL PRINT A LIST OF BIBCODES THAT YOU NEED TO INPUT MANUALLY ON ADS (PAPER FORM)
# TO DOWNLOAD A BIB FILE
# AUTHOR: HARISH VEDANTHAM
#
import sys
import re
import ads
my_token="Noc1OEKJDFiGUBMa0pFvp7SrMxRoPQ6e4veUqseT"
fname = sys.argv[1]
f = open(fname, "r")
lines = f.readlines()
in_header = 1
allcites = []
for i_line in range(len(lines)):
   line = lines[i_line]
   if "begin{document}" in line.replace(" ",""):
      in_header=0
   if in_header:
      continue
   start_inds = [m.start() for m in re.finditer('\\\cite',line)]
   for start_ind in start_inds:
      found=0
      i=0
      while(not found):
         found=line[start_ind+i]=="{"
         i+=1
      citestr_start = start_ind+i
      i=0
      found=0
      while(not found):
         found=line[start_ind+i]=="}"
         i+=1
      citestr_end = start_ind+i-1
      citestr = line[citestr_start:citestr_end].split(",")
      print ("line:%d, %s"%(i_line,citestr))
      allcites+=[x.strip(" ") for x in citestr if len(x)>1]
#
allcites_unique = list(set(allcites))
for cite in allcites_unique:
   print ("%s"%cite)
