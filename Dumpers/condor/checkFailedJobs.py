#!/usr/bin/python
import numpy as n
from ROOT import *
import sys, getopt
from array import array
from optparse import OptionParser
import operator
import os
import glob
import os.path
from os import path

if __name__ == '__main__':


  parser = OptionParser()
  parser.add_option(   "-i", "--inputFile",  dest="inputFile",   default="",      type="string",  help="inputFile" )
  parser.add_option(   "-d", "--dumper",     dest="dumper",      default="False", type="string",  help="dumper"    )

  (options, args) = parser.parse_args()  
  inputFile = options.inputFile
  dumper = 'False'
  if options.dumper!='' and options.dumper!='False':
     dumper = 'True'
  print "inputFile = ", inputFile
  print "dumpper   = ", dumper 

  f = open('arguments_resubmit.txt', 'w')
  with open(inputFile) as f_dump:
     data_dump = f_dump.read()
     lines_dump = data_dump.splitlines() 
     nMissingJobs = 0
     for pos,x in enumerate(lines_dump): 
        if dumper=='False':
          if path.exists(x.split()[1]+'_step2.root'):
            dim = os.path.getsize(x.split()[1]+'_step2.root')
            if dim<5999999: 
              print "Dim:",dim
              nMissingJobs+=1
              f.write(x+'\n')
          else: 
            #print x.split()[1]+'_step2.root'
            nMissingJobs+=1
            f.write(x+'\n')  
        else:
          if path.exists(x.split()[2]):
            dim = os.path.getsize(x.split()[2])
            if dim<14000000: 
              print "Dim:",dim
              nMissingJobs+=1
              f.write(x+'\n')
          else: 
            #print x.split()[1]+'_step2.root'
            nMissingJobs+=1
            f.write(x+'\n')  
        
  f.close()
  print "Missing jobs: ",nMissingJobs
