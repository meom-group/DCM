#!/usr/bin/python
# -*- coding: utf-8 -*-
# class = @Scalability experiment tools@

"""
Created on Sun Jun 25 19:18:08 2017

@author: Polina
"""

import os
import numpy as np
import pandas as pd
import glob
import sys, getopt


def usage(name):
   print ' '
   print 'USAGE: '+name+' -h  -n <script_name> -c  <cores>  -x <nxios_min>'
   print ' '
   print '  PURPOSE:'
   print '     This python script ends up with a bash script (run.sh) holding'
   print '     the ad-hoc command lines for submitting a series of simulations'
   print '     with different domain decomposition, hence total number of cores.'
   print '  '
   print '     It assumes that you have already in the local directory an elementary'
   print '     script launching NEMO and taking 4 online arguments: jpni, jpnj, jpnij, nxios'
   print '     e.g.  ./run_nemo_occigen_scal.sh  20 30 50 3 '
   print '    '
   print '     This script is called from dcmtk_scal_prep which preprocess the processor.layout'
   print '     file produced by MPP_PREP, writing temporary log files (one for a given domain'
   print '     decomposition) named log_<totalcore>'
   print ' '
   print '  OPTIONS:'
   print '     -h : Display this help message'
   print '     -n <script_name> : define nemo scalability script name. Default: ',script
   print '     -c <cores> : gives number of core per compute note. Default :',nc
   print '     -x <cores> : gives the minimum number of cores dedicated to'
   print '           xios_server.exe. Default : ',nxios_min
   print '     -X <cores> : gives the maximum number of cores dedicated to'
   print '           xios_server.exe. Default : ',nxios_max
   print ' '
   sys.exit()

def set_default():
   # set default
   global script
   global nc
   global nxios_min
   global nxios_max
   script = './run_nemo_occigen_scal.sh'
   nc = 28
   nxios_min = 5
   nxios_max = 10

def parse(argv,name):
   global script
   global nc
   global nxios_min
   global nxios_max
   try:
      opts, args = getopt.getopt(argv,"hn:c:x:X:",["help","scal_script=","cores_per_node=","xios_min_cores=","xios_max_cores"])
   except getopt.GetoptError:
      usage(name)
   for opt, arg in opts:
      if opt in ("-h", "--help"):
         usage(name)
      elif opt in ("-n", "--scal_script"):
         script = arg
      elif opt in ("-c", "--cores_per_node"):
         nc = int(arg)
      elif opt in ("-x", "--xios_min_cores"):
         nxios_min = int(arg)
      elif opt in ("-X", "--xios_max_cores"):
         nxios_max = int(arg)

def mkrun():
   # get current dir path
   cwd = os.getcwd()
   # set output path
   outpath = './'

   file_name = [] # log files
   file_name += glob.glob(cwd + '/log_*')
   file_name.sort() # Sort list of file names
   # allocate vars to be written
   jpni = []
   jpnj = []
   proc = []
   xios = []
   # read files
   for i in file_name:
      fl = pd.read_csv(i, delim_whitespace=True,header=None,names=['jpni','jpnj','jpi','jpj','jpij', 'proc','elim','sup'])
      for j in range(len(fl['proc'])):
          jpni.append(fl['jpni'][j])
          jpnj.append(fl['jpnj'][j])
          proc.append(fl['proc'][j])
          print fl['proc'][j], nc, j 
          xi = int(fl['proc'][j]/nc + 1)*nc - fl['proc'][j]
          if ( xi < nxios_min ):
             xi = xi + nc          
          if ( xi > nxios_max ):
             xi = nxios_max
          xios.append(xi)

   # write into file
   outfile = outpath+ 'run.sh'
   # Open file
   outf = []
   out = open((outfile),'w')
   out.write('#!/bin/bash'+'\n')
   for i in range(len(xios)):
       out.write(script+' '+str(jpni[i])+' '+str(jpnj[i])+\
       ' '+str(proc[i])+' '+str(xios[i])+'\n')
   out.close()

if __name__ == "__main__":
   set_default()
   if len(sys.argv) == 1:
      usage(sys.argv[0])
   parse(sys.argv[1:],os.path.basename(sys.argv[0]) )
   mkrun()
