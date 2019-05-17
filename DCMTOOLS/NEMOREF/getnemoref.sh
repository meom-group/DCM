#!/bin/bash

# note that you also need xios2 for compiling the code
#  This NEMO revision was compiled and ran successfully with
# xios rev 1587 of http://forge.ipsl.jussieu.fr/ioserver/svn/XIOS/dev/dev_olga
# to be compiled out of the DCM structure.
#  NEMO fcm files should indicate the xios root directory
#svn co http://forge.ipsl.jussieu.fr/nemo/svn/NEMO/trunk -r 10089 NEMO4
#svn co http://forge.ipsl.jussieu.fr/nemo/svn/NEMO/trunk -r 10374 NEMO4
# as of Jan 29, 2019 :
#svn co -r 10650 https://forge.ipsl.jussieu.fr/nemo/svn/NEMO/releases/release-4.0 NEMO4
# as od May 17, 2019
svn co -r 10992  https://forge.ipsl.jussieu.fr/nemo/svn/NEMO/releases/release-4.0 NEMO4
