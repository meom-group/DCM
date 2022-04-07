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
#svn co -r 10992  https://forge.ipsl.jussieu.fr/nemo/svn/NEMO/releases/release-4.0 NEMO4
# as od May 18, 2019
#svn co -r 10997  https://forge.ipsl.jussieu.fr/nemo/svn/NEMO/releases/release-4.0 NEMO4
# as od May 23, 2019
#svn co -r 11040  https://forge.ipsl.jussieu.fr/nemo/svn/NEMO/releases/release-4.0 NEMO4
# as of June,4  2019
#svn co -r 11075  https://forge.ipsl.jussieu.fr/nemo/svn/NEMO/releases/release-4.0 NEMO4
# as of November, 14 2019
#svn co -r 11902  https://forge.ipsl.jussieu.fr/nemo/svn/NEMO/releases/release-4.0.1 NEMO4
# as of March, 25 2020
#svn co -r 12604  https://forge.ipsl.jussieu.fr/nemo/svn/NEMO/releases/release-4.0.1 NEMO4
# as of March, 27 2020
#svn co -r 12591  https://forge.ipsl.jussieu.fr/nemo/svn/NEMO/releases/r4.0/r4.0.2 NEMO4
# as of February,8 2021
#svn co -r 14325  https://forge.ipsl.jussieu.fr/nemo/svn/NEMO/releases/r4.0/r4.0.5 NEMO4
# as of March,15 2021
#svn co -r 14608  https://forge.ipsl.jussieu.fr/nemo/svn/NEMO/releases/r4.0/r4.0.6 NEMO4
# for 4.2RC
#svn co -r 15299 https://forge.ipsl.jussieu.fr/nemo/svn/NEMO/trunk NEMO4
######  4.2.x (gitlab) starts here ###############
# as of April,7 2022
commit=26911cc471c9316f7a67495d4fd544dce35b758d
git clone git@forge.nemo-ocean.eu:nemo/nemo.git NEMO4
cd NEMO4 ; git checkout $commit
