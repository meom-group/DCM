# Makefile for fortran code in RUNTOOLS/toolkit

FC = gfortran

BINDIR = ../../bin

#VPATH = $(BINDIR)

all: linreg datfinyyyy

linreg: linreg.f90 
	$(FC) linreg.f90 -o linreg

datfinyyyy: datfinyyyy.f90 
	$(FC) datfinyyyy.f90 -o datfinyyyy

clean:
	\rm -f *~ linreg datfinyyy 

install : all
	mv linreg datfinyyyy $(BINDIR)
	cd $(BINDIR) ; ln -sf $(RUNTOOLS)/toolkit/dcmtk_* ./
