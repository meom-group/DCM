NCDF = -I$(NETCDF_INCDIR) -I$(NETCDFF_INCDIR) $(NETCDF_LDFLAGS) $(NETCDFF_LDFLAGS)

#FFLAGS= -O2 $(NCDF) $(NC4)  -CB -fpe0 -ftrapuv -traceback -g  -fp-model precise
#FFLAGS= -O2 $(NCDF) $(NC4)  -CB  -traceback -g
FFLAGS= -O2  $(NCDF) $(NC4) -fp-model precise
F90=ifort
MPF90=mpif90

BINDIR = ../

VPATH = $(BINDIR)




mergefile4.exe: mergefile4.f90 modncfile.o modutil.o
	$(F90) mergefile4.f90 -o $(BINDIR)/mergefile4.exe modncfile.o modutil.o $(FFLAGS)

## Modules

modncfile.o: modncfile.f90 
	$(F90) -c  modncfile.f90  $(FFLAGS)

modutil.o: modutil.f90  modncfile.o
	$(F90) -c  modutil.f90 -o modutil.o  $(FFLAGS)

clean:
	rm -f  $(BINDIR)/*.exe *.mod *.o


	
