# for dev only ... normal use with fcm

NCDF = -I$(NETCDF_INCDIR) -I$(NETCDFF_INCDIR) $(NETCDF_LDFLAGS) $(NETCDFF_LDFLAGS)
F90=ifort

#FFLAGS= -O  $(NCDF) -fp-model precise 
FFLAGS= -O  $(NCDF)  -assume byterecl -convert big_endian -CB -fpe0 -ftrapuv -traceback -g 



bdy_coord_create:  bdy_coord_create.f90
	$(F90) $(FFLAGS) bdy_coord_create.f90 -o bdy_coord_create
