# for dev only ... normal use with fcm

#NCDF = -I$(NETCDF_INCDIR) -I$(NETCDFF_INCDIR) $(NETCDF_LDFLAGS) $(NETCDFF_LDFLAGS)
F90=ifort
MPF90=mpiifort

FFLAGS= -O3  $(NCDF) -fp-model precise   -lnetcdf -lnetcdff
#FFLAGS= -O  $(NCDF)   -assume byterecl -convert big_endian -CB -fpe0 -ftrapuv -traceback -g  -lnetcdf -lnetcdff

EXE = bdy_coord_create.exe bdy_mk_coordinates_from_file.exe  bdy_mk_coordinates_from_file_mpp.exe bdy_mk_mask.f90.exe

all : $(EXE)

bdy_coord_create.exe:  bdy_coord_create.f90
	$(F90) $(FFLAGS) bdy_coord_create.f90 -o bdy_coord_create.exe


bdy_mk_coordinates_from_file.exe:  bdy_mk_coordinates_from_file.f90
	$(F90) $(FFLAGS) bdy_mk_coordinates_from_file.f90 -o bdy_mk_coordinates_from_file.exe

bdy_mk_coordinates_from_file_mpp.exe:  bdy_mk_coordinates_from_file_mpp.f90
	$(MPF90) $(FFLAGS) bdy_mk_coordinates_from_file_mpp.f90 -o bdy_mk_coordinates_from_file_mpp.exe

bdy_mk_mask.f90.exe:  bdy_mk_mask.f90
	$(MPF90) $(FFLAGS)  bdy_mk_mask.f90 -o bdy_mk_mask.exe

clean:
	rm -rf *~ $(EXE)
