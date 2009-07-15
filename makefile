#Makefile for GLOBPE on NEC SX4   f90 (from mrd)

# Specify any desired preprocessor or compiler flags here; -R2 for .L files.

FF = ifort
XFLAGS = -O 
LIBS = -L /tools/netcdf/3.6.0-p1/lib -lnetcdf
INC = -I /tools/netcdf/3.6.0-p1/include

OBJT = terread.o read_ht.o setxyz.o jimcc.o latltoij.o read1km.o read250.o read10km.o amap.o ncdf.o
OBJV = veg.o setxyz.o jimcc.o latltoij.o ncdf.o

OBJr = testr.o SUBR_native_4byte_real.o

terread :$(OBJT)
	$(FF) $(XFLAGS) $(OBJT) $(LIBS) -o terread
veg :$(OBJV)
	$(FF) $(XFLAGS) $(OBJV) $(LIBS) -o veg.80
georead :$(OBJS)
	$(FF) $(XFLAGS) $(OBJS) $(LIBS) -o georead.new
testr : $(OBJr)
	$(FF) $(XFLAGS) $(OBJr) -o testr
read_dem : read_dem.o
	$(FF) $(XFLAGS) read_dem.f -o read_dem

# This section gives the rules for building object modules.

.SUFFIXES:.f90
ecoread.o:
	$(FF) -c $(XFLAGS) $(INC) -assume byterecl ecoread.f90
.f90.o:
	$(FF) -c $(XFLAGS) $(INC) $<
.f.o:
	$(FF) -c $(XFLAGS) $(INC) $<

a.o georead.o jimcc.o jimco.o latltoij.o nterread.o read10km.o read1km.o read250.o read_dem.o read_veg.o setxyz.o terread.o topfilt.o topgencc.o tt.o veg.o veg_read.o :  newmpar.h


