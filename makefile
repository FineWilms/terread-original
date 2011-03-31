#Makefile for GLOBPE on NEC SX4   f90 (from mrd)

# Specify any desired preprocessor or compiler flags here; -R2 for .L files.

FF = ifort
XFLAGS = -O -assume byterecl
LIBS = -L /tools/netcdf/3.6.0-p1/lib -lnetcdf
INC = -I /tools/netcdf/3.6.0-p1/include

OBJT = terread.o read_ht.o read1km.o read250.o read10km.o amap.o rwork.o ncdf.o \
       ccinterp.o latltoij_m.o setxyz_m.o xyzinfo_m.o newmpar_m.o \
       indices_m.o parm_m.o precis_m.o ind_m.o jimco_m.o jimcc_m.o \
       jim_utils.o nfft_m.o
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
.f90.o:
	$(FF) -c $(XFLAGS) $(INC) $<
.f.o:
	$(FF) -c $(XFLAGS) $(INC) $<

terread.o : ccinterp.o rwork.o
read250.o : ccinterp.o rwork.o
read1km.o : ccinterp.o rwork.o
read10km.o : ccinterp.o rwork.o
ccinterp.o : ccinterp.f90 setxyz_m.o xyzinfo_m.o latltoij_m.o newmpar_m.o indices_m.o
latltoij_m.o : latltoij_m.f90 xyzinfo_m.o newmpar_m.o
setxyz_m.o : setxyz_m.f90 newmpar_m.o indices_m.o parm_m.o precis_m.o ind_m.o xyzinfo_m.o jimco_m.o jimcc_m.o 
xyzinfo_m.o : xyzinfo_m.f90 precis_m.o
newmpar_m.o : newmpar_m.f90 
precis_m.o : precis_m.f90
indices_m.o : indices_m.f90
parm_m.o : parm_m.f90 precis_m.o 
ind_m.o : ind_m.f90 newmpar_m.o 
jimcc_m.o : jimcc_m.f90 parm_m.o precis_m.o 
jimco_m.o : jimco_m.f90 precis_m.o jim_utils.o nfft_m.o 
jim_utils.o : jim_utils.f90 precis_m.o 
nfft_m.o : nfft_m.f90 precis_m.o 


