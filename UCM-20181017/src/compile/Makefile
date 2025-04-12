#============================================================
#     Makefile for MUCM
#     original version : 2012-09-07
#     Last update      : 2013-02-16
#     by DOAN QUANG VAN
#============================================================
#FC      = ifort
#FCFLAGS =  -check bounds
#-check nopower -fpe0  -traceback -g  -O0

FC = gfortran
#FCFLAGS = -I${NETCDF}/include
FCFLAGS = -I/usr/local/include

OBJS =  module_sf_ucm_vf.o \
	module_params.o \
	module_vars.o \
	module_urb_ini.o \
	module_phys_funcs.o \
	module_io.o \
	module_setup.o \
	module_initialize.o \
	module_ra_kondo94.o \
	module_ra_driver.o \
	module_sf_slab.o \
	module_sf_ucm.o \
	module_sf_driver.o \
	module_pbl_driver.o \
	module_dyn_urban.o \
	module_dyn_driver.o \
	main.o

EXE = main

$(EXE) : $(OBJS)
	$(FC) -g -o $(EXE) $(FCFLAGS) $(OBJS) -L${NETCDF}/lib/ -lnetcdff #-I${NETCDF}/include

module_sf_ucm_vf.o: ../f90/module_sf_ucm_vf.f90
	$(FC) -c $(FCFLAGS) ../f90/module_sf_ucm_vf.f90

module_params.o: ../f90/module_params.f90
	$(FC) -c $(FCFLAGS) ../f90/module_params.f90

module_vars.o: ../f90/module_vars.f90 \
        ../f90/module_params.f90
	$(FC) -c $(FCFLAGS) ../f90/module_vars.f90

module_urb_ini.o: ../f90/module_urb_ini.f90
	$(FC) -c $(FCFLAGS) ../f90/module_urb_ini.f90

module_phys_funcs.o: ../f90/module_phys_funcs.f90 \
        ../f90/module_params.f90
	$(FC) -c $(FCFLAGS) ../f90/module_phys_funcs.f90

module_io.o: ../f90/module_io.f90 \
	../f90/module_params.f90
	$(FC) -c $(FCFLAGS) ../f90/module_io.f90

module_setup.o: ../f90/module_setup.f90 \
	../f90/module_vars.f90 \
	../f90/module_params.f90 \
	../f90/module_io.f90
	$(FC) -c $(FCFLAGS) ../f90/module_setup.f90

module_initialize.o: ../f90/module_initialize.f90 \
	../f90/module_params.f90 \
	../f90/module_urb_ini.f90
	$(FC) -c $(FCFLAGS) ../f90/module_initialize.f90

module_ra_kondo94.o: ../f90/module_ra_kondo94.f90 \
	../f90/module_params.f90
	$(FC) -c $(FCFLAGS) ../f90/module_ra_kondo94.f90

module_ra_driver.o: ../f90/module_ra_driver.f90 \
	../f90/module_ra_kondo94.f90 \
	../f90/module_params.f90
	$(FC) -c $(FCFLAGS) ../f90/module_ra_driver.f90

module_sf_slab.o: ../f90/module_sf_slab.f90 \
	../f90/module_params.f90 \
	../f90/module_phys_funcs.f90
	$(FC) -c $(FCFLAGS) ../f90/module_sf_slab.f90

module_sf_ucm.o: ../f90/module_sf_ucm.f90 \
	../f90/module_params.f90 \
	../f90/module_phys_funcs.f90
	$(FC) -c $(FCFLAGS) ../f90/module_sf_ucm.f90

module_sf_driver.o: ../f90/module_sf_driver.f90 \
	../f90/module_params.f90 \
	../f90/module_sf_slab.f90 \
	../f90/module_sf_ucm.f90
	$(FC) -c $(FCFLAGS) ../f90/module_sf_driver.f90

module_pbl_driver.o: ../f90/module_pbl_driver.f90 \
	../f90/module_params.f90
	$(FC) -c $(FCFLAGS) ../f90/module_pbl_driver.f90

module_dyn_urban.o: ../f90/module_dyn_urban.f90 \
	../f90/module_params.f90
	$(FC) -c $(FCFLAGS) ../f90/module_dyn_urban.f90

module_dyn_driver.o: ../f90/module_dyn_driver.f90 \
	../f90/module_dyn_urban.f90 \
	../f90/module_params.f90
	$(FC) -c $(FCFLAGS) ../f90/module_dyn_driver.f90

main.o: ../f90/main.f90 \
        ../f90/module_vars.f90 \
        ../f90/module_setup.f90 \
        ../f90/module_initialize.f90 \
        ../f90/module_ra_driver.f90 \
	../f90/module_sf_driver.f90 \
	../f90/module_pbl_driver.f90 \
        ../f90/module_dyn_driver.f90
	$(FC) -c $(FCFLAGS) ../f90/main.f90

clean:
	rm -f *.mod
	rm -f *.o
	rm -f *~
	rm -f ../*~
	rm -f ../../*~
	rm -f ../../main
