#############################################################
F90=mpif90
OPTIONS    = -c -real-size 64 -O3 -openmp 
LOPTIONS   = -O2 -openmp 
##############################################################

objects = \
module_vars.o\
module_multidata.o\
module_mpi.o\
module_vars_pt.o\
module_SEM.o\
module_ibm.o\
ibm.o\
actuator.o\
shapes.o\
fdstag.o\
initial.o\
init_particle.o\
localparameters.o\
alloc_dom.o\
post.o\
flosol.o\
checkdt.o\
bounds.o\
bounds_keps.o\
sipsol.o\
convection.o\
diffusion.o\
newsolv_mg.o\
mgsolver.o\
wall_function.o\
log_law.o\
alloc_pt.o\
MPI_pt.o\
DeltaF_MLS.o\
LPT.o\
timesig.o\
weno.o\
energy.o\
press.o\
roughness_function.o\
rungek.o\
averaging.o\
eddyvis_smag.o\
eddyvis_wale.o\
eddyvis_1eqn.o\
eddyvis_keps.o\
exchange_bc.o\
exchangep.o\
exchangepp.o\
exchangesca.o\
exchange.o\
exchangeu.o\
exchangev.o\
exchangew.o\
exchange_phi.o\
bounds_lsm.o\
lsm.o\
SEM.o\
sediment.o

.SUFFIXES: .for

.for.o:
	$(F90) $(OPTIONS) -o $@ $<

3dFDM.exe: $(objects) 
	$(F90) $(objects) $(LOPTIONS) -o 3dFDM.exe \

clean:
	rm -rfv *.o *.mod 3dFDM_LPT.exe



alloc_dom.o : alloc_dom.for module_mpi.o module_vars.o module_multidata.o 
alloc_pt.o : alloc_pt.for module_vars_pt.o module_multidata.o module_mpi.o module_vars.o 
averaging.o : averaging.for module_vars.o module_multidata.o 
bounds.o : bounds.for ibm.o module_mpi.o module_multidata.o module_vars.o 
bounds_keps.o : bounds_keps.for module_multidata.o module_vars.o 
bounds_lsm.o : bounds_lsm.for module_multidata.o lsm.o module_vars.o 
checkdt.o : checkdt.for module_multidata.o module_mpi.o module_vars.o 
convection.o : convection.for module_multidata.o module_vars.o 
DeltaF_MLS.o : DeltaF_MLS.for module_multidata.o module_vars.o module_mpi.o
diffusion.o : diffusion.for module_multidata.o module_mpi.o module_vars.o 
eddyvis_1eqn.o : eddyvis_1eqn.for module_multidata.o module_vars.o 
eddyvis_smag.o : eddyvis_smag.for module_multidata.o module_vars.o 
eddyvis_wale.o : eddyvis_wale.for module_multidata.o module_vars.o 
eddyvis_keps.o : eddyvis_keps.for module_multidata.o module_vars.o 
energy.o : energy.for module_multidata.o module_mpi.o module_vars.o 
exchange_bc.o : exchange_bc.for module_vars.o module_mpi.o module_multidata.o 
exchange_bcphi.o : exchange_bcphi.for module_vars.o module_mpi.o module_multidata.o 
exchange.o : exchange.for exchange_phi.o module_mpi.o module_multidata.o module_vars.o 
exchangep.o : exchangep.for module_vars.o module_mpi.o module_multidata.o 
exchange_phi.o : exchange_phi.for module_vars.o module_mpi.o module_multidata.o 
exchangepp.o : exchangepp.for module_vars.o module_mpi.o module_multidata.o 
exchangesca.o : exchangesca.for module_vars.o module_mpi.o module_multidata.o 
exchangeu.o : exchangeu.for module_vars.o module_mpi.o module_multidata.o 
exchangev.o : exchangev.for module_vars.o module_mpi.o module_multidata.o 
exchangew.o : exchangew.for module_vars.o module_mpi.o module_multidata.o 
fdstag.o : fdstag.for module_vars.o module_mpi.o 
flosol.o : flosol.for module_vars_pt.o module_multidata.o module_mpi.o module_vars.o 
ibm.o : ibm.for module_ibm.o module_mpi.o module_multidata.o module_vars.o
actuator.o : ibm.for module_ibm.o module_mpi.o module_multidata.o module_vars.o  
initial.o : initial.for module_mpi.o module_multidata.o module_vars.o 
init_particle.o : init_particle.for module_vars_pt.o module_vars.o module_mpi.o module_multidata.o 
localparameters.o : localparameters.for module_multidata.o module_mpi.o module_vars.o 
LPT.o : LPT.for module_vars_pt.o module_vars.o module_mpi.o module_multidata.o 
lsm.o : lsm.for module_mpi.o module_multidata.o module_vars.o 
sediment.o: sediment.for module_multidata.o module_vars.o module_mpi.o
mgsolver.o : mgsolver.for module_multidata.o module_vars.o 
module_mpi.o : module_mpi.for 
module_multidata.o : module_multidata.for 
module_vars.o : module_vars.for 
module_vars_pt.o : module_vars_pt.for 
module_SEM.o     : module_SEM.for
module_ibm.o : module_ibm.for 
MPI_pt.o : MPI_pt.for module_vars_pt.o module_multidata.o module_mpi.o module_vars.o 
newsolv_mg.o : newsolv_mg.for module_multidata.o module_mpi.o module_vars.o 
post.o : post.for module_vars.o module_multidata.o 
press.o : press.for module_multidata.o module_mpi.o module_vars.o 
roughness_function.o : roughness_function.for module_multidata.o module_mpi.o module_vars.o 
rungek.o : rungek.for module_multidata.o module_mpi.o module_vars.o 
shapes.o : shapes.for module_mpi.o ibm.o module_multidata.o module_vars.o 
sipsol.o : sipsol.for module_mpi.o module_multidata.o module_vars.o 
timesig.o : timesig.for module_mpi.o module_vars.o module_multidata.o 
wall_function.o : wall_function.for module_multidata.o module_vars.o 
log_law.o : log_law.for module_multidata.o module_vars.o 
weno.o : weno.for module_multidata.o module_vars.o 
SEM.o : SEM.for module_multidata.o module_vars.o module_SEM.o module_mpi.o
