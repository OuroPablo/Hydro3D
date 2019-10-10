!##########################################################################
        program fdstag
!##########################################################################
        use mpi
        use vars
        implicit none

        call init_parallelisation
     
        call read_mdmap

        call read_control

        call read_infodom

        call alloc_dom
!Mesh generation:
        call localparameters
!Allocation of the main flow variables
        call initial
!LSM: initialisation of the level set field and allocation of variables
	  IF (L_LSM .or. L_LSMbase)  CALL initial_LSM_3D_channel
!Initialisation of the flow field
        call initflowfield
!Rough bed :
        IF (LROUGH)  THEN								
          IF (.not.LRESTART) THEN
          call init_rough
          ELSE
          call rough_restart
         END IF
        END IF
!IBM: reads geom.cin, genetare the geometries, opens files, etc.
        if (LIMB) call imb_initial
!IBM: Assignation of the Lagrangian markers to sub-domains and CPUs
        if (LIMB) call PartLoc_Initial							
!Calculation of the initial flux in the simulation.
        call iniflux
!SCALAR: initial sediment concentration variables
	 if(LSCALAR)   call sediment_init

        if(.not.LRESTART) then
           if (time_averaging) then
              call update_mean
	      if (noise.gt.0.0) call add_noise(noise)
           end if
        end if
!LPT: initialisation of the particles distribution
	  if (LPT) then						
		if (myrank.eq.0) open(unit=202,file='LPT_particles.dat')
		call init_particle	
	  endif

        if ((solver.eq.2).and.(.not.L_LSM)) call coeff

        call MPI_BARRIER (MPI_COMM_WORLD,ierr)
        if(myrank.eq.0) then
           write (numfile,*) '============START ITERATIONS========='
           write (6,*) '============START ITERATIONS========='
        end if
!All fields are initialised, now let's iterate the NS eqs in time:
        call flosol
!End MPI parallelisation
        call end_parallelisation

        end program
!##########################################################################
