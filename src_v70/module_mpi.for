!##########################################################################
        module mpi
!##########################################################################
        implicit none
        include 'mpif.h'

	  SAVE
        integer           :: nprocs     ! nr. of processors
        integer           :: myrank      ! nr of this processor
        integer           :: ierr
        integer           :: MPI_FLT  
        integer           :: status(MPI_STATUS_SIZE)
        real ( kind =8 )  :: wtime
        real ( kind =8 )  :: wtime2
!v70 new tag offsets. Done in collaboration with EPCC !Pablo
!Need to use these ones when using a non-intel compiler, e.g. Cray.
! offset1:  used in most exchange processes
! offset2 and 3: used in exchangepp, exchange_bc and exchange_bc_phi --> when periodic BC are used!

        integer, parameter :: tag_offset1=10000
        integer, parameter :: tag_offset2=500
        integer, parameter :: tag_offset3=250000


! In the older versions this always works:
!        integer, parameter :: tag_offset1=10**5 
!        integer, parameter :: tag_offset2=10**4
!        integer, parameter :: tag_offset3=10**8

        contains
!#########################################################################
        subroutine init_parallelisation
!#########################################################################
        implicit none

!  initialize MPI.
        call MPI_INIT ( ierr )
!  Get the number of processes.
        call MPI_COMM_SIZE ( MPI_COMM_WORLD, nprocs, ierr )
!  Get the individual process iD.
        call MPI_COMM_RANK ( MPI_COMM_WORLD, myrank, ierr )

        MPI_FLT   = MPI_DOUBLE_PRECISION

        if ( myrank == 0 ) then
           wtime = MPI_WTIME ( )
           write ( *, '(a)' ) ' '
           write ( *, '(a)' ) 'HELLO_MPI - Master process:'
           write ( *, '(a)' ) '  FORTRAN90/MPI version'
           write ( *, '(a)' ) ' '
           write ( *, '(a)' ) '  A 3D finite Difference MPI program.'
           write ( *, '(a)' ) ' '
           write ( *, '(a,i8)' ) '  The number of processes is ', nprocs
           write ( *, '(a)' ) ' '
        end if

        end subroutine init_parallelisation
!##########################################################################
        subroutine end_parallelisation
!##########################################################################
        implicit none

        if ( myrank == 0 ) then
           write ( *, '(a)' ) ' '
           write ( *, '(a)' ) '3D_FD_MPI - Master process:'
           write (*,'(a)')'Normal end of execution:"Goodbye, world!".'
           wtime = MPI_WTIME ( ) - wtime
           write ( *, '(a)' ) ' '
           write (*,'(a,g14.6,a)') 'Elapsed wall clock time = ',wtime,
     &'seconds.'
        end if

        call MPI_BARRIER (MPI_COMM_WORLD,ierr)

        call MPI_FINALIZE ( ierr )

        end subroutine end_parallelisation
!##########################################################################

        end module mpi
!##########################################################################
