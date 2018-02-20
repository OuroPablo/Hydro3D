!=======================================================================
!    			LAGRANGIAN PARTICLE TRACKING
!			Bruño Fraga Bugallo
!			Cardiff 2013-2015
!=======================================================================
!######################################################################
      SUBROUTINE MPI_pt
!	Distributes the particles among the processors
!######################################################################
      use vars
      use mpi
      use multidata
	use vars_pt

      implicit none

      integer :: l,o,ii,nxdom,nydom,nzdom,ntot
	real :: lx,ly,lz
	integer,allocatable,dimension(:)::	lpt_proc,lpt_block
	integer,allocatable,dimension(:)::	id_MPI,id_MPI_loc
      double precision, allocatable, dimension(:):: X_MPI,Y_MPI,Z_MPI
      double precision, allocatable, dimension(:):: X_MPI_loc,Y_MPI_loc
      double precision, allocatable, dimension(:):: Z_MPI_loc
      double precision, allocatable, dimension(:):: U_MPI,V_MPI,W_MPI
      double precision, allocatable, dimension(:):: U_MPI_loc,V_MPI_loc
      double precision, allocatable, dimension(:):: W_MPI_loc
      double precision, allocatable, dimension(:):: dp_MPI,dp_MPI_loc

	call MPI_BARRIER (MPI_COMM_WORLD,ierr)
      call MPI_BCAST(np,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

	ntot=np*nprocs

	allocate(X_MPI(ntot),Y_MPI(ntot),Z_MPI(ntot))
	allocate(U_MPI(ntot),V_MPI(ntot),W_MPI(ntot))
	allocate(dp_MPI(ntot),id_MPI(ntot))
	allocate(X_MPI_loc(np),Y_MPI_loc(np),Z_MPI_loc(np))
	allocate(U_MPI_loc(np),V_MPI_loc(np),W_MPI_loc(np))
	allocate(dp_MPI_loc(np),id_MPI_loc(np))
	allocate(lpt_block(np),lpt_proc(np))

      IF (Myrank.eq.0) THEN

	X_MPI = -1.d12

	lx=xen-xst
	ly=yen-yst
	lz=zen-zst

	ptsinproc=0  ; lpt_block=0  ; lpt_proc=0

		do l=1,np

			nxdom=INT((xp_pt(l)-xst-1d-12)/(lx/idom))				!Domain where point belongs in i

			nydom=INT((yp_pt(l)-yst-1d-12)/(ly/jdom))				!Domain where point belongs in j

			nzdom=INT((zp_pt(l)-zst-1d-12)/(lz/kdom))				!Domain where point belongs in k

			lpt_block(l)=idom*jdom*(nzdom)+idom*(nydom)
     &					+nxdom					!Domain where point belongs (0 to num_dom)

			lpt_proc(l)=dom_ad(lpt_block(l))+1 				!Processor where point belongs (1 to num_procs+1)

!			ptsinblk(lpt_block(l)+1)=ptsinblk(lpt_block(l)+1)+1	!array which tells how many points are in each domain

			ptsinproc(lpt_proc(l))=ptsinproc(lpt_proc(l))+1		!array which tells how many points are in each processor


!		write(6,*)lpt_block(l)
!     		write(6,*)xp_pt(l)
!     		write(6,*)(lx/idom)*(nxdom)
!		write(6,*)'======================================='

		enddo

		do o=1,nprocs
!		write(6,*)'1',o,ptsinproc(o)
	  	   ii=1
		   DO l=1,np
			if(lpt_proc(l).eq.(o)) then
				X_MPI(ii+np*(o-1))=xp_pt(l)				!superarrays
				Y_MPI(ii+np*(o-1))=yp_pt(l)
				Z_MPI(ii+np*(o-1))=zp_pt(l)

				U_MPI(ii+np*(o-1))=uop_pt(l)
				V_MPI(ii+np*(o-1))=vop_pt(l)
				W_MPI(ii+np*(o-1))=wop_pt(l)

				dp_MPI(ii+np*(o-1))=dp_pt(l)

				id_MPI(ii+np*(o-1))=lpt_block(l)

				ii=ii+1
			endif
		   ENDDO
		enddo

	ENDIF

      call MPI_BARRIER (MPI_COMM_WORLD,ierr)

        call MPI_BCAST(ptsinproc,nprocs,MPI_INTEGER,0,
     &			MPI_COMM_WORLD,ierr)

        call MPI_SCATTER(ptsinproc,1,MPI_INTEGER,np_loc,1,
     &		MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

!	if (np_loc.eq.0) then
!		np_loc=1
!	endif


!      call MPI_BARRIER (MPI_COMM_WORLD,ierr)

        call MPI_SCATTER(X_MPI,np,MPI_DOUBLE_PRECISION,X_MPI_loc,		!location
     &		np,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_SCATTER(Y_MPI,np,MPI_DOUBLE_PRECISION,Y_MPI_loc,
     &		np,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_SCATTER(Z_MPI,np,MPI_DOUBLE_PRECISION,Z_MPI_loc,
     &		np,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

        call MPI_SCATTER(U_MPI,np,MPI_DOUBLE_PRECISION,U_MPI_loc,		!velocities (of prev tstep)
     &		np,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_SCATTER(V_MPI,np,MPI_DOUBLE_PRECISION,V_MPI_loc,
     &		np,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_SCATTER(W_MPI,np,MPI_DOUBLE_PRECISION,W_MPI_loc,
     &		np,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

        call MPI_SCATTER(dp_MPI,np,MPI_DOUBLE_PRECISION,dp_MPI_loc,	!diameter
     &		np,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

        call MPI_SCATTER(id_MPI,np,MPI_INTEGER,id_MPI_loc,			!block to which they belong
     &		np,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)


	if (np_loc.gt.0) then

      allocate (xp_loc(np_loc),yp_loc(np_loc),zp_loc(np_loc))
      allocate (uop_loc(np_loc),vop_loc(np_loc),wop_loc(np_loc))
	allocate (dp_loc(np_loc))
!	allocate (Fsu(np_loc),Fau(np_loc),Fdu(np_loc),Flu(np_loc))
!	allocate (Fsv(np_loc),Fav(np_loc),Fdv(np_loc),Flv(np_loc))
!	allocate (Faw(np_loc),Fdw(np_loc),Flw(np_loc))
!	allocate (Fgw(np_loc),Fsw(np_loc))
!	allocate (ipux(np_loc),jpvy(np_loc),kpwz(np_loc))
	allocate (id(np_loc))
!	allocate (dh_acumu(np_loc),dh_acumv(np_loc),dh_acumw(np_loc))


		ii=0
		do l=1,np
			if (X_MPI_loc(l).ne.-1.d12) then	
				ii=ii+1		
				xp_loc(l)=X_MPI_loc(l)
				yp_loc(l)=Y_MPI_loc(l)
				zp_loc(l)=Z_MPI_loc(l)

				uop_loc(l)=U_MPI_loc(l)
				vop_loc(l)=V_MPI_loc(l)
				wop_loc(l)=W_MPI_loc(l)

				dp_loc(l)=dp_MPI_loc(l)

				id(l)=id_MPI_loc(l)
			endif
		enddo
		if (ii.ne.np_loc) then
			write(6,*)'MPI ERROR in proc:',myrank
			write(6,*)'np_loc=',np_loc,'=/=',ii
			stop
		endif

	endif

		deallocate(lpt_proc,lpt_block)
		deallocate(X_MPI,Y_MPI,Z_MPI)
		deallocate(U_MPI,V_MPI,W_MPI)
		deallocate(X_MPI_loc,Y_MPI_loc,Z_MPI_loc)
		deallocate(U_MPI_loc,V_MPI_loc,W_MPI_loc)
		deallocate(dp_MPI,dp_MPI_loc)
		deallocate(id_MPI,id_MPI_loc)

      RETURN
      END
