!##########################################################################
        subroutine calmas
!##########################################################################
!
! Calculate divergence free condition
!
        use vars
        use mpi
        use multidata
        implicit none
        integer i,j,k,ib
        double precision fact,fact1,buffer_rmax

        rmax=0.d0
        if(differencing.eq.2) then
        do ib=1,nbp
           fact=0.d0 ; fact1=0.d0
           do  k=dom(ib)%ksp,dom(ib)%kep
              do  i=dom(ib)%isp,dom(ib)%iep
                 do j=dom(ib)%jsp,dom(ib)%jep
                    fact=( 
     &(-dom(ib)%u(i+1,j,k)+27.0*dom(ib)%u(i,j,k)-
     &27.0*dom(ib)%u(i-1,j,k)+dom(ib)%u(i-2,j,k))/(24.0*dom(ib)%dx)+
     &(-dom(ib)%v(i,j+1,k)+27.0*dom(ib)%v(i,j,k)-
     &27.0*dom(ib)%v(i,j-1,k)+dom(ib)%v(i,j-2,k))/(24.0*dom(ib)%dy)+
     &(-dom(ib)%w(i,j,k+1)+27.0*dom(ib)%w(i,j,k)-
     &27.0*dom(ib)%w(i,j,k-1)+dom(ib)%w(i,j,k-2))/(24.0*dom(ib)%dz))
                    dom(ib)%su(i,j,k)=fact/dt

                    fact1=fact*dom(ib)%dx*dom(ib)%dy*dom(ib)%dz
                    resor=max(rmax,abs(fact1))
                    rmax=max(resor,abs(fact1))
                 end do
              end do
           end do
        end do

        else !2nd CDS or WENO:

        do ib=1,nbp
           fact=0.d0 ;  fact1=0.d0
           do  k=dom(ib)%ksp,dom(ib)%kep
              do  i=dom(ib)%isp,dom(ib)%iep
                 do j=dom(ib)%jsp,dom(ib)%jep

                       fact=( 
     & (dom(ib)%u(i,j,k)-dom(ib)%u(i-1,j,k))/dom(ib)%dx+
     & (dom(ib)%v(i,j,k)-dom(ib)%v(i,j-1,k))/dom(ib)%dy+
     & (dom(ib)%w(i,j,k)-dom(ib)%w(i,j,k-1))/dom(ib)%dz)

 			  if (L_LSM) then
                       dom(ib)%sup(i,j,k)=
     & fact*dom(ib)%dx*dom(ib)%dy*dom(ib)%dz/dt

			  else
                       dom(ib)%su(i,j,k)=fact/dt
			  end if

                    fact1=fact*dom(ib)%dx*dom(ib)%dy*dom(ib)%dz
                    resor=max(rmax,abs(fact1))
                    rmax=max(resor,abs(fact1))

                 end do
              end do
           end do
        end do
        end if

        buffer_rmax = rmax

        call MPI_ALLREDUCE(buffer_rmax,rmax,1,MPI_FLT,MPI_MAX,
     &MPI_COMM_WORLD,ierr)

        end subroutine calmas
!##########################################################################
        subroutine calvel
!##########################################################################
!
!Correct velocity field with pressure gradient
!
        use vars
        use mpi
        use multidata
        implicit none
        integer i,j,k,ib

        do ib=1,nbp

           do k=dom(ib)%ksu,dom(ib)%keu
              do i=dom(ib)%isu,dom(ib)%ieu
                 do j=dom(ib)%jsu,dom(ib)%jeu
                   if (L_LSM) then
                     dom(ib)%u(i,j,k)=(dom(ib)%ustar(i,j,k)-dt*alfapr*
     &(dom(ib)%p(i+1,j,k)-dom(ib)%p(i,j,k))/dom(ib)%dx/
     &(0.5*(dom(ib)%dens(i,j,k)+dom(ib)%dens(i+1,j,k))))
		       else
                     dom(ib)%u(i,j,k)=(dom(ib)%ustar(i,j,k)-dt*alfapr*
     & (dom(ib)%p(i+1,j,k)-dom(ib)%p(i,j,k))/dom(ib)%dx)
		       end if
                 end do
              end do
           end do

           do k=dom(ib)%ksv,dom(ib)%kev
              do i=dom(ib)%isv,dom(ib)%iev
                 do j=dom(ib)%jsv,dom(ib)%jev
                   if (L_LSM) then
                     dom(ib)%v(i,j,k)=(dom(ib)%vstar(i,j,k)-dt*alfapr*
     & (dom(ib)%p(i,j+1,k)-dom(ib)%p(i,j,k))/dom(ib)%dy/
     & (0.5*(dom(ib)%dens(i,j,k)+dom(ib)%dens(i,j+1,k))))

			 else
                     dom(ib)%v(i,j,k)=(dom(ib)%vstar(i,j,k)-dt*alfapr*
     & (dom(ib)%p(i,j+1,k)-dom(ib)%p(i,j,k))/dom(ib)%dy)
			 end if
                 end do
              end do
           end do

           do k=dom(ib)%ksw,dom(ib)%kew
              do i=dom(ib)%isw,dom(ib)%iew
                 do j=dom(ib)%jsw,dom(ib)%jew
                   if (L_LSM) then
                     dom(ib)%w(i,j,k)=(dom(ib)%wstar(i,j,k)-dt*alfapr*
     & (dom(ib)%p(i,j,k+1)-dom(ib)%p(i,j,k))/dom(ib)%dz/
     & (0.5*(dom(ib)%dens(i,j,k)+dom(ib)%dens(i,j,k+1))))

			 else
                     dom(ib)%w(i,j,k)=(dom(ib)%wstar(i,j,k)-dt*alfapr*
     & (dom(ib)%p(i,j,k+1)-dom(ib)%p(i,j,k))/dom(ib)%dz)
			 end if
                 end do
             end do
           end do

        end do

        call exchange(1)
        call exchange(2)
        call exchange(3)

        end subroutine calvel
!##########################################################################
        subroutine pbound
!##########################################################################
        use vars
        use multidata
        implicit none
        integer i,j,k,ib,isp,iep,jsp,jep,ksp,kep

        do ib=1,nbp

           isp=dom(ib)%isp; iep=dom(ib)%iep
           jsp=dom(ib)%jsp; jep=dom(ib)%jep
           ksp=dom(ib)%ksp; kep=dom(ib)%kep
!...................... west and east............................
           if (dom(ib)%iprev.lt.0 .and. dom(ib)%bc_west.ne.5) then
              do k=ksp-1,kep+1 
                 do j=jsp-1,jep+1
                    dom(ib)%p(isp-1,j,k)  =dom(ib)%p(isp,j,k)
                 end do
              end do
           end if
           if (dom(ib)%inext.lt.0 .and. dom(ib)%bc_east.ne.5) then
              do k=ksp-1,kep+1
                 do j=jsp-1,jep+1
                    dom(ib)%p(iep+1,j,k)  =dom(ib)%p(iep,j,k)
                 end do
              end do
           end if
!.....................south and north.........................
           if (dom(ib)%jprev.lt.0 .and. dom(ib)%bc_south.ne.5) then
              do k=ksp-1,kep+1
                 do i=isp-1,iep+1
                    dom(ib)%p(i,jsp-1,k)  =dom(ib)%p(i,jsp,k)
                 end do
              end do
           end if
           if (dom(ib)%jnext.lt.0 .and. dom(ib)%bc_north.ne.5) then
              do k=ksp-1,kep+1
                 do i=isp-1,iep+1
                    dom(ib)%p(i,jep+1,k)  =dom(ib)%p(i,jep,k)
                 end do
              end do
           end if
!.....................bottom and top.........................
           if (dom(ib)%kprev.lt.0 .and. dom(ib)%bc_bottom.ne.5) then
              do j=jsp-1,jep+1
                 do i=isp-1,iep+1
                    dom(ib)%p(i,j,ksp-1)  =dom(ib)%p(i,j,ksp)
                 end do
              end do
           end if
           if (dom(ib)%knext.lt.0 .and. dom(ib)%bc_top.ne.5) then
              do j=jsp-1,jep+1
                 do i=isp-1,iep+1
                    dom(ib)%p(i,j,kep+1)  =dom(ib)%p(i,j,kep)
                 end do
              end do
           end if

        end do

        end subroutine pbound
!##########################################################################
        subroutine pressure_forcing
!##########################################################################
        use vars
        use mpi
        use multidata
        implicit none
        integer i,j,k,ib,isp,jspr,jepr,kspr,kepr,ispr,iepr
        double precision qstpp,fakfor,flwsum_loc,flwsum,A_loc,A
        double precision,dimension(2) :: sndbuffer,recbuffer

        MPI_FLT = MPI_DOUBLE_PRECISION

        flwsum_loc = 0.0
        A_loc      = 0.0

        do ib=1,nbp
!           ispr=pl+1; iepr=dom(ib)%ttc_i-pl
           jspr=pl+1; jepr=dom(ib)%ttc_j-pl
           kspr=pl+1; kepr=dom(ib)%ttc_k-pl

           if(dom(ib)%iprev.lt.0) then						!west
              do j=jspr,jepr
                 do k=kspr,kepr
                   isp=dom(ib)%isu-1
			 if (L_LSM) then
			   if (dom(ib)%phi(isp,j,k) .ge. 0.d0) then
                         flwsum_loc=flwsum_loc+
     & dom(ib)%u(isp,j,k)*dom(ib)%dy*dom(ib)%dz*
     & (dom(ib)%dens(isp,j,k)/densl)
                         A_loc=A_loc+dom(ib)%dy*dom(ib)%dz*
     & (dom(ib)%dens(isp,j,k)/densl)
			   end if
			 else if (L_LSMbase) then
			   if (dom(ib)%z(k) .le. length) then
                         flwsum_loc=flwsum_loc+
     & dom(ib)%u(isp,j,k)*dom(ib)%dy*dom(ib)%dz
                         A_loc=A_loc+dom(ib)%dy*dom(ib)%dz
			   end if
			 else
                     flwsum_loc=flwsum_loc+
     & dom(ib)%u(isp,j,k)*dom(ib)%dy*dom(ib)%dz
                     A_loc=A_loc+dom(ib)%dy*dom(ib)%dz
			 end if
                 end do
              end do
           end if
        end do

!.....sum massflux and area over all processors
        sndbuffer(1) = flwsum_loc
        sndbuffer(2) = A_loc

        call MPI_ALLREDUCE(sndbuffer,recbuffer,2,
     &                   MPI_FLT,MPI_SUM,MPI_COMM_WORLD,ierr)

        flwsum = recbuffer(1)
        A      = recbuffer(2)

! --- Calculate and store forcing term --------------------------------
        qstpp = flwsum/A
        fakfor = 0.3
! Whatever the flow is lost, is added back again
! 	  Q(t)=Q(t-1)+0.3?  *(U(0) +U(t-1)-2*U(t))/dt
        forcn=forcn+fakfor*(qzero+qstpn-2.0*qstpp)/dt
        qstpn = qstpp


        if(myrank.eq.0) write (numfile1,'(4F15.6)') 
     &   ctime,forcn,qstpp,flwsum

        end subroutine pressure_forcing
!##########################################################################
        subroutine pressure_1sweep
!##########################################################################
        use vars
        use mpi
        use multidata
        implicit none
        integer :: i,j,k,l,m
        integer :: ib,isp,iep,jsp,jep,ksp,kep
        double precision :: ppref,buffer_ppref
        real ( kind =8 )  :: wtimedum

        MPI_FLT   = MPI_DOUBLE_PRECISION

        nsweep=0
        iter=0
        nsweeps=0

        do ib=1,nbp
           dom(ib)%ap=0.0
           dom(ib)%pp=0.0
           dom(ib)%su=0.0

           isp=dom(ib)%isp; iep=dom(ib)%iep
           jsp=dom(ib)%jsp; jep=dom(ib)%jep
           ksp=dom(ib)%ksp; kep=dom(ib)%kep

           do  k=ksp,kep
              do  i=isp,iep
                 do j=jsp,jep

                    dom(ib)%aw(i,j,k)=1.0/(dom(ib)%dx*dom(ib)%dx)
                    dom(ib)%ae(i,j,k)=1.0/(dom(ib)%dx*dom(ib)%dx)
                    dom(ib)%an(i,j,k)=1.0/(dom(ib)%dy*dom(ib)%dy)
                    dom(ib)%as(i,j,k)=1.0/(dom(ib)%dy*dom(ib)%dy)
                    dom(ib)%at(i,j,k)=1.0/(dom(ib)%dz*dom(ib)%dz)
                    dom(ib)%ab(i,j,k)=1.0/(dom(ib)%dz*dom(ib)%dz)

                    if (dom(ib)%iprev.lt.0 .and. dom(ib)%bc_west.ne.5)
     &  dom(ib)%aw(isp,j,k)=0.0
                    if (dom(ib)%inext.lt.0 .and. dom(ib)%bc_east.ne.5)
     &  dom(ib)%ae(iep,j,k)=0.0
                    if (dom(ib)%jprev.lt.0 .and. dom(ib)%bc_south.ne.5)
     &  dom(ib)%as(i,jsp,k)=0.0
                    if (dom(ib)%jnext.lt.0 .and. dom(ib)%bc_north.ne.5)
     &  dom(ib)%an(i,jep,k)=0.0
                    if (dom(ib)%kprev.lt.0 .and. dom(ib)%bc_bottom.ne.5)
     &  dom(ib)%ab(i,j,ksp)=0.0
                    if (dom(ib)%knext.lt.0 .and. dom(ib)%bc_top.ne.5)
     &  dom(ib)%at(i,j,kep)=0.0

                    dom(ib)%ap(i,j,k) = -1.0*(dom(ib)%aw(i,j,k)+
     &  dom(ib)%ae(i,j,k)+dom(ib)%as(i,j,k)+dom(ib)%an(i,j,k)+
     &  dom(ib)%ab(i,j,k)+dom(ib)%at(i,j,k))
                 end do
              end do
           end do
        end do

        call calvel

        iter=0
 1000   continue
        iter  = iter + 1

        call calmas

        if (rmax.lt.eps.and.iter.gt.1)  goto 3000

        call sipsol(44)

        nsweeps = nsweeps + nsweep

        buffer_ppref=0.d0
 !       do ib=1,nbp
!           if(dom_id(ib).eq.prefdom) then
!              buffer_ppref=dom(ib)%p(ipref,jpref,kpref)+
!     &dom(ib)%pp(ipref,jpref,kpref)
!           end if
!        end do

        call MPI_ALLREDUCE(buffer_ppref,ppref,1,MPI_FLT,MPI_SUM,
     &MPI_COMM_WORLD,ierr)

        do ib=1,nbp
           do  k=dom(ib)%ksp-1,dom(ib)%kep+1
              do  i=dom(ib)%isp-1,dom(ib)%iep+1
                 do j=dom(ib)%jsp-1,dom(ib)%jep+1
                    dom(ib)%p(i,j,k)=(dom(ib)%p(i,j,k)+
     &(dom(ib)%pp(i,j,k)-ppref))
                    dom(ib)%pp(i,j,k)=0.0
                 end do
              end do
           end do
        end do

        call exchange(4)
        call pbound

        call calvel

        if (iter.lt.niter) goto 1000

 3000   continue

        if(rmax.gt.10.d0) then
           if(myrank.eq.0) write(numfile,*),'BIG RMAX!! STOP!!!!!',rmax
           write(6,*),'BIG RMAX!! STOP!!!!!!!!',rmax

           if(myrank.eq.0) then
              open (unit=101, file='final_ctime.dat')
              write (101,'(i8,3F15.6)') 
     &   ntime,ctime,forcn,qstpn
              close(101)
           end if
           call MPI_BARRIER (MPI_COMM_WORLD,ierr)
           call MPI_BARRIER (MPI_COMM_WORLD,ierr)
           stop
        end if

        call boundu
        call boundv
        call boundw

        if (myrank.eq.0) then
           wtimedum = MPI_WTIME ( ) - wtime
           write (6,'(1x,a,i8,a,i8,a,i4,a,i4,a,e13.6,a,e13.6)')
     &      ' myrank:',myrank,' ntime:',ntime,' iters:',iter,
     &      ' sweeps:',nsweeps,' rmax:',rmax
           write (numfile,'(1x,a,i8,a,i8,a,i4,a,i4,a,e13.6,a,e13.6)')
     &      ' myrank:',myrank,' ntime:',ntime,' iters:',iter,
     &      ' sweeps:',nsweeps,' rmax:',rmax
           write(numfile2,'(i8,f15.6,3e20.6)')
     & ntime,wtimedum,rmax,dt,Mdef
        end if

        return
        end subroutine pressure_1sweep
!##########################################################################
