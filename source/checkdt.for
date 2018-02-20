!##########################################################################
        subroutine checkdt
!##########################################################################
        use vars
        use mpi
        use multidata
        implicit none
        integer :: i,j,k,ib
	  double precision :: dxx,dyy,dzz,umax,vmax,wmax,dtmax
	  double precision :: dtmax1,dtvisc,dtvisc1,dtthr
	  double precision :: uc,vc,wc
        double precision :: buffer_umax,buffer_vmax,buffer_wmax
        double precision :: buffer_dtmax,dt1,small
	  double precision :: Cu,Cv,Cw

        umax=0.0
        vmax=0.0
        wmax=0.0
        small=1e-30

        MPI_FLT = MPI_DOUBLE_PRECISION

        do ib=1,nbp
           do i=dom(ib)%isu,dom(ib)%ieu
              do j=dom(ib)%jsu,dom(ib)%jeu
                 do k=dom(ib)%ksu,dom(ib)%keu
                    umax=max(umax,abs(dom(ib)%u(i,j,k)))
                 end do
              end do
           end do
        end do

        buffer_umax=umax
        call MPI_ALLREDUCE (buffer_umax,umax,1,MPI_FLT,MPI_MAX,
     &                   MPI_COMM_WORLD,ierr )

        do ib=1,nbp
           do i=dom(ib)%isv,dom(ib)%iev
              do j=dom(ib)%jsv,dom(ib)%jev
                 do k=dom(ib)%ksv,dom(ib)%kev
                    vmax=max(vmax,abs(dom(ib)%v(i,j,k)))
                 end do
              end do
           end do
        end do

        buffer_vmax=vmax
        call MPI_ALLREDUCE (buffer_vmax,vmax,1,MPI_FLT,MPI_MAX,
     &                   MPI_COMM_WORLD,ierr )

        do ib=1,nbp
           do i=dom(ib)%isw,dom(ib)%iew
              do j=dom(ib)%jsw,dom(ib)%jew
                 do k=dom(ib)%ksw,dom(ib)%kew
                    wmax=max(wmax,abs(dom(ib)%w(i,j,k)))
                 end do
              end do
           end do
        end do

        buffer_wmax=wmax
        call MPI_ALLREDUCE (buffer_wmax,wmax,1,MPI_FLT,MPI_MAX,
     &                   MPI_COMM_WORLD,ierr )


        dtvisc=1e10
        if(SGS) then
           do ib=1,nbp
              dxx=dom(ib)%dx*dom(ib)%dx
              dyy=dom(ib)%dy*dom(ib)%dy
              dzz=dom(ib)%dz*dom(ib)%dz
              do i=dom(ib)%isp,dom(ib)%iep
                 do j=dom(ib)%jsp,dom(ib)%jep
                    do k=dom(ib)%ksp,dom(ib)%kep
                       uc=0.5*(dom(ib)%u(i,j,k)+dom(ib)%u(i-1,j,k))
                       vc=0.5*(dom(ib)%v(i,j,k)+dom(ib)%v(i,j-1,k))
                       wc=0.5*(dom(ib)%w(i,j,k)+dom(ib)%w(i,j,k-1))
                       dtvisc1=1.0/( abs(uc/dom(ib)%dx)+
     & abs(vc/dom(ib)%dy)+abs(wc/dom(ib)%dz)+
     & 2.0*dom(ib)%vis(i,j,k)*(1.0/dxx+1.0/dyy+1.0/dzz)+small)
                       dtvisc=min(dtvisc,dtvisc1)
                       if(LENERGY) then
                          dtthr=0.5*Re*Pr/(1.0/dxx + 1.0/dyy + 1.0/dzz)
                          dtvisc=min(dtvisc,dtthr)
                       end if
                    end do
                 end do
              end do
           end do
        else
           do ib=1,nbp
              dxx=dom(ib)%dx*dom(ib)%dx
              dyy=dom(ib)%dy*dom(ib)%dy
              dzz=dom(ib)%dz*dom(ib)%dz
              dtvisc1=1.0/(1.0/(dxx) + 1.0/(dyy)+ 1.0/(dzz))*Re/2.0 
              dtvisc=min(dtvisc,dtvisc1)
              if(LENERGY) then
                 dtthr=0.5*Re*Pr/(1.0/dxx + 1.0/dyy + 1.0/dzz)
                 dtvisc=min(dtvisc,dtthr)
              end if
           end do
        end if


        dtmax=1e10
        do ib=1,nbp
           dtmax1=min(dom(ib)%dx/umax,dom(ib)%dy/vmax,dom(ib)%dz/wmax)
           dtmax = min(dtmax1,dtmax)
        end do

        dtmax = min(dtvisc,dtmax)

	  wtime_cfl=dt/dtmax			!Pablo 07_2017

        dtmax = safety_factor * dtmax 
        dt=min(dt*1.1,dtmax)

        if (L_LSM)  then
          do ib=1,nbp
            dtvisc=max(mul/densl,mug/densg)*
     & (2.0/(dxx)+2.0/(dyy)+2.0/(dzz))
            Cu=1./((umax/dom(ib)%dx+dtvisc)+sqrt((umax/dom(ib)%dx+
     & dtvisc)**2+4.*abs(grx)/dom(ib)%dx))
            Cv=1./((vmax/dom(ib)%dy+dtvisc)+sqrt((vmax/dom(ib)%dy+
     & dtvisc)**2+4.*abs(gry)/dom(ib)%dy))
            Cw=1./((wmax/dom(ib)%dz+dtvisc)+sqrt((wmax/dom(ib)%dz+
     & dtvisc)**2+4.*abs(grz)/dom(ib)%dz))
            dt = 2.*safety_factor*min(Cu,Cv,Cw)
	    end do
        end if

        buffer_dtmax=dt
        call MPI_ALLREDUCE (buffer_dtmax,dt1,1,MPI_FLT,MPI_MIN,
     &                   MPI_COMM_WORLD,ierr )  

        dt=dt1

	  if (variTS.eq. .true.) dt=dtinit

       dtsum=dtsum+dt

        dtavg=dtsum/(itime-itime_start+1)

        return
        end subroutine
!##########################################################################
