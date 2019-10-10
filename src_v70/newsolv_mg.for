!##########################################################################
        subroutine newsolv_mg
!##########################################################################
        use vars
        use mpi
        use multidata
        implicit none
        integer i,j,k,ib
        double precision :: ppref,buffer_ppref
        real ( kind =8 )  :: wtimedum

        if (L_LSM) then! .or. L_LSMbase)  then
	    call coeff
	  end if

        call calvel
        call calmas

        do ib=1,nbp
           dom(ib)%pp=0.0
        end do

        do iter=1,maxcy

           call mgkcyc

           buffer_ppref=0.d0
!           do ib=1,nbp
!              if(dom_id(ib).eq.prefdom) then
!                 buffer_ppref=dom(ib)%p(ipref,jpref,kpref)+
!     &dom(ib)%pp(ipref,jpref,kpref)
!              end if
!           end do

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
           call calmas

           if (rmax.lt.eps.and.iter.gt.1)  goto 3000

        end do

        if (myrank.eq.0) print*,'not converged!! ',maxcy,rmax
        if (myrank.eq.0) write(numfile,*),'not converged!! ',maxcy,rmax
 3000   continue

        if(rmax.gt.10.0) then
           if(myrank.eq.0) then
		 write(numfile,*),'BIG RMAX!! STOP!!!!!',rmax
             open (unit=101, file='final_ctime.dat')
               write (101,'(i8,3F15.6)') ntime,ctime,forcn,qstpn
             close(101)
           end if
           call MPI_BARRIER (MPI_COMM_WORLD,ierr)
           call MPI_BARRIER (MPI_COMM_WORLD,ierr)
           stop
        end if

        call boundu
        call boundv
        call boundw


        if(myrank.eq.0) then
           wtimedum = MPI_WTIME ( ) - wtime
           write(6,'(1x,a,i8,a,e13.6,a,i8)')
     & ' ntime:',ntime,' rmax:',rmax,' 	  iter',iter
           write(numfile,'(1x,a,i8,a,e13.6,a,i8)')
     & ' ntime:',ntime,' rmax:',rmax,' 	  iter',iter
           write(numfile2,'(i8,f15.6,3e20.6)')
     & ntime,wtimedum,rmax,dt,Mdef
        end if

        return
        end subroutine newsolv_mg
!##########################################################################
        subroutine plotres
!##########################################################################
        use vars
        use mpi
        use multidata
        implicit none
        integer i,j,k,ib,sn
        double precision resid
        character*8 :: chb
        character*25 :: gf,gf2

        if ((mod(itime,n_out).eq.0).and.(itime.gt.itime_start)) then

        do ib=1,nbp
           write(chb,'(i8)') dom_id(ib)
           sn=len(trim(adjustl(chb)))
           chb=repeat('0',(3-sn))//trim(adjustl(chb))
           gf='resid'//trim(adjustl(chb))//'.plt'
           open (unit=88, file=gf)
           write (88,*) 'variables="x","y","z","res"'
           write (88,*)'zone ',
     & ' i=',dom(ib)%iep-dom(ib)%isp+1,', ',
     & ' j=',dom(ib)%jep-dom(ib)%jsp+1,', ',
     & ' k=',dom(ib)%kep-dom(ib)%ksp+1,' f=point'


           gf2='th_resid'//trim(adjustl(chb))//'.plt'
           open (unit=78, file=gf2)
           write (78,*) 'variables="x","y","z","res"'

           do  k=dom(ib)%ksp,dom(ib)%kep
              do j=dom(ib)%jsp,dom(ib)%jep
                 do  i=dom(ib)%isp,dom(ib)%iep
                    resid=( 
     &(dom(ib)%u(i,j,k)-dom(ib)%u(i-1,j,k))*dom(ib)%dy*dom(ib)%dz+
     &(dom(ib)%v(i,j,k)-dom(ib)%v(i,j-1,k))*dom(ib)%dx*dom(ib)%dz+
     &(dom(ib)%w(i,j,k)-dom(ib)%w(i,j,k-1))*dom(ib)%dx*dom(ib)%dy)

                    write (88,88) dom(ib)%xc(i),dom(ib)%yc(j),
     & dom(ib)%zc(k),abs(resid)
                    if(abs(resid).gt.1e-6) write (78,88) dom(ib)%xc(i),
     & dom(ib)%yc(j),dom(ib)%zc(k),abs(resid)
                 end do
              end do
           end do

        close (78)
        close (88)
        end do

        end if

88      format (10e25.8)

        end subroutine plotres
!##########################################################################
