!######################################################################
      module lsm
!######################################################################
	SAVE
      DOUBLE PRECISION,allocatable,dimension(:,:):: f,g,hj_phi,hj_phi1
      DOUBLE PRECISION :: dt_reinit,initial_vol,total_vol
     
      end module lsm
!######################################################################
      subroutine read_lsm_input_data
!######################################################################
      use vars
      use lsm
      use multidata
	use mpi

        open (unit=12, file='in_lsm.cin')

        read (12,*)
        read (12,*) length
        read (12,*) L_LSM
	  read (12,*) L_LSMbase
        read (12,*) ntime_reinit,reldif_LSM
        read (12,*) cfl_lsm
	  read (12,*) LENDS
	  read (12,*) L_anim_phi,L_anim_grd 
	  read (12,*) densl,densg
	  read (12,*) nul,nug
        read (12,*) grx,gry,grz
	  read (12,*) slope

	 close(12)

        mul = nul * densl
        mug = nug * densg

	  if(L_LSM) reinit=.true.

	  if (L_LSM .and. solver.eq.1) then
	   if (myrank.eq.0) then
          print*,'Error: SIP solver not presently compatible with LSM'
	   endif
	   stop
	  endif

	  if (L_LSM .and. differencing.ne.3) then
	   if (myrank.eq.0) then
          print*,'Error: WENO differencing must be used with LSM'
	   endif
	   stop
	  endif

	  if (L_LSMbase .and. L_LSM) then
	   if (myrank.eq.0) then
          print*,'Error: L_LSMbase and L_LSM cannot both be true!'
	   endif
	   stop
	  endif

	  if (L_anim_phi .and. .not.L_LSM) then
	   if (myrank.eq.0) then
          print*,'Error: L_anim_phi cannot be true if L_LSM is false!'
	   endif
	   stop
	  endif

	  if (L_LSM .and. bc_w.ne.1 .or. L_LSMbase .and. bc_w.ne.1) then
	   if (myrank.eq.0) then
          print*,'Warning: LSM might not work well with this inlet BC!'
	   endif
!	   stop
	  endif
	return
	end subroutine
!######################################################################
      subroutine initial_lsm_3d_channel
!######################################################################
      use vars
      use lsm
      use multidata
	use mpi
      implicit none
      integer :: i,j,k,ib,tti,ttj,ttk,sn,is,ie,js,je,ks,ke
      DOUBLE PRECISION :: b,dummy
      character*8  :: chb
      character*31 :: gridfile
      character*80 :: dummyline

	if (myrank.eq.0) then
	  write(*,'(a)') '**********************************************'
	  write(*,'(a)') '*            TWO-PHASE SIMULATION'
	  if (L_LSMbase) then
	    write(*,'(a,f8.3)') '* Base case with rigid lid at z=',length
	  else
	    write(*,'(a,f8.3)') '* Initial water level: z=',length
	  endif
	  write(*,'(a)') '**********************************************'
	endif


 	if (myrank.eq.0) then
          open (unit=1703, file='LSM_phi_mass.dat' )	
	    write(1703,*)'Variables = time , Volume ,  Dvolume'
	endif

      do ib=1,nbp  
        tti=dom(ib)%ttc_i;  ttj=dom(ib)%ttc_j;  ttk=dom(ib)%ttc_k  
        is=dom(ib)%isp; ie=dom(ib)%iep
        js=dom(ib)%jsp; je=dom(ib)%jep
        ks=dom(ib)%ksp; ke=dom(ib)%kep

        allocate (dom(ib)%phi_init(tti,ttj,ttk),
     &dom(ib)%phi_reinit(tti,ttj,ttk),
     &dom(ib)%phi_new(tti,ttj,ttk),dom(ib)%dphi_dx(tti,ttj,ttk),
     &dom(ib)%dphi_dy(tti,ttj,ttk),dom(ib)%dphi_dz(tti,ttj,ttk),
     &dom(ib)%s_phi0(tti,ttj,ttk),dom(ib)%h_phi(tti,ttj,ttk))
! Note: phi and phim are allocated in initial.for
!
! Read in phi field if restarting from previous solution 
!
        if (L_LSM .and. lrestart .or. L_LSMbase .and. lrestart) then    
          write(chb,'(i4)') dom_id(ib)
          sn=len(trim(adjustl(chb)))
          chb=repeat('0',(4-sn))//trim(adjustl(chb))
          gridfile='tecout_phi_'//trim(adjustl(chb))//'.plt' 
          open (unit=703, file=gridfile)
          read (703,*) dummyline
          read (703,*) dummyline
          read (703,*) dummyline
          do k=ks-1,ke 
            do j=js-1,je 
              do i=is-1,ie 
                read (703,73) dummy,dummy,dummy,dom(ib)%phi(i,j,k)
     & ,dom(ib)%phim(i,j,k),dom(ib)%dens(i,j,k),dom(ib)%mu(i,j,k)
              end do
            end do
          end do
 73       format (10e20.8)
          close(703)
          dom(ib)%phi_init = 0.D0
          dom(ib)%phi_new = 0.D0
          dom(ib)%phi_reinit = 0.D0
        else
!
! Initialise uniform phi field, if not restarting
!
          dom(ib)%phi=1.d0
          dom(ib)%phi_init = 0.d0
          dom(ib)%phi_new = 0.d0
          dom(ib)%phi_reinit = 0.d0
      
          if (trim(keyword).eq.'channel') then      ! Channel flow case
            do k=1,ttk                
              do j=1,ttj    
                do i=1,tti
!
! Set initial free surface profile, if not uniform (this is case dependent)
! Initialise length (i.e. distance of free surface from bed)
!========== Richard cube test =========================================
!          if ((dom(ib)%xc(i).ge.0).and.(dom(ib)%xc(i).le.0.035)) then
!            length=0.025
!          else if ((dom(ib)%xc(i).gt.0.035).and.
!     &             (dom(ib)%xc(i).le.0.065)) then
!            length=-0.266666667*dom(ib)%xc(i)+0.0343
!          else if ((dom(ib)%xc(i).gt.0.065).and.
!     &             (dom(ib)%xc(i).le.0.25088)) then 
!            length=0.017
!          end if	
!========== Sibel constriction test ===================================
!          if ((dom(ib)%xc(i).ge.0).and.(dom(ib)%xc(i).le.0.59)) then
!            length=0.076
!          else if ((dom(ib)%xc(i).gt.0.59).and.
!     &             (dom(ib)%xc(i).le.0.885)) then
!            length=-0.213559322*dom(ib)%xc(i)+0.202
!          else if ((dom(ib)%xc(i).gt.0.885).and.
!     &             (dom(ib)%xc(i).le.1.475)) then 
!            length=0.013
!          end if
!======================================================================
!
! Initialise phi (free surface defined by phi=0, phi=-ve above, +ve below)
!
                  if (dom(ib)%zc(k).lt.length)   then
                   dom(ib)%phi(i,j,k) = 1.d0*dabs(dom(ib)%zc(k)-length)
                  else if (dom(ib)%zc(k).gt.length)   then
                   dom(ib)%phi(i,j,k) = -1.d0*dabs(dom(ib)%zc(k)-length)
                  else if (dom(ib)%zc(k).eq.length)  then
                   dom(ib)%phi(i,j,k) = 0.d0
                  end if
                  dom(ib)%s_phi0(i,j,k)=dom(ib)%phi(i,j,k)
                  dom(ib)%phi_reinit(i,j,k)=dom(ib)%phi(i,j,k)
                end do
              end do
            end do

          else if (trim(keyword).eq.'wave') then   ! Solitary wave case

            do k=1,ttk          
              do j=1,ttj    
                do i=1,tti  
                  b=length/(cosh(sqrt(3.*length)/2.*(dom(ib)%xc(i))))**2       
                  if (dom(ib)%zc(k).lt.(b+1))   then
                    dom(ib)%phi(i,j,k) = 1.0*abs(dom(ib)%zc(k)-(b+1))
                  else if (dom(ib)%zc(k).gt.(b+1))   then
                    dom(ib)%phi(i,j,k) = -1.0*abs(dom(ib)%zc(k)-(b+1))
                  else if (dom(ib)%zc(k).eq.(b+1))  then
                    dom(ib)%phi(i,j,k) = 0.0
                  end if
                  dom(ib)%s_phi0(i,j,k)=dom(ib)%phi(i,j,k)
                end do
              end do
            end do

          end if

        end if

        dt_reinit = cfl_lsm*max(dom(ib)%dx,dom(ib)%dy,dom(ib)%dz)

      end do    
!
! Set level set function phi to a signed distance function
!
      if (reinit)  call tvd_rk_reinit
!
! Define density and viscosity above, below and across the surface
!
      call heaviside

      return
      end subroutine initial_lsm_3d_channel
!######################################################################
      subroutine  lsm_3d
!######################################################################
!**********************************************************************
! 3 step runge kutta routine to return convected phi field
!**********************************************************************
      use vars
      use lsm
      use multidata
      use mpi
      implicit none
      integer :: i,j,k,ib
      DOUBLE PRECISION :: uijk,vijk,wijk,h1,h2,h3
!
! First RK step (phi --> phi_new)
!
	
      call dphi_a_v_3d(14)    !Calculate gradients of phi 

      do ib=1,nbp  

        do k=dom(ib)%ksp,dom(ib)%kep
         do i=dom(ib)%isp,dom(ib)%iep 
          do j=dom(ib)%jsp,dom(ib)%jep

           if (i.eq.dom(ib)%isp) then
           uijk=1.0/16.0*(5.0*dom(ib)%u(i-1,j,k)+15.0*dom(ib)%u(i,j,k)-
     &                   5.0*dom(ib)%u(i+1,j,k)+dom(ib)%u(i+2,j,k))
           else if (i.eq.dom(ib)%iep) then
           uijk=1.0/16.0*(dom(ib)%u(i-3,j,k)-5.0*dom(ib)%u(i-2,j,k)+
     &              15.0*dom(ib)%u(i-1,j,k)+5.0*dom(ib)%u(i,j,k))
           else
           uijk=1.0/16.0*(-dom(ib)%u(i-2,j,k)+9.0*dom(ib)%u(i-1,j,k)+
     &                9.0*dom(ib)%u(i,j,k)-dom(ib)%u(i+1,j,k))
           end if

           if (j.eq.dom(ib)%jsp) then
           vijk=1.0/16.0*(5.0*dom(ib)%v(i,j-1,k)+15.0*dom(ib)%v(i,j,k)-
     &                   5.0*dom(ib)%v(i,j+1,k)+dom(ib)%v(i,j+2,k))
           else if (j.eq.dom(ib)%jep) then
           vijk=1.0/16.0*(dom(ib)%v(i,j-3,k)-5.0*dom(ib)%v(i,j-2,k)+
     &              15.0*dom(ib)%v(i,j-1,k)+5.0*dom(ib)%v(i,j,k))
           else
           vijk=1.0/16.0*(-dom(ib)%v(i,j-2,k)+9.0*dom(ib)%v(i,j-1,k)+
     &                9.0*dom(ib)%v(i,j,k)-dom(ib)%v(i,j+1,k))
           end if

           if (k.eq.dom(ib)%ksp) then
           wijk=1.0/16.0*(5.0*dom(ib)%w(i,j,k-1)+15.0*dom(ib)%w(i,j,k)-
     &                   5.0*dom(ib)%w(i,j,k+1)+dom(ib)%w(i,j,k+2))
           else if (k.eq.dom(ib)%kep) then
           wijk=1.0/16.0*(dom(ib)%w(i,j,k-3)-5.0*dom(ib)%w(i,j,k-2)+
     &              15.0*dom(ib)%w(i,j,k-1)+5.0*dom(ib)%w(i,j,k))
           else
           wijk=1.0/16.0*(-dom(ib)%w(i,j,k-2)+9.0*dom(ib)%w(i,j,k-1)+
     &                9.0*dom(ib)%w(i,j,k)-dom(ib)%w(i,j,k+1))
           end if
!
! Convection
! 
           h1 = uijk*dom(ib)%dphi_dx(i,j,k)    
           h2 = vijk*dom(ib)%dphi_dy(i,j,k)   
           h3 = wijk*dom(ib)%dphi_dz(i,j,k)
   
           dom(ib)%phi_new(i,j,k) = dom(ib)%phi(i,j,k)-dt*(h1+h2+h3)

          end do
         end do
        end do

      end do

      call bound_lsm(16)    !(phi_new)
!
! 2nd RK step (phi_new --> pih_init)
!
      call dphi_a_v_3d(16)

      do ib=1,nbp  

        do k=dom(ib)%ksp,dom(ib)%kep
         do i=dom(ib)%isp,dom(ib)%iep 
          do j=dom(ib)%jsp,dom(ib)%jep

           if (i.eq.dom(ib)%isp) then
           uijk=1.0/16.0*(5.0*dom(ib)%u(i-1,j,k)+15.0*dom(ib)%u(i,j,k)-
     &                   5.0*dom(ib)%u(i+1,j,k)+dom(ib)%u(i+2,j,k))
           else if (i.eq.dom(ib)%iep) then
           uijk=1.0/16.0*(dom(ib)%u(i-3,j,k)-5.0*dom(ib)%u(i-2,j,k)+
     &              15.0*dom(ib)%u(i-1,j,k)+5.0*dom(ib)%u(i,j,k))
           else
           uijk=1.0/16.0*(-dom(ib)%u(i-2,j,k)+9.0*dom(ib)%u(i-1,j,k)+
     &                9.0*dom(ib)%u(i,j,k)-dom(ib)%u(i+1,j,k))
           end if

           if (j.eq.dom(ib)%jsp) then
           vijk=1.0/16.0*(5.0*dom(ib)%v(i,j-1,k)+15.0*dom(ib)%v(i,j,k)-
     &                   5.0*dom(ib)%v(i,j+1,k)+dom(ib)%v(i,j+2,k))
           else if (j.eq.dom(ib)%jep) then
           vijk=1.0/16.0*(dom(ib)%v(i,j-3,k)-5.0*dom(ib)%v(i,j-2,k)+
     &              15.0*dom(ib)%v(i,j-1,k)+5.0*dom(ib)%v(i,j,k))
           else
           vijk=1.0/16.0*(-dom(ib)%v(i,j-2,k)+9.0*dom(ib)%v(i,j-1,k)+
     &                9.0*dom(ib)%v(i,j,k)-dom(ib)%v(i,j+1,k))
           end if

           if (k.eq.dom(ib)%ksp) then
           wijk=1.0/16.0*(5.0*dom(ib)%w(i,j,k-1)+15.0*dom(ib)%w(i,j,k)-
     &                   5.0*dom(ib)%w(i,j,k+1)+dom(ib)%w(i,j,k+2))
           else if (k.eq.dom(ib)%kep) then
           wijk=1.0/16.0*(dom(ib)%w(i,j,k-3)-5.0*dom(ib)%w(i,j,k-2)+
     &              15.0*dom(ib)%w(i,j,k-1)+5.0*dom(ib)%w(i,j,k))
           else
           wijk=1.0/16.0*(-dom(ib)%w(i,j,k-2)+9.0*dom(ib)%w(i,j,k-1)+
     &                9.0*dom(ib)%w(i,j,k)-dom(ib)%w(i,j,k+1))
           end if
!
! Convection
! 
           h1 = uijk*dom(ib)%dphi_dx(i,j,k)    
           h2 = vijk*dom(ib)%dphi_dy(i,j,k)   
           h3 = wijk*dom(ib)%dphi_dz(i,j,k)

           dom(ib)%phi_init(i,j,k) = 0.75*dom(ib)%phi(i,j,k) +
     & 0.25*dom(ib)%phi_new(i,j,k) - 0.25*dt*(h1+h2+h3)

          end do
         end do
        end do

      end do

      call bound_lsm(17)   !(phi_init)
!
! 3rd RK step (phi_init --> phi)
!
      call dphi_a_v_3d(17)     

      do ib=1,nbp  

        do k=dom(ib)%ksp,dom(ib)%kep
         do i=dom(ib)%isp,dom(ib)%iep 
          do j=dom(ib)%jsp,dom(ib)%jep

           if (i.eq.dom(ib)%isp) then
           uijk=1.0/16.0*(5.0*dom(ib)%u(i-1,j,k)+15.0*dom(ib)%u(i,j,k)-
     &                   5.0*dom(ib)%u(i+1,j,k)+dom(ib)%u(i+2,j,k))
           else if (i.eq.dom(ib)%iep) then
           uijk=1.0/16.0*(dom(ib)%u(i-3,j,k)-5.0*dom(ib)%u(i-2,j,k)+
     &              15.0*dom(ib)%u(i-1,j,k)+5.0*dom(ib)%u(i,j,k))
           else
           uijk=1.0/16.0*(-dom(ib)%u(i-2,j,k)+9.0*dom(ib)%u(i-1,j,k)+
     &                9.0*dom(ib)%u(i,j,k)-dom(ib)%u(i+1,j,k))
           end if

           if (j.eq.dom(ib)%jsp) then
           vijk=1.0/16.0*(5.0*dom(ib)%v(i,j-1,k)+15.0*dom(ib)%v(i,j,k)-
     &                   5.0*dom(ib)%v(i,j+1,k)+dom(ib)%v(i,j+2,k))
           else if (j.eq.dom(ib)%jep) then
           vijk=1.0/16.0*(dom(ib)%v(i,j-3,k)-5.0*dom(ib)%v(i,j-2,k)+
     &              15.0*dom(ib)%v(i,j-1,k)+5.0*dom(ib)%v(i,j,k))
           else
           vijk=1.0/16.0*(-dom(ib)%v(i,j-2,k)+9.0*dom(ib)%v(i,j-1,k)+
     &                9.0*dom(ib)%v(i,j,k)-dom(ib)%v(i,j+1,k))
           end if

           if (k.eq.dom(ib)%ksp) then
           wijk=1.0/16.0*(5.0*dom(ib)%w(i,j,k-1)+15.0*dom(ib)%w(i,j,k)-
     &                   5.0*dom(ib)%w(i,j,k+1)+dom(ib)%w(i,j,k+2))
           else if (k.eq.dom(ib)%kep) then
           wijk=1.0/16.0*(dom(ib)%w(i,j,k-3)-5.0*dom(ib)%w(i,j,k-2)+
     &              15.0*dom(ib)%w(i,j,k-1)+5.0*dom(ib)%w(i,j,k))
           else
           wijk=1.0/16.0*(-dom(ib)%w(i,j,k-2)+9.0*dom(ib)%w(i,j,k-1)+
     &                9.0*dom(ib)%w(i,j,k)-dom(ib)%w(i,j,k+1))
           end if
!
! Convection
! 
           h1 = uijk*dom(ib)%dphi_dx(i,j,k)    
           h2 = vijk*dom(ib)%dphi_dy(i,j,k)   
           h3 = wijk*dom(ib)%dphi_dz(i,j,k)
     
           dom(ib)%phi(i,j,k) = 1.0/3.0*dom(ib)%phi(i,j,k) +
     & 2.0/3.0*dom(ib)%phi_init(i,j,k) - 2.0/3.0*dt*(h1+h2+h3)   ! We now have updated phi field

!===========================================================================================
! Relaxation zone near East Boundary ===== if x in last L/2, ADJUST ====
!===========================================================================================
! for z=length phi_target = 0
! for z<length phi_target = 1.d0*dabs(dom(ib)%zc(k)-length)
! for z>length phi_target = -1.d0*dabs(dom(ib)%zc(k)-length)

!	X_R_e = (dom(ib)%xc(i)-10.0d0)/2.80d0                                  !== Aristos Christou == X_R = apply relaxation at the very end (L/2)  
!           if (X_R_e.ge.0) then
!	    gamma = 1.0-((exp(X_R_e**3.5)-1.0)/(exp(1.0)-1.0))
!                if (dom(ib)%zc(k).eq.length) then
!                  dom(ib)%phi(i,j,k) = dom(ib)%phi(i,j,k)
!     &  	  *gamma
!                else if (dom(ib)%zc(k).lt.length)   then
!                   dom(ib)%phi(i,j,k) = gamma*dom(ib)%phi(i,j,k)+
!     &            (1.0-gamma)*(1.d0*abs(dom(ib)%zc(k)-length))
!                else if (dom(ib)%zc(k).gt.length)   then
!                   dom(ib)%phi(i,j,k) = gamma*dom(ib)%phi(i,j,k)+
!     &            (1.0-gamma)*(-1.d0*abs(dom(ib)%zc(k)-length))
!                end if
!            end if

!
! Hold level to be held constant at inflow and outflow (if required - may help with stability in inflow-outflow sims)
!
           if (trim(keyword).eq.'channel' .and. lends) then
             if (dom(ib)%iprev.lt.0) then
               if ((i.ge.dom(ib)%isu).and.(i.le.dom(ib)%isu+5)) then   ! Distance phi fixed                 
                 if (dom(ib)%zc(k).lt.length)   then
                   dom(ib)%phi(i,j,k) = 1.0*abs(dom(ib)%zc(k)-length)
                 else if (dom(ib)%zc(k).gt.length)   then
                   dom(ib)%phi(i,j,k) = -1.0*abs(dom(ib)%zc(k)-length)
                 else if (dom(ib)%zc(k).eq.length)  then
                   dom(ib)%phi(i,j,k) = 0.0
                 end if
               end if
             end if
             if  (dom(ib)%inext.lt.0) then
               if ((i.le.dom(ib)%ieu).and.(i.ge.dom(ib)%ieu-5)) then     
                 if (dom(ib)%zc(k).lt.length)   then
                   dom(ib)%phi(i,j,k) = 1.0*abs(dom(ib)%zc(k)-length)
                 else if (dom(ib)%zc(k).gt.length)   then
                   dom(ib)%phi(i,j,k) = -1.0*abs(dom(ib)%zc(k)-length)
                 else if (dom(ib)%zc(k).eq.length)  then
                   dom(ib)%phi(i,j,k) = 0.0
                 end if   
               end if
             end if
           end if

         end do
        end do
       end do

      end do

      call bound_lsm(14)  !(phi)
!
! Reinitialise phi field to a signed distance function
!
      if (reinit .and. mod(itime,1).eq.0) call tvd_rk_reinit     !(phi)
	
      return
      end 
!#######################################################################
      subroutine dphi_a_v_3d(op)     !(fi)            
!#######################################################################
      use vars
      use lsm
      use mpi
      use multidata
      implicit none
      integer :: i,j,k,ifi,op,ib
      DOUBLE PRECISION :: uijk,vijk,wijk,phi_minus_ip12j
      DOUBLE PRECISION :: phi_plus_ip12j,phi_minus_ijp12
      DOUBLE PRECISION :: phi_plus_ijp12

      call hj_weno_dxplus_3d(op)          
      call hj_weno_dyplus_3d(op)      
      call hj_weno_dzplus_3d(op)              
      call hj_weno_dxminus_3d(op)       
      call hj_weno_dyminus_3d(op)    
      call hj_weno_dzminus_3d(op)     
  
      do ib=1,nbp  

        do k=dom(ib)%ksp,dom(ib)%kep
         do i=dom(ib)%isp,dom(ib)%iep 
          do j=dom(ib)%jsp,dom(ib)%jep
         
          if (i.eq.dom(ib)%isp) then
          uijk=1.0/16.0*(5.0*dom(ib)%u(i-1,j,k)+15.0*dom(ib)%u(i,j,k)-
     &                   5.0*dom(ib)%u(i+1,j,k)+dom(ib)%u(i+2,j,k))
          else if (i.eq.dom(ib)%iep) then
          uijk=1.0/16.0*(dom(ib)%u(i-3,j,k)-5.0*dom(ib)%u(i-2,j,k)+
     &              15.0*dom(ib)%u(i-1,j,k)+5.0*dom(ib)%u(i,j,k))
          else
          uijk=1.0/16.0*(-dom(ib)%u(i-2,j,k)+9.0*dom(ib)%u(i-1,j,k)+
     &                9.0*dom(ib)%u(i,j,k)-dom(ib)%u(i+1,j,k))
          end if

          if (j.eq.dom(ib)%jsp) then
          vijk=1.0/16.0*(5.0*dom(ib)%v(i,j-1,k)+15.0*dom(ib)%v(i,j,k)-
     &                   5.0*dom(ib)%v(i,j+1,k)+dom(ib)%v(i,j+2,k))
          else if (j.eq.dom(ib)%jep) then
          vijk=1.0/16.0*(dom(ib)%v(i,j-3,k)-5.0*dom(ib)%v(i,j-2,k)+
     &              15.0*dom(ib)%v(i,j-1,k)+5.0*dom(ib)%v(i,j,k))
          else
          vijk=1.0/16.0*(-dom(ib)%v(i,j-2,k)+9.0*dom(ib)%v(i,j-1,k)+
     &                9.0*dom(ib)%v(i,j,k)-dom(ib)%v(i,j+1,k))
          end if

          if (k.eq.dom(ib)%ksp) then
          wijk=1.0/16.0*(5.0*dom(ib)%w(i,j,k-1)+15.0*dom(ib)%w(i,j,k)-
     &                   5.0*dom(ib)%w(i,j,k+1)+dom(ib)%w(i,j,k+2))
          else if (k.eq.dom(ib)%kep) then
          wijk=1.0/16.0*(dom(ib)%w(i,j,k-3)-5.0*dom(ib)%w(i,j,k-2)+
     &              15.0*dom(ib)%w(i,j,k-1)+5.0*dom(ib)%w(i,j,k))
          else
          wijk=1.0/16.0*(-dom(ib)%w(i,j,k-2)+9.0*dom(ib)%w(i,j,k-1)+
     &                9.0*dom(ib)%w(i,j,k)-dom(ib)%w(i,j,k+1))
          end if

           if (uijk.gt.0.0)   
     &        dom(ib)%dphi_dx(i,j,k) = dom(ib)%dphi_dxminus(i,j,k)
           if (uijk.lt.0.0)      
     &        dom(ib)%dphi_dx(i,j,k) = dom(ib)%dphi_dxplus(i,j,k)
           if (uijk.eq.0.0)
     &        dom(ib)%dphi_dx(i,j,k) = 0.0 

           if (vijk.gt.0.0)   
     &        dom(ib)%dphi_dy(i,j,k) = dom(ib)%dphi_dyminus(i,j,k)
           if (vijk.lt.0.0)      
     &        dom(ib)%dphi_dy(i,j,k) = dom(ib)%dphi_dyplus(i,j,k)
           if (vijk.eq.0.0)
     &        dom(ib)%dphi_dy(i,j,k) = 0.0 

           if (wijk.gt.0.0)   
     &        dom(ib)%dphi_dz(i,j,k) = dom(ib)%dphi_dzminus(i,j,k)
           if (wijk.lt.0.0)      
     &        dom(ib)%dphi_dz(i,j,k) = dom(ib)%dphi_dzplus(i,j,k)
           if (wijk.eq.0.0)
     &        dom(ib)%dphi_dz(i,j,k) = 0.0 

         end do
       end do
      end do      

      end do

      return
      end subroutine dphi_a_v_3d
!######################################################################
      subroutine tvd_rk_reinit
!######################################################################
!**********************************************************************
! 3 step runge kutta routine which returns re-initialised phi field
!**********************************************************************
      use vars
      use lsm
      use multidata
      use mpi
      implicit none
      integer :: i,j,k,it,ib
      logical :: bool
      DOUBLE PRECISION :: max_abs,abs_dphi,abs_phidiff,max_phidiff
      DOUBLE PRECISION :: local_max_abs,local_max_phidiff

      bool=.false.
      it=0

      do while (.not.bool)
!
! 1st step
!     
      call dphi_for_reinit(14)      ! (phi)
      
      do ib=1,nbp  

        do k=dom(ib)%ksp,dom(ib)%kep
         do i=dom(ib)%isp,dom(ib)%iep 
          do j=dom(ib)%jsp,dom(ib)%jep
        
        abs_dphi = dsqrt(dom(ib)%dphi_dx(i,j,k)**2+
     & dom(ib)%dphi_dy(i,j,k)**2+dom(ib)%dphi_dz(i,j,k)**2)

        dom(ib)%s_phi0(i,j,k) = dom(ib)%phi(i,j,k)/
     & dsqrt(dom(ib)%phi(i,j,k)**2+(abs_dphi**2)*dom(ib)%dx**2)

        dom(ib)%phi_reinit(i,j,k) = dom(ib)%phi(i,j,k) +
     & dt_reinit*(dom(ib)%s_phi0(i,j,k)-dom(ib)%s_phi0(i,j,k)*abs_dphi) 

        end do
       end do
      end do
      
      end do
 
      call bound_lsm(15) ! (phi_reinit)    
!
! 2nd step 
!
      call dphi_for_reinit(15) ! (phi_reinit)    

      do ib=1,nbp  

        do k=dom(ib)%ksp,dom(ib)%kep
         do i=dom(ib)%isp,dom(ib)%iep 
          do j=dom(ib)%jsp,dom(ib)%jep

        abs_dphi = dsqrt(dom(ib)%dphi_dx(i,j,k)**2+
     & dom(ib)%dphi_dy(i,j,k)**2+dom(ib)%dphi_dz(i,j,k)**2)

        dom(ib)%s_phi0(i,j,k) = dom(ib)%phi_reinit(i,j,k)/
     & dsqrt(dom(ib)%phi_reinit(i,j,k)**2+(abs_dphi**2)*dom(ib)%dx**2)

        dom(ib)%phi_reinit(i,j,k) = 0.75d0*dom(ib)%phi(i,j,k)+
     &0.25d0*dom(ib)%phi_reinit(i,j,k)+
     &0.25d0*dt_reinit*(dom(ib)%s_phi0(i,j,k)-
     &dom(ib)%s_phi0(i,j,k)*abs_dphi)

          end do
         end do
        end do

      end do

      call bound_lsm(15) ! (phi_reinit)    
!
! 3rd step
!
      call dphi_for_reinit(15) ! (phi_reinit) 
     
      do ib=1,nbp  

        do k=dom(ib)%ksp,dom(ib)%kep
         do i=dom(ib)%isp,dom(ib)%iep 
          do j=dom(ib)%jsp,dom(ib)%jep

        abs_dphi = dsqrt(dom(ib)%dphi_dx(i,j,k)**2+
     & dom(ib)%dphi_dy(i,j,k)**2+dom(ib)%dphi_dz(i,j,k)**2)

        dom(ib)%s_phi0(i,j,k) = dom(ib)%phi_reinit(i,j,k)/
     &dsqrt(dom(ib)%phi_reinit(i,j,k)**2+(abs_dphi**2)*dom(ib)%dx**2)

        dom(ib)%phi_reinit(i,j,k) =1.d0/3.d0*dom(ib)%phi(i,j,k)+
     & 2.d0/3.d0*dom(ib)%phi_reinit(i,j,k)+
     & 2.d0/3.d0*dt_reinit*(dom(ib)%s_phi0(i,j,k)-
     & dom(ib)%s_phi0(i,j,k)*abs_dphi)

          end do
         end do
        end do

      end do

      call bound_lsm(15) ! (phi_reinit)

!
! Check convergence
!
      call dphi_for_reinit(15) ! (phi_reinit)     
          
      max_abs = 0.0
      max_phidiff = 0.0

      do ib=1,nbp  

       do k=dom(ib)%ksp,dom(ib)%kep
        do i=dom(ib)%isp,dom(ib)%iep 
         do j=dom(ib)%jsp,dom(ib)%jep

          abs_dphi = dsqrt(dom(ib)%dphi_dx(i,j,k)**2+
     & dom(ib)%dphi_dy(i,j,k)**2+dom(ib)%dphi_dz(i,j,k)**2)
          max_abs    = max(max_abs,abs_dphi)
          abs_phidiff=dabs(dom(ib)%phi(i,j,k)-dom(ib)%phi_reinit(i,j,k))
          max_phidiff= max(max_phidiff,abs_phidiff)

         end do
        end do
       end do

      end do

      local_max_abs = max_abs

      call mpi_allreduce(local_max_abs,max_abs,1,mpi_flt,mpi_max,
     &mpi_comm_world,ierr)

      local_max_phidiff = max_phidiff

      call mpi_allreduce(local_max_phidiff,max_phidiff,1,mpi_flt,
     &mpi_max,mpi_comm_world,ierr)

      if ((max_phidiff.lt.reldif_lsm).and.(it.ge.1)) then
        bool=.true.
      else
        if (it.ge.ntime_reinit) then
          bool=.true.
        else
          it=it+1
        end if
      end if
   
      do ib=1,nbp
        do k=dom(ib)%ksp,dom(ib)%kep
          do i=dom(ib)%isp,dom(ib)%iep 
            do j=dom(ib)%jsp,dom(ib)%jep  
              dom(ib)%phi(i,j,k) = dom(ib)%phi_reinit(i,j,k)
            end do
          end do
        end do
      end do

      call bound_lsm(14) ! (phi)  

      end do
        
      if (myrank.eq.0) then
        write(*,'(a,f12.6,e15.6,a,i2)')
     & 'normPHI(reinit)',max_abs,max_phidiff,'  its:',it
        write(numfile3,'(i8,f18.8,i8,2f18.8)') ntime,max_abs,it,ctime,dt
      endif

      return
      end subroutine tvd_rk_reinit
!######################################################################
      subroutine dphi_for_reinit(op)  
!######################################################################
      use vars
      use lsm
      use multidata
      implicit none
      integer :: i,j,k,ifi,op,ib
      DOUBLE PRECISION :: s_phi012
      DOUBLE PRECISION, pointer, dimension(:,:,:)::fi
      DOUBLE PRECISION :: lsv,lssig,xm,xp,ym,yp,zm,zp

          call hj_weno_dxplus_3d(op)          
          call hj_weno_dyplus_3d(op)      
          call hj_weno_dzplus_3d(op)              
          call hj_weno_dxminus_3d(op)      
          call hj_weno_dyminus_3d(op)    
          call hj_weno_dzminus_3d(op)      

      do ib=1,nbp

       select case (op)
        case (14) ! (phi)  
          fi => dom(ib)%phi
        case (15) ! (phi_reinit) 
          fi => dom(ib)%phi_reinit
       end select

        do k=dom(ib)%ksp,dom(ib)%kep
         do i=dom(ib)%isp,dom(ib)%iep 
          do j=dom(ib)%jsp,dom(ib)%jep

              lsv = fi(i,j,k)
              lssig = dom(ib)%s_phi0(i,j,k) ! s_phi0=phi initially

          xm=dom(ib)%dphi_dxminus(i,j,k)
          xp=dom(ib)%dphi_dxplus(i,j,k)
          ym=dom(ib)%dphi_dyminus(i,j,k)
          yp=dom(ib)%dphi_dyplus(i,j,k)
          zm=dom(ib)%dphi_dzminus(i,j,k)
          zp=dom(ib)%dphi_dzplus(i,j,k)

      if ((xm*lssig.gt.0.0).and.(xp*lssig.gt.-xm*lssig)) then
            dom(ib)%dphi_dx(i,j,k)=dom(ib)%dphi_dxminus(i,j,k)
      end if

      if((xp*lssig.lt.0.0).and.(xm*lssig.lt.-xp*lssig)) then
            dom(ib)%dphi_dx(i,j,k)=dom(ib)%dphi_dxplus(i,j,k)
      end if 
   
      if ((xp*lssig.gt.0.0).and.(xm*lssig.lt.0.0)) then
            dom(ib)%dphi_dx(i,j,k)=0.0
      end if

      if ((ym*lssig.gt.0.0).and.(yp*lssig.gt.-ym*lssig)) then
            dom(ib)%dphi_dy(i,j,k)=dom(ib)%dphi_dyminus(i,j,k)
      end if

      if ((yp*lssig.lt.0.0).and.(ym*lssig.lt.-yp*lssig)) then
           dom(ib)%dphi_dy(i,j,k)=dom(ib)%dphi_dyplus(i,j,k)
      end if

      if ((yp*lssig.gt.0.0).and.(ym*lssig.lt.0.0)) then
          dom(ib)%dphi_dy(i,j,k) = 0.0
      end if

      if ((zm*lssig.gt.0.0).and.(zp*lssig.gt.-zm*lssig)) then
            dom(ib)%dphi_dz(i,j,k)=dom(ib)%dphi_dzminus(i,j,k)
      end if

      if((zp*lssig.lt.0.0).and.(zm*lssig.lt.-zp*lssig)) then
            dom(ib)%dphi_dz(i,j,k)=dom(ib)%dphi_dzplus(i,j,k)
      end if 
   
      if ((zp*lssig.gt.0.0).and.(zm*lssig.lt.0.0)) then
            dom(ib)%dphi_dz(i,j,k)=0.0
      end if

          end do
         end do
        end do

       end do

      return
      end
!######################################################################
      subroutine hj_weno_dxplus_3d(op)           
!######################################################################
      use vars
      use multidata
      implicit none
      integer :: i,j,k,l,k_star,np,nq,nr,op,ib
      DOUBLE PRECISION, pointer, dimension(:,:,:)::fi
      DOUBLE PRECISION :: q1,q2,q3,c_star,v1,v2,v3,v4,v5,v11,v22
	DOUBLE PRECISION :: v33,v44,v55,s1,s2,s3,a1,a2,a3,w1,w2,w3,e

      do ib=1,nbp

       select case (op)
        case (1)  
          fi => dom(ib)%u
          np = dom(ib)%niul;  nq = dom(ib)%njul;  nr = dom(ib)%nkul
        case (2) 
          fi => dom(ib)%v
          np = dom(ib)%nivl;  nq = dom(ib)%njvl;  nr = dom(ib)%nkvl
        case (3)  
          fi => dom(ib)%w
          np = dom(ib)%niwl;  nq = dom(ib)%njwl;  nr = dom(ib)%nkwl
        case (14) ! (phi)
          fi => dom(ib)%phi
          np = dom(ib)%ttc_i; nq = dom(ib)%ttc_j; nr = dom(ib)%ttc_k
        case (15) ! (phi_reinit) 
          fi => dom(ib)%phi_reinit
          np = dom(ib)%ttc_i; nq = dom(ib)%ttc_j; nr = dom(ib)%ttc_k
        case (16) ! (phi_new)
          fi => dom(ib)%phi_new
          np = dom(ib)%ttc_i; nq = dom(ib)%ttc_j; nr = dom(ib)%ttc_k
        case (17) ! (phi_init)
          fi => dom(ib)%phi_init
          np = dom(ib)%ttc_i; nq = dom(ib)%ttc_j; nr = dom(ib)%ttc_k
       end select

      if (dom(ib)%iprev.lt.0 .and. dom(ib)%bc_west.ne.5) then
       do k=1,nr
        do j=1,nq
            fi(2,j,k)    = fi(3,j,k)
            fi(1,j,k)    = fi(2,j,k)
        end do
       end do 
      else if (dom(ib)%inext.lt.0 .and. dom(ib)%bc_east.ne.5) then
       do k=1,nr
        do j=1,nq
            fi(np-1,j,k) = fi(np-2,j,k)
            fi(np,j,k)   = fi(np-1,j,k)
        end do
       end do 
      end if
!
! Compute 1st order spatial derivative
!
      do k=1,nr
       do i=1,np-1
        do j=1,nq
           dom(ib)%d1(i,j,k) = (fi(i+1,j,k) - fi(i,j,k))/dom(ib)%dx
        end do
       end do
      end do 
!
! Compute 5th order derivative (+ve direction)
!
      dom(ib)%dphi_dxplus = 0.0

      do k=1,nr
       do i=1,np-6
        do j=1,nq

            l=i

           v1 = dom(ib)%d1(l+5,j,k)
           v2 = dom(ib)%d1(l+4,j,k)
           v3 = dom(ib)%d1(l+3,j,k)
           v4 = dom(ib)%d1(l+2,j,k)
           v5 = dom(ib)%d1(l+1,j,k)

           v11 = v1**2
           v22 = v2**2
           v33 = v3**2
           v44 = v4**2
           v55 = v5**2

           e = (1.0e-6) * max(v11,v22,v33,v44,v55) + 1.0d-99

      s1 = (13.0/12.0)*(v1-2.0*v2+v3)**2 + 0.25*(v1-4.0*v2+3.0*v3)**2
      s2 = (13.0/12.0)*(v2-2.0*v3+v4)**2 + 0.25*(v2-v4)**2 
      s3 = (13.0/12.0)*(v3-2.0*v4+v5)**2 + 0.25*(3.0*v3-4.0*v4+v5)**2

           a1 = 0.1/(e+s1)**2
           a2 = 0.6/(e+s2)**2
           a3 = 0.3/(e+s3)**2

           w1 = a1/(a1+a2+a3)
           w2 = a2/(a1+a2+a3)
           w3 = a3/(a1+a2+a3)

      dom(ib)%dphi_dxplus(i+3,j,k) = w1*(v1*(1.0/3.0)-v2*(7.0/6.0)+
     &v3*(11.0/6.0)) + w2*(-v2*(1.0/6.0)+v3*(5.0/6.0)+v4*(1.0/3.0))+
     &w3*(v3*(1.0/3.0)+v4*(5.0/6.0)-v5*(1.0/6.0))

        end do
       end do
      end do

      end do

      return
      end subroutine hj_weno_dxplus_3d
!######################################################################
      subroutine hj_weno_dxminus_3d(op)        
!######################################################################
      use vars
      use multidata
      implicit none
      integer :: i,j,k,l,k_star,np,nq,nr,op,ib
      double precision, pointer, dimension(:,:,:)::fi
      DOUBLE PRECISION :: q1,q2,q3,c_star,v1,v2,v3,v4,v5,v11,v22
	DOUBLE PRECISION :: v33,v44,v55,s1,s2,s3,a1,a2,a3,w1,w2,w3,e

      do ib=1,nbp

       select case (op)
        case (1)  
          fi => dom(ib)%u
          np = dom(ib)%niul;  nq = dom(ib)%njul;  nr = dom(ib)%nkul
        case (2) 
          fi => dom(ib)%v
          np = dom(ib)%nivl;  nq = dom(ib)%njvl;  nr = dom(ib)%nkvl
        case (3)  
          fi => dom(ib)%w
          np = dom(ib)%niwl;  nq = dom(ib)%njwl;  nr = dom(ib)%nkwl
        case (14)  ! (phi)
          fi => dom(ib)%phi
          np = dom(ib)%ttc_i; nq = dom(ib)%ttc_j; nr = dom(ib)%ttc_k
        case (15)  ! (phi_reinit)
          fi => dom(ib)%phi_reinit
          np = dom(ib)%ttc_i; nq = dom(ib)%ttc_j; nr = dom(ib)%ttc_k
        case (16)  ! (phi_new)
          fi => dom(ib)%phi_new
          np = dom(ib)%ttc_i; nq = dom(ib)%ttc_j; nr = dom(ib)%ttc_k
        case (17)  ! (phi_init)
          fi => dom(ib)%phi_init
          np = dom(ib)%ttc_i; nq = dom(ib)%ttc_j; nr = dom(ib)%ttc_k
       end select     

      if (dom(ib)%iprev.lt.0 .and. dom(ib)%bc_west.ne.5) then
       do k=1,nr
        do j=1,nq
            fi(2,j,k)    = fi(3,j,k)
            fi(1,j,k)    = fi(2,j,k)
        end do
       end do 
      else if (dom(ib)%inext.lt.0 .and. dom(ib)%bc_east.ne.5) then
       do k=1,nr
        do j=1,nq
            fi(np-1,j,k) = fi(np-2,j,k)
            fi(np,j,k)   = fi(np-1,j,k)
        end do
       end do 
      end if
!
! Compute 1st order spatial derivative
!
      do k=1,nr
       do i=1,np-1
        do j=1,nq
           dom(ib)%d1(i,j,k) = (fi(i+1,j,k) - fi(i,j,k))/dom(ib)%dx
        end do
       end do
      end do
!
! Compute 5th order derivative (-ve direction)
!
      dom(ib)%dphi_dxminus = 0.0

      do k=1,nr
       do i=1,np-6
        do j=1,nq

           l=i-1

           v1 = dom(ib)%d1(l+1,j,k)
           v2 = dom(ib)%d1(l+2,j,k)
           v3 = dom(ib)%d1(l+3,j,k)
           v4 = dom(ib)%d1(l+4,j,k)
           v5 = dom(ib)%d1(l+5,j,k)

           v11 = v1**2
           v22 = v2**2
           v33 = v3**2
           v44 = v4**2
           v55 = v5**2

           e = (1.0e-6) * max(v11,v22,v33,v44,v55) + 1.0d-99

      s1 = (13.0/12.0)*(v1-2.0*v2+v3)**2 + 0.25*(v1-4.0*v2+3.0*v3)**2
      s2 = (13.0/12.0)*(v2-2.0*v3+v4)**2 + 0.25*(v2-v4)**2
      s3 = (13.0/12.0)*(v3-2.0*v4+v5)**2 + 0.25*(3.0*v3-4.0*v4+v5)**2

           a1 = 0.1/(e+s1)**2
           a2 = 0.6/(e+s2)**2
           a3 = 0.3/(e+s3)**2

           w1 = a1/(a1+a2+a3)
           w2 = a2/(a1+a2+a3)
           w3 = a3/(a1+a2+a3)

      dom(ib)%dphi_dxminus(i+3,j,k) = w1*(v1*(1.0/3.0)-v2*(7.0/6.0)+
     &v3*(11.0/6.0)) + w2*(-v2*(1.0/6.0)+v3*(5.0/6.0)+v4*(1.0/3.0))+
     &w3*(v3*(1.0/3.0)+v4*(5.0/6.0)-v5*(1.0/6.0))

        end do
       end do
      end do

      end do

      return
      end subroutine hj_weno_dxminus_3d
!######################################################################
      subroutine hj_weno_dyplus_3d(op)        
!######################################################################
      use vars
      use multidata
      implicit none
      integer :: i,j,k,l,k_star,np,nq,nr,op,ib
      double precision, pointer, dimension(:,:,:)::fi
      DOUBLE PRECISION :: q1,q2,q3,c_star,v1,v2,v3,v4,v5,v11,v22
	DOUBLE PRECISION :: v33,v44,v55,s1,s2,s3,a1,a2,a3,w1,w2,w3,e

      do ib=1,nbp

       select case (op)
        case (1)  
          fi => dom(ib)%u
          np = dom(ib)%niul;  nq = dom(ib)%njul;  nr = dom(ib)%nkul
        case (2) 
          fi => dom(ib)%v
          np = dom(ib)%nivl;  nq = dom(ib)%njvl;  nr = dom(ib)%nkvl
        case (3)  
          fi => dom(ib)%w
          np = dom(ib)%niwl;  nq = dom(ib)%njwl;  nr = dom(ib)%nkwl
        case (14)  ! (phi)
          fi => dom(ib)%phi
          np = dom(ib)%ttc_i; nq = dom(ib)%ttc_j; nr = dom(ib)%ttc_k
        case (15)  ! (phi_reinit)
          fi => dom(ib)%phi_reinit
          np = dom(ib)%ttc_i; nq = dom(ib)%ttc_j; nr = dom(ib)%ttc_k
        case (16)  ! (phi_new)
          fi => dom(ib)%phi_new
          np = dom(ib)%ttc_i; nq = dom(ib)%ttc_j; nr = dom(ib)%ttc_k
        case (17)  ! (phi_init)
          fi => dom(ib)%phi_init
          np = dom(ib)%ttc_i; nq = dom(ib)%ttc_j; nr = dom(ib)%ttc_k
       end select  

      if (dom(ib)%jprev.lt.0 .and. dom(ib)%bc_south.ne.5) then
       do k=1,nr
        do i=1,np
            fi(i,2,k)    = fi(i,3,k)
            fi(i,1,k)    = fi(i,2,k)
        end do
       end do   
      else if (dom(ib)%jnext.lt.0 .and. dom(ib)%bc_north.ne.5) then 
       do k=1,nr
        do i=1,np
            fi(i,nq-1,k) = fi(i,nq-2,k)
            fi(i,nq,k)   = fi(i,nq-1,k)
        end do
       end do   
      end if
!
! Compute 1st order spatial derivative
!
      do k=1,nr
       do i=1,np
        do j=1,nq-1
           dom(ib)%d1(i,j,k) = (fi(i,j+1,k) - fi(i,j,k))/dom(ib)%dy
        end do
       end do
      end do
!
! Compute 5th order derivative (+ve direction)
!
      dom(ib)%dphi_dyplus = 0.0

      do k=1,nr
       do i=1,np
        do j=1,nq-6
           l=j

           v1 = dom(ib)%d1(i,l+5,k)
           v2 = dom(ib)%d1(i,l+4,k)
           v3 = dom(ib)%d1(i,l+3,k)
           v4 = dom(ib)%d1(i,l+2,k)
           v5 = dom(ib)%d1(i,l+1,k)

           v11 = v1**2
           v22 = v2**2
           v33 = v3**2
           v44 = v4**2
           v55 = v5**2

           e = (1.0e-6) * max(v11,v22,v33,v44,v55) + 1.0d-99

      s1 = (13.0/12.0)*(v1-2.0*v2+v3)**2 + 0.25*(v1-4.0*v2+3.0*v3)**2
      s2 = (13.0/12.0)*(v2-2.0*v3+v4)**2 + 0.25*(v2-v4)**2
      s3 = (13.0/12.0)*(v3-2.0*v4+v5)**2 + 0.25*(3.0*v3-4.0*v4+v5)**2

           a1 = 0.1/(e+s1)**2
           a2 = 0.6/(e+s2)**2
           a3 = 0.3/(e+s3)**2

           w1 = a1/(a1+a2+a3)
           w2 = a2/(a1+a2+a3)
           w3 = a3/(a1+a2+a3)

      dom(ib)%dphi_dyplus(i,j+3,k) = w1*(v1*(1.0/3.0)-v2*(7.0/6.0)+
     &v3*(11.0/6.0)) + w2*(-v2*(1.0/6.0)+v3*(5.0/6.0)+v4*(1.0/3.0))+
     &w3*(v3*(1.0/3.0)+v4*(5.0/6.0)-v5*(1.0/6.0))

        end do
       end do
      end do

      end do

      return
      end subroutine hj_weno_dyplus_3d
!######################################################################
      subroutine hj_weno_dyminus_3d(op)       
!######################################################################
      use vars
      use multidata
      implicit none
      integer :: i,j,k,l,k_star,np,nq,nr,op,ib
      double precision, pointer, dimension(:,:,:)::fi
      DOUBLE PRECISION :: q1,q2,q3,c_star,v1,v2,v3,v4,v5,v11,v22
	DOUBLE PRECISION :: v33,v44,v55,s1,s2,s3,a1,a2,a3,w1,w2,w3,e

      do ib=1,nbp

       select case (op)
        case (1)  
          fi => dom(ib)%u
          np = dom(ib)%niul;  nq = dom(ib)%njul;  nr = dom(ib)%nkul
        case (2) 
          fi => dom(ib)%v
          np = dom(ib)%nivl;  nq = dom(ib)%njvl;  nr = dom(ib)%nkvl
        case (3)  
          fi => dom(ib)%w
          np = dom(ib)%niwl;  nq = dom(ib)%njwl;  nr = dom(ib)%nkwl
        case (14) ! (phi) 
          fi => dom(ib)%phi
          np = dom(ib)%ttc_i; nq = dom(ib)%ttc_j; nr = dom(ib)%ttc_k
        case (15) ! (phi_reinit)
          fi => dom(ib)%phi_reinit
          np = dom(ib)%ttc_i; nq = dom(ib)%ttc_j; nr = dom(ib)%ttc_k
        case (16) ! (phi_new)
          fi => dom(ib)%phi_new
          np = dom(ib)%ttc_i; nq = dom(ib)%ttc_j; nr = dom(ib)%ttc_k
        case (17) ! (phi_init)
          fi => dom(ib)%phi_init
          np = dom(ib)%ttc_i; nq = dom(ib)%ttc_j; nr = dom(ib)%ttc_k
       end select   

      if (dom(ib)%jprev.lt.0 .and. dom(ib)%bc_south.ne.5) then
       do k=1,nr
        do i=1,np
            fi(i,2,k)    = fi(i,3,k)
            fi(i,1,k)    = fi(i,2,k)
        end do
       end do   
      else if (dom(ib)%jnext.lt.0 .and. dom(ib)%bc_north.ne.5) then 
       do k=1,nr
        do i=1,np
            fi(i,nq-1,k) = fi(i,nq-2,k)
            fi(i,nq,k)   = fi(i,nq-1,k)
        end do
       end do   
      end if     
!
! Compute 1st order spatial derivative
!
      do k=1,nr 
       do i=1,np
        do j=1,nq-1
           dom(ib)%d1(i,j,k) = (fi(i,j+1,k) - fi(i,j,k))/dom(ib)%dy
        end do
       end do
      end do
!
! Compute 5th order derivative (-ve direction)
!
      dom(ib)%dphi_dyminus = 0.0

      do k=1,nr 
       do i=1,np
        do j=1,nq-6

           l=j-1

           v1 = dom(ib)%d1(i,l+1,k)
           v2 = dom(ib)%d1(i,l+2,k)
           v3 = dom(ib)%d1(i,l+3,k)
           v4 = dom(ib)%d1(i,l+4,k)
           v5 = dom(ib)%d1(i,l+5,k)

           v11 = v1**2
           v22 = v2**2
           v33 = v3**2
           v44 = v4**2
           v55 = v5**2

           e = (1.0e-6) * max(v11,v22,v33,v44,v55) + 1.0d-99

      s1 = (13.0/12.0)*(v1-2.0*v2+v3)**2 + 0.25*(v1-4.0*v2+3.0*v3)**2
      s2 = (13.0/12.0)*(v2-2.0*v3+v4)**2 + 0.25*(v2-v4)**2
      s3 = (13.0/12.0)*(v3-2.0*v4+v5)**2 + 0.25*(3.0*v3-4.0*v4+v5)**2

           a1 = 0.1/(e+s1)**2
           a2 = 0.6/(e+s2)**2
           a3 = 0.3/(e+s3)**2

           w1 = a1/(a1+a2+a3)
           w2 = a2/(a1+a2+a3)
           w3 = a3/(a1+a2+a3)

      dom(ib)%dphi_dyminus(i,j+3,k) = w1*(v1*(1.0/3.0)-v2*(7.0/6.0)+
     &v3*(11.0/6.0)) + w2*(-v2*(1.0/6.0)+v3*(5.0/6.0)+v4*(1.0/3.0))+
     &w3*(v3*(1.0/3.0)+v4*(5.0/6.0)-v5*(1.0/6.0))

        end do
       end do
      end do

      end do

      return
      end subroutine hj_weno_dyminus_3d
!######################################################################
      subroutine hj_weno_dzplus_3d(op)           
!######################################################################
      use vars
      use multidata
      implicit none
      integer :: i,j,k,l,k_star,np,nq,nr,op,ib
      double precision, pointer, dimension(:,:,:)::fi
      DOUBLE PRECISION :: e,q1,q2,q3,c_star,v1,v2,v3,v4,v5,v11,v22
	DOUBLE PRECISION :: v33,v44,v55,s1,s2,s3,a1,a2,a3,w1,w2,w3

      do ib=1,nbp

       select case (op)
        case (1)  
          fi => dom(ib)%u
          np = dom(ib)%niul;  nq = dom(ib)%njul;  nr = dom(ib)%nkul
        case (2) 
          fi => dom(ib)%v
          np = dom(ib)%nivl;  nq = dom(ib)%njvl;  nr = dom(ib)%nkvl
        case (3)  
          fi => dom(ib)%w
          np = dom(ib)%niwl;  nq = dom(ib)%njwl;  nr = dom(ib)%nkwl
        case (14) ! (phi) 
          fi => dom(ib)%phi
          np = dom(ib)%ttc_i; nq = dom(ib)%ttc_j; nr = dom(ib)%ttc_k
        case (15) ! (phi_reinit)
          fi => dom(ib)%phi_reinit
          np = dom(ib)%ttc_i; nq = dom(ib)%ttc_j; nr = dom(ib)%ttc_k
        case (16) ! (phi_new)
          fi => dom(ib)%phi_new
          np = dom(ib)%ttc_i; nq = dom(ib)%ttc_j; nr = dom(ib)%ttc_k
        case (17) ! (phi_init)
          fi => dom(ib)%phi_init
          np = dom(ib)%ttc_i; nq = dom(ib)%ttc_j; nr = dom(ib)%ttc_k
       end select

      if (dom(ib)%kprev.lt.0 .and. dom(ib)%bc_bottom.ne.5) then
       do i=1,np
        do j=1,nq
            fi(i,j,2)    = fi(i,j,3)
            fi(i,j,1)    = fi(i,j,2)
        end do
       end do    
      else if (dom(ib)%knext.lt.0 .and. dom(ib)%bc_top.ne.5) then 
       do i=1,np
        do j=1,nq
            fi(i,j,nr-1) = fi(i,j,nr-2)
            fi(i,j,nr)   = fi(i,j,nr-1)
        end do
       end do    
      end if
!
! Compute 1st order spatial derivative
!
      do k=1,nr-1
       do i=1,np
        do j=1,nq
           dom(ib)%d1(i,j,k) = (fi(i,j,k+1) - fi(i,j,k))/dom(ib)%dz
        end do
       end do
      end do 
!
! Compute 5th order derivative (-ve direction)
!
      dom(ib)%dphi_dzplus = 0.0

      do k=1,nr-6
       do i=1,np
        do j=1,nq

            l=k

           v1 = dom(ib)%d1(i,j,l+5)
           v2 = dom(ib)%d1(i,j,l+4)
           v3 = dom(ib)%d1(i,j,l+3)
           v4 = dom(ib)%d1(i,j,l+2)
           v5 = dom(ib)%d1(i,j,l+1)

           v11 = v1**2
           v22 = v2**2
           v33 = v3**2
           v44 = v4**2
           v55 = v5**2

           e = (1.0e-6) * max(v11,v22,v33,v44,v55) + 1.0d-99

      s1 = (13.0/12.0)*(v1-2.0*v2+v3)**2 + 0.25*(v1-4.0*v2+3.0*v3)**2
      s2 = (13.0/12.0)*(v2-2.0*v3+v4)**2 + 0.25*(v2-v4)**2 
      s3 = (13.0/12.0)*(v3-2.0*v4+v5)**2 + 0.25*(3.0*v3-4.0*v4+v5)**2

           a1 = 0.1/(e+s1)**2
           a2 = 0.6/(e+s2)**2
           a3 = 0.3/(e+s3)**2

           w1 = a1/(a1+a2+a3)
           w2 = a2/(a1+a2+a3)
           w3 = a3/(a1+a2+a3)

      dom(ib)%dphi_dzplus(i,j,k+3) = w1*(v1*(1.0/3.0)-v2*(7.0/6.0)+
     &v3*(11.0/6.0)) + w2*(-v2*(1.0/6.0)+v3*(5.0/6.0)+v4*(1.0/3.0))+
     &w3*(v3*(1.0/3.0)+v4*(5.0/6.0)-v5*(1.0/6.0))

        end do
       end do
      end do

      end do

      return
      end subroutine hj_weno_dzplus_3d
!######################################################################
      subroutine hj_weno_dzminus_3d(op)        
!######################################################################
      use vars
      use multidata
      implicit none
      integer :: i,j,k,l,k_star,np,nq,nr,op,ib
      double precision, pointer, dimension(:,:,:)::fi
      DOUBLE PRECISION :: q1,q2,q3,c_star,v1,v2,v3,v4,v5,v11,v22
	DOUBLE PRECISION :: v33,v44,v55,s1,s2,s3,a1,a2,a3,w1,w2,w3,e

      do ib=1,nbp

       select case (op)
        case (1)  
          fi => dom(ib)%u
          np = dom(ib)%niul;  nq = dom(ib)%njul;  nr = dom(ib)%nkul
        case (2) 
          fi => dom(ib)%v
          np = dom(ib)%nivl;  nq = dom(ib)%njvl;  nr = dom(ib)%nkvl
        case (3)  
          fi => dom(ib)%w
          np = dom(ib)%niwl;  nq = dom(ib)%njwl;  nr = dom(ib)%nkwl
        case (14)  ! (phi)
          fi => dom(ib)%phi
          np = dom(ib)%ttc_i; nq = dom(ib)%ttc_j; nr = dom(ib)%ttc_k
        case (15)  ! (phi_reinit)
          fi => dom(ib)%phi_reinit
          np = dom(ib)%ttc_i; nq = dom(ib)%ttc_j; nr = dom(ib)%ttc_k
        case (16)  ! (phi_new)
          fi => dom(ib)%phi_new
          np = dom(ib)%ttc_i; nq = dom(ib)%ttc_j; nr = dom(ib)%ttc_k
        case (17)  ! (phi_init)
          fi => dom(ib)%phi_init
          np = dom(ib)%ttc_i; nq = dom(ib)%ttc_j; nr = dom(ib)%ttc_k
       end select   

      if (dom(ib)%kprev.lt.0 .and. dom(ib)%bc_bottom.ne.5) then
       do i=1,np
        do j=1,nq
            fi(i,j,2)    = fi(i,j,3)
            fi(i,j,1)    = fi(i,j,2)
        end do
       end do    
      else if (dom(ib)%knext.lt.0 .and. dom(ib)%bc_top.ne.5) then 
       do i=1,np
        do j=1,nq
            fi(i,j,nr-1) = fi(i,j,nr-2)
            fi(i,j,nr)   = fi(i,j,nr-1)
        end do
       end do    
      end if  
!
! Compute 1st order spatial derivative
!
      do k=1,nr-1
       do i=1,np
        do j=1,nq
           dom(ib)%d1(i,j,k) = (fi(i,j,k+1) - fi(i,j,k))/dom(ib)%dz
        end do
       end do
      end do
!
! Compute 5th order derivative (-ve direction)
!
      dom(ib)%dphi_dzminus = 0.0

      do k=1,nr-6
       do i=1,np
        do j=1,nq

           l=k-1

           v1 = dom(ib)%d1(i,j,l+1)
           v2 = dom(ib)%d1(i,j,l+2)
           v3 = dom(ib)%d1(i,j,l+3)
           v4 = dom(ib)%d1(i,j,l+4)
           v5 = dom(ib)%d1(i,j,l+5)

           v11 = v1**2
           v22 = v2**2
           v33 = v3**2
           v44 = v4**2
           v55 = v5**2

           e = (1.0e-6) * max(v11,v22,v33,v44,v55) + 1.0d-99

      s1 = (13.0/12.0)*(v1-2.0*v2+v3)**2 + 0.25*(v1-4.0*v2+3.0*v3)**2
      s2 = (13.0/12.0)*(v2-2.0*v3+v4)**2 + 0.25*(v2-v4)**2
      s3 = (13.0/12.0)*(v3-2.0*v4+v5)**2 + 0.25*(3.0*v3-4.0*v4+v5)**2

           a1 = 0.1/(e+s1)**2
           a2 = 0.6/(e+s2)**2
           a3 = 0.3/(e+s3)**2

           w1 = a1/(a1+a2+a3)
           w2 = a2/(a1+a2+a3)
           w3 = a3/(a1+a2+a3)

      dom(ib)%dphi_dzminus(i,j,k+3) = w1*(v1*(1.0/3.0)-v2*(7.0/6.0)+
     &v3*(11.0/6.0)) + w2*(-v2*(1.0/6.0)+v3*(5.0/6.0)+v4*(1.0/3.0))+
     &w3*(v3*(1.0/3.0)+v4*(5.0/6.0)-v5*(1.0/6.0))

        end do
       end do
      end do

      end do

      return
      end subroutine hj_weno_dzminus_3d
!#####################################################################
      subroutine heaviside
!######################################################################
      use vars
      use lsm
      use mpi
      use multidata
      implicit none
      integer :: i,j,k,ib,tti,ttj,ttk
      DOUBLE PRECISION :: epsl,n_epsl
      DOUBLE PRECISION, parameter :: pi = 3.14159265359D0
	DOUBLE PRECISION :: vol !!
!	
! Define an infinitely differentiable smoothed heaviside function h_phi
!
      n_epsl = 1.9999999d0		! Transition zone across free surface (2 grid cells width either side)
	vol=0.d0 ;	total_vol=0.d0!!

      do ib=1,nbp  

       epsl = n_epsl*max(dom(ib)%dx,dom(ib)%dy,dom(ib)%dz) 

        do k=dom(ib)%ksp-pl,dom(ib)%kep+pl
         do i=dom(ib)%isp-pl,dom(ib)%iep+pl 
          do j=dom(ib)%jsp-pl,dom(ib)%jep+pl

          if(dom(ib)%phi(i,j,k).lt.(-1.0*epsl))dom(ib)%h_phi(i,j,k)=0.d0  ! h_phi=0 above free surface

          if(dom(ib)%phi(i,j,k).gt.epsl) dom(ib)%h_phi(i,j,k) = 1.d0  	! h_phi=1.0 below free surface

          if(DABS(dom(ib)%phi(i,j,k)).le.epsl) then

          dom(ib)%h_phi(i,j,k)=0.5d0*(1.d0+dom(ib)%phi(i,j,k)/epsl
     & 				+1.d0/pi*DSIN(pi*dom(ib)%phi(i,j,k)/epsl))    

	    vol=vol+dom(ib)%h_phi(i,j,k)*dom(ib)%dx*dom(ib)%dy*dom(ib)%dz   !volume of liquid
          end if

        dom(ib)%dens(i,j,k)=densg+(densl-densg)*dom(ib)%h_phi(i,j,k) !dens = densg above free surface, densl below
        dom(ib)%mu(i,j,k)=mug+(mul-mug)*dom(ib)%h_phi(i,j,k) 	   !mu = mug above free surface, mul below

          end do
         end do
        end do


      end do

	call MPI_ALLREDUCE(vol,total_vol,1,
     &          MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr )  	!!
	
	if(itime.eq.itime_start) initial_vol=total_vol

	if(myrank.eq.0 .and. itime.ne.itime_start)then
     	 WRITE(1703,'(3f15.7)') ctime,total_vol,total_vol/initial_vol
	endif

      call bound_lsm(18) !(dens) 
      call bound_lsm(19) !(mu) 

      return
      end subroutine heaviside
