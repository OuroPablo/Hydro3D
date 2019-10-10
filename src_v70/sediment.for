!=======================================================================
!    			Scalar transport
!			Yan Liu
!			Bruño Fraga Bugallo
!			Cardiff 2015-2016
!=======================================================================
!##########################################################################
        subroutine sediment_init
!##########################################################################

	use vars
	use mpi
	use multidata
	implicit none
	integer :: ib,i,j,k
	double precision :: Vcell

	do ib=1,nbp

!           do k=dom(ib)%ksp,dom(ib)%kep
!           do i=dom(ib)%isp,dom(ib)%iep
!           do j=dom(ib)%jsp,dom(ib)%jep
		if (dom_id(ib).eq.18) then
		Vcell = dom(ib)%dx*dom(ib)%dy*dom(ib)%dz	
!		if (dom(ib)%x(i).ge.0.49.and.dom(ib)%x(i).le.0.51) then
!		if (dom(ib)%y(j).ge.0.49.and.dom(ib)%y(j).le.0.51) then
!		if (dom(ib)%z(k).ge.0.085.and.dom(ib)%z(k).le.0.09) then
		dom(ib)%S(dom(ib)%isp,dom(ib)%jsp,26)=
     & dom(ib)%S(dom(ib)%isp,dom(ib)%jsp,26) + (3.41d-6/Vcell)*dt
!				goto 25
!		endif ; endif ; endif
		endif
!		enddo ; enddo ; enddo
	enddo
!	25 continue
!	do ib=1, nbp
!	dom(ib)%ibfactor = 1

!           do k=dom(ib)%ksp,dom(ib)%kep
!              do i=dom(ib)%isp,dom(ib)%iep
!                 do j=dom(ib)%jsp,dom(ib)%jep

!	if (dom(ib)%zc(k) .le. 4.8) then

!	x_sphere=int(dom(ib)%xc(i)/0.96)*0.96+0.48 
!	y_sphere=int(dom(ib)%yc(j)/0.96)*0.96+0.48
!	z_sphere=int(dom(ib)%zc(k)/0.96)*0.96+0.48
!	distance=((dom(ib)%xc(i)-x_sphere)**2+(dom(ib)%yc(j)-y_sphere)**2
!     &+(dom(ib)%zc(k)-z_sphere)**2)**0.5
!		if(distance .gt. 0.48) then
!		dom(ib)%ibfactor(i,j,k)=1
!		else
!		dom(ib)%ibfactor(i,j,k)=0
!		end if
!	else if (dom(ib)%zc(k) .gt. 4.8) then
!	dom(ib)%ibfactor(i,j,k)=1
!	end if
!		 end do
!	      end do
!	   end do
!	end do

	end

!##########################################################################
        subroutine sediment_4thtest
!##########################################################################
        use vars
        use mpi
        use multidata
        implicit none
        integer :: i,j,k,ib
        double precision :: dxx,dyy,dzz
        double precision :: conv,diff
        double precision :: duSdx,dvSdy,dwSdz,dwsSdz,dSdt
        double precision :: awS,aeS,asS,anS,ab_S,atS,apS
	double precision :: masout
	double precision :: ues,uws,vns,vss,wts,wbs,wsbs
	double precision :: masconv,masdiff,masws
	double precision :: kp,km,ku,kc,kd,b_r,ws,Vcell

        do ib=1,nbp

 	  Vcell = dom(ib)%dx*dom(ib)%dy*dom(ib)%dz	

	  if (dom_id(ib).eq.9) then							!release point for bubble plume
	  	dom(ib)%So(dom(ib)%iep,dom(ib)%jep,30)=
     & 	dom(ib)%So(dom(ib)%iep,dom(ib)%jep,30) + (3.41d-6/Vcell)*dt
	  	dom(ib)%So(dom(ib)%iep,dom(ib)%jep,31)=
     & 	dom(ib)%So(dom(ib)%iep,dom(ib)%jep,31) + (3.41d-6/Vcell)*dt
!		write(6,*)'1',dom(ib)%S(dom(ib)%iep,dom(ib)%jep,29)
!		write(6,*)dom(ib)%S(dom(ib)%iep,dom(ib)%jep,28)
!		write(6,*)Vcell,dt
	  elseif (dom_id(ib).eq.10) then
	  	dom(ib)%So(dom(ib)%isp,dom(ib)%jep,30)=
     & 	dom(ib)%So(dom(ib)%isp,dom(ib)%jep,30) + (3.41d-6/Vcell)*dt
	  	dom(ib)%So(dom(ib)%isp,dom(ib)%jep,31)=
     & 	dom(ib)%So(dom(ib)%isp,dom(ib)%jep,31) + (3.41d-6/Vcell)*dt
!		write(6,*)'2',dom(ib)%S(dom(ib)%isp,dom(ib)%jep,29)
!		write(6,*)dom(ib)%S(dom(ib)%isp,dom(ib)%jep,28)
!		write(6,*)Vcell,dt
	  elseif (dom_id(ib).eq.17) then
	  	dom(ib)%So(dom(ib)%iep,dom(ib)%jsp,30)=
     & 	dom(ib)%So(dom(ib)%iep,dom(ib)%jsp,30) + (3.41d-6/Vcell)*dt
	  	dom(ib)%So(dom(ib)%iep,dom(ib)%jsp,31)=
     & 	dom(ib)%So(dom(ib)%iep,dom(ib)%jsp,31) + (3.41d-6/Vcell)*dt
!		write(6,*)'3',dom(ib)%S(dom(ib)%iep,dom(ib)%jsp,29)
!		write(6,*)dom(ib)%S(dom(ib)%iep,dom(ib)%jsp,28)
!		write(6,*)Vcell,dt
	  elseif (dom_id(ib).eq.18) then
	  	dom(ib)%So(dom(ib)%isp,dom(ib)%jsp,30)=
     & 	dom(ib)%So(dom(ib)%isp,dom(ib)%jsp,30) + (3.41d-6/Vcell)*dt
	  	dom(ib)%So(dom(ib)%isp,dom(ib)%jsp,31)=
     & 	dom(ib)%So(dom(ib)%isp,dom(ib)%jsp,31) + (3.41d-6/Vcell)*dt
!		write(6,*)'4',dom(ib)%S(dom(ib)%isp,dom(ib)%jsp,29)
!		write(6,*)dom(ib)%S(dom(ib)%isp,dom(ib)%jsp,28)
!		write(6,*)Vcell,dt
        endif


	if (itime .eq. itime_start) then
	dom(ib)%sfactor = 1.0
	end if 

        dxx=dom(ib)%dx*dom(ib)%dx
        dyy=dom(ib)%dy*dom(ib)%dy
        dzz=dom(ib)%dz*dom(ib)%dz

           do k=dom(ib)%ksp,dom(ib)%kep
              do i=dom(ib)%isp,dom(ib)%iep
                 do j=dom(ib)%jsp,dom(ib)%jep
!-------Convection
	if(dom(ib)%u(i-1,j,k).gt.0.0) then
	ku=dom(ib)%So(i-2,j,k)
	kc=dom(ib)%So(i-1,j,k)
	kd=dom(ib)%So(i,j,k)
	b_r=max(0.0,
     & min(2.0*((kd-kc)/(kc-ku)),0.75*((kd-kc)/(kc-ku))+0.25,4.0))
	km=(kc+0.5*b_r*(kc-ku))
	else if(dom(ib)%u(i-1,j,k).lt.0.0) then
	ku=dom(ib)%So(i+1,j,k)
	kc=dom(ib)%So(i,j,k)
	kd=dom(ib)%So(i-1,j,k)
	b_r=max(0.0,
     & min(2.0*((kd-kc)/(kc-ku)),0.75*((kd-kc)/(kc-ku))+0.25,4.0))
	km=(kc+0.5*b_r*(kc-ku))
	else
	km=0.5*(dom(ib)%So(i,j,k)+dom(ib)%So(i-1,j,k))
	end if

	if(dom(ib)%u(i,j,k).gt.0.0) then
	ku=dom(ib)%So(i-1,j,k)
	kc=dom(ib)%So(i,j,k)
	kd=dom(ib)%So(i+1,j,k)
	b_r=max(0.0,
     & min(2.0*((kd-kc)/(kc-ku)),0.75*((kd-kc)/(kc-ku))+0.25,4.0))
	kp=(kc+0.5*b_r*(kc-ku))
        else if(dom(ib)%u(i,j,k).lt.0.0) then
           ku=dom(ib)%So(i+2,j,k)
           kc=dom(ib)%So(i+1,j,k)
           kd=dom(ib)%So(i,j,k)
           b_r=max(0.0,
     & min(2.0*((kd-kc)/(kc-ku)),0.75*((kd-kc)/(kc-ku))+0.25,4.0))
           kp=(kc+0.5*b_r*(kc-ku))
        else
           kp=0.5*(dom(ib)%So(i,j,k)+dom(ib)%So(i+1,j,k))
        end if
        duSdx=(dom(ib)%u(i,j,k)*kp-dom(ib)%u(i-1,j,k)*km)/dom(ib)%dx
!------
        if(dom(ib)%v(i,j-1,k).gt.0.0) then
           ku=dom(ib)%So(i,j-2,k)
           kc=dom(ib)%So(i,j-1,k)
           kd=dom(ib)%So(i,j,k)
           b_r=max(0.0,
     & min(2.0*((kd-kc)/(kc-ku)),0.75*((kd-kc)/(kc-ku))+0.25,4.0))
           km=(kc+0.5*b_r*(kc-ku))
        else if(dom(ib)%v(i,j-1,k).lt.0.0) then
           ku=dom(ib)%So(i,j+1,k)
           kc=dom(ib)%So(i,j,k)
           kd=dom(ib)%So(i,j-1,k)
           b_r=max(0.0,
     & min(2.0*((kd-kc)/(kc-ku)),0.75*((kd-kc)/(kc-ku))+0.25,4.0))
           km=(kc+0.5*b_r*(kc-ku))
        else
           km=0.5*(dom(ib)%So(i,j,k)+dom(ib)%So(i,j-1,k))
        end if

        if(dom(ib)%v(i,j,k).gt.0.0) then
           ku=dom(ib)%So(i,j-1,k)
           kc=dom(ib)%So(i,j,k)
           kd=dom(ib)%So(i,j+1,k)
           b_r=max(0.0,
     & min(2.0*((kd-kc)/(kc-ku)),0.75*((kd-kc)/(kc-ku))+0.25,4.0))
           kp=(kc+0.5*b_r*(kc-ku))

        else if(dom(ib)%v(i,j,k).lt.0.0) then
           ku=dom(ib)%So(i,j+2,k)
           kc=dom(ib)%So(i,j+1,k)
           kd=dom(ib)%So(i,j,k)
           b_r=max(0.0,
     & min(2.0*((kd-kc)/(kc-ku)),0.75*((kd-kc)/(kc-ku))+0.25,4.0))
           kp=(kc+0.5*b_r*(kc-ku))
        else
           kp=0.5*(dom(ib)%So(i,j,k)+dom(ib)%So(i,j+1,k))
        end if
        dvSdy=(dom(ib)%v(i,j,k)*kp-dom(ib)%v(i,j-1,k)*km)/dom(ib)%dy
!-------
        if(dom(ib)%w(i,j,k-1).gt.0.0) then
           ku=dom(ib)%So(i,j,k-2)
           kc=dom(ib)%So(i,j,k-1)
           kd=dom(ib)%So(i,j,k)
           b_r=max(0.0,
     & min(2.0*((kd-kc)/(kc-ku)),0.75*((kd-kc)/(kc-ku))+0.25,4.0))
           km=(kc+0.5*b_r*(kc-ku))

        else if(dom(ib)%w(i,j,k-1).lt.0.0) then
           ku=dom(ib)%So(i,j,k+1)
           kc=dom(ib)%So(i,j,k)
           kd=dom(ib)%So(i,j,k-1)
           b_r=max(0.0,
     & min(2.0*((kd-kc)/(kc-ku)),0.75*((kd-kc)/(kc-ku))+0.25,4.0))
           km=(kc+0.5*b_r*(kc-ku))
        else
           km=0.5*(dom(ib)%So(i,j,k)+dom(ib)%So(i,j,k-1))
        end if
        if(dom(ib)%w(i,j,k).gt.0.0) then
           ku=dom(ib)%So(i,j,k-1)
           kc=dom(ib)%So(i,j,k)
           kd=dom(ib)%So(i,j,k+1)
           b_r=max(0.0,
     & min(2.0*((kd-kc)/(kc-ku)),0.75*((kd-kc)/(kc-ku))+0.25,4.0))
           kp=(kc+0.5*b_r*(kc-ku))

        else if(dom(ib)%w(i,j,k).lt.0.0) then
           ku=dom(ib)%So(i,j,k+2)
           kc=dom(ib)%So(i,j,k+1)
           kd=dom(ib)%So(i,j,k)
           b_r=max(0.0,
     & min(2.0*((kd-kc)/(kc-ku)),0.75*((kd-kc)/(kc-ku))+0.25,4.0))
           kp=(kc+0.5*b_r*(kc-ku))
        else
           kp=0.5*(dom(ib)%So(i,j,k)+dom(ib)%So(i,j,k+1))
        end if
        dwSdz=(dom(ib)%w(i,j,k)*kp-dom(ib)%w(i,j,k-1)*km)/dom(ib)%dz

        conv=(duSdx+dvSdy+dwSdz)

!-------Diffusion
                 awS=-dom(ib)%vis(i,j,k)/(dxx*Sc_t)
                 aeS=-dom(ib)%vis(i,j,k)/(dxx*Sc_t)
                 anS=-dom(ib)%vis(i,j,k)/(dyy*Sc_t)
                 asS=-dom(ib)%vis(i,j,k)/(dyy*Sc_t)
                 atS=-dom(ib)%vis(i,j,k)/(dzz*Sc_t)
                 ab_S=-dom(ib)%vis(i,j,k)/(dzz*Sc_t)

                 apS = -1.0*(awS+aeS+asS+anS+ab_S+atS)

        	diff=(apS*dom(ib)%So(i,j,k)+
     & anS*dom(ib)%So(i,j+1,k) + asS*dom(ib)%So(i,j-1,k)+
     & aeS*dom(ib)%So(i+1,j,k) + awS*dom(ib)%So(i-1,j,k)+
     & atS*dom(ib)%So(i,j,k+1) + ab_S*dom(ib)%So(i,j,k-1))

		dom(ib)%S(i,j,k)=dom(ib)%So(i,j,k)-dt*(conv+diff)

!	if (dom(ib)%S(i,j,k) .lt. 0.0) then
!	write (81,*) dom(ib)%S(i,j,k), dom(ib)%So(i,j,k)
!	write (81,*) i,j,k
!	write (81,*) conv,diff
!	write (81,*) duSdx,dvSdy,dwSdz
!	end if

	end do
	end do
	end do

        end do

	do ib =1, nbp
           do k=dom(ib)%ksp,dom(ib)%kep
              do i=dom(ib)%isp,dom(ib)%iep
                 do j=dom(ib)%jsp,dom(ib)%jep
	if (dom(ib)%S(i,j,k) .gt. 100) then
	write(6,*)'ERROR: scalar too big'
	stop
	end if
	end do;end do;end do;end do

        call exchange(8)

        call boundS

        return
        end subroutine sediment_4thtest
!##########################################################################
        subroutine boundS
!##########################################################################
	use mpi
        use vars
        use multidata
        implicit none
        integer :: i,j,k,ib,ni,nj,nk,ly
        integer :: is,ie,js,je,ks,ke
	double precision :: absz,absy

        if (PERIODIC) call exchange_bc(8,pl_ex)

        do ly=0,pl_ex

        do ib=1,nbp
           ni=dom(ib)%ttc_i; nj=dom(ib)%ttc_j; nk=dom(ib)%ttc_k
           is=dom(ib)%isp; ie=dom(ib)%iep
           js=dom(ib)%jsp; je=dom(ib)%jep
           ks=dom(ib)%ksp; ke=dom(ib)%kep
	
! Boundary Conditions for S
!..............................................................................
!=== West ===>   ..  4=wall  ..    1=Inflow
!..............................................................................
        if (dom(ib)%iprev.lt.0) then
           if (dom(ib)%Tbc_west.eq.4) then
              do k=ks-1,ke+1; do j=js-1,je+1
                 dom(ib)%S(is-1-ly,j,k)= dom(ib)%S(is+ly,j,k)
              end do; end do

           else if (dom(ib)%Tbc_west.eq.1) then					!CHANGE
              do k=ks-1,ke+1; do j=js-1,je+1
                 dom(ib)%S(is-1-ly,j,k)= 1.0
              end do; end do
	   else if (dom(ib)%Tbc_west.eq. 11) then
	      do k=ks-1,ke+1; do j=js-1,je+1
	      absz=abs(dom(ib)%zc(k)-5.22)
	      absy=abs(dom(ib)%yc(j)-4.32)
	      if (absz .le. 0.5*dom(ib)%dz) then 
	      	if (absy .le. 0.5*dom(ib)%dy) then
		   dom(ib)%S(is-1-ly,j,k)= 1.0
	      	else 
	           dom(ib)%S(is-1-ly,j,k)= 0.0
	      	end if
	      else 
	         dom(ib)%S(is-1-ly,j,k)= 0.0
	      end if	
	      end do; end do
           end if
        end if
!...............................................................................
!=== East ===>   ..  4=wall  ..    2=Outflow
!...............................................................................
        if (dom(ib)%inext.lt.0) then
           if (dom(ib)%Tbc_east.eq.4) then
              do k=ks-1,ke+1; do j=js-1,je+1
                 dom(ib)%S(ie+1+ly,j,k)= dom(ib)%S(ie-ly,j,k)	
              end do; end do

           else if (dom(ib)%Tbc_east.eq.2) then
              do k=ks-1,ke+1; do j=js-1,je+1
                 dom(ib)%S(ie+1+ly,j,k)= dom(ib)%S(ie-ly,j,k)
              end do; end do
           end if
        end if
!...............................................................................
!=== South ===>  ..  4=wall  ..      
!...............................................................................
        if (dom(ib)%jprev.lt.0) then
           if (dom(ib)%Tbc_south.eq.4) then 
              do k=ks-1,ke+1; do i=is-1,ie+1
                 dom(ib)%S(i,js-1-ly,k)= dom(ib)%S(i,js+ly,k)	
              end do; end do

           end if
        end if
!.............................................................................
!=== North ===>  ..  4=wall  ..    
!.............................................................................
        if (dom(ib)%jnext.lt.0) then
           if (dom(ib)%Tbc_north.eq.4) then
              do k=ks-1,ke+1; do i=is-1,ie+1
                 dom(ib)%S(i,je+1+ly,k) = dom(ib)%S(i,je-ly,k) 	
              end do; end do

           end if
        end if
!...............................................................................
!=== Bottom ===>  ..  6=Net deposition  ..   7=Erosion
!...............................................................................
        if (dom(ib)%kprev.lt.0) then
           if (dom(ib)%Tbc_bottom.eq.6) then
              do j=js-1,je+1; do i=is-1,ie+1
                 dom(ib)%S(i,j,ks-1-ly)= 0.0	
              end do; end do

           else if (dom(ib)%Tbc_bottom.eq.7) then
              do j=js-1,je+1; do i=is-1,ie+1
                 dom(ib)%S(i,j,ks-1-ly)= 1.5 * dom(ib)%S(i,j,ks+ly)	
              end do; end do
           end if
        end if
!.............................................................................
!=== Top ===>  ..  8=free surface
!.............................................................................
        if (dom(ib)%knext.lt.0) then
           if (dom(ib)%Tbc_top.eq.8) then
              do j=js-1,je+1; do i=is-1,ie+1
                 dom(ib)%S(i,j,ke+1+ly) = 0.0	
              end do; end do

           end if
        end if

!==============================================================================
        end do

	end do ! ly

        end subroutine boundS
!#############################################################################
