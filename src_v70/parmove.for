!##########################################################################
        module module_parmove
!##########################################################################
	  double precision a_p,b_p,c_p,sigy_XZ,umax,vmax,wmax
	  integer :: n_p,n_pnum,XZunit
	  double precision x_p(300),y_p(300),z_p(300)
        double precision u_p(300),v_p(300),w_p(300)
        double precision u_p0(300),v_p0(300),w_p0(300)
        double precision phi_xp(300),phi_yp(300),phi_zp(300)
        double precision omga_xp(300),omga_yp(300),omga_zp(300)
        double precision Fdx_p(300),Fdy_p(300),Fdz_p(300)
        double precision fx_p(300),fy_p(300),fz_p(300)
        double precision Fl_p(300),Fm_p(300),ff_p(300),FB_p(300)
	  double precision,allocatable :: dh_loc(:,:,:),I_nr(:,:,:)
	  double precision,allocatable :: J_nr(:,:,:),K_nr(:,:,:)
	  double precision,allocatable :: u_int(:),v_int(:),w_int(:)
        end module

!##########################################################################
        subroutine particle_ini
!##########################################################################
       use module_parmove
       use mpi
       use vars
       use multidata
       implicit none 
       integer ::i,sn,l
       CHARACTER*8  :: chb
       CHARACTER*25 :: gf

! read 'particle.cin' to get the initial information of particles
       open(unit=1,file='particle.cin')
       read (1,*) n_p
       read (1,*) a_p,b_p,c_p
       do i=1,n_p
       read (1,*) x_p(i),y_p(i),z_p(i),
     & u_p(i),v_p(i),w_p(i)
     & ,phi_yp(i),omga_yp(i)
       end do

!open forces file to write
       open(unit=5555,file='debug.dat')
       write (5555,*)
     & 'variables=t,x,y,z,u,v,w,phi,o,Fdx,
     &Fdy,Fdz,Fl,Fm,G,ff,FB,fx,fy,fz'

! open output file 'loc_p*.dat' to prepare writing       
        if (myrank.eq.0) then 
        do l=1,n_p
        n_pnum=l*500+38
        write(chb,'(i8)') l
        sn=len(trim(adjustl(chb)))
        chb=repeat('0',(4-sn))//trim(adjustl(chb))
        gf='loc_p'//trim(adjustl(chb))//'.dat'
        open(unit=n_pnum,file=gf)
        write(n_pnum,*) 'title=', 'loc_p'
        !write(n_pnum,*) 'variables="t","x","y","z","u","v","w"
    ! & "phx","phy","phz","wx","wy","wz"'
        !write(n_pnum,88) 0.0,x_p(l),y_p(l),z_p(l),
    ! & u_p(l),v_p(l),w_p(l),phi_xp(l),phi_yp(l),phi_zp(l),
    ! & omga_xp(l),omga_yp(l),omga_zp(l)
        write(n_pnum,*) 'variables="t","x","y","z","u","v","w",
     & "phy","wy"'
        write(n_pnum,88) 0.0,x_p(l),y_p(l),z_p(l),u_p(l),v_p(l),w_p(l)
     & ,phi_yp(l),omga_yp(l)
       end do
       end if

88      format (15e19.8)
       end
        
!##########################################################################
        subroutine particle_move
!##########################################################################
       use module_parmove
	use mpi
	use vars
	use multidata
        implicit none
	integer :: ib,ni,nj,nk,ks,ke,i,j,k,l,nxdom,nydom,nzdom
	integer :: is,ie,js,je,strlen,sn,ndoms,nway,nl
	double precision ::lx,ly,lz,m_p,G_p,a1,a2,a3,zd,pi,dstard,thetcr
      double precision ::k1,du,dv,dw,REp,Cd,Cl,Fdx,Fdy,Fdz,Fl,ff,Fd,tt
      double precision ::ad_p,bd_p,cd_p,du1,dw1,Fdx1,Fdz1,u_ref,d_vel
      double precision ::b_a,c_a,dvdx,dwdx,dudy,dwdy,dudz,dvdz,kk,u_s
      double precision ::omga_x,omga_y,omga_z,S_xy,S_xz,S_yz,Fm,Rew,Cw
      double precision ::Ix,Iy,Iz,FB,nxl,dh,dist,l_dist
      integer,allocatable,dimension(:)::locdom_p
!      double precision,allocatable,dimension(:)::fx_p,fy_p,fz_p
      CHARACTER*8  :: chb
      CHARACTER*25 :: gf
      
	if (itime.eq.itime_start) then
		allocate(dh_loc(n_p,27,3)) ; dh_loc=0.d0	!Interpolation function	
		allocate(I_nr(n_p,1,3)) ; I_nr=0.d0		!1-U-grid ;2-V-grid; 3-W-grid
		allocate(J_nr(n_p,1,3)) ; J_nr=0.d0		!
		allocate(K_nr(n_p,1,3)) ; K_nr=0.d0		!
	      allocate(u_int(n_p),v_int(n_p),w_int(n_p))
	endif	

      allocate (locdom_p(n_p))
!      allocate (fx_p(n_p),fy_p(n_p),fz_p(n_p))

! calculating some parameters           
           u_ref=0.316 
           u_s=0.024
           pi=4.D0*DATAN(1.D0)
           ad_p=a_p*10.0
           bd_p=b_p*10.0
           cd_p=c_p*10.0
           G_p=9.81*1.65*pi*a_p*b_p*c_p/6.0*(1.0/u_ref*10.0)**2.0
           m_p=1.0/6.0*2.65*pi*ad_p*bd_p*cd_p
           b_a=b_p/a_p
           c_a=c_p/a_p
           zd=0.01735/4.0	!0.5*cd_p 			!CHANGED PABLO
            kk=9.0/(1.0+1.65+0.5)/pi**0.5/
     & ((ad_p+bd_p+cd_p)/3.0*u_s/u_ref*Re)**0.5
           lx=xen-xst
           ly=yen-yst
           lz=zen-zst
           ndoms=idom*jdom*kdom
                   
! loop in domains to find particles in which domains
	do ib=1,nbp
	   ni=dom(ib)%ttc_i; nj=dom(ib)%ttc_j; nk=dom(ib)%ttc_k
           ks=dom(ib)%ksp; ke=dom(ib)%kep
           js=dom(ib)%jsp; je=dom(ib)%jep
           is=dom(ib)%isp; ie=dom(ib)%iep 

  ! loop in particles  
       do l=1,n_p
           ff=0.0
           nxdom=INT((x_p(l)-xst)/(lx/idom))
           nydom=INT((y_p(l)-yst)/(ly/jdom))
           nzdom=INT((z_p(l)-zst)/(lz/kdom))
           locdom_p(l)=idom*jdom*nzdom+idom*nydom+nxdom

! record particles velocity at the beginning of this time loop
         u_p0(l)=u_p(l)
         v_p0(l)=v_p(l)
         w_p0(l)=w_p(l)

! find the closest uvw for particles and calculate slip velocity 
       if(locdom_p(l).eq.dom_id(ib)) then     !checking in which domain 

	 nway=2
	if(nway.eq.1) then
!NEW WAY TO FIND THE CLOSEST EULERIAN CELL
! u      
        do k=ks-1,ke 
           if ((z_p(l).ge.dom(ib)%z(k-1)).and.
     &         (z_p(l).lt.dom(ib)%z(k))) then  
         do j=js-1,je 
           if ((y_p(l).ge.dom(ib)%y(j-1)).and.
     &         (y_p(l).lt.dom(ib)%y(j))) then
          do i=is-1,ie
           if ((x_p(l).ge.dom(ib)%xc(i)).and.
     &         (x_p(l).lt.dom(ib)%xc(i+1))) then
         
           if(dom(ib)%u(i,j,k).lt.0.0) then   !x,y,z axis
             du=0.0-u_p(l)
           else
             du=dom(ib)%u(i,j,k)-u_p(l) 
           end if 
           end if
          end do
           end if
         end do
           end if 
        end do  
! v
       do k=ks-1,ke 
           if ((z_p(l).ge.dom(ib)%z(k-1)).and.
     &         (z_p(l).lt.dom(ib)%z(k))) then 
         do j=js-1,je 
           if ((y_p(l).ge.dom(ib)%yc(j)).and.
     &         (y_p(l).lt.dom(ib)%yc(j+1))) then
          do i=is-1,ie
           if ((x_p(l).ge.dom(ib)%x(i-1)).and.
     &         (x_p(l).lt.dom(ib)%x(i))) then

           dv=dom(ib)%v(i,j,k)-v_p(l)

           end if
          end do
           end if
         end do
           end if 
        end do          
 ! w
       do k=ks-1,ke 
           if ((z_p(l).ge.dom(ib)%zc(k)).and.
     &         (z_p(l).lt.dom(ib)%zc(k+1))) then  
         do j=js-1,je 
           if ((y_p(l).ge.dom(ib)%y(j-1)).and.
     &         (y_p(l).lt.dom(ib)%y(j))) then
          do i=is-1,ie
           if ((x_p(l).ge.dom(ib)%x(i-1)).and.
     &         (x_p(l).lt.dom(ib)%x(i))) then
        
           dw=dom(ib)%w(i,j,k)-w_p(l)

           end if
          end do
           end if
         end do
           end if 
        end do 
	else if (nway.eq.2) then

!u_p,v_p,w_p set to 0 here?

	  nxl=1.49999d0		!phi 3 function kernel width
	  nl=0			!number of neighbours
	  dist=100000.d0 
	  u_int(L)=0.d0 ; v_int(L)=0.d0 ; w_int(L)=0.d0
         DO I = is-1,ie 
       IF (dom(ib)%x(i) .gt.(x_p(L)+nxl*dom(ib)%dx) .or.
     &     dom(ib)%x(i) .lt.(x_p(L)-nxl*dom(ib)%dx)) THEN
           DO J = js-1,je 
       IF (dom(ib)%yc(j).gt.(y_p(L)+nxl*dom(ib)%dy) .or.
     &     dom(ib)%yc(j).lt.(y_p(L)-nxl*dom(ib)%dy)) THEN
            DO K = ks-1,ke 
       IF (dom(ib)%zc(k).gt.(z_p(L)+nxl*dom(ib)%dz) .or.
     &     dom(ib)%zc(k).lt.(z_p(L)-nxl*dom(ib)%dz)) THEN 
	    nl=nl+1									!number of neighbour
         dh_loc(L,nl,1)= dh(dom(ib)%dx,dom(ib)%dy,dom(ib)%dz,
     &  dom(ib)%X(I),dom(ib)%YC(J),dom(ib)%ZC(K),x_p(L),y_p(L),z_p(L),5)  !INTERPOLATION function

         u_int(L)=u_int(L)+dom(ib)%u(I,J,K)*dh_loc(L,nl,1)     	!interpolated velocity at particle's location

	   l_dist=SQRT((dom(ib)%x(i)-x_p(L))**2+(dom(ib)%yc(i)-y_p(L))**2
     &        +(dom(ib)%zc(i)-z_p(L))**2)					!Distance particle and fluid cell

	   if(l_dist.lt.dist) then
		I_nr(L,1,1)=I ; J_nr(L,1,1)=J ; K_nr(L,1,1)=K   	!Closest neighbour of particle l in U-grid
		dist=l_dist
	   endif 
		if (nl.ge.27)  GOTO 700				!maximum number of neighbours possible
 	    ENDIF
           END DO
	   ENDIF 
          END DO
	  ENDIF
         END DO        
700	continue
! v
	  nl=0 ;  dist=100000.d0 
         DO I = is-1,ie 
       IF (dom(ib)%xc(i) .gt.(x_p(L)+nxl*dom(ib)%dx) .or.
     &     dom(ib)%xc(i) .lt.(x_p(L)-nxl*dom(ib)%dx)) THEN
           DO J = js-1,je 
       IF (dom(ib)%y(j) .gt.(y_p(L)+nxl*dom(ib)%dy) .or.
     &     dom(ib)%y(j) .lt.(y_p(L)-nxl*dom(ib)%dy)) THEN
            DO K = ks-1,ke 
       IF (dom(ib)%zc(k).gt.(z_p(L)+nxl*dom(ib)%dz) .or.
     &     dom(ib)%zc(k).lt.(z_p(L)-nxl*dom(ib)%dz)) THEN 
	    nl=nl+1
        dh_loc(L,nl,2)= dh(dom(ib)%dx,dom(ib)%dy,dom(ib)%dz,
     &  dom(ib)%XC(I),dom(ib)%Y(J),dom(ib)%ZC(K),x_p(L),y_p(L),z_p(L),5)  				!INTERPOLATION function
         v_int(L)=v_int(L)+dom(ib)%v(I,J,K)*dh_loc(L,nl,2)     
	   l_dist=SQRT((dom(ib)%xc(i)-x_p(L))**2+(dom(ib)%y(i)-y_p(L))**2
     &         +(dom(ib)%zc(i)-z_p(L))**2)				!Distance particle and fluid cell
	   if(l_dist.lt.dist) then
		I_nr(L,1,2)=I ; J_nr(L,1,2)=J ; K_nr(L,1,2)=K   	!Closest neighbour of particle l
		dist=l_dist
	   endif 

		if (nl.ge.27)  GOTO 701
 	    ENDIF
           END DO
	   ENDIF 
          END DO
	  ENDIF
         END DO        
701	continue
! w
	  nl=0 ;  dist=100000.d0 
         DO I = is-1,ie 
       IF (dom(ib)%xc(i) .gt.(x_p(L)+nxl*dom(ib)%dx) .or.
     &     dom(ib)%xc(i) .lt.(x_p(L)-nxl*dom(ib)%dx)) THEN
           DO J = js-1,je 
       IF (dom(ib)%yc(j).gt.(y_p(L)+nxl*dom(ib)%dy) .or.
     &     dom(ib)%yc(j).lt.(y_p(L)-nxl*dom(ib)%dy)) THEN
            DO K = ks-1,ke 
       IF (dom(ib)%z(k) .gt.(z_p(L)+nxl*dom(ib)%dz) .or.
     &     dom(ib)%z(k) .lt.(z_p(L)-nxl*dom(ib)%dz)) THEN 
	    nl=nl+1
        dh_loc(L,nl,3)= dh(dom(ib)%dx,dom(ib)%dy,dom(ib)%dz,
     & dom(ib)%XC(I),dom(ib)%Y(J),dom(ib)%ZC(K) ,x_p(L),y_p(L),z_p(L),5)  				!INTERPOLATION function
         w_int(L)=w_int(L)+dom(ib)%w(I,J,K)*dh_loc(L,nl,3)     
	   l_dist=SQRT((dom(ib)%xc(i)-x_p(L))**2+(dom(ib)%yc(i)-y_p(L))**2
     &         +(dom(ib)%z(i)-z_p(L))**2)					!Distance particle and fluid cell
	   if(l_dist.lt.dist) then
		I_nr(L,1,3)=I ; J_nr(L,1,3)=J ; K_nr(L,1,3)=K   	!Closest neighbour of particle l
		dist=l_dist
	   endif 

		if (nl.ge.27)  GOTO 702
 	    ENDIF
           END DO
	   ENDIF 
          END DO
	  ENDIF
         END DO        
702	continue

!Calculate gradients from the closest cell
!U-grid (x,yc,zc) : I_nr(L,1,1),J_nr(L,1,1),K_nr(L,1,1)
!V-grid (xc,y,zc) : I_nr(L,1,2),J_nr(L,1,2),K_nr(L,1,2)
!W-grid (xc,yc,z) : I_nr(L,1,3),J_nr(L,1,3),K_nr(L,1,3)
!Calculate the slip velocity
          if(u_int(l).lt.0.0) then   !x,y,z axis
           du=0.0-u_p(l)
          else
           du=u_int(l)-u_p(l) 
          end if 
           dv=v_int(l)-v_p(l)
           dw=w_int(l)-w_p(l)

	endif

!particle-frame coordination        
         du1=du*cos(phi_yp(l))-dw*sin(phi_yp(l))   !x1,y1,z1 axis
         dw1=dw*cos(phi_yp(l))+du*sin(phi_yp(l))
        
! calculate drag force
         d_vel=(du**2.0+dv**2.0+dw**2.0)**0.5
         REp=d_vel*((a_p+b_p+c_p)/3.0) 
         ! REp=(du**2.0+dw**2.0)**0.5*((a_p+b_p+c_p)/3.0)
     &   /10.0**(-6.0)*u_ref
         Cd=24.0/REp*(1.0+0.15*REp**0.681)+0.407/(1.0+8710.0/REp)  
         
         Fdx1=d_vel*du1*Cd*pi*bd_p*cd_p/8.0
         Fdy=d_vel*dv*Cd*pi*ad_p*cd_p/8.0
         Fdz1=d_vel*dw1*Cd*pi*ad_p*bd_p/8.0
         
          !  Fdx=abs(du)*du*Cd*pi*bd_p*cd_p/8.0  
          Fdx=Fdx1*cos(phi_yp(l))+Fdz1*sin(phi_yp(l))
          ! Fdz=abs(dw)*dw*Cd*pi*ad_p*bd_p/8.0  
          Fdz=Fdz1*cos(phi_yp(l))-Fdx1*sin(phi_yp(l))
         
! calculate lift force 
          Cl=0.3                     !can change 0.5
!          Fl=((dom(ib)%u(i+1,j+1,k+2)-u_p(l))**2.0		!OLD
!     & -(dom(ib)%u(i+1,j+1,k+1)-u_p(l))**2.0)
!     & *Cl*pi*bd_p/8.0      !*ad_p
!     & *(ad_p*abs(cos(phi_yp(l)))+cd_p*abs(sin(phi_yp(l))))
        do k=ks-1,ke 
           if ((z_p(l).ge.dom(ib)%zc(k)).and.
     &     (z_p(l).lt.dom(ib)%zc(k+1))) then 
         do j=js-1,je 
           if ((y_p(l).ge.dom(ib)%y(j-1)).and.
     &     (y_p(l).lt.dom(ib)%y(j))) then
          do i=is-1,ie
           if ((x_p(l).ge.dom(ib)%xc(i)).and.
     &     (x_p(l).lt.dom(ib)%xc(i+1))) then
!use here : I_nr(L,1,3),J_nr(L,1,3),K_nr(L,1,3)+1 ??
            Fl=((dom(ib)%u(i,j,k+1)-u_p(l))**2.0
     & -(dom(ib)%u(i,j,k)-u_p(l))**2.0)
     & *Cl*pi*bd_p/8.0      !*ad_p
     & *(ad_p*abs(cos(phi_yp(l)))+cd_p*abs(sin(phi_yp(l))))

           end if
          end do
           end if
         end do
           end if 
        end do 

         if (Fl.lt.0.0) Fl=0.0 
         Fm=pi*ad_p*bd_p*cd_p/8.0*d_vel*omga_yp(l)
         if (Fm.lt.0.0) Fm=0.0 
         fz_p(l)=Fdz+Fl-G_p+Fm  

! consider friction when particles lie on bed
         if (z_p(l).le.zd)then
            k1=0.6                         !can change
            ff=k1*fz_p(l)*(-1.0)
            ff=max(0.0,ff)
            Fd=(Fdx**2+Fdy**2)**0.5 
          if (Fd.eq.0)then
            fx_p(l)=0.0
            fy_p(l)=0.0
          else
            if (u_p(l).eq.0)then  
            !fx_p(l)=Fdx-ff
            fx_p(l)=Fdx-ff*abs(Fdx/Fd)
            fx_p(l)=max(0.0,fx_p(l))
            else
            !fx_p(l)=Fdx-ff*u_p(l)/abs(u_p(l))
            fx_p(l)=Fdx-ff*abs(Fdx/Fd)*u_p(l)/abs(u_p(l))
            end if
            
            if (v_p(l).eq.0)then
            fy_p(l)=Fdy-ff*abs(Fdy/Fd)
            fy_p(l)=max(0.0,fy_p(l))
            else
            fy_p(l)=Fdy-ff*abs(Fdy/Fd)*v_p(l)/abs(v_p(l))
           
            end if 
          end if 
            else
            fx_p(l)=Fdx
            fy_p(l)=Fdy 
          end if
          
!calculate rotation parameters of flow 			
!use here "3" (I_nr(L,1,3),J_nr(L,1,3),K_nr(L,1,3))
        do k=ks-1,ke
           if ((z_p(l).ge.dom(ib)%z(k-1)).and.
     &     (z_p(l).lt.dom(ib)%z(k))) then 
         do j=js-1,je 
           if ((y_p(l).ge.dom(ib)%yc(j)).and.
     &     (y_p(l).lt.dom(ib)%yc(j+1))) then
          do i=is-1,ie
           if ((x_p(l).ge.dom(ib)%xc(i)).and.
     &     (x_p(l).lt.dom(ib)%xc(i+1))) then
     
             dudy=(dom(ib)%u(i,j+1,k)-dom(ib)%u(i,j,k))/dom(ib)%dy
             dvdx=(dom(ib)%v(i+1,j,k)-dom(ib)%v(i,j,k))/dom(ib)%dx
           end if
          end do
           end if
         end do
           end if 
        end do 
!use here "2"
       do k=ks-1,ke
           if ((z_p(l).ge.dom(ib)%zc(k)).and.
     &     (z_p(l).lt.dom(ib)%zc(k+1))) then 
         do j=js-1,je
           if ((y_p(l).ge.dom(ib)%y(j-1)).and.
     &     (y_p(l).lt.dom(ib)%y(j))) then 
          do i=is-1,ie
           if ((x_p(l).ge.dom(ib)%xc(i)).and.
     &     (x_p(l).lt.dom(ib)%xc(i+1))) then
                  
             dudz=(dom(ib)%u(i,j,k+1)-dom(ib)%u(i,j,k))/dom(ib)%dz
             dwdx=(dom(ib)%w(i+1,j,k)-dom(ib)%w(i,j,k))/dom(ib)%dx
           end if
          end do
           end if
         end do
           end if 
        end do 
!use here "1"
      do k=ks-1,ke
           if ((z_p(l).ge.dom(ib)%zc(k)).and.
     &     (z_p(l).lt.dom(ib)%zc(k+1))) then
         do j=js-1,je 
           if ((y_p(l).ge.dom(ib)%yc(j)).and.
     &     (y_p(l).lt.dom(ib)%yc(j+1))) then
          do i=is-1,ie
           if ((x_p(l).ge.dom(ib)%x(i-1)).and.
     &     (x_p(l).lt.dom(ib)%x(i))) then
      
              dvdz=(dom(ib)%v(i,j,k+1)-dom(ib)%v(i,j,k))/dom(ib)%dz
              dwdy=(dom(ib)%w(i,j+1,k)-dom(ib)%w(i,j,k))/dom(ib)%dy
            end if
          end do
           end if
         end do
           end if 
        end do
         
	 
       omga_y=0.5*(dudz-dwdx)
       !omga_z=0.5*(dudy-dvdx)
       !omga_x=0.5*(dvdz-dwdy)
       !S_xy=0.5*(dudy+dvdx)
       S_xz=0.5*(dudz+dwdx)
       !S_yz=0.5*(dvdz+dwdy)


! N-S Equation to calculate particles translational movement
       ! velocity     
        u_p(l)=u_p(l)+fx_p(l)*dt/(m_p+m_p/2.65*0.5)  !0.5 added mass
        u_p(l)=max(u_p(l),0.0)
        v_p(l)=v_p(l)+fy_p(l)*dt/(m_p+m_p/2.65*0.5)
        w_p(l)=w_p(l)+fz_p(l)*dt/(m_p+m_p/2.65*0.5)
     & +dt*(2.0*kk*(dt**0.5)*w_p(l)/dt)/(1.0+2.0*kk*(dt**0.5))  !Basset force 
         
	if(u_p(l).gt.umax)write(6,*)'Particle goes too fast'


        FB=(2.0*kk*(dt**0.5)*w_p(l)/dt)/(1.0+2.0*kk*(dt**0.5))
     & *(m_p+m_p/2.65*0.5) 
        ! location    
        a1=x_p(l)+u_p(l)*dt
        if (a1.ge.lx) then
        x_p(l)=a1-lx
      !   else if((a1.lt.dom(ib)%x(is-1)).and.
      !& (nxdom.eq.0)) then
      !   x_p(l)=a1+lx
        else
        x_p(l)=a1
        end if
        
            a2=y_p(l)+v_p(l)*dt 
            if (a2.ge.ly) then
            y_p(l)=ly-0.0001
            v_p(l)=0.0
            else if(a2.le.0.0) then
            y_p(l)=0.0001
            v_p(l)=0.0
             else
            y_p(l)=a2 
             end if
         
        a3=z_p(l)+w_p(l)*dt  
        if(a3.gt.zd) then
         if(a3.lt.lz) then 
         z_p(l)=a3
         else 
         z_p(l)=lz-0.0001
         w_p(l)=0.0
         end if
        else 
          z_p(l)=zd
          w_p(l)=0.0
        end if
        
! calculate particles rotational movement        
       Rew=abs(omga_y-omga_yp(l))	!NOW WE HAVE "abs" BUT NOT IN THE OLD CODE
     &  *((ad_p+bd_p+cd_p)/6.0)**2.0*Re     
      if(Rew.ge.0.0.and.Rew.le.10.0)then
       Cw=16.0*pi/(Rew+1.d-5)
      elseif(Rew.gt.10.0.and.Rew.le.1000.0) then
       Cw=16.0*pi/Rew*(1.0+Rew**2.0/288184.0)
      elseif(Rew.gt.1000.0.and.Rew.le.4000.0) then
       Cw=16.0*pi/Rew/3.0+6.64/Rew**0.5
      elseif(Rew.gt.4000.0) then
       Cw=0.10	!Check this number from XXXX paper
      else
        write(6,*) 'Rew is out of range',Rew,Cw
      end if

        Ix=1.0/5.0*m_p*(bd_p**2.0+cd_p**2.0)
        Iy=1.0/5.0*m_p*(ad_p**2.0+cd_p**2.0)
        Iz=1.0/5.0*m_p*(ad_p**2.0+bd_p**2.0)
       omga_yp(l)=omga_yp(l)+dt/Iy*Cw/2.0*
     & ((ad_p+bd_p+cd_p)/6.0)**5.0*(omga_y-omga_yp(l))
     & *abs(omga_y-omga_yp(l))  
       
        !omga_xp(l)=omga_xp(l)+dt*15.0/Re/bd_p/cd_p/2.65
    ! & *((b_a**2-c_a**2)/(b_a**2+c_a**2)*S_yz+(omga_x-omga_xp(l)))
!        omga_yp(l)=omga_yp(l)+dt*15.0/Re/ad_p/cd_p/2.65
!     & *((c_a**2-1.0)/(1.0+c_a**2)*S_xz+(omga_y-omga_yp(l)))	!IN THE OLD CODE

!         if (omga_yp(l).lt.0.0)then
!          omga_yp(l)=0.0					!VALID IN THE OLD CODE
!         end if

        !omga_zp(l)=omga_zp(l)+dt*15.0/Re/ad_p/bd_p/2.65
    ! & *((1.0-b_a**2)/(1.0+b_a**2)*S_yz+(omga_z-omga_zp(l)))
         
         !phi_xp(l)=phi_xp(l)+omga_xp(l)*dt
         phi_yp(l)=phi_yp(l)+omga_yp(l)*dt
         !phi_zp(l)=phi_zp(l)+omga_zp(l)*dt

!----------preparing to write forces in particles 176
        if (l.eq.176) then  
	   Fdx_p(l)=Fdx ; Fdy_p(l)=Fdy ; Fdz_p(l)=Fdz
         Fl_p(l)=Fl   ; Fm_p(l)=Fm   ; ff_p(l)=ff    ; FB_p(l)=FB
        end if

       end if	
       end do	!n_p
       end do	!ib
     
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)

!send calculated particles' new infomation to all other domains       
       do l=1,n_p
	  call MPI_Bcast(x_p(l),1,MPI_FLT,locdom_p(l),MPI_COMM_WORLD,ierr)
        call MPI_Bcast(y_p(l),1,MPI_FLT,locdom_p(l),MPI_COMM_WORLD,ierr)
        call MPI_Bcast(z_p(l),1,MPI_FLT,locdom_p(l),MPI_COMM_WORLD,ierr)
        call MPI_Bcast(u_p(l),1,MPI_FLT,locdom_p(l),MPI_COMM_WORLD,ierr)
        call MPI_Bcast(v_p(l),1,MPI_FLT,locdom_p(l),MPI_COMM_WORLD,ierr)
        call MPI_Bcast(w_p(l),1,MPI_FLT,locdom_p(l),MPI_COMM_WORLD,ierr)
       call MPI_Bcast(phi_yp(l),1,MPI_FLT,
     &  locdom_p(l),MPI_COMM_WORLD,ierr)
        !call MPI_Bcast(phi_zp(l),1,MPI_FLT,locdom_p(l),MPI_COMM_WORLD,ierr)
        !call MPI_Bcast(omga_xp(l),1,MPI_FLT,locdom_p(l),MPI_COMM_WORLD,ierr)
       call MPI_Bcast(omga_yp(l),1,MPI_FLT,
     &  locdom_p(l),MPI_COMM_WORLD,ierr)
       ! call MPI_Bcast(omga_zp(l),1,MPI_FLT,locdom_p(l),MPI_COMM_WORLD,ierr)
        if (l.eq.176) then
      call MPI_Bcast(Fdx_p(l),1,MPI_FLT,locdom_p(l),MPI_COMM_WORLD,ierr)
      call MPI_Bcast(Fdy_p(l),1,MPI_FLT,locdom_p(l),MPI_COMM_WORLD,ierr)
      call MPI_Bcast(Fdz_p(l),1,MPI_FLT,locdom_p(l),MPI_COMM_WORLD,ierr)
	 call MPI_Bcast(fx_p(l),1,MPI_FLT,locdom_p(l),MPI_COMM_WORLD,ierr)
       call MPI_Bcast(fy_p(l),1,MPI_FLT,locdom_p(l),MPI_COMM_WORLD,ierr)
       call MPI_Bcast(fz_p(l),1,MPI_FLT,locdom_p(l),MPI_COMM_WORLD,ierr)
	 call MPI_Bcast(Fl_p(l),1,MPI_FLT,locdom_p(l),MPI_COMM_WORLD,ierr)
       call MPI_Bcast(Fm_p(l),1,MPI_FLT,locdom_p(l),MPI_COMM_WORLD,ierr)
       call MPI_Bcast(ff_p(l),1,MPI_FLT,locdom_p(l),MPI_COMM_WORLD,ierr)
       call MPI_Bcast(FB_p(l),1,MPI_FLT,locdom_p(l),MPI_COMM_WORLD,ierr)
        end if 
       end do

! write particles every step of move in output files        
       do ib=1, nbp
        do l=1,n_p
        n_pnum=l*500+38
        tt=(itime)*dtavg    !-30000
       ! if (mod(itime,10).eq.0) then 
        if (myrank.eq.0) then
!-----------------write forces in particles 176
        if (l.eq.176) then
        write(5555,88) itime*1.0,x_p(l),y_p(l),z_p(l),
     & u_p(l),v_p(l),w_p(l),phi_yp(l),omga_yp(l),
     & Fdx_p(l),Fdy_p(l),Fdz_p(l),Fl_p(l),Fm_p(l),
     & G_p,ff_p(l),FB_p(l),fx_p(l),fy_p(l),fz_p(l)               !!!!!!!
        end if !l=176

        !write(n_pnum,88) tt,x_p(l),y_p(l),z_p(l),
    ! & u_p(l),v_p(l),w_p(l),phi_xp(l),phi_yp(l),phi_zp(l),
    ! & omga_xp(l),omga_yp(l),omga_zp(l)
        write(n_pnum,88) tt,x_p(l),y_p(l),z_p(l),u_p(l),v_p(l),w_p(l)
     & ,phi_yp(l),omga_yp(l)
      
        end if
       ! end if
       end do
      end do

88      format (15e19.8)
	end
!##########################################################################
        subroutine particle_feedback
!##########################################################################
        use module_parmove
        use vars
        use mpi
        use multidata
        implicit none
        integer :: i,j,k,ib,l,nl,nway
        double precision :: visc,diff,fre,pi,m_p,dist,l_dist,nxl,dh

         pi=4.D0*DATAN(1.D0)
         m_p=1.0/6.0*2.65*pi*a_p*b_p*c_p*1000.0   !particles mass
	   nxl=1.49999d0		!phi 3 function kernel width
        do ib=1,nbp

!U-grid (x,yc,zc) : I_nr(L,1,1),J_nr(L,1,1),K_nr(L,1,1)
!V-grid (xc,y,zc) : I_nr(L,1,2),J_nr(L,1,2),K_nr(L,1,2)
!W-grid (xc,yc,z) : I_nr(L,1,3),J_nr(L,1,3),K_nr(L,1,3)
!------loop in particles        
         do l=1,n_p    
	nway=2
	if (nway.eq.1) then
!-------feedback in velocity u 
         do k=dom(ib)%ksu-1,dom(ib)%keu
           if ((z_p(l).ge.dom(ib)%z(k-1)).and.
     &     (z_p(l).lt.dom(ib)%z(k))) then
         do j=dom(ib)%jsu-1,dom(ib)%jeu 
           if ((y_p(l).ge.dom(ib)%y(j-1)).and.
     &     (y_p(l).lt.dom(ib)%y(j))) then
         do i=dom(ib)%isu-1,dom(ib)%ieu
           if((x_p(l).ge.dom(ib)%xc(i)).and.
     &     (x_p(l).lt.dom(ib)%xc(i+1))) then

          fre=-m_p*(u_p(l)-u_p0(l))/dt/dom(ib)%dx/dom(ib)%dy/dom(ib)%dz
	!IF(l.eq.20) write(6,*)'before',dom(ib)%ustar(i,j,k)
          dom(ib)%ustar(i,j,k)=dom(ib)%ustar(i,j,k)+dt*fre
	!IF(l.eq.20) write(6,*)'after ',dom(ib)%ustar(i,j,k)
          end if
         end do
          end if
         end do
          end if
         end do
!-------feedback in velocity v
         do k=dom(ib)%ksv-1,dom(ib)%kev
           if ((z_p(l).ge.dom(ib)%z(k-1)).and.
     &     (z_p(l).lt.dom(ib)%z(k))) then
         do j=dom(ib)%jsv-1,dom(ib)%jev
           if ((y_p(l).ge.dom(ib)%yc(j)).and.
     &     (y_p(l).lt.dom(ib)%yc(j+1))) then
         do i=dom(ib)%isv-1,dom(ib)%iev
           if ((x_p(l).ge.dom(ib)%x(i-1)).and.
     &     (x_p(l).lt.dom(ib)%x(i))) then       
           fre=-m_p*(v_p(l)-v_p0(l))/dt/dom(ib)%dx/dom(ib)%dy/dom(ib)%dz
           dom(ib)%vstar(i,j,k)=dom(ib)%vstar(i,j,k)+dt*fre
          end if
         end do
          end if
         end do
          end if
         end do
!-------feedback in velocity w
         do k=dom(ib)%ksw-1,dom(ib)%kew
           if((z_p(l).ge.dom(ib)%zc(k)).and.
     &     (z_p(l).lt.dom(ib)%zc(k+1))) then
         do j=dom(ib)%jsw-1,dom(ib)%jew
           if ((y_p(l).ge.dom(ib)%y(j-1)).and.
     &     (y_p(l).lt.dom(ib)%y(j))) then
         do i=dom(ib)%isw-1,dom(ib)%iew
           if ((x_p(l).ge.dom(ib)%x(i-1)).and.
     &     (x_p(l).lt.dom(ib)%x(i))) then
           fre=-m_p*(w_p(l)-w_p0(l))/dt/dom(ib)%dx/dom(ib)%dy/dom(ib)%dz
           dom(ib)%wstar(i,j,k)=dom(ib)%wstar(i,j,k)+dt*fre
          end if
         end do
          end if
         end do
          end if
         end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	else if(nway.eq.2) then
	  nl=0			!number of neighbours
        fre=-m_p*(u_p(l)-u_p0(l))/(dt*dom(ib)%dx*dom(ib)%dy*dom(ib)%dz)
        DO K = dom(ib)%ksu-1,dom(ib)%keu
       IF (dom(ib)%z(k).gt.(z_p(L)+nxl*dom(ib)%dz) .or.
     &     dom(ib)%z(k).lt.(z_p(L)-nxl*dom(ib)%dz)) THEN 
         DO J = dom(ib)%jsu-1,dom(ib)%jeu 
       IF (dom(ib)%y(j).gt.(y_p(L)+nxl*dom(ib)%dy) .or.
     &     dom(ib)%y(j).lt.(y_p(L)-nxl*dom(ib)%dy)) THEN
         DO I = dom(ib)%isu-1,dom(ib)%ieu
       IF (dom(ib)%xc(i).gt.(x_p(L)+nxl*dom(ib)%dx) .or.
     &     dom(ib)%xc(i).lt.(x_p(L)-nxl*dom(ib)%dx)) THEN
	    nl=nl+1									!number of neighbour
        dh_loc(L,nl,1)= dh(dom(ib)%dx,dom(ib)%dy,dom(ib)%dz,
     &  dom(ib)%XC(I),dom(ib)%Y(J),dom(ib)%Z(K),x_p(L),y_p(L),z_p(L),5)  !INTERPOLATION function
        dom(ib)%ustar(i,j,k)=dom(ib)%ustar(i,j,k)+dt*fre*dh_loc(L,nl,1)
		if (nl.ge.27)  GOTO 700				!maximum number of neighbours possible
 	    ENDIF
           END DO
	   ENDIF 
          END DO
	  ENDIF
         END DO        
700	continue
	  nl=0			!number of neighbours
        fre=-m_p*(v_p(l)-v_p0(l))/(dt*dom(ib)%dx*dom(ib)%dy*dom(ib)%dz)
        DO K = dom(ib)%ksv-1,dom(ib)%kev
       IF (dom(ib)%z(k).gt.(z_p(L)+nxl*dom(ib)%dz) .or.
     &     dom(ib)%z(k).lt.(z_p(L)-nxl*dom(ib)%dz)) THEN 
         DO J = dom(ib)%jsv-1,dom(ib)%jev
       IF (dom(ib)%yc(j).gt.(y_p(L)+nxl*dom(ib)%dy) .or.
     &     dom(ib)%yc(j).lt.(y_p(L)-nxl*dom(ib)%dy)) THEN
         DO I = dom(ib)%isv-1,dom(ib)%iev
       IF (dom(ib)%x(i).gt.(x_p(L)+nxl*dom(ib)%dx) .or.
     &     dom(ib)%x(i).lt.(x_p(L)-nxl*dom(ib)%dx)) THEN
	    nl=nl+1									!number of neighbour
        dh_loc(L,nl,2)= dh(dom(ib)%dx,dom(ib)%dy,dom(ib)%dz,
     &  dom(ib)%X(I),dom(ib)%YC(J),dom(ib)%Z(K),x_p(L),y_p(L),z_p(L),5)  !INTERPOLATION function
        dom(ib)%vstar(i,j,k)=dom(ib)%vstar(i,j,k)+dt*fre*dh_loc(L,nl,2)
		if (nl.ge.27)  GOTO 701				!maximum number of neighbours possible
 	    ENDIF
           END DO
	   ENDIF 
          END DO
	  ENDIF
         END DO        
701	continue
	  nl=0			!number of neighbours
        fre=-m_p*(w_p(l)-w_p0(l))/(dt*dom(ib)%dx*dom(ib)%dy*dom(ib)%dz)
        DO K = dom(ib)%ksw-1,dom(ib)%kew
       IF (dom(ib)%zc(k).gt.(z_p(L)+nxl*dom(ib)%dz) .or.
     &     dom(ib)%zc(k).lt.(z_p(L)-nxl*dom(ib)%dz)) THEN 
         DO J = dom(ib)%jsw-1,dom(ib)%jew
       IF (dom(ib)%y(j).gt.(y_p(L)+nxl*dom(ib)%dy) .or.
     &     dom(ib)%y(j).lt.(y_p(L)-nxl*dom(ib)%dy)) THEN
         DO I = dom(ib)%isw-1,dom(ib)%iew
       IF (dom(ib)%x(i).gt.(x_p(L)+nxl*dom(ib)%dx) .or.
     &     dom(ib)%x(i).lt.(x_p(L)-nxl*dom(ib)%dx)) THEN
	    nl=nl+1									!number of neighbour
        dh_loc(L,nl,3)= dh(dom(ib)%dx,dom(ib)%dy,dom(ib)%dz,
     &  dom(ib)%X(I),dom(ib)%Y(J),dom(ib)%ZC(K),x_p(L),y_p(L),z_p(L),5)  !INTERPOLATION function
        dom(ib)%wstar(i,j,k)=dom(ib)%wstar(i,j,k)+dt*fre*dh_loc(L,nl,3)
		if (nl.ge.27)  GOTO 702				!maximum number of neighbours possible
 	    ENDIF
           END DO
	   ENDIF 
          END DO
	  ENDIF
         END DO        
702	continue

	endif

        end do   !l=1,n_p 
        end do !ib=1,nbp
        return
        end subroutine particle_feedback
