!##########################################################################
        subroutine eddyv_wale
!##########################################################################
        use vars
        use multidata
        implicit none
        integer :: i,j,k
        integer :: ib,is,ie,js,je,ks,ke
        double precision :: dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz
        double precision :: vr_a,vr_b,vr_c,vr_d
        double precision :: s12,s13,s23
        double precision :: s11d,s22d,s33d,s12d,s13d,s23d
        double precision :: ss,sdsd,eqnA,denom
        double precision :: rrey,h1,h2,h3,rh123
        double precision :: ufv_c,ufv_n1,ufv_n2
        double precision :: vfv_c,vfv_n1,vfv_n2
        double precision :: wfv_c,wfv_n1,wfv_n2
        double precision :: cw,delta_grid,l_s
        double precision :: dx,dy,dz
        integer :: sn
        character*8 :: chb1
        character*25 :: gf


        rrey=1.0/Re
        cw=0.46
    
        do ib=1,nbp

        delta_grid=(dom(ib)%dx*dom(ib)%dy*dom(ib)%dz)**(1.0/3.0)

        is=dom(ib)%isp; ie=dom(ib)%iep
        js=dom(ib)%jsp; je=dom(ib)%jep
        ks=dom(ib)%ksp; ke=dom(ib)%kep

        do k=ks,ke
           do j=js,je
              do i=is,ie

!====================================================
        vr_a = 0.25*( dom(ib)%u(i,j,k)   + dom(ib)%u(i-1,j,k) +
     &	 	      dom(ib)%u(i,j+1,k) + dom(ib)%u(i-1,j+1,k) )
        vr_b = 0.25*( dom(ib)%u(i,j,k)   + dom(ib)%u(i-1,j,k) +
     &	 	      dom(ib)%u(i,j-1,k) + dom(ib)%u(i-1,j-1,k) )

        vr_c = 0.25*( dom(ib)%u(i,j,k)   + dom(ib)%u(i-1,j,k) +
     &	 	      dom(ib)%u(i,j,k+1) + dom(ib)%u(i-1,j,k+1) )
        vr_d = 0.25*( dom(ib)%u(i,j,k)   + dom(ib)%u(i-1,j,k) +
     &	 	      dom(ib)%u(i,j,k-1) + dom(ib)%u(i-1,j,k-1) )

        dudx = ( dom(ib)%u(i,j,k) - dom(ib)%u(i-1,j,k) )/dom(ib)%dx
        dudy = ( vr_a - vr_b )/dom(ib)%dy
        dudz = ( vr_c - vr_d )/dom(ib)%dz
!====================================================
        vr_a = 0.25*( dom(ib)%v(i,j,k)   + dom(ib)%v(i,j-1,k) +
     &	 	      dom(ib)%v(i+1,j,k) + dom(ib)%v(i+1,j-1,k) )
        vr_b = 0.25*( dom(ib)%v(i,j,k)   + dom(ib)%v(i,j-1,k) +
     &	 	      dom(ib)%v(i-1,j,k) + dom(ib)%v(i-1,j-1,k) )

        vr_c = 0.25*( dom(ib)%v(i,j,k)   + dom(ib)%v(i,j-1,k) +
     &	 	      dom(ib)%v(i,j,k+1) + dom(ib)%v(i,j-1,k+1) )
        vr_d = 0.25*( dom(ib)%v(i,j,k)   + dom(ib)%v(i,j-1,k) +
     &	 	      dom(ib)%v(i,j,k-1) + dom(ib)%v(i,j-1,k-1) )

        dvdy = ( dom(ib)%v(i,j,k) - dom(ib)%v(i,j-1,k) )/dom(ib)%dy
        dvdx = ( vr_a - vr_b )/dom(ib)%dx
        dvdz = ( vr_c - vr_d )/dom(ib)%dz
!====================================================
        vr_a = 0.25*( dom(ib)%w(i,j,k)   + dom(ib)%w(i,j,k-1) +
     &	 	      dom(ib)%w(i+1,j,k) + dom(ib)%w(i+1,j,k-1) )
        vr_b = 0.25*( dom(ib)%w(i,j,k)   + dom(ib)%w(i,j,k-1) +
     &	 	      dom(ib)%w(i-1,j,k) + dom(ib)%w(i-1,j,k-1) )

        vr_c = 0.25*( dom(ib)%w(i,j,k)   + dom(ib)%w(i,j,k-1) +
     &	 	      dom(ib)%w(i,j+1,k) + dom(ib)%w(i,j+1,k-1) )
        vr_d = 0.25*( dom(ib)%w(i,j,k)   + dom(ib)%w(i,j,k-1) +
     &	 	      dom(ib)%w(i,j-1,k) + dom(ib)%w(i,j-1,k-1) )

        dwdz = ( dom(ib)%w(i,j,k) - dom(ib)%w(i,j,k-1) )/dom(ib)%dz
        dwdx = ( vr_a - vr_b )/dom(ib)%dx
        dwdy = ( vr_c - vr_d )/dom(ib)%dy


!==========================================================================
! ..... WALL BOUNDARIES
!==========================================================================
        if (i.eq.dom(ib)%isp) then         
           if (dom(ib)%iprev.lt.0) then
		if (dom(ib)%bc_west.ge.61 .or. dom(ib)%bc_west.eq.4) then

              h1=dom(ib)%dx; h2=2.0*dom(ib)%dx; h3=dom(ib)%dx
              rh123=1.0/(h1*h2*h3)
           
              ufv_c =0.5*( dom(ib)%u(i,j,k)  +dom(ib)%u(i-1,j,k  ) )
              ufv_n1=0.5*( dom(ib)%u(i+1,j,k)+dom(ib)%u(i,j,k) )
              ufv_n2=0.5*( dom(ib)%u(i+2,j,k)+dom(ib)%u(i+1,j,k) )

              vfv_c =0.5*( dom(ib)%v(i,j,k)  +dom(ib)%v(i,j-1,k) )
              vfv_n1=0.5*( dom(ib)%v(i+1,j,k)+dom(ib)%v(i+1,j-1,k) )
              vfv_n2=0.5*( dom(ib)%v(i+2,j,k)+dom(ib)%v(i+2,j-1,k) )

              wfv_c =0.5*( dom(ib)%w(i,j,k)  +dom(ib)%w(i,j,k-1) )
              wfv_n1=0.5*( dom(ib)%w(i+1,j,k)+dom(ib)%w(i+1,j,k-1) )
              wfv_n2=0.5*( dom(ib)%w(i+2,j,k)+dom(ib)%w(i+2,j,k-1) )

              dudx=(h2*h2*(ufv_n1-ufv_c)+h1*h1*(ufv_c-ufv_n2))*rh123
              dvdx=(h2*h2*(vfv_n1-vfv_c)+h1*h1*(vfv_c-vfv_n2))*rh123
              dwdx=(h2*h2*(wfv_n1-wfv_c)+h1*h1*(wfv_c-wfv_n2))*rh123
		endif
           end if
        end if

        if (i.eq.dom(ib)%iep) then         
           if (dom(ib)%inext.lt.0) then
		if (dom(ib)%bc_east.ge.61 .or. dom(ib)%bc_east.eq.4) then

              h1=dom(ib)%dx; h2=2.0*dom(ib)%dx; h3=dom(ib)%dx
              rh123=-1.0/(h1*h2*h3)
           
              ufv_c =0.5*( dom(ib)%u(i,j,k)  +dom(ib)%u(i-1,j,k  ) )
              ufv_n1=0.5*( dom(ib)%u(i-1,j,k)+dom(ib)%u(i-2,j,k) )
              ufv_n2=0.5*( dom(ib)%u(i-2,j,k)+dom(ib)%u(i-3,j,k) )

              vfv_c =0.5*( dom(ib)%v(i,j,k)  +dom(ib)%v(i,j-1,k) )
              vfv_n1=0.5*( dom(ib)%v(i-1,j,k)+dom(ib)%v(i-1,j-1,k) )
              vfv_n2=0.5*( dom(ib)%v(i-2,j,k)+dom(ib)%v(i-2,j-1,k) )

              wfv_c =0.5*( dom(ib)%w(i,j,k)  +dom(ib)%w(i,j,k-1) )
              wfv_n1=0.5*( dom(ib)%w(i-1,j,k)+dom(ib)%w(i-1,j,k-1) )
              wfv_n2=0.5*( dom(ib)%w(i-2,j,k)+dom(ib)%w(i-2,j,k-1) )

              dudx=(h2*h2*(ufv_n1-ufv_c)+h1*h1*(ufv_c-ufv_n2))*rh123
              dvdx=(h2*h2*(vfv_n1-vfv_c)+h1*h1*(vfv_c-vfv_n2))*rh123
              dwdx=(h2*h2*(wfv_n1-wfv_c)+h1*h1*(wfv_c-wfv_n2))*rh123
		endif
           end if
        end if

        if (j.eq.dom(ib)%jsp) then         
           if (dom(ib)%jprev.lt.0) then
		if (dom(ib)%bc_south.ge.61 .or. dom(ib)%bc_south.eq.4) then

              h1=dom(ib)%dy; h2=2.0*dom(ib)%dy; h3=dom(ib)%dy
              rh123=1.0/(h1*h2*h3)
           
              ufv_c =0.5*( dom(ib)%u(i,j,k)  +dom(ib)%u(i-1,j,k  ) )
              ufv_n1=0.5*( dom(ib)%u(i,j+1,k)+dom(ib)%u(i-1,j+1,k) )
              ufv_n2=0.5*( dom(ib)%u(i,j+2,k)+dom(ib)%u(i-1,j+2,k) )

              vfv_c =0.5*( dom(ib)%v(i,j,k)  +dom(ib)%v(i,j-1,k  ) )
              vfv_n1=0.5*( dom(ib)%v(i,j+1,k)+dom(ib)%v(i,j,k) )
              vfv_n2=0.5*( dom(ib)%v(i,j+2,k)+dom(ib)%v(i,j+1,k) )

              wfv_c =0.5*( dom(ib)%w(i,j,k)  +dom(ib)%w(i,j,k-1) )
              wfv_n1=0.5*( dom(ib)%w(i,j+1,k)+dom(ib)%w(i,j+1,k-1) )
              wfv_n2=0.5*( dom(ib)%w(i,j+2,k)+dom(ib)%w(i,j+2,k-1) )

              dudy=(h2*h2*(ufv_n1-ufv_c)+h1*h1*(ufv_c-ufv_n2))*rh123
              dvdy=(h2*h2*(vfv_n1-vfv_c)+h1*h1*(vfv_c-vfv_n2))*rh123
              dwdy=(h2*h2*(wfv_n1-wfv_c)+h1*h1*(wfv_c-wfv_n2))*rh123
		endif
           end if
        end if

        if (j.eq.dom(ib)%jep) then         
           if (dom(ib)%jnext.lt.0) then
		if (dom(ib)%bc_north.ge.61 .or. dom(ib)%bc_north.eq.4) then

              h1=dom(ib)%dy; h2=2.0*dom(ib)%dy; h3=dom(ib)%dy
              rh123=-1.0/(h1*h2*h3)
           
              ufv_c =0.5*( dom(ib)%u(i,j,k)  +dom(ib)%u(i-1,j,k  ) )
              ufv_n1=0.5*( dom(ib)%u(i,j-1,k)+dom(ib)%u(i-1,j-1,k) )
              ufv_n2=0.5*( dom(ib)%u(i,j-2,k)+dom(ib)%u(i-1,j-2,k) )

              vfv_c =0.5*( dom(ib)%v(i,j,k)  +dom(ib)%v(i,j-1,k  ) )
              vfv_n1=0.5*( dom(ib)%v(i,j-1,k)+dom(ib)%v(i,j-2,k) )
              vfv_n2=0.5*( dom(ib)%v(i,j-2,k)+dom(ib)%v(i,j-3,k) )

              wfv_c =0.5*( dom(ib)%w(i,j,k)  +dom(ib)%w(i,j,k-1) )
              wfv_n1=0.5*( dom(ib)%w(i,j-1,k)+dom(ib)%w(i,j-1,k-1) )
              wfv_n2=0.5*( dom(ib)%w(i,j-2,k)+dom(ib)%w(i,j-2,k-1) )

              dudy=(h2*h2*(ufv_n1-ufv_c)+h1*h1*(ufv_c-ufv_n2))*rh123
              dvdy=(h2*h2*(vfv_n1-vfv_c)+h1*h1*(vfv_c-vfv_n2))*rh123
              dwdy=(h2*h2*(wfv_n1-wfv_c)+h1*h1*(wfv_c-wfv_n2))*rh123
		endif
           end if
        end if

        if (k.eq.dom(ib)%ksp) then         
           if (dom(ib)%kprev.lt.0) then 
		if (dom(ib)%bc_bottom.ge.61 .or. dom(ib)%bc_bottom.eq.4) then

              h1=dom(ib)%dz; h2=2.0*dom(ib)%dz; h3=dom(ib)%dz
              rh123=1.0/(h1*h2*h3)
           
              ufv_c =0.5*( dom(ib)%u(i,j,k)  +dom(ib)%u(i-1,j,k) )
              ufv_n1=0.5*( dom(ib)%u(i,j,k+1)+dom(ib)%u(i-1,j,k+1) )
              ufv_n2=0.5*( dom(ib)%u(i,j,k+2)+dom(ib)%u(i-1,j,k+2) )

              vfv_c =0.5*( dom(ib)%v(i,j,k)  +dom(ib)%v(i,j-1,k) )
              vfv_n1=0.5*( dom(ib)%v(i,j,k+1)+dom(ib)%v(i,j-1,k+1) )
              vfv_n2=0.5*( dom(ib)%v(i,j,k+2)+dom(ib)%v(i,j-1,k+2) )

              wfv_c =0.5*( dom(ib)%w(i,j,k)  +dom(ib)%w(i,j,k-1) )
              wfv_n1=0.5*( dom(ib)%w(i,j,k+1)+dom(ib)%w(i,j,k) )
              wfv_n2=0.5*( dom(ib)%w(i,j,k+2)+dom(ib)%w(i,j,k+1) )

              dudz=(h2*h2*(ufv_n1-ufv_c)+h1*h1*(ufv_c-ufv_n2))*rh123
              dvdz=(h2*h2*(vfv_n1-vfv_c)+h1*h1*(vfv_c-vfv_n2))*rh123
              dwdz=(h2*h2*(wfv_n1-wfv_c)+h1*h1*(wfv_c-wfv_n2))*rh123
		endif
           end if
        end if

        if (k.eq.dom(ib)%kep) then         
           if (dom(ib)%knext.lt.0) then
		if (dom(ib)%bc_top.ge.61 .or. dom(ib)%bc_top.eq.4) then

              h1=dom(ib)%dz; h2=2.0*dom(ib)%dz; h3=dom(ib)%dz
              rh123=-1.0/(h1*h2*h3)
           
              ufv_c =0.5*( dom(ib)%u(i,j,k)  +dom(ib)%u(i-1,j,k) )
              ufv_n1=0.5*( dom(ib)%u(i,j,k-1)+dom(ib)%u(i-1,j,k-1) )
              ufv_n2=0.5*( dom(ib)%u(i,j,k-2)+dom(ib)%u(i-1,j,k-2) )

              vfv_c =0.5*( dom(ib)%v(i,j,k)  +dom(ib)%v(i,j-1,k) )
              vfv_n1=0.5*( dom(ib)%v(i,j,k-1)+dom(ib)%v(i,j-1,k-1) )
              vfv_n2=0.5*( dom(ib)%v(i,j,k-2)+dom(ib)%v(i,j-1,k-2) )

              wfv_c =0.5*( dom(ib)%w(i,j,k)  +dom(ib)%w(i,j,k-1) )
              wfv_n1=0.5*( dom(ib)%w(i,j,k-1)+dom(ib)%w(i,j,k-2) )
              wfv_n2=0.5*( dom(ib)%w(i,j,k-2)+dom(ib)%w(i,j,k-3) )

              dudz=(h2*h2*(ufv_n1-ufv_c)+h1*h1*(ufv_c-ufv_n2))*rh123
              dvdz=(h2*h2*(vfv_n1-vfv_c)+h1*h1*(vfv_c-vfv_n2))*rh123
              dwdz=(h2*h2*(wfv_n1-wfv_c)+h1*h1*(wfv_c-wfv_n2))*rh123
		endif
           end if
        end if

!==========================================================================
! ..... EDDY VISCOSITY CALCULATION
!==========================================================================
        l_s  = (cw * delta_grid)**2.0

        s12 = 0.5 * (dudy + dvdx)

        s13 = 0.5 * (dudz + dwdx)

        s23 = 0.5 * (dvdz + dwdy)

        ss  = ( dudx*dudx + dvdy*dvdy   + dwdz*dwdz  +
     &      2.0*s12*s12   + 2.0*s13*s13 + 2.0*s23*s23 )


        eqnA=dudx*dudx+dvdy*dvdy+dwdz*dwdz+
     & 2.0*dudy*dvdx+2.0*dudz*dwdx+2.0*dvdz*dwdy

        s11d = dudx*dudx+dudy*dvdx+dudz*dwdx-eqnA/3.0

        s22d = dvdx*dudy+dvdy*dvdy+dvdz*dwdy-eqnA/3.0

        s33d = dwdx*dudz+dwdy*dvdz+dwdz*dwdz-eqnA/3.0

        s12d =0.5*( dudx*dudy+dudy*dvdy+dudz*dwdy+
     &              dvdx*dudx+dvdy*dvdx+dvdz*dwdx)

        s13d =0.5*( dudx*dudz+dudy*dvdz+dudz*dwdz+
     &              dwdx*dudx+dwdy*dvdx+dwdz*dwdx)

        s23d =0.5*( dvdx*dudz+dvdy*dvdz+dvdz*dwdz+
     &              dwdx*dudy+dwdy*dvdy+dwdz*dwdy)

        sdsd = ( s11d*s11d +     s22d*s22d +     s33d*s33d  +
     &       2.0*s12d*s12d + 2.0*s13d*s13d + 2.0*s23d*s23d  )

        if (L_LSM) rrey=dom(ib)%mu(i,j,k)/dom(ib)%dens(i,j,k)

        if(sdsd.ne.0.0) then
           denom=(ss**2.5+sdsd**1.25)
           if (denom.eq.0.0) then
			!print*,'error, denominator is zero!'
           		dom(ib)%vis(i,j,k) = rrey
	     endif
           dom(ib)%vis(i,j,k) = rrey + l_s * sdsd**1.5 / denom
        else
           dom(ib)%vis(i,j,k) = rrey 
        end if

              end do
           end do
        end do

        end do

!        call exchange(7)

        do ib=1,nbp

        is=dom(ib)%isp; ie=dom(ib)%iep
        js=dom(ib)%jsp; je=dom(ib)%jep
        ks=dom(ib)%ksp; ke=dom(ib)%kep

           dx=dom(ib)%dx
           dy=dom(ib)%dy
           dz=dom(ib)%dz

!..........................................................................
!=== West ===> ..  4=wall  ..   1=Inlet
!..........................................................................
        if (dom(ib)%iprev.lt.0) then
           if (dom(ib)%bc_west.eq.4) then
              do k=ks-1,ke+1
                 do j=js-1,je+1
                    dom(ib)%vis(is,j,k)= rrey
                    dom(ib)%vis(is-1,j,k)= dom(ib)%vis(is,j,k)
                 end do
              end do
           else if (dom(ib)%bc_west.eq.1) then
              do k=ks-1,ke+1
                 do j=js-1,je+1
                    dom(ib)%vis(is-1,j,k)= dom(ib)%vis(is,j,k)
                 end do
              end do
           else if (dom(ib)%bc_west.ge.63) then
	     call wall_function(1,dom(ib)%bc_west)
              do k=ks-1,ke+1
                 do j=js-1,je+1
                    !dom(ib)%vis(is,j,k)= dom(ib)%tauww2(j,k)*dx		!m2/s
                    dom(ib)%vis(is-1,j,k)= dom(ib)%vis(is,j,k)
                 end do
              end do
           else if (dom(ib)%bc_west.eq.61.or.dom(ib)%bc_west.eq.62) then
	     call log_law(1,dom(ib)%bc_west)
              do k=ks-1,ke+1
                 do j=js-1,je+1
                    !dom(ib)%vis(is,j,k)= dom(ib)%tauww2(j,k)*dx		!m2/s
                    dom(ib)%vis(is-1,j,k)= dom(ib)%vis(is,j,k)
                 end do
              end do
           else
              do k=ks-1,ke+1
                 do j=js-1,je+1
                    dom(ib)%vis(is-1,j,k)= dom(ib)%vis(is,j,k)
                 end do
              end do
           end if
        else
           do k=ks-1,ke+1
              do j=js-1,je+1
                 dom(ib)%vis(is-1,j,k)= dom(ib)%vis(is,j,k)
              end do
           end do
        end if
!..........................................................................
!=== East ===> ..  4=wall  ..   2=Outflow
!..........................................................................
        if (dom(ib)%inext.lt.0) then
           if (dom(ib)%bc_east.eq.4) then
              do k=ks-1,ke+1
                 do j=js-1,je+1
                    dom(ib)%vis(ie,j,k)= rrey
                    dom(ib)%vis(ie+1,j,k)= dom(ib)%vis(ie,j,k)
                 end do
              end do
           else if (dom(ib)%bc_east.eq.2) then
              do k=ks-1,ke+1
                 do j=js-1,je+1
                    dom(ib)%vis(ie+1,j,k)= dom(ib)%vis(ie,j,k)
                 end do
              end do
           else if (dom(ib)%bc_east.ge.63) then
	     call wall_function(2,dom(ib)%bc_east)
              do k=ks-1,ke+1
                 do j=js-1,je+1
                    !dom(ib)%vis(ie,j,k)= dom(ib)%tauwe2(j,k)*dx		!m2/s
                    dom(ib)%vis(ie+1,j,k)= dom(ib)%vis(ie,j,k)
                 end do
              end do
           else if (dom(ib)%bc_east.eq.61.or.dom(ib)%bc_east.eq.62) then
	     call log_law(2,dom(ib)%bc_east)
              do k=ks-1,ke+1
                 do j=js-1,je+1
                    !dom(ib)%vis(ie,j,k)= dom(ib)%tauwe2(j,k)*dx		!m2/s
                    dom(ib)%vis(ie+1,j,k)= dom(ib)%vis(ie,j,k)
                 end do
              end do
           else
              do k=ks-1,ke+1
                 do j=js-1,je+1
                    dom(ib)%vis(ie+1,j,k)= dom(ib)%vis(ie,j,k)
                 end do
              end do
           end if
        else
           do k=ks-1,ke+1
              do j=js-1,je+1
                 dom(ib)%vis(ie+1,j,k)= dom(ib)%vis(ie,j,k)
              end do
           end do
        end if
!..........................................................................
!=== South ===> ..   4=wall ..   3=Symmetry
!..........................................................................
        if (dom(ib)%jprev.lt.0) then
           if (dom(ib)%bc_south.eq.4) then
              do k=ks-1,ke+1
                 do i=is-1,ie+1
                    dom(ib)%vis(i,js,k)= rrey
                    dom(ib)%vis(i,js-1,k)= dom(ib)%vis(i,js,k)
                 end do
              end do
           else if (dom(ib)%bc_south.ge.63) then
	     call wall_function(3,dom(ib)%bc_south)
              do k=ks-1,ke+1
                 do i=is-1,ie+1
                   ! dom(ib)%vis(i,js,k)= dom(ib)%tauws2(i,k)*dy		!m2/s
                    dom(ib)%vis(i,js-1,k)= dom(ib)%vis(i,js,k)
                 end do
              end do
           else if(dom(ib)%bc_south.eq.61.or.dom(ib)%bc_south.eq.62)then
	     call log_law(3,dom(ib)%bc_south)
              do k=ks-1,ke+1
                 do i=is-1,ie+1
                    !dom(ib)%vis(i,js,k)= dom(ib)%tauws2(i,k)*dy		!m2/s
                    dom(ib)%vis(i,js-1,k)= dom(ib)%vis(i,js,k)
                 end do
              end do
           else
              do k=ks-1,ke+1
                 do i=is-1,ie+1
                    dom(ib)%vis(i,js-1,k)= dom(ib)%vis(i,js,k)
                 end do
              end do
           end if
        else
           do k=ks-1,ke+1
              do i=is-1,ie+1
                 dom(ib)%vis(i,js-1,k)= dom(ib)%vis(i,js,k)
              end do
           end do
        end if
!..........................................................................
!=== North ===>  ..   4=wall ..   44=moving wall ..  3=Symmetry
!..........................................................................
        if (dom(ib)%jnext.lt.0) then
           if (dom(ib)%bc_north.eq.4) then
              do k=ks-1,ke+1
                 do i=is-1,ie+1
                    dom(ib)%vis(i,je,k) = rrey
                    dom(ib)%vis(i,je+1,k) =  dom(ib)%vis(i,je,k)
                 end do
              end do
           else if (dom(ib)%bc_north.eq.3) then
              do k=ks-1,ke+1
                 do i=is-1,ie+1
                    dom(ib)%vis(i,je+1,k) =  dom(ib)%vis(i,je,k)
                 end do
              end do
           else if (dom(ib)%bc_north.ge.63) then
	     call wall_function(4,dom(ib)%bc_north)
              do k=ks-1,ke+1
                 do i=is-1,ie+1
                   ! dom(ib)%vis(i,je,k)= dom(ib)%tauwn2(i,k)*dy		!m2/s
                    dom(ib)%vis(i,je+1,k) =  dom(ib)%vis(i,je,k)
                 end do
              end do
           else if(dom(ib)%bc_north.eq.61.or.dom(ib)%bc_north.eq.62)then
	     call log_law(4,dom(ib)%bc_north)
              do k=ks-1,ke+1
                 do i=is-1,ie+1
                  !  dom(ib)%vis(i,je,k)= dom(ib)%tauwn2(i,k)*dy		!m2/s
                    dom(ib)%vis(i,je+1,k) =  dom(ib)%vis(i,je,k)
                 end do
              end do
           else
              do k=ks-1,ke+1
                 do i=is-1,ie+1
                    dom(ib)%vis(i,je+1,k) =  dom(ib)%vis(i,je,k)
                 end do
              end do
           end if
        else
           do k=ks-1,ke+1
              do i=is-1,ie+1
                 dom(ib)%vis(i,je+1,k) =  dom(ib)%vis(i,je,k)
              end do
           end do
        end if
!..........................................................................
!=== Bottom ===> ..   4=wall ..   3=Symmetry
!..........................................................................
        if (dom(ib)%kprev.lt.0) then
           if (dom(ib)%bc_bottom.eq.4) then 
              do j=js-1,je+1
                 do i=is-1,ie+1
                    dom(ib)%vis(i,j,ks)= rrey
                    dom(ib)%vis(i,j,ks-1)= dom(ib)%vis(i,j,ks)
                 end do
              end do
           else if (dom(ib)%bc_bottom.eq.3) then
              do j=js-1,je+1
                 do i=is-1,ie+1
                    dom(ib)%vis(i,j,ks-1)= dom(ib)%vis(i,j,ks)
                 end do
              end do
           else if (dom(ib)%bc_bottom.ge.63) then
	     call wall_function(5,dom(ib)%bc_bottom)
              do j=js-1,je+1
                 do i=is-1,ie+1
                   ! dom(ib)%vis(i,j,ks)= dom(ib)%tauwb2(i,j)*dz
                    dom(ib)%vis(i,j,ks-1)= dom(ib)%vis(i,j,ks)
                 end do
              end do
           else if(dom(ib)%bc_bottom.eq.61.or.
     &					dom(ib)%bc_bottom.eq.62)then
	     call log_law(5,dom(ib)%bc_bottom)
              do j=js-1,je+1
                 do i=is-1,ie+1
                  !  dom(ib)%vis(i,j,ks)= dom(ib)%tauwb2(i,j)*dz
                    dom(ib)%vis(i,j,ks-1)= dom(ib)%vis(i,j,ks)
                 end do
              end do
           else
              do j=js-1,je+1
                 do i=is-1,ie+1
                    dom(ib)%vis(i,j,ks-1)= dom(ib)%vis(i,j,ks)
                 end do
              end do
           end if
        else
           do j=js-1,je+1
              do i=is-1,ie+1
                 dom(ib)%vis(i,j,ks-1)= dom(ib)%vis(i,j,ks)
              end do
           end do
        end if
!..........................................................................
!=== Top ===>  ..   4=wall ..     3=Symmetry
!..........................................................................
        if (dom(ib)%knext.lt.0) then
           if (dom(ib)%bc_top.eq.4) then             
              do j=js-1,je+1
                 do i=is-1,ie+1
                    dom(ib)%vis(i,j,ke)   = rrey
                    dom(ib)%vis(i,j,ke+1)   = dom(ib)%vis(i,j,ke)
                 end do
              end do
           else if (dom(ib)%bc_top.eq.3) then
              do j=js-1,je+1
                 do i=is-1,ie+1
                    dom(ib)%vis(i,j,ke+1)   = dom(ib)%vis(i,j,ke)
                 end do
              end do
           else if (dom(ib)%bc_top.ge.63) then
	     call wall_function(6,dom(ib)%bc_top)
              do j=js-1,je+1
                 do i=is-1,ie+1
                   ! dom(ib)%vis(i,j,ke)= dom(ib)%tauwt2(i,j)*dz
                    dom(ib)%vis(i,j,ke+1)   = dom(ib)%vis(i,j,ke)
                 end do
              end do
           else if(dom(ib)%bc_top.eq.61.or.
     &					dom(ib)%bc_top.eq.62)then
	     call log_law(6,dom(ib)%bc_top)
              do j=js-1,je+1
                 do i=is-1,ie+1
                   ! dom(ib)%vis(i,j,ke)   = dom(ib)%tauwt2(i,j)*dz
                    dom(ib)%vis(i,j,ke+1)   = dom(ib)%vis(i,j,ke)
                 end do
              end do
           else
              do j=js-1,je+1
                 do i=is-1,ie+1
                    dom(ib)%vis(i,j,ke+1)   = dom(ib)%vis(i,j,ke)
                 end do
              end do
           end if
        else
           do j=js-1,je+1
              do i=is-1,ie+1
                 dom(ib)%vis(i,j,ke+1)   = dom(ib)%vis(i,j,ke)
              end do
           end do
        end if

        end do

        return
        end subroutine eddyv_wale
!##########################################################################
