!##########################################################################
        subroutine eddyv_smag
!##########################################################################
        use vars
        use multidata
        implicit none
        integer :: i,j,k,chk_wallboun
        integer :: ib,is,ie,js,je,ks,ke
        double precision dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz
        double precision vr_a,vr_b,vr_c,vr_d
        double precision s12,s13,s23,sbet,rrey,utauw,he1
        double precision h1,h2,h3,rh123,dist
        double precision ufv_c,ufv_n1,ufv_n2
        double precision vfv_c,vfv_n1,vfv_n2
        double precision wfv_c,wfv_n1,wfv_n2
        double precision dx,dy,dz
        double precision cs1,delta_grid,l_s,dnmin,yplus,ratio,damp
	  logical vandriest

        vandriest=.true.
        rrey=1.0/Re
        cs1=0.1

        call tauw_noslip								!friction velocity for no slip boundaries

        do ib=1,nbp

        delta_grid=(dom(ib)%dx*dom(ib)%dy*dom(ib)%dz)**(1.0/3.0)

        is=dom(ib)%isp; ie=dom(ib)%iep
        js=dom(ib)%jsp; je=dom(ib)%jep
        ks=dom(ib)%ksp; ke=dom(ib)%kep

        do k=ks,ke
           do j=js,je
              do i=is,ie

              chk_wallboun=0

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
	     end if
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
	     end if
           end if
        end if

!==========================================================================
! ..... DAMPING
!==========================================================================
        utauw = 1.e3
        dnmin = 1e10

        !if (dom(ib)%bc_west.eq.4) then
        if (dom(ib)%bc_west.ge.63) then
		call wall_function(1,dom(ib)%bc_west)
        elseif (dom(ib)%bc_west.eq.61.or.dom(ib)%bc_west.eq.62) then
		call log_law(1,dom(ib)%bc_west)
        end if
           dist=dom(ib)%xc(i)-dom(ib)%xsl
           if(dist.le.dnmin) then
              dnmin = dist
              if (dom(ib)%iprev.lt.0) then
                 chk_wallboun=1
                 he1   = sqrt(abs( dom(ib)%tauww(j,k) ))
                 utauw = min (he1,utauw)
              end if
           end if


        if (dom(ib)%bc_east.ge.63) then
		call wall_function(2,dom(ib)%bc_east)
        elseif (dom(ib)%bc_east.eq.61.or.dom(ib)%bc_east.eq.62) then
		call log_law(2,dom(ib)%bc_east)
        end if
           dist=dom(ib)%xel-dom(ib)%xc(i)
           if(dist.le.dnmin) then
              dnmin = dist
              if (dom(ib)%inext.lt.0) then
                 chk_wallboun=1
                 he1   = sqrt(abs( dom(ib)%tauwe(j,k) ))
                 utauw = min (he1,utauw)
              end if
           end if

        if (dom(ib)%bc_south.ge.63) then
		call wall_function(3,dom(ib)%bc_south)
        elseif (dom(ib)%bc_south.eq.61.or.dom(ib)%bc_south.eq.62) then
		call log_law(3,dom(ib)%bc_south)
        end if
           dist=dom(ib)%yc(j)-dom(ib)%ysl
           if(dist.le.dnmin) then
              dnmin = dist
              if (dom(ib)%jprev.lt.0) then
                 chk_wallboun=1
                 he1   = sqrt(abs( dom(ib)%tauws(i,k) ))
                 utauw = min (he1,utauw)
              end if
           end if

        if (dom(ib)%bc_north.ge.63) then
		call wall_function(4,dom(ib)%bc_north)
        elseif (dom(ib)%bc_north.eq.61.or.dom(ib)%bc_north.eq.62) then
		call log_law(4,dom(ib)%bc_north)
        end if
           dist=dom(ib)%yel-dom(ib)%yc(j)
           if(dist.le.dnmin) then
              dnmin = dist
              if (dom(ib)%jnext.lt.0) then
                 chk_wallboun=1
                 he1   = sqrt(abs( dom(ib)%tauwn(i,k) ))
                 utauw = min (he1,utauw)
              end if
           end if

        if (dom(ib)%bc_bottom.ge.63) then
		call wall_function(5,dom(ib)%bc_bottom)
        elseif (dom(ib)%bc_bottom.eq.61.or.dom(ib)%bc_bottom.eq.62) then
		call log_law(5,dom(ib)%bc_bottom)
        end if
           dist=dom(ib)%zc(k)-dom(ib)%zsl
           if(dist.le.dnmin) then
              dnmin = dist
              if (dom(ib)%kprev.lt.0) then
                 chk_wallboun=1
                 he1   = sqrt(abs( dom(ib)%tauwb(i,j) ))
                 utauw = min (he1,utauw)
              end if
           end if

        if (dom(ib)%bc_top.ge.63) then
		call wall_function(6,dom(ib)%bc_top)
        elseif (dom(ib)%bc_top.eq.61.or.dom(ib)%bc_top.eq.62) then
		call log_law(6,dom(ib)%bc_top)
        end if
           dist=dom(ib)%zel-dom(ib)%zc(k)
           if(dist.le.dnmin) then
              dnmin = dist
              if (dom(ib)%knext.lt.0) then
                 chk_wallboun=1
                 he1   = sqrt(abs( dom(ib)%tauwt(i,j) ))
                 utauw = min (he1,utauw)
              end if
           end if


        damp= 1.0

        if (vandriest .and. chk_wallboun.eq.1) then
           if(abs(dnmin-1e10).lt.0.001) print*,'errorrrrr1'
           yplus = dnmin * Re * utauw
           ratio = min (yplus/25.0,100.0)
           ratio = ratio ** 3

           if (ratio.lt.12.) then
              damp = sqrt (1. - exp (-ratio) )
           endif
        endif

!==========================================================================
! ..... EDDY VISCOSITY CALCULATION
!==========================================================================
        l_s  = (cs1 * delta_grid * damp)**2.0

        s12 = 0.5 * (dudy + dvdx)

        s13 = 0.5 * (dudz + dwdx)

        s23 = 0.5 * (dvdz + dwdy)

        sbet = SQRT (2. * ( dudx*dudx + dvdy*dvdy + dwdz*dwdz +
     &               2. * ( s12*s12   + s13*s13   + s23*s23 ) )  )

        if (L_LSM) rrey=dom(ib)%mu(i,j,k)/dom(ib)%dens(i,j,k)

        dom(ib)%vis(i,j,k) = ( rrey + l_s * sbet )

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
        end subroutine eddyv_smag
!##########################################################################
        subroutine tauw_noslip
!##########################################################################
        use vars
        use multidata
        implicit none
        integer i,j,k,ib
        double precision delta,n_x,n_y,n_z,vnor,vtan
        double precision uc,vc,wc,rrey,small


        rrey=1.0/Re
        small = 1.e-30

        do ib=1,nbp

        dom(ib)%tauww=0.0; dom(ib)%tauwe=0.0
        dom(ib)%tauwn=0.0; dom(ib)%tauws=0.0
        dom(ib)%tauwt=0.0; dom(ib)%tauwb=0.0

        if (dom(ib)%iprev.lt.0) then
           if (dom(ib)%bc_west.eq.4) then						!if it's 6 it was calculated before
              delta=0.5*dom(ib)%dx
              n_x=1.0
              n_y=0.0
              n_z=0.0
              i=dom(ib)%isp
              do k=dom(ib)%ksp,dom(ib)%kep
                 do j=dom(ib)%jsp,dom(ib)%jep
                    uc=0.5*(dom(ib)%u(i,j,k)+dom(ib)%u(i-1,j,k))
                    vc=0.5*(dom(ib)%v(i,j,k)+dom(ib)%v(i,j-1,k))
                    wc=0.5*(dom(ib)%w(i,j,k)+dom(ib)%w(i,j,k-1))
                    vnor=n_x*uc+n_y*vc+n_z*wc
                    vtan=sqrt(abs(uc*uc+vc*vc+wc*wc-vnor*vnor+small))
                    dom(ib)%tauww(j,k) = rrey * vtan / delta
                 end do
              end do
           end if
        end if

        if (dom(ib)%inext.lt.0) then
           if (dom(ib)%bc_east.eq.4) then
              delta=0.5*dom(ib)%dx
              n_x=-1.0
              n_y=0.0
              n_z=0.0
              i=dom(ib)%iep
              do k=dom(ib)%ksp,dom(ib)%kep
                 do j=dom(ib)%jsp,dom(ib)%jep
                    uc=0.5*(dom(ib)%u(i,j,k)+dom(ib)%u(i-1,j,k))
                    vc=0.5*(dom(ib)%v(i,j,k)+dom(ib)%v(i,j-1,k))
                    wc=0.5*(dom(ib)%w(i,j,k)+dom(ib)%w(i,j,k-1))
                    vnor=n_x*uc+n_y*vc+n_z*wc
                    vtan=sqrt(abs(uc*uc+vc*vc+wc*wc-vnor*vnor+small))
                    dom(ib)%tauwe(j,k) = rrey * vtan / delta
                 end do
              end do
           end if
        end if

        if (dom(ib)%jprev.lt.0) then
           if (dom(ib)%bc_south.eq.4) then
              delta=0.5*dom(ib)%dy
              n_x=0.0
              n_y=1.0
              n_z=0.0
              j=dom(ib)%jsp
              do k=dom(ib)%ksp,dom(ib)%kep
                 do i=dom(ib)%isp,dom(ib)%iep
                    uc=0.5*(dom(ib)%u(i,j,k)+dom(ib)%u(i-1,j,k))
                    vc=0.5*(dom(ib)%v(i,j,k)+dom(ib)%v(i,j-1,k))
                    wc=0.5*(dom(ib)%w(i,j,k)+dom(ib)%w(i,j,k-1))
                    vnor=n_x*uc+n_y*vc+n_z*wc
                    vtan=sqrt(abs(uc*uc+vc*vc+wc*wc-vnor*vnor+small))
                    dom(ib)%tauws(i,k) = rrey * vtan / delta
                 end do
              end do
           end if
        end if

        if (dom(ib)%jnext.lt.0) then
           if (dom(ib)%bc_north.eq.4) then
              delta=0.5*dom(ib)%dy
              n_x=0.0
              n_y=-1.0
              n_z=0.0
              j=dom(ib)%jep
              do k=dom(ib)%ksp,dom(ib)%kep
                 do i=dom(ib)%isp,dom(ib)%iep
                    uc=0.5*(dom(ib)%u(i,j,k)+dom(ib)%u(i-1,j,k))
                    vc=0.5*(dom(ib)%v(i,j,k)+dom(ib)%v(i,j-1,k))
                    wc=0.5*(dom(ib)%w(i,j,k)+dom(ib)%w(i,j,k-1))
                    vnor=n_x*uc+n_y*vc+n_z*wc
                    vtan=sqrt(abs(uc*uc+vc*vc+wc*wc-vnor*vnor+small))
                    dom(ib)%tauwn(i,k) = rrey * vtan / delta
                 end do
              end do
           end if
        end if

        if (dom(ib)%kprev.lt.0) then
           if (dom(ib)%bc_bottom.eq.4) then
              delta=0.5*dom(ib)%dz
              n_x=0.0
              n_y=0.0
              n_z=1.0
              k=dom(ib)%ksp
              do j=dom(ib)%jsp,dom(ib)%jep
                 do i=dom(ib)%isp,dom(ib)%iep
                    uc=0.5*(dom(ib)%u(i,j,k)+dom(ib)%u(i-1,j,k))
                    vc=0.5*(dom(ib)%v(i,j,k)+dom(ib)%v(i,j-1,k))
                    wc=0.5*(dom(ib)%w(i,j,k)+dom(ib)%w(i,j,k-1))
                    vnor=n_x*uc+n_y*vc+n_z*wc
                    vtan=sqrt(abs(uc*uc+vc*vc+wc*wc-vnor*vnor+small))
                    dom(ib)%tauwb(i,j) = rrey * vtan / delta
                 end do
              end do
           end if
        end if

        if (dom(ib)%knext.lt.0) then
           if (dom(ib)%bc_top.eq.4) then
              delta=0.5*dom(ib)%dz
              n_x=0.0
              n_y=0.0
              n_z=-1.0
              k=dom(ib)%kep
              do j=dom(ib)%jsp,dom(ib)%jep
                 do i=dom(ib)%isp,dom(ib)%iep
                    uc=0.5*(dom(ib)%u(i,j,k)+dom(ib)%u(i-1,j,k))
                    vc=0.5*(dom(ib)%v(i,j,k)+dom(ib)%v(i,j-1,k))
                    wc=0.5*(dom(ib)%w(i,j,k)+dom(ib)%w(i,j,k-1))
                    vnor=n_x*uc+n_y*vc+n_z*wc
                    vtan=sqrt(abs(uc*uc+vc*vc+wc*wc-vnor*vnor+small))
                    dom(ib)%tauwt(i,j) = rrey * vtan / delta
                 end do
              end do
           end if
        end if

        end do

        return
        end subroutine tauw_noslip
!##########################################################################
