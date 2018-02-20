!##########################################################################
      subroutine wall_function(bound,cond)
!			Bruño Fraga Bugallo
!			Cardiff 2014
!			werner/wengle type boundary conditions
!##########################################################################
        use vars
        use multidata
        implicit none
        integer i,j,k,ib,bound,cond
        double precision delta,n_x,n_y,n_z,vnor,vtan,dvtan,sub
        double precision uc,vc,wc,rrey,small,dycell,rycell,vtankr
        double precision tausub,taupow
	  double precision aaa,bbb,const1,const2,const3,const4

        rrey=1.0/Re
        small = 1.e-30

	  SELECT CASE (cond)

	  CASE (63)
!.....specify constants for 1/6 power law ..............................
      data aaa,     bbb
     &   / 8.3, 0.1666666666/

	  CASE (64)
!.....specify constants for 1/7 power law ..............................
      data aaa,     bbb
     &   / 8.3, 0.1428571429/

	  CASE (65)
!.....specify constants for 1/8 power law ..............................
      data aaa,     bbb
     &   / 8.3, 0.125/
	
	end select

!.....constant factors .................................................

      const1 = 0.5 * (1. - bbb) * aaa ** ((1. + bbb) / (1. - bbb))
      const2 = (1. + bbb) / aaa
      const3 = aaa ** (2. / (1. - bbb))
      const4 = 2. / (1. + bbb)

        do ib=1,nbp

	  SELECT CASE (bound)

	  CASE (1) 
              delta=0.5*dom(ib)%dx
              n_x=1.0
              n_y=0.0
              n_z=0.0
              i=dom(ib)%isp
              do k=dom(ib)%ksp-1,dom(ib)%kep+1
                 do j=dom(ib)%jsp-1,dom(ib)%jep+1
                    uc=0.5*(dom(ib)%u(i,j,k)+dom(ib)%u(i-1,j,k))
                    vc=0.5*(dom(ib)%v(i,j,k)+dom(ib)%v(i,j-1,k))
                    wc=0.5*(dom(ib)%w(i,j,k)+dom(ib)%w(i,j,k-1))
                    vnor=n_x*uc+n_y*vc+n_z*wc
                    vtan=sqrt(abs(uc*uc+vc*vc+wc*wc-vnor*vnor+small))

            	  dycell = 2. * delta
            	  rycell = 1. / dycell
            	  vtankr = 0.5 * rrey * rycell * const3
            	  dvtan  = vtankr - vtan
            	  sub    = MAX (SIGN(1.,dvtan),0.)

            	  tausub   = rrey * vtan / delta
            	  taupow   = ( const1 * (rrey * rycell)**(1.+bbb) +
     &                 ( const2 * (rrey * rycell)**bbb) * vtan)
     &                    ** const4
		 	  dom(ib)%tauww(j,k)=(sub*tausub+(1.-sub)*taupow)	!tau_1
			  dom(ib)%tauww2(j,k)=dom(ib)%tauww(j,k)/vtan		!needs to be multiplied by a velocity component to provide tau_1j
                 end do
              end do


	     CASE (2)
              delta=0.5*dom(ib)%dx
              n_x=-1.0
              n_y=0.0
              n_z=0.0
              i=dom(ib)%iep
              do k=dom(ib)%ksp-1,dom(ib)%kep+1
                 do j=dom(ib)%jsp-1,dom(ib)%jep+1
                    uc=0.5*(dom(ib)%u(i,j,k)+dom(ib)%u(i-1,j,k))
                    vc=0.5*(dom(ib)%v(i,j,k)+dom(ib)%v(i,j-1,k))
                    wc=0.5*(dom(ib)%w(i,j,k)+dom(ib)%w(i,j,k-1))
                    vnor=n_x*uc+n_y*vc+n_z*wc
                    vtan=sqrt(abs(uc*uc+vc*vc+wc*wc-vnor*vnor+small))

            	  dycell = 2. * delta
            	  rycell = 1. / dycell
            	  vtankr = 0.5 * rrey * rycell * const3
            	  dvtan  = vtankr - vtan
            	  sub    = MAX (SIGN(1.,dvtan),0.)

            	  tausub   = rrey * vtan / delta
            	  taupow   = ( const1 * (rrey * rycell)**(1.+bbb) +
     &                 ( const2 * (rrey * rycell)**bbb) * vtan)
     &                    ** const4
		 	  dom(ib)%tauwe(j,k)=(sub*tausub+(1.-sub)*taupow)	!units m2/s2
			  dom(ib)%tauwe2(j,k)=dom(ib)%tauwe(j,k)/vtan		!units m/s
                 end do
              end do

           CASE(3)
              delta=0.5*dom(ib)%dy
              n_x=0.0
              n_y=1.0
              n_z=0.0
              j=dom(ib)%jsp

              do k=dom(ib)%ksp-1,dom(ib)%kep+1
                 do i=dom(ib)%isp-1,dom(ib)%iep+1
                    uc=0.5*(dom(ib)%u(i,j,k)+dom(ib)%u(i-1,j,k))
                    vc=0.5*(dom(ib)%v(i,j,k)+dom(ib)%v(i,j-1,k))
                    wc=0.5*(dom(ib)%w(i,j,k)+dom(ib)%w(i,j,k-1))
                    vnor=n_x*uc+n_y*vc+n_z*wc
                    vtan=sqrt(abs(uc*uc+vc*vc+wc*wc-vnor*vnor+small))

            	  dycell = 2. * delta
            	  rycell = 1. / dycell
            	  vtankr = 0.5 * rrey * rycell * const3
            	  dvtan  = vtankr - vtan
            	  sub    = MAX (SIGN(1.,dvtan),0.)

            	  tausub   = rrey * vtan / delta
            	  taupow   = ( const1 * (rrey * rycell)**(1.+bbb) +
     &                 ( const2 * (rrey * rycell)**bbb) * vtan)
     &                    ** const4
		 	  dom(ib)%tauws(i,k)=(sub*tausub+(1.-sub)*taupow)
			  dom(ib)%tauws2(i,k)=dom(ib)%tauws(i,k)/vtan
                 end do
              end do


           CASE (4)
              delta=0.5*dom(ib)%dy
              n_x=0.0
              n_y=-1.0
              n_z=0.0
              j=dom(ib)%jep
              do k=dom(ib)%ksp-1,dom(ib)%kep+1
                 do i=dom(ib)%isp-1,dom(ib)%iep+1
                    uc=0.5*(dom(ib)%u(i,j,k)+dom(ib)%u(i-1,j,k))
                    vc=0.5*(dom(ib)%v(i,j,k)+dom(ib)%v(i,j-1,k))
                    wc=0.5*(dom(ib)%w(i,j,k)+dom(ib)%w(i,j,k-1))
                    vnor=n_x*uc+n_y*vc+n_z*wc
                    vtan=sqrt(abs(uc*uc+vc*vc+wc*wc-vnor*vnor+small))

            	  dycell = 2. * delta
            	  rycell = 1. / dycell
            	  vtankr = 0.5 * rrey * rycell * const3
            	  dvtan  = vtankr - vtan
            	  sub    = MAX (SIGN(1.,dvtan),0.)

            	  tausub   = rrey * vtan / delta
            	  taupow   = ( const1 * (rrey * rycell)**(1.+bbb) +
     &                 ( const2 * (rrey * rycell)**bbb) * vtan)
     &                    ** const4
		 	  dom(ib)%tauwn(i,k)=(sub*tausub+(1.-sub)*taupow)
			  dom(ib)%tauwn2(i,k)=dom(ib)%tauwn(i,k)/vtan
                 end do
              end do

           CASE (5)

              delta=0.5*dom(ib)%dz
              n_x=0.0
              n_y=0.0
              n_z=1.0
              k=dom(ib)%ksp
              do j=dom(ib)%jsp-1,dom(ib)%jep+1
                 do i=dom(ib)%isp-1,dom(ib)%iep+1
                    uc=0.5*(dom(ib)%u(i,j,k)+dom(ib)%u(i-1,j,k))
                    vc=0.5*(dom(ib)%v(i,j,k)+dom(ib)%v(i,j-1,k))
                    wc=0.5*(dom(ib)%w(i,j,k)+dom(ib)%w(i,j,k-1))
			  vnor=n_x*uc+n_y*vc+n_z*wc
                    vtan=sqrt(abs(uc*uc+vc*vc+wc*wc-vnor*vnor+small))

            	  dycell = 2. * delta
            	  rycell = 1. / dycell
            	  vtankr = 0.5 * rrey * rycell * const3
            	  dvtan  = vtankr - vtan
            	  sub    = MAX (SIGN(1.,dvtan),0.)

            	  tausub   = rrey * vtan / delta
            	  taupow   = ( const1 * (rrey * rycell)**(1.+bbb) +
     &                 ( const2 * (rrey * rycell)**bbb) * vtan)
     &                    ** const4
		 	  dom(ib)%tauwb(i,j)=(sub*tausub+(1.-sub)*taupow)
			  dom(ib)%tauwb2(i,j)=dom(ib)%tauwb(i,j)/vtan
                 end do
              end do


           CASE (6)
              delta=0.5*dom(ib)%dz
              n_x=0.0
              n_y=0.0
              n_z=-1.0
              k=dom(ib)%kep
              do j=dom(ib)%jsp-1,dom(ib)%jep+1
                 do i=dom(ib)%isp-1,dom(ib)%iep+1
                    uc=0.5*(dom(ib)%u(i,j,k)+dom(ib)%u(i-1,j,k))
                    vc=0.5*(dom(ib)%v(i,j,k)+dom(ib)%v(i,j-1,k))
                    wc=0.5*(dom(ib)%w(i,j,k)+dom(ib)%w(i,j,k-1))
                    vnor=n_x*uc+n_y*vc+n_z*wc
                    vtan=sqrt(abs(uc*uc+vc*vc+wc*wc-vnor*vnor+small))

            	  dycell = 2. * delta
            	  rycell = 1. / dycell
            	  vtankr = 0.5 * rrey * rycell * const3
            	  dvtan  = vtankr - vtan
            	  sub    = MAX (SIGN(1.,dvtan),0.)

            	  tausub   = rrey * vtan / delta
            	  taupow   = ( const1 * (rrey * rycell)**(1.+bbb) +
     &                 ( const2 * (rrey * rycell)**bbb) * vtan)
     &                    ** const4
		 	  dom(ib)%tauwt(i,j)=(sub*tausub+(1.-sub)*taupow)
			  dom(ib)%tauwt2(i,j)=dom(ib)%tauwt(i,j)/vtan
                 end do
              end do

	  end select
        end do

        return
        end subroutine wall_function
!##########################################################################
