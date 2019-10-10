!##########################################################################
        subroutine diffusion
!##########################################################################
        use vars
        use mpi
        use multidata
        implicit none
        integer :: i,j,k,ib
        integer :: is,ie,js,je,ks,ke
        double precision :: dxx,dyy,dzz

! ------ Compute coefficients for diffusion terms

        if(diff_sch.eq.1) then 
           fac=dt
        else if(diff_sch.eq.2) then 
           fac=dt/2.0
        end if


        do ib=1,nbp
           dom(ib)%ap = 0.0
           dom(ib)%su = 0.0
        end do

        call boundu

        do ib=1,nbp
        dxx=dom(ib)%dx*dom(ib)%dx
        dyy=dom(ib)%dy*dom(ib)%dy
        dzz=dom(ib)%dz*dom(ib)%dz

        is=dom(ib)%isu; ie=dom(ib)%ieu
        js=dom(ib)%jsu; je=dom(ib)%jeu
        ks=dom(ib)%ksu; ke=dom(ib)%keu

        do k=ks,ke
           do i=is,ie
              do j=js,je

                 dom(ib)%aw(i,j,k)=-fac*dom(ib)%vis(i,j,k)/dxx
                 dom(ib)%ae(i,j,k)=-fac*dom(ib)%vis(i,j,k)/dxx
                 dom(ib)%an(i,j,k)=-fac*dom(ib)%vis(i,j,k)/dyy
                 dom(ib)%as(i,j,k)=-fac*dom(ib)%vis(i,j,k)/dyy
                 dom(ib)%at(i,j,k)=-fac*dom(ib)%vis(i,j,k)/dzz
                 dom(ib)%ab(i,j,k)=-fac*dom(ib)%vis(i,j,k)/dzz

                 if (dom(ib)%iprev.lt.0)  dom(ib)%aw(is,j,k)=0.0
                 if (dom(ib)%inext.lt.0)  dom(ib)%ae(ie,j,k)=0.0
                 if (dom(ib)%jprev.lt.0)  dom(ib)%as(i,js,k)=0.0
                 if (dom(ib)%jnext.lt.0)  dom(ib)%an(i,je,k)=0.0
                 if (dom(ib)%kprev.lt.0)  dom(ib)%ab(i,j,ks)=0.0
                 if (dom(ib)%knext.lt.0)  dom(ib)%at(i,j,ke)=0.0

                 dom(ib)%ap(i,j,k) = dom(ib)%ap(i,j,k)-1.0*(
     &  dom(ib)%aw(i,j,k)+dom(ib)%ae(i,j,k)+
     &  dom(ib)%as(i,j,k)+dom(ib)%an(i,j,k)+
     &  dom(ib)%ab(i,j,k)+dom(ib)%at(i,j,k))+1.0

                 if(diff_sch.eq.1) then 
        dom(ib)%su(i,j,k)=dom(ib)%ustar(i,j,k)+dom(ib)%su(i,j,k)
                 else if(diff_sch.eq.2) then 
                    dom(ib)%su(i,j,k)=((1.0-
     &dom(ib)%ap(i,j,k))*dom(ib)%u(i,j,k)-
     &dom(ib)%aw(i,j,k)*dom(ib)%u(i-1,j,k)-
     &dom(ib)%ae(i,j,k)*dom(ib)%u(i+1,j,k)-
     &dom(ib)%as(i,j,k)*dom(ib)%u(i,j-1,k)-
     &dom(ib)%an(i,j,k)*dom(ib)%u(i,j+1,k)-
     &dom(ib)%ab(i,j,k)*dom(ib)%u(i,j,k-1)-
     &dom(ib)%at(i,j,k)*dom(ib)%u(i,j,k+1)+
     &dom(ib)%ustar(i,j,k))+dom(ib)%su(i,j,k)
                 end if

                 if (pressureforce) dom(ib)%su(i,j,k)=
     & dom(ib)%su(i,j,k)+forcn*dt

                 if (LENERGY) dom(ib)%su(i,j,k)=
     & dom(ib)%su(i,j,k)+dt*grx*(1.d0-0.5*beta*
     & (dom(ib)%T(i+1,j,k)+dom(ib)%T(i,j,k)) )

                 dom(ib)%ustar(i,j,k)=0.0
              end do
           end do
        end do
        end do

        call sipsol(11)

        do ib=1,nbp
           dom(ib)%ap = 0.0
           dom(ib)%su = 0.0
        end do

        call boundv

        do ib=1,nbp
        dxx=dom(ib)%dx*dom(ib)%dx
        dyy=dom(ib)%dy*dom(ib)%dy
        dzz=dom(ib)%dz*dom(ib)%dz

        is=dom(ib)%isv; ie=dom(ib)%iev
        js=dom(ib)%jsv; je=dom(ib)%jev
        ks=dom(ib)%ksv; ke=dom(ib)%kev

        do k=ks,ke
           do i=is,ie
              do j=js,je

                 dom(ib)%aw(i,j,k)=-fac*dom(ib)%vis(i,j,k)/dxx
                 dom(ib)%ae(i,j,k)=-fac*dom(ib)%vis(i,j,k)/dxx
                 dom(ib)%an(i,j,k)=-fac*dom(ib)%vis(i,j,k)/dyy
                 dom(ib)%as(i,j,k)=-fac*dom(ib)%vis(i,j,k)/dyy
                 dom(ib)%at(i,j,k)=-fac*dom(ib)%vis(i,j,k)/dzz
                 dom(ib)%ab(i,j,k)=-fac*dom(ib)%vis(i,j,k)/dzz

                 if (dom(ib)%iprev.lt.0)  dom(ib)%aw(is,j,k)=0.0
                 if (dom(ib)%inext.lt.0)  dom(ib)%ae(ie,j,k)=0.0
                 if (dom(ib)%jprev.lt.0)  dom(ib)%as(i,js,k)=0.0
                 if (dom(ib)%jnext.lt.0)  dom(ib)%an(i,je,k)=0.0
                 if (dom(ib)%kprev.lt.0)  dom(ib)%ab(i,j,ks)=0.0
                 if (dom(ib)%knext.lt.0)  dom(ib)%at(i,j,ke)=0.0

                 dom(ib)%ap(i,j,k) = dom(ib)%ap(i,j,k)-1.0*(
     &  dom(ib)%aw(i,j,k)+dom(ib)%ae(i,j,k)+
     &  dom(ib)%as(i,j,k)+dom(ib)%an(i,j,k)+
     &  dom(ib)%ab(i,j,k)+dom(ib)%at(i,j,k))+1.0

                 if(diff_sch.eq.1) then 
        dom(ib)%su(i,j,k)=dom(ib)%vstar(i,j,k)+dom(ib)%su(i,j,k)
                 else if(diff_sch.eq.2) then 
                    dom(ib)%su(i,j,k)=((1.0-
     &dom(ib)%ap(i,j,k))*dom(ib)%v(i,j,k)-
     &dom(ib)%aw(i,j,k)*dom(ib)%v(i-1,j,k)-
     &dom(ib)%ae(i,j,k)*dom(ib)%v(i+1,j,k)-
     &dom(ib)%as(i,j,k)*dom(ib)%v(i,j-1,k)-
     &dom(ib)%an(i,j,k)*dom(ib)%v(i,j+1,k)-
     &dom(ib)%ab(i,j,k)*dom(ib)%v(i,j,k-1)-
     &dom(ib)%at(i,j,k)*dom(ib)%v(i,j,k+1)+
     &dom(ib)%vstar(i,j,k))+dom(ib)%su(i,j,k)
                 end if

                 if (LENERGY) dom(ib)%su(i,j,k)=
     & dom(ib)%su(i,j,k)+dt*gry*(1.d0-0.5*beta*
     & (dom(ib)%T(i,j+1,k)+dom(ib)%T(i,j,k)) )

                 dom(ib)%vstar(i,j,k)=0.0
              end do
           end do
        end do
        end do

        call sipsol(22)

        do ib=1,nbp
           dom(ib)%ap = 0.0
           dom(ib)%su = 0.0
        end do

        call boundw

        do ib=1,nbp
        dxx=dom(ib)%dx*dom(ib)%dx
        dyy=dom(ib)%dy*dom(ib)%dy
        dzz=dom(ib)%dz*dom(ib)%dz

        is=dom(ib)%isw; ie=dom(ib)%iew
        js=dom(ib)%jsw; je=dom(ib)%jew
        ks=dom(ib)%ksw; ke=dom(ib)%kew

        do k=ks,ke
           do i=is,ie
              do j=js,je

                 dom(ib)%aw(i,j,k)=-fac*dom(ib)%vis(i,j,k)/dxx
                 dom(ib)%ae(i,j,k)=-fac*dom(ib)%vis(i,j,k)/dxx
                 dom(ib)%an(i,j,k)=-fac*dom(ib)%vis(i,j,k)/dyy
                 dom(ib)%as(i,j,k)=-fac*dom(ib)%vis(i,j,k)/dyy
                 dom(ib)%at(i,j,k)=-fac*dom(ib)%vis(i,j,k)/dzz
                 dom(ib)%ab(i,j,k)=-fac*dom(ib)%vis(i,j,k)/dzz

                 if (dom(ib)%iprev.lt.0)  dom(ib)%aw(is,j,k)=0.0
                 if (dom(ib)%inext.lt.0)  dom(ib)%ae(ie,j,k)=0.0
                 if (dom(ib)%jprev.lt.0)  dom(ib)%as(i,js,k)=0.0
                 if (dom(ib)%jnext.lt.0)  dom(ib)%an(i,je,k)=0.0
                 if (dom(ib)%kprev.lt.0)  dom(ib)%ab(i,j,ks)=0.0
                 if (dom(ib)%knext.lt.0)  dom(ib)%at(i,j,ke)=0.0

                 dom(ib)%ap(i,j,k) = dom(ib)%ap(i,j,k)-1.0*(
     &  dom(ib)%aw(i,j,k)+dom(ib)%ae(i,j,k)+
     &  dom(ib)%as(i,j,k)+dom(ib)%an(i,j,k)+
     &  dom(ib)%ab(i,j,k)+dom(ib)%at(i,j,k))+1.0

                 if(diff_sch.eq.1) then 
        dom(ib)%su(i,j,k)=dom(ib)%wstar(i,j,k)+dom(ib)%su(i,j,k)
                 else if(diff_sch.eq.2) then 
                    dom(ib)%su(i,j,k)=((1.0-
     &dom(ib)%ap(i,j,k))*dom(ib)%w(i,j,k)-
     &dom(ib)%aw(i,j,k)*dom(ib)%w(i-1,j,k)-
     &dom(ib)%ae(i,j,k)*dom(ib)%w(i+1,j,k)-
     &dom(ib)%as(i,j,k)*dom(ib)%w(i,j-1,k)-
     &dom(ib)%an(i,j,k)*dom(ib)%w(i,j+1,k)-
     &dom(ib)%ab(i,j,k)*dom(ib)%w(i,j,k-1)-
     &dom(ib)%at(i,j,k)*dom(ib)%w(i,j,k+1)+
     &dom(ib)%wstar(i,j,k))+dom(ib)%su(i,j,k)
                 end if

                 if (LENERGY) dom(ib)%su(i,j,k)=
     & dom(ib)%su(i,j,k)+dt*grz*(1.d0-0.5*beta*
     & (dom(ib)%T(i,j,k+1)+dom(ib)%T(i,j,k)) )

                 dom(ib)%wstar(i,j,k)=0.0
              end do
           end do
        end do
        end do

        call sipsol(33)

        if (LROUGH) call rough_velocity

        return
        end subroutine diffusion
!##########################################################################
