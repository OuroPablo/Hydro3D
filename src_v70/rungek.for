!##########################################################################
        subroutine rungek_conv4th
!##########################################################################
        use vars
        use mpi
        use multidata
        implicit none
        integer :: i,j,k,ib,ii,jj,kk
        double precision :: du2dx,dv2dy,dw2dz
        double precision :: duvdx,duvdy,duwdx,duwdz,dvwdy,dvwdz
        double precision :: sp2,sp1,sn,sm1,sm2

        call boundu

        do ib=1,nbp

           do k=dom(ib)%ksu,dom(ib)%keu
              do i=dom(ib)%isu,dom(ib)%ieu
                 do j=dom(ib)%jsu,dom(ib)%jeu

        ii=i+2; sp2=((-dom(ib)%u(ii+1,j,k)+9.0*dom(ib)%u(ii,j,k)+
     & 9.0*dom(ib)%u(ii-1,j,k)-dom(ib)%u(ii-2,j,k))/16.0)**2
        ii=i+1; sp1=((-dom(ib)%u(ii+1,j,k)+9.0*dom(ib)%u(ii,j,k)+
     & 9.0*dom(ib)%u(ii-1,j,k)-dom(ib)%u(ii-2,j,k))/16.0)**2
        ii=i;   sn=((-dom(ib)%u(ii+1,j,k)+9.0*dom(ib)%u(ii,j,k)+
     & 9.0*dom(ib)%u(ii-1,j,k)-dom(ib)%u(ii-2,j,k))/16.0)**2
        ii=i-1; sm1=((-dom(ib)%u(ii+1,j,k)+9.0*dom(ib)%u(ii,j,k)+
     & 9.0*dom(ib)%u(ii-1,j,k)-dom(ib)%u(ii-2,j,k))/16.0)**2
                 du2dx=(-sp2+27.0*sp1-27.0*sn+sm1)/(24.0*dom(ib)%dx)

        jj=j+1; sp1=((-dom(ib)%u(i,jj+2,k)+9.0*dom(ib)%u(i,jj+1,k)+
     & 9.0*dom(ib)%u(i,jj,k)-dom(ib)%u(i,jj-1,k))/16.0)*
     &((-dom(ib)%v(i+2,jj,k)+9.0*dom(ib)%v(i+1,jj,k)+
     & 9.0*dom(ib)%v(i,jj,k)-dom(ib)%v(i-1,jj,k))/16.0)
        jj=j;   sn=((-dom(ib)%u(i,jj+2,k)+9.0*dom(ib)%u(i,jj+1,k)+
     & 9.0*dom(ib)%u(i,jj,k)-dom(ib)%u(i,jj-1,k))/16.0)*
     &((-dom(ib)%v(i+2,jj,k)+9.0*dom(ib)%v(i+1,jj,k)+
     & 9.0*dom(ib)%v(i,jj,k)-dom(ib)%v(i-1,jj,k))/16.0)
        jj=j-1; sm1=((-dom(ib)%u(i,jj+2,k)+9.0*dom(ib)%u(i,jj+1,k)+
     & 9.0*dom(ib)%u(i,jj,k)-dom(ib)%u(i,jj-1,k))/16.0)*
     &((-dom(ib)%v(i+2,jj,k)+9.0*dom(ib)%v(i+1,jj,k)+
     & 9.0*dom(ib)%v(i,jj,k)-dom(ib)%v(i-1,jj,k))/16.0)
        jj=j-2; sm2=((-dom(ib)%u(i,jj+2,k)+9.0*dom(ib)%u(i,jj+1,k)+
     & 9.0*dom(ib)%u(i,jj,k)-dom(ib)%u(i,jj-1,k))/16.0)*
     &((-dom(ib)%v(i+2,jj,k)+9.0*dom(ib)%v(i+1,jj,k)+
     & 9.0*dom(ib)%v(i,jj,k)-dom(ib)%v(i-1,jj,k))/16.0)
                 duvdy=(-sp1+27.0*sn-27.0*sm1+sm2)/(24.0*dom(ib)%dy)

        kk=k+1; sp1=((-dom(ib)%u(i,j,kk+2)+9.0*dom(ib)%u(i,j,kk+1)+
     & 9.0*dom(ib)%u(i,j,kk)-dom(ib)%u(i,j,kk-1))/16.0)*
     &((-dom(ib)%w(i+2,j,kk)+9.0*dom(ib)%w(i+1,j,kk)+
     & 9.0*dom(ib)%w(i,j,kk)-dom(ib)%w(i-1,j,kk))/16.0)
        kk=k;   sn=((-dom(ib)%u(i,j,kk+2)+9.0*dom(ib)%u(i,j,kk+1)+
     & 9.0*dom(ib)%u(i,j,kk)-dom(ib)%u(i,j,kk-1))/16.0)*
     &((-dom(ib)%w(i+2,j,kk)+9.0*dom(ib)%w(i+1,j,kk)+
     & 9.0*dom(ib)%w(i,j,kk)-dom(ib)%w(i-1,j,kk))/16.0)
        kk=k-1; sm1=((-dom(ib)%u(i,j,kk+2)+9.0*dom(ib)%u(i,j,kk+1)+
     & 9.0*dom(ib)%u(i,j,kk)-dom(ib)%u(i,j,kk-1))/16.0)*
     &((-dom(ib)%w(i+2,j,kk)+9.0*dom(ib)%w(i+1,j,kk)+
     & 9.0*dom(ib)%w(i,j,kk)-dom(ib)%w(i-1,j,kk))/16.0)
        kk=k-2; sm2=((-dom(ib)%u(i,j,kk+2)+9.0*dom(ib)%u(i,j,kk+1)+
     & 9.0*dom(ib)%u(i,j,kk)-dom(ib)%u(i,j,kk-1))/16.0)*
     &((-dom(ib)%w(i+2,j,kk)+9.0*dom(ib)%w(i+1,j,kk)+
     & 9.0*dom(ib)%w(i,j,kk)-dom(ib)%w(i-1,j,kk))/16.0)
                 duwdz=(-sp1+27.0*sn-27.0*sm1+sm2)/(24.0*dom(ib)%dz)

                 dom(ib)%ustar(i,j,k)=(dom(ib)%uoo(i,j,k)-
     & dt*alfapr*(du2dx+duvdy+duwdz))

                 end do
              end do
           end do
        end do

        call boundv

        do ib=1,nbp

           do k=dom(ib)%ksv,dom(ib)%kev
              do i=dom(ib)%isv,dom(ib)%iev
                 do j=dom(ib)%jsv,dom(ib)%jev

        jj=j+2; sp2=((-dom(ib)%v(i,jj+1,k)+9.0*dom(ib)%v(i,jj,k)+
     & 9.0*dom(ib)%v(i,jj-1,k)-dom(ib)%v(i,jj-2,k))/16.0)**2
        jj=j+1; sp1=((-dom(ib)%v(i,jj+1,k)+9.0*dom(ib)%v(i,jj,k)+
     & 9.0*dom(ib)%v(i,jj-1,k)-dom(ib)%v(i,jj-2,k))/16.0)**2
        jj=j  ; sn=((-dom(ib)%v(i,jj+1,k)+9.0*dom(ib)%v(i,jj,k)+
     & 9.0*dom(ib)%v(i,jj-1,k)-dom(ib)%v(i,jj-2,k))/16.0)**2
        jj=j-1; sm1=((-dom(ib)%v(i,jj+1,k)+9.0*dom(ib)%v(i,jj,k)+
     & 9.0*dom(ib)%v(i,jj-1,k)-dom(ib)%v(i,jj-2,k))/16.0)**2
                 dv2dy=(-sp2+27.0*sp1-27.0*sn+sm1)/(24.0*dom(ib)%dy)

        ii=i+1; sp1=((-dom(ib)%v(ii+2,j,k)+9.0*dom(ib)%v(ii+1,j,k)+
     & 9.0*dom(ib)%v(ii,j,k)-dom(ib)%v(ii-1,j,k))/16.0)*
     &((-dom(ib)%u(ii,j+2,k)+9.0*dom(ib)%u(ii,j+1,k)+
     & 9.0*dom(ib)%u(ii,j,k)-dom(ib)%u(ii,j-1,k))/16.0)
        ii=i;   sn=((-dom(ib)%v(ii+2,j,k)+9.0*dom(ib)%v(ii+1,j,k)+
     & 9.0*dom(ib)%v(ii,j,k)-dom(ib)%v(ii-1,j,k))/16.0)*
     &((-dom(ib)%u(ii,j+2,k)+9.0*dom(ib)%u(ii,j+1,k)+
     & 9.0*dom(ib)%u(ii,j,k)-dom(ib)%u(ii,j-1,k))/16.0)
        ii=i-1; sm1=((-dom(ib)%v(ii+2,j,k)+9.0*dom(ib)%v(ii+1,j,k)+
     & 9.0*dom(ib)%v(ii,j,k)-dom(ib)%v(ii-1,j,k))/16.0)*
     &((-dom(ib)%u(ii,j+2,k)+9.0*dom(ib)%u(ii,j+1,k)+
     & 9.0*dom(ib)%u(ii,j,k)-dom(ib)%u(ii,j-1,k))/16.0)
        ii=i-2; sm2=((-dom(ib)%v(ii+2,j,k)+9.0*dom(ib)%v(ii+1,j,k)+
     & 9.0*dom(ib)%v(ii,j,k)-dom(ib)%v(ii-1,j,k))/16.0)*
     &((-dom(ib)%u(ii,j+2,k)+9.0*dom(ib)%u(ii,j+1,k)+
     & 9.0*dom(ib)%u(ii,j,k)-dom(ib)%u(ii,j-1,k))/16.0)
                 duvdx=(-sp1+27.0*sn-27.0*sm1+sm2)/(24.0*dom(ib)%dx)

        kk=k+1; sp1=((-dom(ib)%v(i,j,kk+2)+9.0*dom(ib)%v(i,j,kk+1)+
     & 9.0*dom(ib)%v(i,j,kk)-dom(ib)%v(i,j,kk-1))/16.0)*
     &((-dom(ib)%w(i,j+2,kk)+9.0*dom(ib)%w(i,j+1,kk)+
     & 9.0*dom(ib)%w(i,j,kk)-dom(ib)%w(i,j-1,kk))/16.0)
        kk=k;   sn=((-dom(ib)%v(i,j,kk+2)+9.0*dom(ib)%v(i,j,kk+1)+
     & 9.0*dom(ib)%v(i,j,kk)-dom(ib)%v(i,j,kk-1))/16.0)*
     &((-dom(ib)%w(i,j+2,kk)+9.0*dom(ib)%w(i,j+1,kk)+
     & 9.0*dom(ib)%w(i,j,kk)-dom(ib)%w(i,j-1,kk))/16.0)
        kk=k-1; sm1=((-dom(ib)%v(i,j,kk+2)+9.0*dom(ib)%v(i,j,kk+1)+
     & 9.0*dom(ib)%v(i,j,kk)-dom(ib)%v(i,j,kk-1))/16.0)*
     &((-dom(ib)%w(i,j+2,kk)+9.0*dom(ib)%w(i,j+1,kk)+
     & 9.0*dom(ib)%w(i,j,kk)-dom(ib)%w(i,j-1,kk))/16.0)
        kk=k-2; sm2=((-dom(ib)%v(i,j,kk+2)+9.0*dom(ib)%v(i,j,kk+1)+
     & 9.0*dom(ib)%v(i,j,kk)-dom(ib)%v(i,j,kk-1))/16.0)*
     &((-dom(ib)%w(i,j+2,kk)+9.0*dom(ib)%w(i,j+1,kk)+
     & 9.0*dom(ib)%w(i,j,kk)-dom(ib)%w(i,j-1,kk))/16.0)
                 dvwdz=(-sp1+27.0*sn-27.0*sm1+sm2)/(24.0*dom(ib)%dz)

                 dom(ib)%vstar(i,j,k)=(dom(ib)%voo(i,j,k)-
     & dt*alfapr*(duvdx+dv2dy+dvwdz))

                 end do
              end do
           end do
        end do


        call boundw

        do ib=1,nbp

           do k=dom(ib)%ksw,dom(ib)%kew
              do i=dom(ib)%isw,dom(ib)%iew
                 do j=dom(ib)%jsw,dom(ib)%jew

        kk=k+2; sp2=((-dom(ib)%w(i,j,kk+1)+9.0*dom(ib)%w(i,j,kk)+
     & 9.0*dom(ib)%w(i,j,kk-1)-dom(ib)%w(i,j,kk-2))/16.0)**2
        kk=k+1; sp1=((-dom(ib)%w(i,j,kk+1)+9.0*dom(ib)%w(i,j,kk)+
     & 9.0*dom(ib)%w(i,j,kk-1)-dom(ib)%w(i,j,kk-2))/16.0)**2
        kk=k  ; sn=((-dom(ib)%w(i,j,kk+1)+9.0*dom(ib)%w(i,j,kk)+
     & 9.0*dom(ib)%w(i,j,kk-1)-dom(ib)%w(i,j,kk-2))/16.0)**2
        kk=k-1; sm1=((-dom(ib)%w(i,j,kk+1)+9.0*dom(ib)%w(i,j,kk)+
     & 9.0*dom(ib)%w(i,j,kk-1)-dom(ib)%w(i,j,kk-2))/16.0)**2
                 dw2dz=(-sp2+27.0*sp1-27.0*sn+sm1)/(24.0*dom(ib)%dz)

        ii=i+1; sp1=((-dom(ib)%w(ii+2,j,k)+9.0*dom(ib)%w(ii+1,j,k)+
     & 9.0*dom(ib)%w(ii,j,k)-dom(ib)%w(ii-1,j,k))/16.0)*
     &((-dom(ib)%u(ii,j,k+2)+9.0*dom(ib)%u(ii,j,k+1)+
     & 9.0*dom(ib)%u(ii,j,k)-dom(ib)%u(ii,j,k-1))/16.0)
        ii=i;   sn=((-dom(ib)%w(ii+2,j,k)+9.0*dom(ib)%w(ii+1,j,k)+
     & 9.0*dom(ib)%w(ii,j,k)-dom(ib)%w(ii-1,j,k))/16.0)*
     &((-dom(ib)%u(ii,j,k+2)+9.0*dom(ib)%u(ii,j,k+1)+
     & 9.0*dom(ib)%u(ii,j,k)-dom(ib)%u(ii,j,k-1))/16.0)
        ii=i-1; sm1=((-dom(ib)%w(ii+2,j,k)+9.0*dom(ib)%w(ii+1,j,k)+
     & 9.0*dom(ib)%w(ii,j,k)-dom(ib)%w(ii-1,j,k))/16.0)*
     &((-dom(ib)%u(ii,j,k+2)+9.0*dom(ib)%u(ii,j,k+1)+
     & 9.0*dom(ib)%u(ii,j,k)-dom(ib)%u(ii,j,k-1))/16.0)
        ii=i-2; sm2=((-dom(ib)%w(ii+2,j,k)+9.0*dom(ib)%w(ii+1,j,k)+
     & 9.0*dom(ib)%w(ii,j,k)-dom(ib)%w(ii-1,j,k))/16.0)*
     &((-dom(ib)%u(ii,j,k+2)+9.0*dom(ib)%u(ii,j,k+1)+
     & 9.0*dom(ib)%u(ii,j,k)-dom(ib)%u(ii,j,k-1))/16.0)
                 duwdx=(-sp1+27.0*sn-27.0*sm1+sm2)/(24.0*dom(ib)%dx)

        jj=j+1; sp1=((-dom(ib)%w(i,jj+2,k)+9.0*dom(ib)%w(i,jj+1,k)+
     & 9.0*dom(ib)%w(i,jj,k)-dom(ib)%w(i,jj-1,k))/16.0)*
     &((-dom(ib)%v(i,jj,k+2)+9.0*dom(ib)%v(i,jj,k+1)+
     & 9.0*dom(ib)%v(i,jj,k)-dom(ib)%v(i,jj,k-1))/16.0)
        jj=j;   sn=((-dom(ib)%w(i,jj+2,k)+9.0*dom(ib)%w(i,jj+1,k)+
     & 9.0*dom(ib)%w(i,jj,k)-dom(ib)%w(i,jj-1,k))/16.0)*
     &((-dom(ib)%v(i,jj,k+2)+9.0*dom(ib)%v(i,jj,k+1)+
     & 9.0*dom(ib)%v(i,jj,k)-dom(ib)%v(i,jj,k-1))/16.0)
        jj=j-1; sm1=((-dom(ib)%w(i,jj+2,k)+9.0*dom(ib)%w(i,jj+1,k)+
     & 9.0*dom(ib)%w(i,jj,k)-dom(ib)%w(i,jj-1,k))/16.0)*
     &((-dom(ib)%v(i,jj,k+2)+9.0*dom(ib)%v(i,jj,k+1)+
     & 9.0*dom(ib)%v(i,jj,k)-dom(ib)%v(i,jj,k-1))/16.0)
        jj=j-2; sm2=((-dom(ib)%w(i,jj+2,k)+9.0*dom(ib)%w(i,jj+1,k)+
     & 9.0*dom(ib)%w(i,jj,k)-dom(ib)%w(i,jj-1,k))/16.0)*
     &((-dom(ib)%v(i,jj,k+2)+9.0*dom(ib)%v(i,jj,k+1)+
     & 9.0*dom(ib)%v(i,jj,k)-dom(ib)%v(i,jj,k-1))/16.0)
                 dvwdy=(-sp1+27.0*sn-27.0*sm1+sm2)/(24.0*dom(ib)%dy)

                 dom(ib)%wstar(i,j,k)=(dom(ib)%woo(i,j,k)-
     & dt*alfapr*(duwdx+dvwdy+dw2dz))  

                 end do
              end do
           end do

        end do


        return
        end subroutine rungek_conv4th
!##########################################################################
        subroutine rungek_diff4th
!##########################################################################
        use vars
        use mpi
        use multidata
        implicit none
        integer :: i,j,k,ib
        integer :: is,ie,js,je,ks,ke
        double precision :: dxx,dyy,dzz,visc,diff

        fac=1.0

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
                 visc=dom(ib)%vis(i,j,k)
 
                 diff=visc*(
     &((-dom(ib)%u(i+2,j,k)+16.0*dom(ib)%u(i+1,j,k)
     &  -30.0*dom(ib)%u(i,j,k)+16.0*dom(ib)%u(i-1,j,k)
     &  -dom(ib)%u(i-2,j,k))/(12.0*dxx))+
     &((-dom(ib)%u(i,j+2,k)+16.0*dom(ib)%u(i,j+1,k)
     &  -30.0*dom(ib)%u(i,j,k)+16.0*dom(ib)%u(i,j-1,k)
     &  -dom(ib)%u(i,j-2,k))/(12.0*dyy))+
     &((-dom(ib)%u(i,j,k+2)+16.0*dom(ib)%u(i,j,k+1)
     &  -30.0*dom(ib)%u(i,j,k)+16.0*dom(ib)%u(i,j,k-1)
     &  -dom(ib)%u(i,j,k-2))/(12.0*dzz)))

                 dom(ib)%ustar(i,j,k)=(dom(ib)%ustar(i,j,k)+
     & dt*alfapr*diff) 

                 if (pressureforce) dom(ib)%ustar(i,j,k)=
     & dom(ib)%ustar(i,j,k)+dt*alfapr*forcn

              end do
           end do
        end do


        is=dom(ib)%isv; ie=dom(ib)%iev
        js=dom(ib)%jsv; je=dom(ib)%jev
        ks=dom(ib)%ksv; ke=dom(ib)%kev

        do k=ks,ke
           do i=is,ie
              do j=js,je
                 visc=dom(ib)%vis(i,j+1,k)

                 diff=visc*(
     &((-dom(ib)%v(i+2,j,k)+16.0*dom(ib)%v(i+1,j,k)
     &  -30.0*dom(ib)%v(i,j,k)+16.0*dom(ib)%v(i-1,j,k)
     &  -dom(ib)%v(i-2,j,k))/(12.0*dxx))+
     &((-dom(ib)%v(i,j+2,k)+16.0*dom(ib)%v(i,j+1,k)
     &  -30.0*dom(ib)%v(i,j,k)+16.0*dom(ib)%v(i,j-1,k)
     &  -dom(ib)%v(i,j-2,k))/(12.0*dyy))+
     &((-dom(ib)%v(i,j,k+2)+16.0*dom(ib)%v(i,j,k+1)
     &  -30.0*dom(ib)%v(i,j,k)+16.0*dom(ib)%v(i,j,k-1)
     &  -dom(ib)%v(i,j,k-2))/(12.0*dzz)))
 
                 dom(ib)%vstar(i,j,k)=(dom(ib)%vstar(i,j,k)+
     & dt*alfapr*diff)

                 if (pressureforce_y) dom(ib)%vstar(i,j,k)=
     & dom(ib)%vstar(i,j,k)+dt*alfapr*forcn_y

              end do
           end do
        end do

        is=dom(ib)%isw; ie=dom(ib)%iew
        js=dom(ib)%jsw; je=dom(ib)%jew
        ks=dom(ib)%ksw; ke=dom(ib)%kew

        do k=ks,ke
           do i=is,ie
              do j=js,je
                 visc=dom(ib)%vis(i,j,k+1)

                 diff=visc*(
     &((-dom(ib)%w(i+2,j,k)+16.0*dom(ib)%w(i+1,j,k)
     &  -30.0*dom(ib)%w(i,j,k)+16.0*dom(ib)%w(i-1,j,k)
     &  -dom(ib)%w(i-2,j,k))/(12.0*dxx))+
     &((-dom(ib)%w(i,j+2,k)+16.0*dom(ib)%w(i,j+1,k)
     &  -30.0*dom(ib)%w(i,j,k)+16.0*dom(ib)%w(i,j-1,k)
     &  -dom(ib)%w(i,j-2,k))/(12.0*dyy))+
     &((-dom(ib)%w(i,j,k+2)+16.0*dom(ib)%w(i,j,k+1)
     &  -30.0*dom(ib)%w(i,j,k)+16.0*dom(ib)%w(i,j,k-1)
     &  -dom(ib)%w(i,j,k-2))/(12.0*dzz)))
 
                 dom(ib)%wstar(i,j,k)=(dom(ib)%wstar(i,j,k)+
     & dt*alfapr*diff) 

                 if (pressureforce_z) dom(ib)%wstar(i,j,k)=
     & dom(ib)%wstar(i,j,k)+dt*alfapr*forcn_z

              end do
           end do
        end do

        end do

        if (LENERGY) call mom_buo
        if (LROUGH) call rough_velocity
        return
        end subroutine rungek_diff4th
!##########################################################################
        subroutine rungek_conv2nd
!##########################################################################
        use vars
        use mpi
        use multidata
        implicit none
        integer :: i,j,k,ib
        double precision :: du2dx,dv2dy,dw2dz
        double precision :: up12,um12,vp12,vm12,wp12,wm12
        double precision :: duvdx,duvdy,duwdx,duwdz,dvwdy,dvwdz

        call boundu

        do ib=1,nbp

           do k=dom(ib)%ksu,dom(ib)%keu
              do i=dom(ib)%isu,dom(ib)%ieu
                 do j=dom(ib)%jsu,dom(ib)%jeu

                 up12=0.5*(dom(ib)%u(i,j,k)+dom(ib)%u(i+1,j,k))
                 um12=0.5*(dom(ib)%u(i,j,k)+dom(ib)%u(i-1,j,k))
                 du2dx=(up12**2-um12**2)/dom(ib)%dx
                 duvdy=0.25*( (dom(ib)%v(i,j,k)+dom(ib)%v(i+1,j,k))*
     &(dom(ib)%u(i,j,k)+dom(ib)%u(i,j+1,k)) -
     &(dom(ib)%v(i,j-1,k)+dom(ib)%v(i+1,j-1,k))*
     &(dom(ib)%u(i,j,k)+dom(ib)%u(i,j-1,k)) )/dom(ib)%dy
                 duwdz=0.25*( (dom(ib)%w(i,j,k)+dom(ib)%w(i+1,j,k))*
     &(dom(ib)%u(i,j,k)+dom(ib)%u(i,j,k+1)) -
     &(dom(ib)%w(i,j,k-1)+dom(ib)%w(i+1,j,k-1))*
     &(dom(ib)%u(i,j,k)+dom(ib)%u(i,j,k-1)) )/dom(ib)%dz

                 dom(ib)%ustar(i,j,k)=(dom(ib)%uoo(i,j,k)-
     & dt*alfapr*(du2dx+duvdy+duwdz))

                 end do
              end do
           end do
        end do

        call boundv

        do ib=1,nbp

           do k=dom(ib)%ksv,dom(ib)%kev
              do i=dom(ib)%isv,dom(ib)%iev
                 do j=dom(ib)%jsv,dom(ib)%jev

                 vp12=0.5*(dom(ib)%v(i,j,k)+dom(ib)%v(i,j+1,k))
                 vm12=0.5*(dom(ib)%v(i,j,k)+dom(ib)%v(i,j-1,k))
                 dv2dy=(vp12**2-vm12**2)/dom(ib)%dy
                 duvdx=0.25*( (dom(ib)%u(i,j,k)+dom(ib)%u(i,j+1,k))*
     &(dom(ib)%v(i,j,k)+dom(ib)%v(i+1,j,k)) -
     &(dom(ib)%u(i-1,j,k)+dom(ib)%u(i-1,j+1,k))*
     &(dom(ib)%v(i,j,k)+dom(ib)%v(i-1,j,k)) )/dom(ib)%dx
                 dvwdz=0.25*( (dom(ib)%w(i,j,k)+dom(ib)%w(i,j+1,k))*
     &(dom(ib)%v(i,j,k)+dom(ib)%v(i,j,k+1)) -
     &(dom(ib)%w(i,j,k-1)+dom(ib)%w(i,j+1,k-1))*
     &(dom(ib)%v(i,j,k)+dom(ib)%v(i,j,k-1)) )/dom(ib)%dz

                 dom(ib)%vstar(i,j,k)=(dom(ib)%voo(i,j,k)-
     & dt*alfapr*(duvdx+dv2dy+dvwdz))

                 end do
              end do
           end do
        end do


        call boundw

        do ib=1,nbp

           do k=dom(ib)%ksw,dom(ib)%kew
              do i=dom(ib)%isw,dom(ib)%iew
                 do j=dom(ib)%jsw,dom(ib)%jew

                 wp12=0.5*(dom(ib)%w(i,j,k)+dom(ib)%w(i,j,k+1))
                 wm12=0.5*(dom(ib)%w(i,j,k)+dom(ib)%w(i,j,k-1))
                 dw2dz=(wp12**2-wm12**2)/dom(ib)%dz
                 duwdx=0.25*( (dom(ib)%u(i,j,k)+dom(ib)%u(i,j,k+1))*
     &(dom(ib)%w(i,j,k)+dom(ib)%w(i+1,j,k)) -
     &(dom(ib)%u(i-1,j,k)+dom(ib)%u(i-1,j,k+1))*
     &(dom(ib)%w(i,j,k)+dom(ib)%w(i-1,j,k)) )/dom(ib)%dx
                 dvwdy=0.25*( (dom(ib)%v(i,j,k)+dom(ib)%v(i,j,k+1))*
     &(dom(ib)%w(i,j,k)+dom(ib)%w(i,j+1,k)) -
     &(dom(ib)%v(i,j-1,k)+dom(ib)%v(i,j-1,k+1))*
     &(dom(ib)%w(i,j,k)+dom(ib)%w(i,j-1,k)) )/dom(ib)%dy

                 dom(ib)%wstar(i,j,k)=(dom(ib)%woo(i,j,k)-
     & dt*alfapr*(duwdx+dvwdy+dw2dz))  

                 end do
              end do
           end do

        end do

        return
        end subroutine rungek_conv2nd
!##########################################################################
        subroutine rungek_convWENO
!##########################################################################
        use vars
        use mpi
        use multidata
        implicit none
        integer :: i,j,k,ib
        double precision :: du2dx,dv2dy,dw2dz
        double precision :: duvdx,duvdy,duwdx,duwdz,dvwdy,dvwdz
        double precision :: uijk,vijk,wijk,dudx,dudy,dudz
        double precision :: dvdx,dvdy,dvdz,dwdx,dwdy,dwdz

        call boundu
        call HJ_WENO_dx(1)
        call HJ_WENO_dy(1)
        call HJ_WENO_dz(1)

        do ib=1,nbp

           do k=dom(ib)%ksu,dom(ib)%keu
              do i=dom(ib)%isu,dom(ib)%ieu
                 do j=dom(ib)%jsu,dom(ib)%jeu

                    if (dom(ib)%u(i,j,k).gt.0.0) then
                       dudx = dom(ib)%dphi_dxminus(i,j,k)
                    else if (dom(ib)%u(i,j,k).lt.0.0) then
                       dudx = dom(ib)%dphi_dxplus(i,j,k)
                    else  
                       dudx = 0.0  
                    end if             

                    vijk=0.25*(dom(ib)%v(i,j,k)+dom(ib)%v(i+1,j,k)+
     &     dom(ib)%v(i,j-1,k)+dom(ib)%v(i+1,j-1,k))
                    if (vijk.gt.0.0) then
                       dudy = dom(ib)%dphi_dyminus(i,j,k)
                    else if (vijk.lt.0.0) then  
                       dudy = dom(ib)%dphi_dyplus(i,j,k)
                    else  
                       dudy = 0.0   
                    end if     

                    wijk=0.25*(dom(ib)%w(i,j,k)+dom(ib)%w(i+1,j,k)+
     &     dom(ib)%w(i,j,k-1)+dom(ib)%w(i+1,j,k-1))
                    if (wijk.gt.0.0) then
                       dudz = dom(ib)%dphi_dzminus(i,j,k)
                    else if (wijk.lt.0.0) then 
                       dudz = dom(ib)%dphi_dzplus(i,j,k)
                    else 
                       dudz = 0.0    
                    end if

                    du2dx=dom(ib)%u(i,j,k)*dudx
                    duvdy=vijk*dudy
                    duwdz=wijk*dudz

                 dom(ib)%ustar(i,j,k)=(dom(ib)%uoo(i,j,k)-
     & dt*alfapr*(du2dx+duvdy+duwdz))

                 end do
              end do
           end do
        end do

        call boundv
        call HJ_WENO_dx(2)
        call HJ_WENO_dy(2)
        call HJ_WENO_dz(2)

        do ib=1,nbp

           do k=dom(ib)%ksv,dom(ib)%kev
              do i=dom(ib)%isv,dom(ib)%iev
                 do j=dom(ib)%jsv,dom(ib)%jev

                    uijk=0.25*(dom(ib)%u(i,j,k)+dom(ib)%u(i,j+1,k)+
     &     dom(ib)%u(i-1,j,k)+dom(ib)%u(i-1,j+1,k))
                    if (uijk.gt.0.0) then
                       dvdx = dom(ib)%dphi_dxminus(i,j,k)
                    else if (uijk.lt.0.0) then
                       dvdx = dom(ib)%dphi_dxplus(i,j,k)
                    else  
                       dvdx = 0.0  
                    end if             

                    if (dom(ib)%v(i,j,k).gt.0.0) then
                       dvdy = dom(ib)%dphi_dyminus(i,j,k)
                    else if (dom(ib)%v(i,j,k).lt.0.0) then  
                       dvdy = dom(ib)%dphi_dyplus(i,j,k)
                    else  
                       dvdy = 0.0   
                    end if     

                    wijk=0.25*(dom(ib)%w(i,j,k)+dom(ib)%w(i,j+1,k)+
     &     dom(ib)%w(i,j,k-1)+dom(ib)%w(i,j+1,k-1))
                    if (wijk.gt.0.0) then
                       dvdz = dom(ib)%dphi_dzminus(i,j,k)
                    else if (wijk.lt.0.0) then 
                       dvdz = dom(ib)%dphi_dzplus(i,j,k)
                    else 
                       dvdz = 0.0    
                    end if

                    duvdx=uijk*dvdx
                    dv2dy=dom(ib)%v(i,j,k)*dvdy
                    dvwdz=wijk*dvdz

                 dom(ib)%vstar(i,j,k)=(dom(ib)%voo(i,j,k)-
     & dt*alfapr*(duvdx+dv2dy+dvwdz))

                 end do
              end do
           end do
        end do


        call boundw
        call HJ_WENO_dx(3)
        call HJ_WENO_dy(3)
        call HJ_WENO_dz(3)

        do ib=1,nbp

           do k=dom(ib)%ksw,dom(ib)%kew
              do i=dom(ib)%isw,dom(ib)%iew
                 do j=dom(ib)%jsw,dom(ib)%jew

                    uijk=0.25*(dom(ib)%u(i,j,k)+dom(ib)%u(i,j,k+1)+
     &     dom(ib)%u(i-1,j,k)+dom(ib)%u(i-1,j,k+1))
                    if (uijk.gt.0.0) then
                       dwdx = dom(ib)%dphi_dxminus(i,j,k)
                    else if (uijk.lt.0.0) then
                       dwdx = dom(ib)%dphi_dxplus(i,j,k)
                    else  
                       dwdx = 0.0  
                    end if             

                    vijk=0.25*(dom(ib)%v(i,j,k)+dom(ib)%v(i,j,k+1)+
     &     dom(ib)%v(i,j-1,k)+dom(ib)%v(i,j-1,k+1))
                    if (vijk.gt.0.0) then
                       dwdy = dom(ib)%dphi_dyminus(i,j,k)
                    else if (vijk.lt.0.0) then  
                       dwdy = dom(ib)%dphi_dyplus(i,j,k)
                    else  
                       dwdy = 0.0   
                    end if     

                    if (dom(ib)%w(i,j,k).gt.0.0) then
                       dwdz = dom(ib)%dphi_dzminus(i,j,k)
                    else if (dom(ib)%w(i,j,k).lt.0.0) then 
                       dwdz = dom(ib)%dphi_dzplus(i,j,k)
                    else 
                       dwdz = 0.0    
                    end if

                    duwdx=uijk*dwdx
                    dvwdy=vijk*dwdy
                    dw2dz=dom(ib)%w(i,j,k)*dwdz

                 dom(ib)%wstar(i,j,k)=(dom(ib)%woo(i,j,k)-
     & dt*alfapr*(duwdx+dvwdy+dw2dz))  

                 end do
              end do
           end do

        end do

        return
        end subroutine rungek_convWENO
!##########################################################################
        subroutine rungek_diff2nd
!##########################################################################
        use vars
        use mpi
        use multidata
        implicit none
        integer :: i,j,k,ib
        integer :: is,ie,js,je,ks,ke
        double precision :: dxx,dyy,dzz,visc,diff

        fac=1.d0

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
                 visc=dom(ib)%vis(i+1,j,k)
 
                 dom(ib)%aw(i,j,k)=-visc/dxx
                 dom(ib)%ae(i,j,k)=-visc/dxx
                 dom(ib)%an(i,j,k)=-visc/dyy
                 dom(ib)%as(i,j,k)=-visc/dyy
                 dom(ib)%at(i,j,k)=-visc/dzz
                 dom(ib)%ab(i,j,k)=-visc/dzz

                 dom(ib)%ap(i,j,k) = -1.0*(
     &dom(ib)%aw(i,j,k)+dom(ib)%ae(i,j,k)+
     &dom(ib)%as(i,j,k)+dom(ib)%an(i,j,k)+
     &dom(ib)%ab(i,j,k)+dom(ib)%at(i,j,k))

                 diff=-(dom(ib)%ap(i,j,k)*dom(ib)%u(i,j,k)+
     & dom(ib)%an(i,j,k)*dom(ib)%u(i,j+1,k) +
     & dom(ib)%as(i,j,k)*dom(ib)%u(i,j-1,k)+
     & dom(ib)%ae(i,j,k)*dom(ib)%u(i+1,j,k) +
     & dom(ib)%aw(i,j,k)*dom(ib)%u(i-1,j,k)+
     & dom(ib)%at(i,j,k)*dom(ib)%u(i,j,k+1) +
     & dom(ib)%ab(i,j,k)*dom(ib)%u(i,j,k-1))

          dom(ib)%ustar(i,j,k)=dom(ib)%ustar(i,j,k)+dt*alfapr*diff 

        if (L_LSM) then! .or. L_LSMbase) then
          dom(ib)%ustar(i,j,k)=dom(ib)%ustar(i,j,k)+
     & dt*alfapr*(grx-grz*sin(atan(slope)))
        end if

        if (pressureforce) then
	    if (L_LSM) then
	      if (dom(ib)%phi(i,j,k) .ge. 0.0) then
 	        dom(ib)%ustar(i,j,k)=dom(ib)%ustar(i,j,k)+dt*alfapr*forcn
	      endif
	    else if (L_LSMbase) then
		if (dom(ib)%zc(k).le.length) then
		dom(ib)%ustar(i,j,k)=dom(ib)%ustar(i,j,k)+dt*alfapr*forcn
		endif
	    else
 	      dom(ib)%ustar(i,j,k)=dom(ib)%ustar(i,j,k)+dt*alfapr*forcn
	    endif
	  endif

              end do
           end do
        end do

        is=dom(ib)%isv; ie=dom(ib)%iev
        js=dom(ib)%jsv; je=dom(ib)%jev
        ks=dom(ib)%ksv; ke=dom(ib)%kev

        do k=ks,ke
           do i=is,ie
              do j=js,je
                 visc=dom(ib)%vis(i,j+1,k)

                 dom(ib)%aw(i,j,k)=-visc/dxx
                 dom(ib)%ae(i,j,k)=-visc/dxx
                 dom(ib)%an(i,j,k)=-visc/dyy
                 dom(ib)%as(i,j,k)=-visc/dyy
                 dom(ib)%at(i,j,k)=-visc/dzz
                 dom(ib)%ab(i,j,k)=-visc/dzz

                 dom(ib)%ap(i,j,k) = -1.0*(
     &dom(ib)%aw(i,j,k)+dom(ib)%ae(i,j,k)+
     &dom(ib)%as(i,j,k)+dom(ib)%an(i,j,k)+
     &dom(ib)%ab(i,j,k)+dom(ib)%at(i,j,k))

                 diff=-(dom(ib)%ap(i,j,k)*dom(ib)%v(i,j,k)+
     & dom(ib)%an(i,j,k)*dom(ib)%v(i,j+1,k) +
     & dom(ib)%as(i,j,k)*dom(ib)%v(i,j-1,k)+
     & dom(ib)%ae(i,j,k)*dom(ib)%v(i+1,j,k) +
     & dom(ib)%aw(i,j,k)*dom(ib)%v(i-1,j,k)+
     & dom(ib)%at(i,j,k)*dom(ib)%v(i,j,k+1) +
     & dom(ib)%ab(i,j,k)*dom(ib)%v(i,j,k-1))
 
           dom(ib)%vstar(i,j,k)=dom(ib)%vstar(i,j,k)+dt*alfapr*diff

        if (L_LSM) then! .or. L_LSMbase) then
          dom(ib)%vstar(i,j,k)=dom(ib)%vstar(i,j,k)+dt*alfapr*gry
        end if

        if (pressureforce_y) then				!Pablo 03/2018
	    if (L_LSM) then
	      if (dom(ib)%phi(i,j,k) .ge. 0.0) then
 	       dom(ib)%vstar(i,j,k)=dom(ib)%vstar(i,j,k)+dt*alfapr*forcn_y
	      endif
	    else if (L_LSMbase) then
		if (dom(ib)%zc(k).le.length) then
 	       dom(ib)%vstar(i,j,k)=dom(ib)%vstar(i,j,k)+dt*alfapr*forcn_y
		endif
	    else
 	       dom(ib)%vstar(i,j,k)=dom(ib)%vstar(i,j,k)+dt*alfapr*forcn_y
	    endif
	  endif

              end do
           end do
        end do

        is=dom(ib)%isw; ie=dom(ib)%iew
        js=dom(ib)%jsw; je=dom(ib)%jew
        ks=dom(ib)%ksw; ke=dom(ib)%kew

        do k=ks,ke
           do i=is,ie
              do j=js,je
                 visc=dom(ib)%vis(i,j,k+1)

                 dom(ib)%aw(i,j,k)=-visc/dxx
                 dom(ib)%ae(i,j,k)=-visc/dxx
                 dom(ib)%an(i,j,k)=-visc/dyy
                 dom(ib)%as(i,j,k)=-visc/dyy
                 dom(ib)%at(i,j,k)=-visc/dzz
                 dom(ib)%ab(i,j,k)=-visc/dzz

                 dom(ib)%ap(i,j,k) = -1.0*(
     &dom(ib)%aw(i,j,k)+dom(ib)%ae(i,j,k)+
     &dom(ib)%as(i,j,k)+dom(ib)%an(i,j,k)+
     &dom(ib)%ab(i,j,k)+dom(ib)%at(i,j,k))

                 diff=-(dom(ib)%ap(i,j,k)*dom(ib)%w(i,j,k)+
     & dom(ib)%an(i,j,k)*dom(ib)%w(i,j+1,k) +
     & dom(ib)%as(i,j,k)*dom(ib)%w(i,j-1,k)+
     & dom(ib)%ae(i,j,k)*dom(ib)%w(i+1,j,k) +
     & dom(ib)%aw(i,j,k)*dom(ib)%w(i-1,j,k)+
     & dom(ib)%at(i,j,k)*dom(ib)%w(i,j,k+1) +
     & dom(ib)%ab(i,j,k)*dom(ib)%w(i,j,k-1))
 
           dom(ib)%wstar(i,j,k)=dom(ib)%wstar(i,j,k)+dt*alfapr*diff

        if (L_LSM) then
!===========================================================
!Artificial damping zone Absorbing waves
!===========================================================
!	   X_R_e = (dom(ib)%xc(i)-3.2d0)/0.8d0                              !sponge layer zone,artificial viscosity (2L)
	    
!            if (X_R_e.ge.0.0)  then
!	     gamma = (exp(X_R_e**2.0)-1.0)/(exp(1.0)-1.0)
!	     dom(ib)%wstar(i,j,k)=dom(ib)%wstar(i,j,k)+dt*alfapr*grz*
!     & cos(atan(slope))-dt*alfapr*(100.0+0.0*abs(dom(ib)%wstar(i,j,k)))
!     & *gamma*dom(ib)%wstar(i,j,k)

      dom(ib)%wstar(i,j,k)=dom(ib)%wstar(i,j,k)+dt*alfapr*grz*
     & cos(atan(slope))
        end if

        if (pressureforce_z) then
	    if (L_LSM) then
	      if (dom(ib)%phi(i,j,k) .ge. 0.0) then
 	       dom(ib)%wstar(i,j,k)=dom(ib)%wstar(i,j,k)+dt*alfapr*forcn_z
	      endif
	    else if (L_LSMbase) then
		if (dom(ib)%zc(k).le.length) then
 	       dom(ib)%wstar(i,j,k)=dom(ib)%wstar(i,j,k)+dt*alfapr*forcn_z
		endif
	    else
 	       dom(ib)%wstar(i,j,k)=dom(ib)%wstar(i,j,k)+dt*alfapr*forcn_z
	    endif
	  endif

              end do
           end do
        end do

        end do

        if (LENERGY) call mom_buo
        if (LROUGH) call rough_velocity

        return
        end subroutine rungek_diff2nd
!##########################################################################
