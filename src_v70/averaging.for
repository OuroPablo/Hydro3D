!##########################################################################
        subroutine update_mean
!##########################################################################
!
!     calculate averaged values of u, v, w, (p),
!                                  u'u', v'v', w'w',
!                                  u'v', u'w', v'w'
!
!     installed options (for constant grid spacing in aver. directions):
!
!     only time averaging
        use multidata
        use vars
        implicit none
!        double precision    :: facp1,facm1,facp2,facm2
        double precision :: ufuf,vfvf,wfwf,ufvf,ufwf,vfwf
        double precision :: ucf,vcf,wcf,pfpf,TfTf,SfSf,SfUf,SfVf,SfWf
        double precision :: dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz
	  double precision :: uf_ip1,uf_im1,uf_jp1,uf_jm1,uf_kp1,uf_km1
	  double precision :: vf_ip1,vf_im1,vf_jp1,vf_jm1,vf_kp1,vf_km1
	  double precision :: wf_ip1,wf_im1,wf_jp1,wf_jm1,wf_kp1,wf_km1
        integer :: i,j,k,ib

        do ib=1,nbp

!.....For first order moments
        if (ctime.ge.t_start_averaging1) then
          do k=1,dom(ib)%ttc_k
            do j=1,dom(ib)%ttc_j
              do i=1,dom(ib)%ttc_i
	          if (L_LSM) then
                  if (dom(ib)%phi(i,j,k).ge.0.0) then
		        dom(ib)%ntav1(i,j,k)=dom(ib)%ntav1(i,j,k)+1
                    dom(ib)%facp1(i,j,k)=1./dom(ib)%ntav1(i,j,k)
        		  dom(ib)%facm1(i,j,k)=1.-dom(ib)%facp1(i,j,k)
			end if
		    else
			dom(ib)%ntav1(i,j,k)=dom(ib)%ntav1(i,j,k)+1
                  dom(ib)%facp1(i,j,k)=1./dom(ib)%ntav1(i,j,k)
        		dom(ib)%facm1(i,j,k)=1.-dom(ib)%facp1(i,j,k)
		    end if
              end do
            end do
          end do
        endif

!.....For second order moments
        if (ctime.ge.t_start_averaging2) then
          do k=1,dom(ib)%ttc_k
            do j=1,dom(ib)%ttc_j
              do i=1,dom(ib)%ttc_i
	          if (L_LSM) then
			if (dom(ib)%phi(i,j,k).ge.0.0) then
		        dom(ib)%ntav2(i,j,k)=dom(ib)%ntav2(i,j,k)+1
                    dom(ib)%facp2(i,j,k)=1./dom(ib)%ntav2(i,j,k)
        		  dom(ib)%facm2(i,j,k)=1.-dom(ib)%facp2(i,j,k)
			end if
		    else
		      dom(ib)%ntav2(i,j,k)=dom(ib)%ntav2(i,j,k)+1
                  dom(ib)%facp2(i,j,k)=1./dom(ib)%ntav2(i,j,k)
        		dom(ib)%facm2(i,j,k)=1.-dom(ib)%facp2(i,j,k)
		    end if
              end do
            end do
          end do
        endif

!        do ib=1,nbp

        do k=1,dom(ib)%ttc_k
           do j=1,dom(ib)%ttc_j
              do i=1,dom(ib)%ttc_i
		    if (L_LSM) then
	            if (dom(ib)%phi(i,j,k).ge.0.0) then
                    dom(ib)%um(i,j,k)=dom(ib)%facm1(i,j,k)*
     & dom(ib)%um(i,j,k)+dom(ib)%facp1(i,j,k)*dom(ib)%u(i,j,k)
                    ufuf=(dom(ib)%u(i,j,k)-dom(ib)%um(i,j,k))*
     & (dom(ib)%u(i,j,k)-dom(ib)%um(i,j,k))
                    dom(ib)%uum(i,j,k)=dom(ib)%facm2(i,j,k)*
     & dom(ib)%uum(i,j,k)+dom(ib)%facp2(i,j,k) * ufuf
	            end if
		    else
                    dom(ib)%um(i,j,k)=dom(ib)%facm1(i,j,k)*
     & dom(ib)%um(i,j,k)+dom(ib)%facp1(i,j,k)*dom(ib)%u(i,j,k)
                    ufuf = (dom(ib)%u(i,j,k)-dom(ib)%um(i,j,k))*
     & (dom(ib)%u(i,j,k)-dom(ib)%um(i,j,k))
                    dom(ib)%uum(i,j,k)=dom(ib)%facm2(i,j,k)*
     & dom(ib)%uum(i,j,k)+dom(ib)%facp2(i,j,k)*ufuf
		    end if
              end do
           end do
        end do

        do k=1,dom(ib)%ttc_k
           do j=1,dom(ib)%ttc_j
              do i=1,dom(ib)%ttc_i
		    if (L_LSM) then
	            if (dom(ib)%phi(i,j,k).ge.0.0) then
                    dom(ib)%vm(i,j,k)=dom(ib)%facm1(i,j,k)*
     & dom(ib)%vm(i,j,k)+dom(ib)%facp1(i,j,k)*dom(ib)%v(i,j,k)
                    vfvf=(dom(ib)%v(i,j,k)-dom(ib)%vm(i,j,k))*
     & (dom(ib)%v(i,j,k)-dom(ib)%vm(i,j,k))
                    dom(ib)%vvm(i,j,k)=dom(ib)%facm2(i,j,k)*
     & dom(ib)%vvm(i,j,k)+dom(ib)%facp2(i,j,k)*vfvf
	            end if
		    else
                    dom(ib)%vm(i,j,k)=dom(ib)%facm1(i,j,k)*
     & dom(ib)%vm(i,j,k)+dom(ib)%facp1(i,j,k)*dom(ib)%v(i,j,k)
                    vfvf=(dom(ib)%v(i,j,k)-dom(ib)%vm(i,j,k))*
     & (dom(ib)%v(i,j,k)-dom(ib)%vm(i,j,k))
                    dom(ib)%vvm(i,j,k)=dom(ib)%facm2(i,j,k)*
     & dom(ib)%vvm(i,j,k)+dom(ib)%facp2(i,j,k)*vfvf
		    end if
              end do
           end do
        end do

        do k=1,dom(ib)%ttc_k
           do j=1,dom(ib)%ttc_j
              do i=1,dom(ib)%ttc_i
		    if (L_LSM) then
	            if (dom(ib)%phi(i,j,k).ge.0.0) then
                    dom(ib)%wm(i,j,k)=dom(ib)%facm1(i,j,k)*
     & dom(ib)%wm(i,j,k)+dom(ib)%facp1(i,j,k)*dom(ib)%w(i,j,k)
                    wfwf=(dom(ib)%w(i,j,k)-dom(ib)%wm(i,j,k))*
     & (dom(ib)%w(i,j,k)-dom(ib)%wm(i,j,k))
                    dom(ib)%wwm(i,j,k)=dom(ib)%facm2(i,j,k)*
     & dom(ib)%wwm(i,j,k)+dom(ib)%facp2(i,j,k)*wfwf
	            end if
		    else
                    dom(ib)%wm(i,j,k)=dom(ib)%facm1(i,j,k)*
     & dom(ib)%wm(i,j,k)+dom(ib)%facp1(i,j,k)*dom(ib)%w(i,j,k)
                    wfwf=(dom(ib)%w(i,j,k)-dom(ib)%wm(i,j,k))*
     & (dom(ib)%w(i,j,k)-dom(ib)%wm(i,j,k))
                    dom(ib)%wwm(i,j,k)=dom(ib)%facm2(i,j,k)*
     & dom(ib)%wwm(i,j,k)+dom(ib)%facp2(i,j,k)*wfwf
		    end if
              end do
           end do
        end do

        do k=1,dom(ib)%ttc_k
           do j=1,dom(ib)%ttc_j
              do i=1,dom(ib)%ttc_i
		    if (L_LSM) then
			dom(ib)%phim(i,j,k)=dom(ib)%facm1(i,j,k)*
     & dom(ib)%phim(i,j,k)+dom(ib)%facp1(i,j,k)*dom(ib)%phi(i,j,k)
	            if (dom(ib)%phi(i,j,k).ge.0.0) then
                  dom(ib)%pm(i,j,k)=dom(ib)%facm1(i,j,k)*
     & dom(ib)%pm(i,j,k)+dom(ib)%facp1(i,j,k)*dom(ib)%p(i,j,k)
                  pfpf=(dom(ib)%p(i,j,k)-dom(ib)%pm(i,j,k))*
     & (dom(ib)%p(i,j,k)-dom(ib)%pm(i,j,k))
                  dom(ib)%ppm(i,j,k)=dom(ib)%facm2(i,j,k)*
     & dom(ib)%ppm(i,j,k)+dom(ib)%facp2(i,j,k) * pfpf
                  dom(ib)%Tm(i,j,k)=dom(ib)%facm1(i,j,k)*
     & dom(ib)%Tm(i,j,k)+dom(ib)%facp1(i,j,k)*dom(ib)%T(i,j,k)
                  TfTf=(dom(ib)%T(i,j,k)-dom(ib)%Tm(i,j,k))*
     & (dom(ib)%T(i,j,k)-dom(ib)%Tm(i,j,k))
                  dom(ib)%Ttm(i,j,k)=dom(ib)%facm2(i,j,k)*
     & dom(ib)%Ttm(i,j,k)+dom(ib)%facp2(i,j,k)*TfTf
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
	IF(LSCALAR.EQ..TRUE.) THEN
                  dom(ib)%Sm(i,j,k)=dom(ib)%facm1(i,j,k)*
     & dom(ib)%Sm(i,j,k)+dom(ib)%facp1(i,j,k)*dom(ib)%S(i,j,k)
                  SfSf=(dom(ib)%S(i,j,k)-dom(ib)%Sm(i,j,k))*
     & (dom(ib)%S(i,j,k)-dom(ib)%Sm(i,j,k))
                  dom(ib)%Stm(i,j,k)=dom(ib)%facm2(i,j,k)*
     & dom(ib)%Stm(i,j,k)+dom(ib)%facp2(i,j,k)*SfSf
                  SfUf=(dom(ib)%S(i,j,k)-dom(ib)%Sm(i,j,k))*
     & (dom(ib)%U(i,j,k)-dom(ib)%Um(i,j,k))
                  dom(ib)%SUtm(i,j,k)=dom(ib)%facm2(i,j,k)*
     & dom(ib)%SUtm(i,j,k)+dom(ib)%facp2(i,j,k)*SfUf

                  SfVf=(dom(ib)%S(i,j,k)-dom(ib)%Sm(i,j,k))*
     & (dom(ib)%V(i,j,k)-dom(ib)%Vm(i,j,k))
                  dom(ib)%SVtm(i,j,k)=dom(ib)%facm2(i,j,k)*
     & dom(ib)%SVtm(i,j,k)+dom(ib)%facp2(i,j,k)*SfVf

                  SfWf=(dom(ib)%S(i,j,k)-dom(ib)%Sm(i,j,k))*
     & (dom(ib)%W(i,j,k)-dom(ib)%Wm(i,j,k))
                  dom(ib)%SWtm(i,j,k)=dom(ib)%facm2(i,j,k)*
     & dom(ib)%SWtm(i,j,k)+dom(ib)%facp2(i,j,k)*SfWf
	ENDIF
	            end if 
 		    else              
                  dom(ib)%pm(i,j,k)=dom(ib)%facm1(i,j,k)*
     & dom(ib)%pm(i,j,k)+dom(ib)%facp1(i,j,k)*dom(ib)%p(i,j,k)
                  pfpf=(dom(ib)%p(i,j,k)-dom(ib)%pm(i,j,k))*
     & (dom(ib)%p(i,j,k)-dom(ib)%pm(i,j,k))
                  dom(ib)%ppm(i,j,k)=dom(ib)%facm2(i,j,k)*
     & dom(ib)%ppm(i,j,k)+dom(ib)%facp2(i,j,k)*pfpf
                  dom(ib)%Tm(i,j,k)=dom(ib)%facm1(i,j,k)*
     & dom(ib)%Tm(i,j,k)+dom(ib)%facp1(i,j,k)*dom(ib)%T(i,j,k)
                  TfTf=(dom(ib)%T(i,j,k)-dom(ib)%Tm(i,j,k))*
     & (dom(ib)%T(i,j,k)-dom(ib)%Tm(i,j,k))
                  dom(ib)%Ttm(i,j,k)=dom(ib)%facm2(i,j,k)*
     & dom(ib)%Ttm(i,j,k)+dom(ib)%facp2(i,j,k)*TfTf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
	IF(LSCALAR.EQ..TRUE.) THEN
                  dom(ib)%Sm(i,j,k)=dom(ib)%facm1(i,j,k)*
     & dom(ib)%Sm(i,j,k)+dom(ib)%facp1(i,j,k)*dom(ib)%S(i,j,k)
                  SfSf=(dom(ib)%S(i,j,k)-dom(ib)%Sm(i,j,k))*
     & (dom(ib)%S(i,j,k)-dom(ib)%Sm(i,j,k))
                  dom(ib)%Stm(i,j,k)=dom(ib)%facm2(i,j,k)*
     & dom(ib)%Stm(i,j,k)+dom(ib)%facp2(i,j,k)*SfSf

                  SfUf=(dom(ib)%S(i,j,k)-dom(ib)%Sm(i,j,k))*
     & (dom(ib)%U(i,j,k)-dom(ib)%Um(i,j,k))
                  dom(ib)%SUtm(i,j,k)=dom(ib)%facm2(i,j,k)*
     & dom(ib)%SUtm(i,j,k)+dom(ib)%facp2(i,j,k)*SfUf

                  SfVf=(dom(ib)%S(i,j,k)-dom(ib)%Sm(i,j,k))*
     & (dom(ib)%V(i,j,k)-dom(ib)%Vm(i,j,k))
                  dom(ib)%SVtm(i,j,k)=dom(ib)%facm2(i,j,k)*
     & dom(ib)%SVtm(i,j,k)+dom(ib)%facp2(i,j,k)*SfVf

                  SfWf=(dom(ib)%S(i,j,k)-dom(ib)%Sm(i,j,k))*
     & (dom(ib)%W(i,j,k)-dom(ib)%Wm(i,j,k))
                  dom(ib)%SWtm(i,j,k)=dom(ib)%facm2(i,j,k)*
     & dom(ib)%SWtm(i,j,k)+dom(ib)%facp2(i,j,k)*SfWf
	ENDIF
	          end if
              end do
           end do
        end do

        do k=2,dom(ib)%ttc_k-1
           do j=2,dom(ib)%ttc_j-1
              do i=2,dom(ib)%ttc_i-1
           dom(ib)%vism(i,j,k)=dom(ib)%facm1(i,j,k)*dom(ib)%vism(i,j,k)+
     &  dom(ib)%facp1(i,j,k)*dom(ib)%vis(i,j,k)		
              end do
           end do
        end do

        do k=2,dom(ib)%ttc_k-1
           do j=2,dom(ib)%ttc_j-1
              do i=2,dom(ib)%ttc_i-1
	uf_ip1 = dom(ib)%um(i+1,j,k)-dom(ib)%u(i+1,j,k)
	uf_im1 = dom(ib)%um(i-1,j,k)-dom(ib)%u(i-1,j,k)
	uf_jp1 = dom(ib)%um(i,j+1,k)-dom(ib)%u(i,j+1,k)
	uf_jm1 = dom(ib)%um(i,j-1,k)-dom(ib)%u(i,j-1,k)
	uf_kp1 = dom(ib)%um(i,j,k+1)-dom(ib)%u(i,j,k+1)
	uf_km1 = dom(ib)%um(i,j,k-1)-dom(ib)%u(i,j,k-1)

	vf_ip1 = dom(ib)%vm(i+1,j,k)-dom(ib)%v(i+1,j,k)
	vf_im1 = dom(ib)%vm(i-1,j,k)-dom(ib)%v(i-1,j,k)
	vf_jp1 = dom(ib)%vm(i,j+1,k)-dom(ib)%v(i,j+1,k)
	vf_jm1 = dom(ib)%vm(i,j-1,k)-dom(ib)%v(i,j-1,k)
	vf_kp1 = dom(ib)%vm(i,j,k+1)-dom(ib)%v(i,j,k+1)
	vf_km1 = dom(ib)%vm(i,j,k-1)-dom(ib)%v(i,j,k-1)

	wf_ip1 = dom(ib)%wm(i+1,j,k)-dom(ib)%w(i+1,j,k)
	wf_im1 = dom(ib)%wm(i-1,j,k)-dom(ib)%w(i-1,j,k)
	wf_jp1 = dom(ib)%wm(i,j+1,k)-dom(ib)%w(i,j+1,k)
	wf_jm1 = dom(ib)%wm(i,j-1,k)-dom(ib)%w(i,j-1,k)
	wf_kp1 = dom(ib)%wm(i,j,k+1)-dom(ib)%w(i,j,k+1)
	wf_km1 = dom(ib)%wm(i,j,k-1)-dom(ib)%w(i,j,k-1)

		  dudx = (uf_ip1-uf_im1)/(2*dom(ib)%dx)
		  dudy = (uf_jp1-uf_jm1)/(2*dom(ib)%dy)
		  dudz = (uf_kp1-uf_km1)/(2*dom(ib)%dz)

		  dvdx = (vf_ip1-vf_im1)/(2*dom(ib)%dx)
		  dvdy = (vf_jp1-vf_jm1)/(2*dom(ib)%dy)
		  dvdz = (vf_kp1-vf_km1)/(2*dom(ib)%dz)

		  dwdx = (wf_ip1-wf_im1)/(2*dom(ib)%dx)
		  dwdy = (wf_jp1-wf_jm1)/(2*dom(ib)%dy)
		  dwdz = (wf_kp1-wf_km1)/(2*dom(ib)%dz)

        dom(ib)%eps(i,j,k) = dom(ib)%vis(i,j,k) *
     &  (2.d0*dudx**2+dvdx**2+dwdx**2
     &   +dudy**2+2.d0*dvdy**2+dwdy**2
     &   +dudz**2+dvdz**2+2.d0*dwdz**2
     &   +2.d0*dudy*dvdx+2.d0*dudz*dwdx+2.d0*dvdz*dwdy)

           dom(ib)%epsm(i,j,k)=dom(ib)%facm1(i,j,k)*dom(ib)%epsm(i,j,k)+
     &  dom(ib)%facp1(i,j,k)*dom(ib)%eps(i,j,k)					
              end do
           end do
        end do

        do k=2,dom(ib)%ttc_k 
           do j=2,dom(ib)%ttc_j 
              do i=2,dom(ib)%ttc_i
		    if (L_LSM) then 
	            if (dom(ib)%phi(i,j,k).ge.0.0) then
                    ucf=0.5*((dom(ib)%u(i-1,j,k)-dom(ib)%um(i-1,j,k))+
     & (dom(ib)%u(i,j,k)-dom(ib)%um(i,j,k)))
                    vcf=0.5*((dom(ib)%v(i,j-1,k)-dom(ib)%vm(i,j-1,k))+
     & (dom(ib)%v(i,j,k)-dom(ib)%vm(i,j,k)))
                    wcf=0.5*((dom(ib)%w(i,j,k-1)-dom(ib)%wm(i,j,k-1))+
     & (dom(ib)%w(i,j,k)-dom(ib)%wm(i,j,k)))
                    ufvf = ucf * vcf
                    ufwf = ucf * wcf
                    vfwf = vcf * wcf
                    dom(ib)%uvm(i,j,k)=dom(ib)%facm2(i,j,k)*
     & dom(ib)%uvm(i,j,k)+dom(ib)%facp2(i,j,k)*ufvf
                    dom(ib)%uwm(i,j,k)=dom(ib)%facm2(i,j,k)*
     & dom(ib)%uwm(i,j,k)+dom(ib)%facp2(i,j,k)*ufwf
                    dom(ib)%vwm(i,j,k)=dom(ib)%facm2(i,j,k)*
     & dom(ib)%vwm(i,j,k)+dom(ib)%facp2(i,j,k)*vfwf
	            end if 
	          else
                  ucf=0.5*((dom(ib)%u(i-1,j,k)-dom(ib)%um(i-1,j,k))+
     & (dom(ib)%u(i,j,k)-dom(ib)%um(i,j,k)))
                  vcf=0.5*((dom(ib)%v(i,j-1,k)-dom(ib)%vm(i,j-1,k))+
     & (dom(ib)%v(i,j,k)-dom(ib)%vm(i,j,k)))
                  wcf=0.5*((dom(ib)%w(i,j,k-1)-dom(ib)%wm(i,j,k-1))+
     & (dom(ib)%w(i,j,k)-dom(ib)%wm(i,j,k)))
                  ufvf = ucf * vcf
                  ufwf = ucf * wcf
                  vfwf = vcf * wcf
                  dom(ib)%uvm(i,j,k)=dom(ib)%facm2(i,j,k)*
     & dom(ib)%uvm(i,j,k)+dom(ib)%facp2(i,j,k)*ufvf
                  dom(ib)%uwm(i,j,k)=dom(ib)%facm2(i,j,k)*
     & dom(ib)%uwm(i,j,k)+dom(ib)%facp2(i,j,k)*ufwf
                  dom(ib)%vwm(i,j,k)=dom(ib)%facm2(i,j,k)*
     & dom(ib)%vwm(i,j,k)+dom(ib)%facp2(i,j,k)*vfwf
	          end if
              end do
           end do
        end do

        end do

        end subroutine update_mean
!##########################################################################
        subroutine add_noise(fnoise)
!##########################################################################
        use vars
        use multidata
        implicit none
        double precision    :: fnoise
        double precision    :: random_number_normal
        integer :: i,j,k,ib

!......add some gausian noise to flowfield

        call RANDOM_SEED

        do ib=1,nbp

        do k=1,dom(ib)%ttc_k
           do j=1,dom(ib)%ttc_j
              do i=1,dom(ib)%ttc_i
                 if (L_LSMbase) then
		      if (dom(ib)%zc(k).le.length) then
			 dom(ib)%u(i,j,k) = dom(ib)%u(i,j,k) +
     & random_number_normal(0.0,fnoise)
	 	      endif
		     else if (L_LSM) then
			if (dom(ib)%phi(i,j,k) .ge. 0.0) then
                  dom(ib)%u(i,j,k) = dom(ib)%u(i,j,k) +
     & random_number_normal(0.0,fnoise)
			endif
		     else
                  dom(ib)%u(i,j,k) = dom(ib)%u(i,j,k) +
     & random_number_normal(0.0,fnoise)
		     endif
              end do
           end do
        end do

        do k=1,dom(ib)%ttc_k
           do j=1,dom(ib)%ttc_j
              do i=1,dom(ib)%ttc_i
                 if (L_LSMbase) then
		      if (dom(ib)%zc(k).le.length) then
                   dom(ib)%v(i,j,k) = dom(ib)%v(i,j,k) +
     & random_number_normal(0.0,fnoise)
	 	      endif
		     else if (L_LSM) then
			if (dom(ib)%phi(i,j,k) .ge. 0.0) then
                    dom(ib)%v(i,j,k) = dom(ib)%v(i,j,k) +
     & random_number_normal(0.0,fnoise)
			endif
		     else
                  dom(ib)%v(i,j,k) = dom(ib)%v(i,j,k) +
     & random_number_normal(0.0,fnoise)
                 endif
              end do
           end do
        end do

        do k=1,dom(ib)%ttc_k
           do j=1,dom(ib)%ttc_j
              do i=1,dom(ib)%ttc_i
                 if (L_LSMbase) then
		      if (dom(ib)%zc(k).le.length) then
                   dom(ib)%w(i,j,k) = dom(ib)%w(i,j,k) +
     & random_number_normal(0.0,fnoise)
	 	      endif
		     else if (L_LSM) then
			if (dom(ib)%phi(i,j,k) .ge. 0.0) then
                    dom(ib)%w(i,j,k) = dom(ib)%w(i,j,k) +
     & random_number_normal(0.0,fnoise)
			endif
		     else
                  dom(ib)%w(i,j,k) = dom(ib)%w(i,j,k) +
     & random_number_normal(0.0,fnoise)
                 endif
              end do
           end do
        end do

        end do

        end subroutine add_noise
!##########################################################################
        function random_number_normal(mean,sigma) result( fn_val )
!##########################################################################
!       Generate random numbers
!       with a normal distribution with given mean and standard deviaton.
!
!       Generate a random normal deviate using the polar method.
!       Reference: Marsaglia,G. & Bray,T.A. 'A convenient method for generating
!                  normal variables', Siam Rev., vol.6, 260-264, 1964.
!       (source from internet)

        implicit none
        double precision       :: fn_val,mean,sigma
        double precision       :: ull, sumall
        double precision, save :: vll, sln
        logical, save   :: second = .false.
        double precision, parameter :: one = 1.0, vsmall = tiny( one )

        if (second) then

!...... If second, use the second random number generated on last call
           second = .false.
           fn_val = vll*sln

        else
!...... First call; generate a pair of random normals
           second = .true.
           do
              call random_number(ull)
              call random_number(vll)

              ull = scale( ull, 1 ) - one
              vll = scale( vll, 1 ) - one

!.........vsmall added to prevent LOG(zero) / zero
              sumall = ull*ull + vll*vll + vsmall
              if(sumall < one) exit
           end do

           sln = sqrt(- scale( log(sumall), 1 ) / sumall)
           fn_val = ull*sln
        end if

!.....set mean and standart deviation
        fn_val = fn_val * sigma + mean

        return

        end function random_number_normal
!##########################################################################
