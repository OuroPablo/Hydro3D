!##########################################################################
        subroutine alloc_dom
!##########################################################################
        use multidata
        use vars
        use mpi
        implicit none
        integer i,j,k,ib,dir
        integer sync_dir,cpu_next,cpu_prev,rdivmy,rdivng,my_cor
        integer ng_p,ng_n

        if(rdivmax.gt.1) then

!===============================================================
        do sync_dir = 1,3
        do ib=1,nbp

           if (sync_dir.eq.1)  then
              cpu_prev=dom(ib)%iprev
              cpu_next=dom(ib)%inext
           else if (sync_dir.eq.2)  then
              cpu_prev=dom(ib)%jprev
              cpu_next=dom(ib)%jnext
           else
              cpu_prev=dom(ib)%kprev
              cpu_next=dom(ib)%knext
           end if

           if (cpu_prev.ge.0) then
              rdivmy=rdiv(dom_id(ib))
              rdivng=rdiv(cpu_prev)

              if (rdivmy.ne.rdivng) then
                 if (rdivmy.gt.rdivng) then
                    if (rdivmy .ne. 2*rdivng) then
                       print*, '====ERROR1===>  wrong rdiv',dom_id(ib)
                       stop
                    end if
                 else if (rdivmy.lt.rdivng) then
                    if (2*rdivmy .ne. rdivng) then
                       print*, '====ERROR2===>  wrong rdiv',dom_id(ib)
                       stop
                    end if
                 else
                    print*, '===ERROR==='
                 end if
              end if
           end if

           if (cpu_next.ge.0) then
              rdivmy=rdiv(dom_id(ib))
              rdivng=rdiv(cpu_next)

              if (rdivmy.ne.rdivng) then
                 if (rdivmy.gt.rdivng) then
                    if (rdivmy .ne. 2*rdivng) then
                       print*, '====ERROR3===>  wrong rdiv',dom_id(ib)
                       stop
                    end if
                 else if (rdivmy.lt.rdivng) then
                    if (2*rdivmy .ne. rdivng) then
                       print*, '====ERROR4===>  wrong rdiv',dom_id(ib)
                       stop
                    end if
                 else
                    print*, '===ERROR==='
                 end if
              end if
           end if


           if (sync_dir.eq.1)  then

!======================================================================
!******* Corner-Neighbors
!======================================================================
           do dir=1,4

           select case (dir)
              Case (1)  
                 ng_p=dom(ib)%corprev1
                 ng_n=dom(ib)%cornext1
              Case (2)
                 ng_p=dom(ib)%corprev2
                 ng_n=dom(ib)%cornext2
              Case (3)  
                 ng_p=dom(ib)%corprev3
                 ng_n=dom(ib)%cornext3
              Case (4)
                 ng_p=dom(ib)%corprev4
                 ng_n=dom(ib)%cornext4
           end select

           if (ng_p.ge.0) then
              rdivmy=rdiv(dom_id(ib))
              rdivng=rdiv(ng_p)

              if (rdivmy.ne.rdivng) then
                 if (rdivmy.gt.rdivng) then
                    if (rdivmy .ne. 2*rdivng) then
        print*, '====ERRORcor-p===>  wrong rdiv',dom_id(ib),dir
                       stop
                    end if
                 else !if (rdivmy.lt.rdivng) then
                    if (2*rdivmy .ne. rdivng) then
        print*, '====ERRORcor-p===>  wrong rdiv',dom_id(ib),dir
                       stop
                    end if
                 end if
              end if
           end if
           if (ng_n.ge.0) then
              rdivmy=rdiv(dom_id(ib))
              rdivng=rdiv(ng_n)

              if (rdivmy.ne.rdivng) then
                 if (rdivmy.gt.rdivng) then
                    if (rdivmy .ne. 2*rdivng) then
        print*, '====ERRORcor-n===>  wrong rdiv',dom_id(ib),dir
                       stop
                    end if
                 else !if (rdivmy.lt.rdivng) then
                    if (2*rdivmy .ne. rdivng) then
        print*, '====ERRORcor-n===>  wrong rdiv',dom_id(ib),dir
                       stop
                    end if
                 end if
              end if
           end if

           end do
!======================================================================
!******* Edge-Neighbors
!======================================================================
           do dir=1,6

           select case (dir)
              Case (1)  
                 ng_p=dom(ib)%edgprev1
                 ng_n=dom(ib)%edgnext1
              Case (2)
                 ng_p=dom(ib)%edgprev2
                 ng_n=dom(ib)%edgnext2
              Case (3)  
                 ng_p=dom(ib)%edgprev3
                 ng_n=dom(ib)%edgnext3
              Case (4)
                 ng_p=dom(ib)%edgprev4
                 ng_n=dom(ib)%edgnext4
              Case (5)  
                 ng_p=dom(ib)%edgprev5
                 ng_n=dom(ib)%edgnext5
              Case (6)
                 ng_p=dom(ib)%edgprev6
                 ng_n=dom(ib)%edgnext6
           end select

           if (ng_p.ge.0) then
              rdivmy=rdiv(dom_id(ib))
              rdivng=rdiv(ng_p)

              if (rdivmy.ne.rdivng) then
                 if (rdivmy.gt.rdivng) then
                    if (rdivmy .ne. 2*rdivng) then
        print*, '====ERRORedg-p===>  wrong rdiv',dom_id(ib),dir
                       stop
                    end if
                 else !if (rdivmy.lt.rdivng) then
                    if (2*rdivmy .ne. rdivng) then
        print*, '====ERRORedg-p===>  wrong rdiv',dom_id(ib),dir
                       stop
                    end if
                 end if
              end if
           end if
           if (ng_n.ge.0) then
              rdivmy=rdiv(dom_id(ib))
              rdivng=rdiv(ng_n)

              if (rdivmy.ne.rdivng) then
                 if (rdivmy.gt.rdivng) then
                    if (rdivmy .ne. 2*rdivng) then
        print*, '====ERRORedg-n===>  wrong rdiv',dom_id(ib),dir
                       stop
                    end if
                 else !if (rdivmy.lt.rdivng) then
                    if (2*rdivmy .ne. rdivng) then
        print*, '====ERRORedg-n===>  wrong rdiv',dom_id(ib),dir
                       stop
                    end if
                 end if
              end if
           end if

           end do
!======================================================================
!======================================================================
           end if !if (sync_dir.eq.1)  then
        end do
        end do
!===============================================================

        do ib=1,nbp

        if(dom(ib)%bc_west.eq.5 .or. dom(ib)%bc_east.eq.5) then
!=== West ===> 
           if (dom(ib)%iprev.lt.0) then
              if (dom(ib)%inext.ge.0) then
                 my_cor=dom(ib)%per_ip
                 if(rdiv(dom_id(ib)).ne.rdiv(my_cor)) then
                    print*,'==ERROR==>wrong rdiv for perBCw',dom_id(ib)
                    stop
                 end if
              end if
           end if
!=== East ===> 
           if (dom(ib)%inext.lt.0) then
              if (dom(ib)%iprev.ge.0) then
                 my_cor=dom(ib)%per_in
                 if(rdiv(dom_id(ib)).ne.rdiv(my_cor)) then
                    print*,'==ERROR==>wrong rdiv for perBCe',dom_id(ib)
                    stop
                 end if
              end if
           end if
        end if

        if(dom(ib)%bc_south.eq.5 .or. dom(ib)%bc_north.eq.5) then
!=== South ===> 
           if (dom(ib)%jprev.lt.0) then
              if (dom(ib)%jnext.ge.0) then
                 my_cor=dom(ib)%per_jp  
                 if(rdiv(dom_id(ib)).ne.rdiv(my_cor)) then
                    print*,'==ERROR==>wrong rdiv for perBCs',dom_id(ib)
                    stop
                 end if
              end if
           end if
!=== North ===> 
           if (dom(ib)%jnext.lt.0) then
              if (dom(ib)%jprev.ge.0) then
                 my_cor=dom(ib)%per_jn
                 if(rdiv(dom_id(ib)).ne.rdiv(my_cor)) then
                    print*,'==ERROR==>wrong rdiv for perBCn',dom_id(ib)
                    stop
                 end if
              end if
           end if
        end if

        if(dom(ib)%bc_bottom.eq.5 .or. dom(ib)%bc_top.eq.5) then
!=== Bottom ===> 
           if (dom(ib)%kprev.lt.0) then
              if (dom(ib)%knext.ge.0) then
                 my_cor=dom(ib)%per_kp
                 if(rdiv(dom_id(ib)).ne.rdiv(my_cor)) then
                    print*,'==ERROR==>wrong rdiv for perBCb',dom_id(ib)
                    stop
                 end if
              end if
           end if
!=== Top ===> 
           if (dom(ib)%knext.lt.0) then
              if (dom(ib)%kprev.ge.0) then
                 my_cor=dom(ib)%per_kn
                 if(rdiv(dom_id(ib)).ne.rdiv(my_cor)) then
                    print*,'==ERROR==>wrong rdiv for perBCt',dom_id(ib)
                    stop
                 end if
              end if
           end if
        end if

        end do
!===============================================================
        end if

!        call datainfo

        return
        end 
!##########################################################################
        subroutine read_mdmap
!##########################################################################
        use multidata
        use vars
        use mpi
        implicit none
        integer i,j,k,ib
        integer nptemp,myranktemp,nbtemp,cpu_no
        integer,allocatable,dimension(:) :: domtemp,buf_domindid
        character*80 :: dummyline

	if (myrank.eq.0) then
        numfile=1001
        open (unit=numfile, file='output.dat')
        end if

        open (unit=12, file='mdmap.cin')

        read (12,*) num_domains !number of domains
        read (12,*) nptemp !number of processors

        if (nprocs .ne. nptemp) then
           print*, '=====ERROR====='
           print*, 'number of cpus do not match map file'
           stop
        end if

        if (num_domains .gt. 9999) then
           print*, '=====ERROR====='
           print*, 'number of domains are exceeding',
     & ' the limit in exchange subroutine'
           stop
        end if

        allocate (dom_id(num_domains),domtemp(num_domains))
        allocate (dom_indid(0:num_domains-1),dom_ad(0:num_domains-1))
        allocate (buf_domindid(0:num_domains-1))
       	allocate (imbinblk(num_domains))	!Pablo

        dom_ad=-1
        read (12,*) dummyline

        nbpmax=0
        do i=1,nprocs
           read(12,*) myranktemp,nbtemp,domtemp(1:nbtemp)
           nbpmax=max(nbpmax,nbtemp)

           do ib=1,nbtemp
              dom_ad(domtemp(ib))=myranktemp
           end do

           if(myranktemp.eq.myrank) then
              nbp=nbtemp !number of domains for this processor
              dom_id(1:nbp)=domtemp(1:nbp)

              do ib=0,num_domains-1
                 dom_indid(ib)=-1
              end do
              do ib=1,nbp
                 dom_indid(dom_id(ib))=ib
              end do

           end if
        end do

        read (12,*) dummyline
        close(12)


        buf_domindid = dom_indid

        call MPI_ALLREDUCE(buf_domindid,dom_indid,num_domains,
     &MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ierr)

        if(myrank.eq.0) then
           do i=0,num_domains-1
              if(dom_ad(i).eq.-1) then
                 print*,'unidentified domain in mdmap no:',i
                 print*,'ERROR! check mdmap.cin'
                 stop
              end if
              print*, 'dom:',i,' cpu:',dom_ad(i),' ib:',dom_indid(i) 
        write(numfile,*) 'dom:',i,' cpu:',dom_ad(i),' ib:',dom_indid(i)
           end do
        end if

        allocate (dom(nbp))

        call MPI_BARRIER (MPI_COMM_WORLD,ierr)

        return
        end subroutine read_mdmap
!##########################################################################
        subroutine read_infodom
!##########################################################################
        use multidata
        use vars
        use mpi
        implicit none

        integer :: i,j,k,ib,ndoms,say,ndo
        integer :: sync_dir,cpu_next,cpu_prev
 !       double precision::xs2,xe2,ys2,ye2,zs2,ze2
        character*80 :: dummyline

        open (unit=12, file='infodom.cin')

        read (12,*) ndoms !number of domains
        read (12,*) dummyline

        if (ndoms .ne. num_domains) then
           print*, '=====ERROR====='
           print*, 'number of domains does not match infodom file'
           stop
        end if

        allocate (xcor(0:ndoms-1,2),ycor(0:ndoms-1,2),zcor(0:ndoms-1,2))
        allocate (rdiv(0:ndoms-1))

        rdivmax=1
        do i=0,ndoms-1
           read(12,*) ndo,rdiv(ndo),xcor(ndo,1),xcor(ndo,2),ycor(ndo,1),
     & ycor(ndo,2),zcor(ndo,1),zcor(ndo,2)
           if(rdiv(ndo).gt.rdivmax) rdivmax=rdiv(ndo)
        end do

        read (12,*) dummyline								
	  read (12,*) idom
	  read (12,*) jdom
	  read (12,*) kdom

        close(12)

        do ib=1,nbp
           dom(ib)%xsl=xcor(dom_id(ib),1)
           dom(ib)%xel=xcor(dom_id(ib),2)

           dom(ib)%ysl=ycor(dom_id(ib),1)
           dom(ib)%yel=ycor(dom_id(ib),2)

           dom(ib)%zsl=zcor(dom_id(ib),1)
           dom(ib)%zel=zcor(dom_id(ib),2)
        end do

        xst=xcor(0,1); xen=xcor(0,2)
        yst=ycor(0,1); yen=ycor(0,2)
        zst=zcor(0,1); zen=zcor(0,2)

        do i=1,ndoms-1
           xst=min(xst,xcor(i,1))
           xen=max(xen,xcor(i,2))

           yst=min(yst,ycor(i,1))
           yen=max(yen,ycor(i,2))

           zst=min(zst,zcor(i,1))
           zen=max(zen,zcor(i,2))
        end do

        prefdom=-1
        do i=0,ndoms-1
           if(xcor(i,1).eq.xst .and. ycor(i,1).eq.yst
     & .and. zcor(i,1).eq.zst) prefdom=i
        end do

        if(myrank.eq.0) print*,'pressure reference domain is ',prefdom
 
        if(prefdom.lt.0) then
!           print*,'ERROR about prefdom!!!!!!' 
           prefdom=0
!           stop
        end if

        do ib=1,nbp

           dom(ib)%iprev= -1
           dom(ib)%inext= -1
           dom(ib)%jprev= -1
           dom(ib)%jnext= -1
           dom(ib)%kprev= -1
           dom(ib)%knext= -1

!================== I-PREVIOUS NEIGHBOR ===================================
           say=0; 
           do i=0,ndoms-1
              if(i.ne.dom_id(ib)) then
                 if(abs(dom(ib)%xsl-xcor(i,2)).lt.1E-08) then

                    if(abs(dom(ib)%ysl-ycor(i,1)).lt.1E-08
     & .and. abs(dom(ib)%yel-ycor(i,2)).lt.1E-08
     & .and. abs(dom(ib)%zsl-zcor(i,1)).lt.1E-08
     & .and. abs(dom(ib)%zel-zcor(i,2)).lt.1E-08) then
                       say=say+1
                       dom(ib)%iprev=i
                    end if 

                 end if
              end if
           end do

           if(say.gt.1) then
              print*,'dom#:',dom_id(ib),' has more than one',
     & ' previous neighbor in x direction'
              stop
           end if

!================== I-NEXT NEIGHBOR ===================================
           say=0; 
           do i=0,ndoms-1
              if(i.ne.dom_id(ib)) then
                 if(abs(dom(ib)%xel-xcor(i,1)).lt.1E-08) then

                    if(abs(dom(ib)%ysl-ycor(i,1)).lt.1E-08
     & .and. abs(dom(ib)%yel-ycor(i,2)).lt.1E-08
     & .and. abs(dom(ib)%zsl-zcor(i,1)).lt.1E-08
     & .and. abs(dom(ib)%zel-zcor(i,2)).lt.1E-08) then
                       say=say+1
                       dom(ib)%inext=i
                    end if 

                 end if
              end if
           end do

           if(say.gt.1) then
              print*,'dom#:',dom_id(ib),' has more than one',
     & ' next neighbor in x direction'
              stop
           end if


!================== J-PREVIOUS NEIGHBOR ===================================
           say=0; 
           do i=0,ndoms-1
              if(i.ne.dom_id(ib)) then
                 if(abs(dom(ib)%ysl-ycor(i,2)).lt.1E-08) then

                    if(abs(dom(ib)%xsl-xcor(i,1)).lt.1E-08
     & .and. abs(dom(ib)%xel-xcor(i,2)).lt.1E-08
     & .and. abs(dom(ib)%zsl-zcor(i,1)).lt.1E-08
     & .and. abs(dom(ib)%zel-zcor(i,2)).lt.1E-08) then
                       say=say+1
                       dom(ib)%jprev=i
                    end if 

                 end if
              end if
           end do

           if(say.gt.1) then
              print*,'dom#:',dom_id(ib),' has more than one',
     & ' previous neighbor in y direction'
              stop
           end if

!================== J-NEXT NEIGHBOR ===================================
           say=0; 
           do i=0,ndoms-1
              if(i.ne.dom_id(ib)) then
                 if(abs(dom(ib)%yel-ycor(i,1)).lt.1E-08) then

                    if(abs(dom(ib)%xsl-xcor(i,1)).lt.1E-08
     & .and. abs(dom(ib)%xel-xcor(i,2)).lt.1E-08
     & .and. abs(dom(ib)%zsl-zcor(i,1)).lt.1E-08
     & .and. abs(dom(ib)%zel-zcor(i,2)).lt.1E-08) then
                       say=say+1
                       dom(ib)%jnext=i
                    end if 

                 end if
              end if
           end do

           if(say.gt.1) then
              print*,'dom#:',dom_id(ib),' has more than one',
     & ' next neighbor in y direction'
              stop
           end if


!================== K-PREVIOUS NEIGHBOR ===================================
           say=0; 
           do i=0,ndoms-1
              if(i.ne.dom_id(ib)) then
                 if(abs(dom(ib)%zsl-zcor(i,2)).lt.1E-08) then

                    if(abs(dom(ib)%xsl-xcor(i,1)).lt.1E-08
     & .and. abs(dom(ib)%xel-xcor(i,2)).lt.1E-08
     & .and. abs(dom(ib)%ysl-ycor(i,1)).lt.1E-08
     & .and. abs(dom(ib)%yel-ycor(i,2)).lt.1E-08) then
                       say=say+1
                       dom(ib)%kprev=i
                    end if 

                 end if
              end if
           end do

           if(say.gt.1) then
              print*,'dom#:',dom_id(ib),' has more than one',
     & ' previous neighbor in z direction'
              stop
           end if

!================== K-NEXT NEIGHBOR ===================================
           say=0; 
           do i=0,ndoms-1
              if(i.ne.dom_id(ib)) then
                 if(abs(dom(ib)%zel-zcor(i,1)).lt.1E-08) then

                    if(abs(dom(ib)%xsl-xcor(i,1)).lt.1E-08
     & .and. abs(dom(ib)%xel-xcor(i,2)).lt.1E-08
     & .and. abs(dom(ib)%ysl-ycor(i,1)).lt.1E-08
     & .and. abs(dom(ib)%yel-ycor(i,2)).lt.1E-08) then
                       say=say+1
                       dom(ib)%knext=i
                    end if 

                 end if
              end if
           end do

           if(say.gt.1) then
              print*,'dom#:',dom_id(ib),' has more than one',
     & ' next neighbor in z direction'
              stop
           end if

!==========================================================================

!        write (6,5501) dom_id(ib),dom(ib)%iprev,dom(ib)%inext,
!     & dom(ib)%jprev,dom(ib)%jnext,dom(ib)%kprev,dom(ib)%knext

        end do

!**************************************************************************
!============= SET CORNER NEIGHBORS
!**************************************************************************
        do ib=1,nbp

           dom(ib)%corprev1= -1
           dom(ib)%corprev2= -1
           dom(ib)%corprev3= -1
           dom(ib)%corprev4= -1
           dom(ib)%cornext1= -1
           dom(ib)%cornext2= -1
           dom(ib)%cornext3= -1
           dom(ib)%cornext4= -1

!==================CORNER Previous #1 NEIGHBOR ==================
           say=0; 
           do i=0,ndoms-1
              if(i.ne.dom_id(ib)) then
                 if(abs(dom(ib)%xsl-xcor(i,2)).lt.1E-08 .and.
     &              abs(dom(ib)%ysl-ycor(i,2)).lt.1E-08 .and.
     &              abs(dom(ib)%zsl-zcor(i,2)).lt.1E-08) then
                    dom(ib)%corprev1=i; say=say+1
                 end if 
              end if
           end do
           if(say.gt.1) then
              print*,dom_id(ib),' has more than one cr-ng'; stop
           end if
!==================CORNER Previous #2 NEIGHBOR ==================
           say=0; 
           do i=0,ndoms-1
              if(i.ne.dom_id(ib)) then
                 if(abs(dom(ib)%xsl-xcor(i,2)).lt.1E-08 .and.
     &              abs(dom(ib)%yel-ycor(i,1)).lt.1E-08 .and.
     &              abs(dom(ib)%zsl-zcor(i,2)).lt.1E-08) then
                    dom(ib)%corprev2=i; say=say+1
                 end if 
              end if
           end do
           if(say.gt.1) then
              print*,dom_id(ib),' has more than one cr-ng'; stop
           end if
!==================CORNER Previous #3 NEIGHBOR ==================
           say=0; 
           do i=0,ndoms-1
              if(i.ne.dom_id(ib)) then
                 if(abs(dom(ib)%xel-xcor(i,1)).lt.1E-08 .and.
     &              abs(dom(ib)%yel-ycor(i,1)).lt.1E-08 .and.
     &              abs(dom(ib)%zsl-zcor(i,2)).lt.1E-08) then
                    dom(ib)%corprev3=i; say=say+1
                 end if 
              end if
           end do
           if(say.gt.1) then
              print*,dom_id(ib),' has more than one cr-ng'; stop
           end if
!==================CORNER Previous #4 NEIGHBOR ==================
           say=0; 
           do i=0,ndoms-1
              if(i.ne.dom_id(ib)) then
                 if(abs(dom(ib)%xel-xcor(i,1)).lt.1E-08 .and.
     &              abs(dom(ib)%ysl-ycor(i,2)).lt.1E-08 .and.
     &              abs(dom(ib)%zsl-zcor(i,2)).lt.1E-08) then
                    dom(ib)%corprev4=i; say=say+1
                 end if 
              end if
           end do
           if(say.gt.1) then
              print*,dom_id(ib),' has more than one cr-ng'; stop
           end if
!==================CORNER Next #1 NEIGHBOR ==================
           say=0; 
           do i=0,ndoms-1
              if(i.ne.dom_id(ib)) then
                 if(abs(dom(ib)%xel-xcor(i,1)).lt.1E-08 .and.
     &              abs(dom(ib)%yel-ycor(i,1)).lt.1E-08 .and.
     &              abs(dom(ib)%zel-zcor(i,1)).lt.1E-08) then
                    dom(ib)%cornext1=i; say=say+1
                 end if 
              end if
           end do
           if(say.gt.1) then
              print*,dom_id(ib),' has more than one cr-ng'; stop
           end if
!==================CORNER Next #2 NEIGHBOR ==================
           say=0; 
           do i=0,ndoms-1
              if(i.ne.dom_id(ib)) then
                 if(abs(dom(ib)%xel-xcor(i,1)).lt.1E-08 .and.
     &              abs(dom(ib)%ysl-ycor(i,2)).lt.1E-08 .and.
     &              abs(dom(ib)%zel-zcor(i,1)).lt.1E-08) then
                    dom(ib)%cornext2=i; say=say+1
                 end if 
              end if
           end do
           if(say.gt.1) then
              print*,dom_id(ib),' has more than one cr-ng'; stop
           end if
!==================CORNER Next #3 NEIGHBOR ==================
           say=0; 
           do i=0,ndoms-1
              if(i.ne.dom_id(ib)) then
                 if(abs(dom(ib)%xsl-xcor(i,2)).lt.1E-08 .and.
     &              abs(dom(ib)%ysl-ycor(i,2)).lt.1E-08 .and.
     &              abs(dom(ib)%zel-zcor(i,1)).lt.1E-08) then
                    dom(ib)%cornext3=i; say=say+1
                 end if 
              end if
           end do
           if(say.gt.1) then
              print*,dom_id(ib),' has more than one cr-ng'; stop
           end if
!==================CORNER Next #4 NEIGHBOR ==================
           say=0; 
           do i=0,ndoms-1
              if(i.ne.dom_id(ib)) then
                 if(abs(dom(ib)%xsl-xcor(i,2)).lt.1E-08 .and.
     &              abs(dom(ib)%yel-ycor(i,1)).lt.1E-08 .and.
     &              abs(dom(ib)%zel-zcor(i,1)).lt.1E-08) then
                    dom(ib)%cornext4=i; say=say+1
                 end if 
              end if
           end do
           if(say.gt.1) then
              print*,dom_id(ib),' has more than one cr-ng'; stop
           end if


!        write (6,5601) dom_id(ib),dom(ib)%corprev1,dom(ib)%corprev2,
!     & dom(ib)%corprev3,dom(ib)%corprev4,dom(ib)%cornext1,
!     & dom(ib)%cornext2,dom(ib)%cornext3,dom(ib)%cornext4
        end do

!**************************************************************************
!============= SET EDGE NEIGHBORS
!**************************************************************************
        do ib=1,nbp

           dom(ib)%edgprev1= -1
           dom(ib)%edgprev2= -1
           dom(ib)%edgprev3= -1
           dom(ib)%edgprev4= -1
           dom(ib)%edgprev5= -1
           dom(ib)%edgprev6= -1
           dom(ib)%edgnext1= -1
           dom(ib)%edgnext2= -1
           dom(ib)%edgnext3= -1
           dom(ib)%edgnext4= -1
           dom(ib)%edgnext5= -1
           dom(ib)%edgnext6= -1

!==================EDGE Previous #1 NEIGHBOR ==================
           say=0; 
           do i=0,ndoms-1
              if(i.ne.dom_id(ib)) then
                 if(abs(dom(ib)%xsl-xcor(i,2)).lt.1E-08 .and.
     &              abs(dom(ib)%zsl-zcor(i,2)).lt.1E-08) then
                 if(abs(dom(ib)%ysl-ycor(i,1)).lt.1E-08
     & .and.        abs(dom(ib)%yel-ycor(i,2)).lt.1E-08) then
                       dom(ib)%edgprev1=i; say=say+1
                    end if 
                 end if 
              end if
           end do
           if(say.gt.1) then
              print*,dom_id(ib),' has more than one edg-ng'; stop
           end if
!==================EDGE Previous #2 NEIGHBOR ==================
           say=0; 
           do i=0,ndoms-1
              if(i.ne.dom_id(ib)) then
                 if(abs(dom(ib)%yel-ycor(i,1)).lt.1E-08 .and.
     &              abs(dom(ib)%zsl-zcor(i,2)).lt.1E-08) then
                 if(abs(dom(ib)%xsl-xcor(i,1)).lt.1E-08
     & .and.        abs(dom(ib)%xel-xcor(i,2)).lt.1E-08) then
                       dom(ib)%edgprev2=i; say=say+1
                    end if 
                 end if 
              end if
           end do
           if(say.gt.1) then
              print*,dom_id(ib),' has more than one edg-ng'; stop
           end if
!==================EDGE Previous #3 NEIGHBOR ==================
           say=0; 
           do i=0,ndoms-1
              if(i.ne.dom_id(ib)) then
                 if(abs(dom(ib)%xel-xcor(i,1)).lt.1E-08 .and.
     &              abs(dom(ib)%zsl-zcor(i,2)).lt.1E-08) then
                 if(abs(dom(ib)%ysl-ycor(i,1)).lt.1E-08
     & .and.        abs(dom(ib)%yel-ycor(i,2)).lt.1E-08) then
                       dom(ib)%edgprev3=i; say=say+1
                    end if 
                 end if 
              end if
           end do
           if(say.gt.1) then
              print*,dom_id(ib),' has more than one edg-ng'; stop
           end if
!==================EDGE Previous #4 NEIGHBOR ==================
           say=0; 
           do i=0,ndoms-1
              if(i.ne.dom_id(ib)) then
                 if(abs(dom(ib)%ysl-ycor(i,2)).lt.1E-08 .and.
     &              abs(dom(ib)%zsl-zcor(i,2)).lt.1E-08) then
                 if(abs(dom(ib)%xsl-xcor(i,1)).lt.1E-08
     & .and.        abs(dom(ib)%xel-xcor(i,2)).lt.1E-08) then
                       dom(ib)%edgprev4=i; say=say+1
                    end if 
                 end if 
              end if
           end do
           if(say.gt.1) then
              print*,dom_id(ib),' has more than one edg-ng'; stop
           end if
!==================EDGE Previous #5 NEIGHBOR ==================
           say=0; 
           do i=0,ndoms-1
              if(i.ne.dom_id(ib)) then
                 if(abs(dom(ib)%xsl-xcor(i,2)).lt.1E-08 .and.
     &              abs(dom(ib)%ysl-ycor(i,2)).lt.1E-08) then
                 if(abs(dom(ib)%zsl-zcor(i,1)).lt.1E-08
     & .and.        abs(dom(ib)%zel-zcor(i,2)).lt.1E-08) then
                       dom(ib)%edgprev5=i; say=say+1
                    end if 
                 end if 
              end if
           end do
           if(say.gt.1) then
              print*,dom_id(ib),' has more than one edg-ng'; stop
           end if
!==================EDGE Previous #6 NEIGHBOR ==================
           say=0; 
           do i=0,ndoms-1
              if(i.ne.dom_id(ib)) then
                 if(abs(dom(ib)%xsl-xcor(i,2)).lt.1E-08 .and.
     &              abs(dom(ib)%yel-ycor(i,1)).lt.1E-08) then
                 if(abs(dom(ib)%zsl-zcor(i,1)).lt.1E-08
     & .and.        abs(dom(ib)%zel-zcor(i,2)).lt.1E-08) then
                       dom(ib)%edgprev6=i; say=say+1
                    end if 
                 end if 
              end if
           end do
           if(say.gt.1) then
              print*,dom_id(ib),' has more than one edg-ng'; stop
           end if
!==================EDGE Next #1 NEIGHBOR ==================
           say=0; 
           do i=0,ndoms-1
              if(i.ne.dom_id(ib)) then
                 if(abs(dom(ib)%xel-xcor(i,1)).lt.1E-08 .and.
     &              abs(dom(ib)%zel-zcor(i,1)).lt.1E-08) then
                 if(abs(dom(ib)%ysl-ycor(i,1)).lt.1E-08
     & .and.        abs(dom(ib)%yel-ycor(i,2)).lt.1E-08) then
                       dom(ib)%edgnext1=i; say=say+1
                    end if 
                 end if 
              end if
           end do
           if(say.gt.1) then
              print*,dom_id(ib),' has more than one edg-ng'; stop
           end if
!==================EDGE Next #2 NEIGHBOR ==================
           say=0; 
           do i=0,ndoms-1
              if(i.ne.dom_id(ib)) then
                 if(abs(dom(ib)%ysl-ycor(i,2)).lt.1E-08 .and.
     &              abs(dom(ib)%zel-zcor(i,1)).lt.1E-08) then
                 if(abs(dom(ib)%xsl-xcor(i,1)).lt.1E-08
     & .and.        abs(dom(ib)%xel-xcor(i,2)).lt.1E-08) then
                       dom(ib)%edgnext2=i; say=say+1
                    end if 
                 end if 
              end if
           end do
           if(say.gt.1) then
              print*,dom_id(ib),' has more than one edg-ng'; stop
           end if
!==================EDGE Next #3 NEIGHBOR ==================
           say=0; 
           do i=0,ndoms-1
              if(i.ne.dom_id(ib)) then
                 if(abs(dom(ib)%xsl-xcor(i,2)).lt.1E-08 .and.
     &              abs(dom(ib)%zel-zcor(i,1)).lt.1E-08) then
                 if(abs(dom(ib)%ysl-ycor(i,1)).lt.1E-08
     & .and.        abs(dom(ib)%yel-ycor(i,2)).lt.1E-08) then
                       dom(ib)%edgnext3=i; say=say+1
                    end if 
                 end if 
              end if
           end do
           if(say.gt.1) then
              print*,dom_id(ib),' has more than one edg-ng'; stop
           end if
!==================EDGE Next #4 NEIGHBOR ==================
           say=0; 
           do i=0,ndoms-1
              if(i.ne.dom_id(ib)) then
                 if(abs(dom(ib)%yel-ycor(i,1)).lt.1E-08 .and.
     &              abs(dom(ib)%zel-zcor(i,1)).lt.1E-08) then
                 if(abs(dom(ib)%xsl-xcor(i,1)).lt.1E-08
     & .and.        abs(dom(ib)%xel-xcor(i,2)).lt.1E-08) then
                       dom(ib)%edgnext4=i; say=say+1
                    end if 
                 end if 
              end if
           end do
           if(say.gt.1) then
              print*,dom_id(ib),' has more than one edg-ng'; stop
           end if
!==================EDGE Next #5 NEIGHBOR ==================
           say=0; 
           do i=0,ndoms-1
              if(i.ne.dom_id(ib)) then
                 if(abs(dom(ib)%xel-xcor(i,1)).lt.1E-08 .and.
     &              abs(dom(ib)%yel-ycor(i,1)).lt.1E-08) then
                 if(abs(dom(ib)%zsl-zcor(i,1)).lt.1E-08
     & .and.        abs(dom(ib)%zel-zcor(i,2)).lt.1E-08) then
                       dom(ib)%edgnext5=i; say=say+1
                    end if 
                 end if 
              end if
           end do
           if(say.gt.1) then
              print*,dom_id(ib),' has more than one edg-ng'; stop
           end if
!==================EDGE Next #6 NEIGHBOR ==================
           say=0; 
           do i=0,ndoms-1
              if(i.ne.dom_id(ib)) then
                 if(abs(dom(ib)%xel-xcor(i,1)).lt.1E-08 .and.
     &              abs(dom(ib)%ysl-ycor(i,2)).lt.1E-08) then
                 if(abs(dom(ib)%zsl-zcor(i,1)).lt.1E-08
     & .and.        abs(dom(ib)%zel-zcor(i,2)).lt.1E-08) then
                       dom(ib)%edgnext6=i; say=say+1
                    end if 
                 end if 
              end if
           end do
           if(say.gt.1) then
              print*,dom_id(ib),' has more than one edg-ng'; stop
           end if

!        write (6,5602) dom_id(ib),dom(ib)%edgprev1,dom(ib)%edgprev2,
!     & dom(ib)%edgprev3,dom(ib)%edgprev4,dom(ib)%edgprev5,
!     & dom(ib)%edgprev6,dom(ib)%edgnext1,dom(ib)%edgnext2,
!     & dom(ib)%edgnext3,dom(ib)%edgnext4,dom(ib)%edgnext5,
!     & dom(ib)%edgnext6

        end do
!**************************************************************************
        call MPI_BARRIER (MPI_COMM_WORLD,ierr)
!**************************************************************************
        do ib=1,nbp

           dom(ib)%per_ip= -1
           dom(ib)%per_in= -1
           dom(ib)%per_jp= -1
           dom(ib)%per_jn= -1
           dom(ib)%per_kp= -1
           dom(ib)%per_kn= -1

!================== I-PREVIOUS NEIGHBOR ===================================
        if(dom(ib)%bc_west.eq.5 .and. dom(ib)%iprev.lt.0) then

           say=0; 
           do i=0,ndoms-1
              if(i.ne.dom_id(ib)) then
                 if(abs(xen-xcor(i,2)).lt.1E-08) then

                    if(abs(dom(ib)%ysl-ycor(i,1)).lt.1E-08
     & .and. abs(dom(ib)%yel-ycor(i,2)).lt.1E-08
     & .and. abs(dom(ib)%zsl-zcor(i,1)).lt.1E-08
     & .and. abs(dom(ib)%zel-zcor(i,2)).lt.1E-08) then
                       say=say+1
                       dom(ib)%per_ip=i
                    end if 

                 end if
              end if
           end do

           if(say.gt.1) then
              print*,'dom#:',dom_id(ib),' has more than one',
     & ' periodic previous neighbor in x direction'
              stop
           end if

!           write (6,5502) dom_id(ib),dom(ib)%per_ip
        end if
!================== I-NEXT NEIGHBOR ===================================
        if(dom(ib)%bc_east.eq.5 .and. dom(ib)%inext.lt.0) then

           say=0; 
           do i=0,ndoms-1
              if(i.ne.dom_id(ib)) then
                 if(abs(xst-xcor(i,1)).lt.1E-08) then

                    if(abs(dom(ib)%ysl-ycor(i,1)).lt.1E-08
     & .and. abs(dom(ib)%yel-ycor(i,2)).lt.1E-08
     & .and. abs(dom(ib)%zsl-zcor(i,1)).lt.1E-08
     & .and. abs(dom(ib)%zel-zcor(i,2)).lt.1E-08) then
                       say=say+1
                       dom(ib)%per_in=i
                    end if 

                 end if
              end if
           end do

           if(say.gt.1) then
              print*,'dom#:',dom_id(ib),' has more than one',
     & ' periodic next neighbor in x direction'
              stop
           end if

!           write (6,5503) dom_id(ib),dom(ib)%per_in
        end if
!================== J-PREVIOUS NEIGHBOR ===================================
        if(dom(ib)%bc_south.eq.5 .and. dom(ib)%jprev.lt.0) then

           say=0; 
           do i=0,ndoms-1
              if(i.ne.dom_id(ib)) then
                 if(abs(yen-ycor(i,2)).lt.1E-08) then

                    if(abs(dom(ib)%xsl-xcor(i,1)).lt.1E-08
     & .and. abs(dom(ib)%xel-xcor(i,2)).lt.1E-08
     & .and. abs(dom(ib)%zsl-zcor(i,1)).lt.1E-08
     & .and. abs(dom(ib)%zel-zcor(i,2)).lt.1E-08) then
                       say=say+1
                       dom(ib)%per_jp=i
                    end if 

                 end if
              end if
           end do

           if(say.gt.1) then
              print*,'dom#:',dom_id(ib),' has more than one',
     & ' periodic previous neighbor in y direction'
              stop
           end if

!           write (6,5504) dom_id(ib),dom(ib)%per_jp
        end if
!================== J-NEXT NEIGHBOR ===================================
        if(dom(ib)%bc_north.eq.5 .and. dom(ib)%jnext.lt.0) then

           say=0; 
           do i=0,ndoms-1
              if(i.ne.dom_id(ib)) then
                 if(abs(yst-ycor(i,1)).lt.1E-08) then

                    if(abs(dom(ib)%xsl-xcor(i,1)).lt.1E-08
     & .and. abs(dom(ib)%xel-xcor(i,2)).lt.1E-08
     & .and. abs(dom(ib)%zsl-zcor(i,1)).lt.1E-08
     & .and. abs(dom(ib)%zel-zcor(i,2)).lt.1E-08) then
                       say=say+1
                       dom(ib)%per_jn=i
                    end if 

                 end if
              end if
           end do

           if(say.gt.1) then
              print*,'dom#:',dom_id(ib),' has more than one',
     & ' periodic next neighbor in y direction'
              stop
           end if

!           write (6,5505) dom_id(ib),dom(ib)%per_jn
        end if
!================== K-PREVIOUS NEIGHBOR ===================================
        if(dom(ib)%bc_bottom.eq.5 .and. dom(ib)%kprev.lt.0) then

           say=0; 
           do i=0,ndoms-1
              if(i.ne.dom_id(ib)) then
                 if(abs(zen-zcor(i,2)).lt.1E-08) then

                    if(abs(dom(ib)%xsl-xcor(i,1)).lt.1E-08
     & .and. abs(dom(ib)%xel-xcor(i,2)).lt.1E-08
     & .and. abs(dom(ib)%ysl-ycor(i,1)).lt.1E-08
     & .and. abs(dom(ib)%yel-ycor(i,2)).lt.1E-08) then
                       say=say+1
                       dom(ib)%per_kp=i
                    end if 

                 end if
              end if
           end do

           if(say.gt.1) then
              print*,'dom#:',dom_id(ib),' has more than one',
     & ' periodic previous neighbor in z direction'
              stop
           end if

!           write (6,5506) dom_id(ib),dom(ib)%per_kp
        end if
!================== K-NEXT NEIGHBOR ===================================
        if(dom(ib)%bc_top.eq.5 .and. dom(ib)%knext.lt.0) then

           say=0; 
           do i=0,ndoms-1
              if(i.ne.dom_id(ib)) then
                 if(abs(zst-zcor(i,1)).lt.1E-08) then

                    if(abs(dom(ib)%xsl-xcor(i,1)).lt.1E-08
     & .and. abs(dom(ib)%xel-xcor(i,2)).lt.1E-08
     & .and. abs(dom(ib)%ysl-ycor(i,1)).lt.1E-08
     & .and. abs(dom(ib)%yel-ycor(i,2)).lt.1E-08) then
                       say=say+1
                       dom(ib)%per_kn=i
                    end if 

                 end if
              end if
           end do

           if(say.gt.1) then
              print*,'dom#:',dom_id(ib),' has more than one',
     & ' periodic next neighbor in z direction'
              stop
           end if

!           write (6,5507) dom_id(ib),dom(ib)%per_kn
        end if
!==========================================================================
        end do


!**************************************************************************
!**************************************************************************

           PERIODIC=.false.
           do ib=1,nbp
              if(dom(ib)%bc_west.eq.5 .or. dom(ib)%bc_east.eq.5 .or.
     & dom(ib)%bc_south.eq.5 .or. dom(ib)%bc_north.eq.5 .or.
     & dom(ib)%bc_bottom.eq.5 .or. dom(ib)%bc_top.eq.5) then
                 PERIODIC=.true.
              end if
           end do

!**************************************************************************
        call MPI_BARRIER (MPI_COMM_WORLD,ierr)
        call MPI_BARRIER (MPI_COMM_WORLD,ierr)
!**************************************************************************


!**************************************************************************
        do ib=1,nbp

           dom(ib)%coarse_ng =.false.
           dom(ib)%fine_ng   =.false.
 
           do sync_dir = 1,3

              if (sync_dir.eq.1)  then
                 cpu_prev=dom(ib)%iprev
                 cpu_next=dom(ib)%inext
              else if (sync_dir.eq.2)  then
                 cpu_prev=dom(ib)%jprev
                 cpu_next=dom(ib)%jnext
              else if (sync_dir.eq.3)  then
                 cpu_prev=dom(ib)%kprev
                 cpu_next=dom(ib)%knext
              end if


              if (cpu_prev.ge.0) then
                 if(rdiv(dom_id(ib)).gt.rdiv(cpu_prev)) then
                    dom(ib)%coarse_ng =.true.
                 else if(rdiv(dom_id(ib)).lt.rdiv(cpu_prev)) then
                    dom(ib)%fine_ng =.true.
                 end if
              end if

              if (cpu_next.ge.0) then
                 if(rdiv(dom_id(ib)).gt.rdiv(cpu_next)) then
                    dom(ib)%coarse_ng =.true.
                 else if(rdiv(dom_id(ib)).lt.rdiv(cpu_next)) then
                    dom(ib)%fine_ng =.true.
                 end if
              end if

           end do


           if (dom(ib)%corprev1.ge.0) then
              if(rdiv(dom_id(ib)).gt.rdiv(dom(ib)%corprev1)) then
                 dom(ib)%coarse_ng =.true.
              else if(rdiv(dom_id(ib)).lt.rdiv(dom(ib)%corprev1)) then
                 dom(ib)%fine_ng =.true.
              end if
           end if
           if (dom(ib)%corprev2.ge.0) then
              if(rdiv(dom_id(ib)).gt.rdiv(dom(ib)%corprev2)) then
                 dom(ib)%coarse_ng =.true.
              else if(rdiv(dom_id(ib)).lt.rdiv(dom(ib)%corprev2)) then
                 dom(ib)%fine_ng =.true.
              end if
           end if
           if (dom(ib)%corprev3.ge.0) then
              if(rdiv(dom_id(ib)).gt.rdiv(dom(ib)%corprev3)) then
                 dom(ib)%coarse_ng =.true.
              else if(rdiv(dom_id(ib)).lt.rdiv(dom(ib)%corprev3)) then
                 dom(ib)%fine_ng =.true.
              end if
           end if
           if (dom(ib)%corprev4.ge.0) then
              if(rdiv(dom_id(ib)).gt.rdiv(dom(ib)%corprev4)) then
                 dom(ib)%coarse_ng =.true.
              else if(rdiv(dom_id(ib)).lt.rdiv(dom(ib)%corprev4)) then
                 dom(ib)%fine_ng =.true.
              end if
           end if
           if (dom(ib)%cornext1.ge.0) then
              if(rdiv(dom_id(ib)).gt.rdiv(dom(ib)%cornext1)) then
                 dom(ib)%coarse_ng =.true.
              else if(rdiv(dom_id(ib)).lt.rdiv(dom(ib)%cornext1)) then
                 dom(ib)%fine_ng =.true.
              end if
           end if
           if (dom(ib)%cornext2.ge.0) then
              if(rdiv(dom_id(ib)).gt.rdiv(dom(ib)%cornext2)) then
                 dom(ib)%coarse_ng =.true.
              else if(rdiv(dom_id(ib)).lt.rdiv(dom(ib)%cornext2)) then
                 dom(ib)%fine_ng =.true.
              end if
           end if
           if (dom(ib)%cornext3.ge.0) then
              if(rdiv(dom_id(ib)).gt.rdiv(dom(ib)%cornext3)) then
                 dom(ib)%coarse_ng =.true.
              else if(rdiv(dom_id(ib)).lt.rdiv(dom(ib)%cornext3)) then
                 dom(ib)%fine_ng =.true.
              end if
           end if
           if (dom(ib)%cornext4.ge.0) then
              if(rdiv(dom_id(ib)).gt.rdiv(dom(ib)%cornext4)) then
                 dom(ib)%coarse_ng =.true.
              else if(rdiv(dom_id(ib)).lt.rdiv(dom(ib)%cornext4)) then
                 dom(ib)%fine_ng =.true.
              end if
           end if

           if (dom(ib)%edgprev1.ge.0) then
              if(rdiv(dom_id(ib)).gt.rdiv(dom(ib)%edgprev1)) then
                 dom(ib)%coarse_ng =.true.
              else if(rdiv(dom_id(ib)).lt.rdiv(dom(ib)%edgprev1)) then
                 dom(ib)%fine_ng =.true.
              end if
           end if
           if (dom(ib)%edgprev2.ge.0) then
              if(rdiv(dom_id(ib)).gt.rdiv(dom(ib)%edgprev2)) then
                 dom(ib)%coarse_ng =.true.
              else if(rdiv(dom_id(ib)).lt.rdiv(dom(ib)%edgprev2)) then
                 dom(ib)%fine_ng =.true.
              end if
           end if
           if (dom(ib)%edgprev3.ge.0) then
              if(rdiv(dom_id(ib)).gt.rdiv(dom(ib)%edgprev3)) then
                 dom(ib)%coarse_ng =.true.
              else if(rdiv(dom_id(ib)).lt.rdiv(dom(ib)%edgprev3)) then
                 dom(ib)%fine_ng =.true.
              end if
           end if
           if (dom(ib)%edgprev4.ge.0) then
              if(rdiv(dom_id(ib)).gt.rdiv(dom(ib)%edgprev4)) then
                 dom(ib)%coarse_ng =.true.
              else if(rdiv(dom_id(ib)).lt.rdiv(dom(ib)%edgprev4)) then
                 dom(ib)%fine_ng =.true.
              end if
           end if
           if (dom(ib)%edgprev5.ge.0) then
              if(rdiv(dom_id(ib)).gt.rdiv(dom(ib)%edgprev5)) then
                 dom(ib)%coarse_ng =.true.
              else if(rdiv(dom_id(ib)).lt.rdiv(dom(ib)%edgprev5)) then
                 dom(ib)%fine_ng =.true.
              end if
           end if
           if (dom(ib)%edgprev6.ge.0) then
              if(rdiv(dom_id(ib)).gt.rdiv(dom(ib)%edgprev6)) then
                 dom(ib)%coarse_ng =.true.
              else if(rdiv(dom_id(ib)).lt.rdiv(dom(ib)%edgprev6)) then
                 dom(ib)%fine_ng =.true.
              end if
           end if
           if (dom(ib)%edgnext1.ge.0) then
              if(rdiv(dom_id(ib)).gt.rdiv(dom(ib)%edgnext1)) then
                 dom(ib)%coarse_ng =.true.
              else if(rdiv(dom_id(ib)).lt.rdiv(dom(ib)%edgnext1)) then
                 dom(ib)%fine_ng =.true.
              end if
           end if
           if (dom(ib)%edgnext2.ge.0) then
              if(rdiv(dom_id(ib)).gt.rdiv(dom(ib)%edgnext2)) then
                 dom(ib)%coarse_ng =.true.
              else if(rdiv(dom_id(ib)).lt.rdiv(dom(ib)%edgnext2)) then
                 dom(ib)%fine_ng =.true.
              end if
           end if
           if (dom(ib)%edgnext3.ge.0) then
              if(rdiv(dom_id(ib)).gt.rdiv(dom(ib)%edgnext3)) then
                 dom(ib)%coarse_ng =.true.
              else if(rdiv(dom_id(ib)).lt.rdiv(dom(ib)%edgnext3)) then
                 dom(ib)%fine_ng =.true.
              end if
           end if
           if (dom(ib)%edgnext4.ge.0) then
              if(rdiv(dom_id(ib)).gt.rdiv(dom(ib)%edgnext4)) then
                 dom(ib)%coarse_ng =.true.
              else if(rdiv(dom_id(ib)).lt.rdiv(dom(ib)%edgnext4)) then
                 dom(ib)%fine_ng =.true.
              end if
           end if
           if (dom(ib)%edgnext5.ge.0) then
              if(rdiv(dom_id(ib)).gt.rdiv(dom(ib)%edgnext5)) then
                 dom(ib)%coarse_ng =.true.
              else if(rdiv(dom_id(ib)).lt.rdiv(dom(ib)%edgnext5)) then
                 dom(ib)%fine_ng =.true.
              end if
           end if
           if (dom(ib)%edgnext6.ge.0) then
              if(rdiv(dom_id(ib)).gt.rdiv(dom(ib)%edgnext6)) then
                 dom(ib)%coarse_ng =.true.
              else if(rdiv(dom_id(ib)).lt.rdiv(dom(ib)%edgnext6)) then
                 dom(ib)%fine_ng =.true.
              end if
           end if
        end do

!**************************************************************************

        if (myrank.eq.0) then
           write (6,*) '===== end of connectivity information ===== '
           write (6,*) ' '
           write (6,*) ' '
           write (numfile,*) '=== end of connectivity information === '
           write (numfile,*) ' '
           write (numfile,*) ' '
        end if


5501  format(/1x,'neighbors of dom#:',i5,' ip:',i5,' in:',i5,
     & ' jp:', i5,' jn:',i5,' kp:',i5,' kn:',i5)
5601  format(/1x,'cor-neighbors of dom#:',i5,' pr1:',i5,' pr2:',i5,
     & ' pr3:', i5,' pr4:',i5,' nx1:',i5,' nx2:',i5,
     & ' nx3:', i5,' nx4:',i5)
5602  format(/1x,'edge-neighbors of dom#:',i5,' pr1:',i5,' pr2:',i5,
     & ' pr3:', i5,' pr4:',i5,' pr5:',i5,' pr6:',i5,
     & ' nx1:', i5,' nx2:',i5,' nx3:',i5,' nx4:',i5,
     & ' nx5:', i5,' nx6:',i5)

5502  format(/1x,'periodic ng of dom#:',i5,' i-prev periodic ng:',i5)
5503  format(/1x,'periodic ng of dom#:',i5,' i-next periodic ng:',i5)
5504  format(/1x,'periodic ng of dom#:',i5,' j-prev periodic ng:',i5)
5505  format(/1x,'periodic ng of dom#:',i5,' j-next periodic ng:',i5)
5506  format(/1x,'periodic ng of dom#:',i5,' k-prev periodic ng:',i5)
5507  format(/1x,'periodic ng of dom#:',i5,' k-next periodic ng:',i5)

!	xst=max(xst,xs2)
!	xen=min(xen,xe2)
!	yst=max(yst,ys2)
!	yen=min(yen,ye2)
!	zst=max(zst,zs2)
!	zen=min(zen,ze2)

        return
        end subroutine read_infodom
!##########################################################################
        subroutine datainfo
!##########################################################################
        use multidata
        use mpi
        implicit none
        integer :: ib,aa


        aa=110+myrank

        do ib=1,nbp
           write(aa,*) 'dom_id:',dom_id(ib),' ib:',ib
           
           write(aa,*) 'iprev',dom(ib)%iprev
           if(dom(ib)%iprev.ge.0) then
              write(aa,*) ' neig_add:',dom_ad(dom(ib)%iprev)
              write(aa,*) ' neig_ib:',dom_indid(dom(ib)%iprev)
           end if
           write(aa,*) '	'
           write(aa,*) 'inext',dom(ib)%inext
           if(dom(ib)%inext.ge.0) then
              write(aa,*) ' neig_add:',dom_ad(dom(ib)%inext)
              write(aa,*) ' neig_ib:',dom_indid(dom(ib)%inext)
           end if
           write(aa,*) '	'
           write(aa,*) 'jprev',dom(ib)%jprev
           if(dom(ib)%jprev.ge.0) then
              write(aa,*) ' neig_add:',dom_ad(dom(ib)%jprev)
              write(aa,*) ' neig_ib:',dom_indid(dom(ib)%jprev)
           end if
           write(aa,*) '	'
           write(aa,*) 'jnext',dom(ib)%jnext
           if(dom(ib)%jnext.ge.0) then
              write(aa,*) ' neig_add:',dom_ad(dom(ib)%jnext)
              write(aa,*) ' neig_ib:',dom_indid(dom(ib)%jnext)
           end if
           write(aa,*) '	'
           write(aa,*) 'kprev',dom(ib)%kprev
           if(dom(ib)%kprev.ge.0) then
              write(aa,*) ' neig_add:',dom_ad(dom(ib)%kprev)
              write(aa,*) ' neig_ib:',dom_indid(dom(ib)%kprev)
           end if
           write(aa,*) '	'
           write(aa,*) 'knext',dom(ib)%knext
           if(dom(ib)%knext.ge.0) then
              write(aa,*) ' neig_add:',dom_ad(dom(ib)%knext)
              write(aa,*) ' neig_ib:',dom_indid(dom(ib)%knext)
           end if
           write(aa,*) '*************'  
        end do


        return
        end subroutine datainfo
!##########################################################################
