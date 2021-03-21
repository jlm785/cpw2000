  subroutine GetStringFromSpeciesAndL(s, l, m, str)
  integer s, l, m
  
   
  character*5 str
  character*2 s_str, l_str, m_str

  return

  
  if(s .eq.-1) then
    write(s_str,'(a2)') "sa"
  else
    write(s_str,'("s",i1)') s
  endif
  
  if(l .eq.-1) then
    write(l_str,'(a2)') "la"
  else
    write(l_str,'("l",i1)') l
  endif

  if(m .eq.-10) then
    write(m_str,'(a2)') "ma"
  else
    if (m .gt. 0) then
      write(m_str,'("m+",i1)') m
    else
      write(m_str,'("m-",i1)') m    
    endif
    
  endif

  str = s_str // "_" // l_str // "_" // "_" // m_str
  
  end subroutine 
  
       program band

!      reads a file with band information and then calls out_band_xmgrace

!      adapted 2 June 2019 jlm.
!      modified for questions 18 June 2019 clr.
!      copyright  Jose Luis Martins/Carlos Loia Reis/INESC-MN

       implicit none
       integer, parameter                 :: REAL64 = selected_real_kind(12)

!      main variables

       character(len=60)                  ::  filename                   !  file to be read
       character(len=65)                  ::  filenamexmgr               !  file to be written
       integer                            ::  io                         !  tape number 

       character(len=50)                  ::  title                      !  title for plots
       character(len=140)                 ::  subtitle                   !  subtitle for plots
       integer                            ::  nstyle                     !  choice of plot style

       integer                            ::  neig                       !  number of bands
       integer                            ::  nrk                        !  number of k-vectors

       real(REAL64)                       ::  eref                       !  reference energy for plot
       integer                            ::  nocc                       !  number of occupied states (different color)

       integer                            ::  nvert                      !  number of vertical lines in plot
       integer                            ::  nlines                     !  number of lines in reciprocal space

       integer                            ::  nbaslcao                   !  number of atomic orbitals 

!      allocatable variables

       real(REAL64), allocatable          ::  xk(:)                      !  x coordinate of k-point in plot
       real(REAL64), allocatable          ::  e_of_k(:,:)                !  band energies of k-point in plot
       real(REAL64), allocatable          ::  pkn(:,:)                   !  size of k-point in plot
       real(REAL64), allocatable          ::  pkn_fld(:,:)                !  size of k-point in plot

       real(REAL64), allocatable          ::  xcvert(:)                  !  x coordinate of vertical line

       logical, allocatable               ::  ljump(:)                   !  indicates if the new line contains a jump from the preceeding
       integer, allocatable               ::  nkstep(:)                  !  number of steps in line
       character(len=6), allocatable      ::  label(:)                   !  label of symmetry k-points
       real(REAL64), allocatable          ::  xklab(:)                   !  x coordinate of label

       integer, allocatable               ::  infolcao(:,:)              !  information about the original atomic orbital.  (type of atom, atom of that type, n,l,m)
       real(REAL64), allocatable          ::  basxpsi(:,:,:)             !  |<bas|psi>|^2 for each k 


       integer                            :: ntype
       character(len=2), allocatable      :: nameat(:)

!      local variables

       integer          ::  ichoice

!      constants

       real(REAL64), parameter  :: ZERO = 0.0_REAL64, UM = 1.0_REAL64

!      counters

       integer      ::   n, n10, j, k, jj
       
       character(len=1)                    ::  yesno
       integer                             :: ioreplay = 335
       
       integer                             :: itype_in
       integer                             :: itask, ios
       integer                             :: l_in, m_in, iatom, natom_in
       integer, allocatable                :: at_list(:)
       
       integer                             :: degen_l(0:2)

       real(REAL64) :: xscale
       character*5 :: info_str
 
       real(REAL64) :: sum_pkn
       
       integer, allocatable :: iselect(:)
       integer              :: use_iselect

      degen_l(0) =1.0_REAL64
      degen_l(1) =1.0_REAL64
      degen_l(2) =1.0_REAL64

       open(unit = ioreplay, file = 'replay_band.dat',                   &
     &                  status = 'UNKNOWN', form = 'FORMATTED')
       


       write(6,*)
       write(6,*)  '  Which file do you want to read?'
       write(6,*)
       write(6,*)  '  1:  BAND.DAT (use this one)'
       write(6,*)  '  2:  BAND_SO.DAT (not implemented)'
       write(6,*)  '  3:  Choose name'
       write(6,*)
       write(6,*)  '  Enter choice?'
       write(6,*)

       read(5,*)  ichoice
       write(ioreplay,*) ichoice ,'   Which file do you want to read ?'

       if(ichoice == 1) then
         filename = 'BAND.DAT'
         filenamexmgr = 'band_fld.agr'
       elseif(ichoice == 2) then
         filename = 'BAND_SO.DAT'
         filenamexmgr = 'band_fld_so.agr'
       else
         write(6,*)  '  Enter filename (max 60 char)'
         write(6,*)
         read(5,*)  filename
         filenamexmgr = filename//'.agr'
       endif
       
         
       io = 15

       open(unit=io,file=filename,form='formatted')

       read(io,'(a50)') title
       read(io,'(a140)') subtitle
       read(io,'(i10)') nstyle

       read(io,'(2i10)') neig,nrk
       

       allocate(xk(nrk))
       allocate(e_of_k(neig,nrk))
       allocate(pkn(nrk,neig))
       allocate(pkn_fld(nrk,neig))

       do n = 1,nrk
       do j = 1,neig
!         pkn(n,j) = 2.0 - (j*1.9)/(neig*1.0)
       enddo
       enddo

       n10 = neig/10
       do n = 1,nrk
         read(io,'(f14.6)') xk(n)
         if(n10 > 0) then
           do j = 1,n10
             read(io,'(5x,10f14.6)') (e_of_k(jj,n),jj=10*(j-1)+1,10*j)
           enddo
         endif
         if(neig > 10*n10) then
           read(io,'(5x,10f14.6)') (e_of_k(jj,n),jj=10*n10+1,neig)
         endif
       enddo

       read(io,'(5x,f14.6,i8)') eref,nocc

       read(io,'(2i10)') nvert,nlines

       allocate(xcvert(nvert))
       allocate(ljump(nlines))
       allocate(nkstep(nlines))
       allocate(label(nvert+nlines))
       allocate(xklab(nvert+nlines))

       n10 = nvert/10
       if(n10 > 0) then
         do j = 1,n10
           read(io,'(5x,10f14.6)') (xcvert(jj),jj=10*(j-1)+1,10*j)
         enddo
       endif
       if(nvert > 10*n10) then
         read(io,'(5x,10f14.6)') (xcvert(jj),jj=10*n10+1,nvert)
       endif
 
       n10 = nlines/10
       if(n10 > 0) then
         do j = 1,n10
           read(io,'(5x,10(l5,i5))') (ljump(jj),nkstep(jj),              &
     &                                     jj=10*(j-1)+1,10*j)
         enddo
       endif
       if(nlines > 10*n10) then
         read(io,'(5x,1010(l5,i5))') (ljump(jj),nkstep(jj),              &
     &                                     jj=10*n10+1,nlines)
       endif
  
       n10 = (nvert+nlines)/10
       if(n10 > 0) then
         do j = 1,n10
           read(io,'(5x,10(2x,a6,f14.6))') (label(jj),xklab(jj),         &
     &                                     jj=10*(j-1)+1,10*j)
         enddo
       endif
       if(nvert+nlines > 10*n10) then
         read(io,'(5x,1010(2x,a6,f14.6))') (label(jj),xklab(jj),         &
     &                                     jj=10*n10+1,nvert+nlines)
       endif

       read(io,'(5x,i8)') nbaslcao

       allocate(infolcao(5,nbaslcao))
       allocate(basxpsi(nbaslcao,neig,nrk))

       do n = 1,nbaslcao
         read(io,'(5x,5i8,5x,"type, atom, n, l ,m")')                    &
     &                                (infolcao(j,n),j=1,5)
       enddo

       n10 = nbaslcao/10
       do k = 1,nrk
         do n = 1,neig
          if(n10 > 0) then
             do j = 1,n10
               read(io,'(5x,10f14.6)') (basxpsi(jj,n,k),                 &
     &                                      jj=10*(j-1)+1,10*j)
             enddo
           endif
           if(nbaslcao > 10*n10) then
               read(io,'(5x,10f14.6)') (basxpsi(jj,n,k),                 &
     &                                      jj=10*n10+1,nbaslcao)
           endif
         enddo
       enddo
       
       !! extra unfolding info
       do k = 1,nrk
         do n = 1,neig
            read(io,*) pkn_fld(k,n)
          enddo
       enddo
       !!

       !! extra atomic info
       
         read(io,*) ntype
         
         allocate(nameat(ntype))
            
         do n = 1,ntype
          read(io,'(a2)') nameat(n)
         enddo

       close(unit=io)
      
       open (unit = 64, file ="band_info.txt", form ="formatted")
       
       write(64,*) neig
       write(64,*) nbaslcao
       
       do j= 1, nbaslcao
        write(64,'(3x,6i5)') j,infolcao(1,j), infolcao(2,j), infolcao(3,j), infolcao(4,j), infolcao(5,j)        
       enddo

       write(64,'(i5)') ntype
            
        do n = 1,ntype
          write(64,'(a2)') nameat(n)
        enddo
       
       close(unit=64)
       
       write(*,'("   Total number of  eigenvalues is ", i5)'),           neig
       write(*,'("   Total number of localized orbitals is    ", i5)'),  nbaslcao
       write(*, *) '  this system has the following infolcao table:'
       write(*,*)
       write(*,'(3x, 3a6, a2, 3x, a2, 3x, a2)') "iorb ", "itype", "iatom", "n", "l", "m"      
       write(*,*)       
       do j= 1, nbaslcao
        write(*,'(3x,6i5)') j,infolcao(1,j), infolcao(2,j), infolcao(3,j), infolcao(4,j), infolcao(5,j)        
       enddo
       write(*,*)       

      write(*,*) '  do you want to proceed?'
      read(5,*,iostat=ios) yesno
      write(ioreplay,*) yesno ,'   do you want to proceed?'
                  
      if(yesno == 'y' .or. yesno == 'Y') then
      else
        stop        
      endif

      write(*,*) '  do you want to give iselect table?'
      read(5,*,iostat=ios) yesno
      write(ioreplay,*) yesno ,'   do you want to give iselect table?'
      
      if(yesno == 'y' .or. yesno == 'Y') then

       allocate(iselect(nbaslcao))
       do j= 1, nbaslcao
        read(5,*) iselect(j)       
       enddo       
       use_iselect = 1
       
          do k = 1,nrk
          do n = 1,neig
            pkn(k,n) = ZERO
            do j = 1,nbaslcao         
              if( iselect(j) == 1) then
                pkn(k,n) = pkn(k,n) + basxpsi(j,n,k)
              endif
            enddo
          enddo
          write(6,'(i5,200f10.4)') k,(pkn(k,n),n=1,neig)
          enddo
          
      goto 667
       
      else      
       use_iselect =0
      endif
      
      


      !!! Decision Tree
      
      write(*,*) '  do you want to see all atoms?'
      read(5,*,iostat=ios) yesno
      write(ioreplay,*) yesno ,'   do you want to see all atoms?'
      
      if(yesno == 'y' .or. yesno == 'Y') then

          write(*,*) '  Enter l quantum number  (integer beteween 0 and 2) (-1 for all) '
          read(5,*,iostat=ios) l_in
          write(ioreplay,*) l_in,'   l quantum number (-1 for all)'

          write(*,*) '  Enter m quantum number  (integer beteween -l and l) (-10 for all) '
          read(5,*,iostat=ios) m_in
          write(ioreplay,*) m_in,'   m quantum number (-10 for all)'

          
          call GetStringFromSpeciesAndL(-1, l_in, m_in, info_str)
          
          !!! 1          
          !!! PKN Maniplulation 
          do k = 1,nrk
          do n = 1,neig
            pkn(k,n) = ZERO
            do j = 1,nbaslcao         
              if( (infolcao(4,j) == l_in .or. l_in == -1) .and. (infolcao(5,j) == m_in .or. m_in == -10)) then
                pkn(k,n) = pkn(k,n) + basxpsi(j,n,k)
              endif
            enddo
          enddo
          write(6,'(i5,200f10.4)') k,(pkn(k,n),n=1,neig)
          enddo
          !!!
          
                
      else ! No I dont want to see all atoms
        
        write(*,*) '  Enter atomic species itype (integer).'
        write(*,*) '  see block Chemical_Species_Label in cpw.in '
        read(5,*,iostat=ios) itype_in
        write(ioreplay,*) itype_in,'   atomic species'
        
        write(*,*) ' Do you want to see all atoms of this type? (yes/no)' 
        read(5,*,iostat=ios) yesno        
        write(ioreplay,*) yesno,'   Do you want to see all atoms of this type?'        
        
        if(yesno == 'y' .or. yesno == 'Y') then
        
          write(*,*) '  Enter l quantum number  (integer beteween 0 and 2) (-1 for all)  '
          read(5,*,iostat=ios) l_in
          write(ioreplay,*) l_in,'   l quantum number (-1 for all)'
          
          write(*,*) '  Enter m quantum number  (integer beteween -l and l) (-10 for all) '
          read(5,*,iostat=ios) m_in
          write(ioreplay,*) m_in,'   m quantum number (-10 for all)'
                    
          call GetStringFromSpeciesAndL(itype_in, l_in, m_in, info_str)
          
          !!!! 2
          !!! PKN Maniplulation 
          do k = 1,nrk
          do n = 1,neig
            pkn(k,n) = ZERO
            do j = 1,nbaslcao         
              if((infolcao(4,j) == l_in .or. l_in == -1) .and. (infolcao(5,j) == m_in .or. m_in == -10) .and. infolcao(1,j) == itype_in) then
                pkn(k,n) = pkn(k,n) + basxpsi(j,n,k)
              endif
            enddo
          enddo
          write(6,'(i5,200f10.4)') k,(pkn(k,n),n=1,neig)
          enddo
          !!!
                  
        else ! No I dont want to see all atoms of this type
        
          write(*,*) '  how many atoms? '
          read(5,*,iostat=ios) natom_in
          write(ioreplay,*) natom_in,'   how many atoms?'
          allocate (at_list(natom_in))
          
          write(*,*) '  enter atom list (1 index per line) '          
          do iatom = 1, natom_in
            read(5,*,iostat=ios) at_list(iatom)
            write(ioreplay,*) at_list(iatom), '   iatom'
          enddo
        
          write(*,*) '  Enter l quantum number (integer beteween 0 and 2) (-1 for all)   '
          read(5,*,iostat=ios) l_in
          write(ioreplay,*) l_in,'   l quantum number (integer beteween 0 and 2) (-1 for all) '
          
          write(*,*) '  Enter m quantum number  (integer beteween -l and l) (-10 for all) '
          read(5,*,iostat=ios) m_in
          write(ioreplay,*) m_in,'   m quantum number (-10 for all)'
          
          write(*,*) '  Enter scale  (real number to multilply weights)'
          write(*,*) '  use 1.0 if you dont want to use this feature'
          read(5,*,iostat=ios) xscale
          write(ioreplay,*) xscale,'   Enter scale'
          
          !!!! 3
          !!! PKN Maniplulation 
          do k = 1,nrk
          do n = 1,neig
            pkn(k,n) = ZERO
            do j = 1,nbaslcao
!            write(*,*) 'infolcao', infolcao(2,j)
              do iatom = 1, natom_in
                if((infolcao(4,j) == l_in .or. l_in == -1) .and. (infolcao(5,j) == m_in .or. m_in == -10) .and. infolcao(1,j) == itype_in .and. infolcao(2,j) == at_list(iatom) ) then
                  
!                  if(k==1 .and. n==1) then
!                  write(*,*) 'using atom #',  at_list(iatom)
!                  write(*,*) 'itype #',  infolcao(1,j)
!                  write(*,*) 'iatom #',  infolcao(2,j)
!                  write(*,*) 'n     #',  infolcao(3,j)
!                  write(*,*) 'l     #',  infolcao(4,j)
!                  write(*,*) 'm     #',  infolcao(5,j)
!                  endif
                  
                  pkn(k,n) = pkn(k,n) + basxpsi(j,n,k)*xscale
                endif
              enddo
            enddo
          enddo
!          write(6,'(i5,200f10.4)') k,(pkn(k,n),n=1,neig)
          enddo
          !!!
!          stop
        
        endif
      
      endif
      
!      stop

667 continue

       do k = 1,nrk
          sum_pkn = 0.0
            do n = 1,neig
              sum_pkn = sum_pkn + pkn(k,n)
            enddo

            do n = 1,neig
              pkn(k,n) = 2.0*pkn(k,n)
            enddo


            write(*,'("sum for kpt#",i8,f12.5)') k, sum_pkn 
       enddo



       do k = 1,nrk
         do n = 1,neig
!            pkn(k,n) = pkn(k,n) * 1.5

            pkn(k,n) = pkn(k,n)* pkn_fld(k,n)
          enddo
       enddo
       
       write(subtitle,*) ""

       call out_band_fold_xmgrace(filenamexmgr,io,                       &
     &        title,subtitle,nstyle,                                     &
     &        pkn,neig,nrk,xk,e_of_k,eref,nocc,                          &
     &        nvert,xcvert,nlines,ljump,nkstep,label,xklab)
     
       call system("cp " // filenamexmgr // "band_fld_"//info_str//".agr");
       
       call system("cat Default.agr band_fld.agr    > band_fld_all.agr");
       call system("cat Default.agr band_fld_so.agr > band_fld_so_all.agr");

       deallocate(xk)
       deallocate(e_of_k)

       deallocate(xcvert)
       deallocate(ljump)
       deallocate(nkstep)
       deallocate(label)
       deallocate(xklab)

       stop
       end 
