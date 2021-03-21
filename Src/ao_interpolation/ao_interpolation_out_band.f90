!------------------------------------------------------------!
! This file is distributed as part of the cpw2000 code and   !
! under the terms of the GNU General Public License. See the !
! file `LICENSE' in the root directory of the cpw2000        !
! distribution, or http://www.gnu.org/copyleft/gpl.txt       !
!                                                            !
! The webpage of the cpw2000 code is not yet written         !
!                                                            !
! The cpw2000 code is hosted on GitHub:                      !
!                                                            !
! https://github.com/jlm785/cpw2000                          !
!------------------------------------------------------------!

!>     This subroutine calculates the band structure along a path
!>     using atomic orbital interpolation

       subroutine ao_interpolation_out_band(title,subtitle,              &
     &        noiData,ztot,adot)
       
!      Written by Carlos Loia Reis before 2017, based on previous code,
!      Modified, documentation, May 2020. JLM

!      copyright  Carlos Loia Reis /INESC-MN.

!      version 4.98


       use NonOrthoInterp

       implicit none

       integer, parameter                 :: REAL64 = selected_real_kind(12)

!      input

       character(len=50), intent(in)      ::  title                      !<  title for plots
       character(len=140), intent(in)     ::  subtitle                   !<  subtitle for plots

       type(noiData_t)                    :: noiData                     !<  see NonOrthoInterp
       
       real(REAL64), intent(in)           ::  ztot                       !<  total charge density (electrons/cell)
       real(REAL64), intent(in)           ::  adot(3,3)                  !<  metric in direct space
       
!      allocatable arrays for Brillouin zone path

       integer                            ::  nlines                     !  number of lines in reciprocal space
       integer, allocatable               ::  nkstep(:)                  !  number of steps in line
       logical, allocatable               ::  ljump(:)                   !  indicates if the new line contains a jump from the preceeding
       integer                            ::  nvert                      !  number of vertical lines in plot
       real(REAL64), allocatable          ::  xcvert(:)                  !  x coordinate of vertical line
       real(REAL64), allocatable          ::  xk(:)                      !  x coordinate of k-point in plot
       real(REAL64), allocatable          ::  rk(:,:)                    !  x coordinate of k-point in plot
       real(REAL64), allocatable          ::  e_of_k(:,:)                !  band energies of k-point in plot
       real(REAL64), allocatable          ::  e_of_k_so(:,:)             !  spin-orbit band energies of k-point in plot
       character(len=6), allocatable      ::  label(:)                   !  label of symmetry k-points
       real(REAL64), allocatable          ::  xklab(:)                   !  x coordinate of label

!      local variables


       integer           ::  neig           !  number of eigenvectors required (maybe modified on output)
       real(REAL64)      ::  rkpt(3)        !  j-th component in lattice coordinates of the k-point

       real(REAL64)      ::  eref           !  reference energy for plot
       integer           ::  nocc           !  number of occupied states (different color)
       integer           ::  nstyle         !  choice of plot style

       integer           ::  irk
       integer           ::  iotape
       integer           ::  nrk2
       logical           ::  lfile

!      counters

       integer    ::  j, n

       real(REAL64), allocatable :: ev_interp(:)
   

       iotape = 13
       call out_band_circuit_size('BAND_LINES.DAT',iotape,1,adot,        &
     &                  neig,nrk2,nlines,nvert)

       allocate(xk(nrk2))
       allocate(rk(3,nrk2))
       allocate(xcvert(nvert))
       allocate(ljump(nlines))
       allocate(nkstep(nlines))
       allocate(label(nvert+nlines))
       allocate(xklab(nvert+nlines))

       call out_band_get_circuit('BAND_LINES.DAT',iotape,1,adot,         &
     &                  xk,rk,xcvert,ljump,nkstep,label,xklab,           &
     &                  neig,nrk2,nlines,nvert)


       allocate(e_of_k(neig,nrk2))
       allocate(e_of_k_so(2*neig,nrk2))
      
!       mxdbnd = neig
       
       allocate(ev_interp(noiData%nband))

       do irk=1,nrk2

!        loop over k-points

         do j=1,3
           rkpt(j) = rk(j,irk)
         enddo
         
         write(*,'(i5,3f8.5)') irk, rkpt(1),rkpt(2),rkpt(3)
         
         call NonOrthoInterpRun(noiData,rkpt,ev_interp)

         if(noiData%lso==1) then
          do j=1,2*neig
            e_of_k_so(j,irk) = ev_interp(j)
          enddo
         else
          do j=1,neig
            e_of_k(j,irk) = ev_interp(j)
          enddo         
         
         endif
       enddo

       iotape = 15
       nstyle = 1

!      writes the output files for xmgrace and gnuplot
          
         if(noiData%lso==1) then       
            n = min(nint(ztot + 0.01),2*neig)
            eref = e_of_k_so(n,1)
            do irk = 1,nrk2
            do j=1,n
              if(e_of_k_so(j,irk) > eref) eref = e_of_k_so(j,irk)
            enddo
            enddo
       nocc = n
!       do irk = 1,nrk2
!         if(e_of_k(n+1,irk) < eref) then
!           nocc = 0
!
!           exit
!
!         endif
!       enddo

       call out_band_xmgrace('band_interp_so.agr',iotape,                &
     &    title,subtitle,nstyle,                                         &
     &    2*neig,nrk2,xk,e_of_k_so,eref,nocc,                            &
     &    nvert,xcvert,nlines,ljump,nkstep,label,xklab)

       call out_band_dots_xmgrace('dots_band_interp_so.agr',iotape,      &
     &    nint(ztot),2*neig,nrk2,xk,e_of_k_so,      &
     &    eref,nvert,xcvert,nlines,ljump,nkstep,label,xklab)
     
         else         
            n = min(nint(0.5*ztot + 0.01),neig)
            eref = e_of_k(n,1)
            do irk = 1,nrk2
            do j=1,n
              if(e_of_k(j,irk) > eref) eref = e_of_k(j,irk)
            enddo
            enddo

            nocc = n
!       do irk = 1,nrk2
!         if(e_of_k(n+1,irk) < eref) then
!           nocc = 0
!
!           exit
!
!         endif
!       enddo

        call out_band_xmgrace('band_interp.agr',iotape,                  &
     &        title,subtitle,nstyle,                                     &
     &        neig,nrk2,xk,e_of_k,eref,nocc,                             &
     &        nvert,xcvert,nlines,ljump,nkstep,label,xklab)
     
        call out_band_dots_xmgrace('dots_band_interp.agr',iotape,        &
     &        nint(ztot)/2,neig,nrk2,xk,e_of_k,                          &
     &        eref,nvert,xcvert,nlines,ljump,nkstep,label,xklab)
             
         endif
         
!      cheat code         
!      loop over k-points
       lfile = .false.
        open(unit=144,file='BAND_MESH.DAT',form='formatted', err=10)
       lfile= .true.
10      continue

        if(lfile .and. (command_argument_count() >=1) ) then
!        if(lfile .and. (IARGC()>=1)) then
          write(*,*) 'BAND_MESH.DAT detected !!!'
          write(*,*) 'begining band structure calculation with old style files.'
          open(127,file="BAND_ENERGY.DAT",form="formatted")        
          read(144,*) nrk2,neig
          if (noiData%lso==1) then
            write(127,*) '#',nrk2,2*nint(ztot)
          else
            write(127,*) '#',nrk2,nint(ztot)
          endif
              
          write(127,*) '#   Energies in eV, k at end of line'
                  
          do irk=1,nrk2
              read(144,*) rkpt(1),rkpt(2),rkpt(3)
              call NonOrthoInterpRun(noiData,rkpt,ev_interp)
              write(*,'(i5,3f8.3)') irk, rkpt(1),rkpt(2),rkpt(3)
              if (noiData%lso==1) then
                WRITE(127,'(200F22.8)') (irk-1)*1.0D0, (ev_interp(J)*27.212,J=1,2*nint(ztot)),rkpt(1),rkpt(2),rkpt(3)
              else
                WRITE(127,'(200F22.8)') (irk-1)*1.0D0, (ev_interp(J)*27.212,J=1,nint(ztot)),rkpt(1),rkpt(2),rkpt(3)
              endif       
          enddo
          close(unit=127)
          close(unit=144)
        endif

!      Cleanup

       deallocate(nkstep)
       deallocate(ljump)
       deallocate(xcvert)
       deallocate(xk)
       deallocate(rk)
       deallocate(e_of_k)
       deallocate(e_of_k_so)
       deallocate(label)
       deallocate(xklab)

       return
       end subroutine ao_interpolation_out_band
