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

!>     This subroutine calculates the band structure path from
!>     dat on filename (default BAND_LINES.DAT) or invents default.

       subroutine out_band_get_circuit(filename, iotape, ninterp, adot,  &
     &                  xk, rk, xcvert, ljump, nkstep, label, xklab,     &
     &                  neig, nrk, nlines, nvert)

!      Written, April 10, 2014. jlm
!      Modified, ninterp, 11 June 2020. JLM
!      Modified, dx < EPS, 19 August 2020. JLM
!      copyright  Jose Luis Martins/INESC-MN

!      version 4.98

       implicit none

       integer, parameter          :: REAL64 = selected_real_kind(12)


!      input

       character(len=*), intent(in)       ::  filename                   !<  file to be read
       integer, intent(in)                ::  iotape                     !<  tape number 
       real(REAL64), intent(in)           ::  adot(3,3)                  !<  metric in direct space

       integer, intent(in)                ::  ninterp                    !<  density of points is reduced by ninterp.  In normal cases ninterp = 1.

       integer, intent(in)                ::  neig                       !<  number of eigenvectors required
       integer, intent(in)                ::  nrk                        !<  number of k-points in the circuit
       integer, intent(in)                ::  nlines                     !<  number of lines in reciprocal space
       integer, intent(in)                ::  nvert                      !<  number of vertical lines

!      output

       real(REAL64), intent(out)          ::  xk(nrk)                    !<  x coordinate of k-point in plot
       real(REAL64), intent(out)          ::  rk(3,nrk)                  !<  k-point to be plotted in lattice coordinates
       real(REAL64), intent(out)          ::  xcvert(nvert)              !<  x coordinate of vertical line
       logical, intent(out)               ::  ljump(nlines)              !<  indicates if the new line contains a jump from the preceeding
       integer, intent(out)               ::  nkstep(nlines)             !<  number of steps in line
       character(len=6), intent(out)      ::  label(nvert+nlines)        !<  label of symmetry k-points
       real(REAL64), intent(out)          ::  xklab(nvert+nlines)        !<  x coordinate of label

!      local allocatable arrays

       real(REAL64), allocatable          ::  rkbegin(:,:)               !  begin of line in reciprocal space (lattice coordinates)
       real(REAL64), allocatable          ::  rkend(:,:)                 !  end of line in reciprocal space (lattice coordinates)
       real(REAL64), allocatable          ::  rkdist(:)                  !  length of line in reciprocal space (lattice coordinates)

       character(len=6), allocatable      ::  labbeg(:)                  !  label of symmetry k-points begin of line
       character(len=6), allocatable      ::  labend(:)                  !  label of symmetry k-points end of line
       character(len=6), allocatable      ::  lablines(:)                !  label of symmetry lines

!      local variables 

       integer           ::  neigloc,nrkloc,nlinesloc,nvertloc
       real(REAL64)      ::  vcell, bdot(3,3)
       real(REAL64)      ::  dist, xjump
       real(REAL64)      ::  dx                       !  target average distance between k-points, if it is negative or zero or non-existent use other data

       logical           ::  ldens
       integer           ::  ioerr, ioerr3

       character(len=120) ::  fline 

!      constants

       real(REAL64), parameter  :: ZERO = 0.0_REAL64
       real(REAL64), parameter  :: EPS = 0.00000001_REAL64

!      counters

       integer    ::  i, j, n
       integer    ::  irk, jrk


       call adot_to_bdot(adot,vcell,bdot)

!      allocations

       allocate(rkbegin(3,nlines))
       allocate(rkend(3,nlines))
       allocate(rkdist(nlines))
       allocate(labbeg(nlines))
       allocate(labend(nlines))
       allocate(lablines(nlines))
       
!      opens file

       open(unit=iotape,file=filename, status='old',iostat=ioerr,        &
     &      form = 'formatted')

       if(ioerr == 0) then

!        reads the first line which may have different formats

         read(iotape,'(a120)') fline

         ldens = .TRUE.

         read(fline,*,iostat=ioerr3) nlinesloc,neigloc,dx

         if(ioerr3 /= 0) then

           ldens = .FALSE.

           read(fline,*) nlinesloc,neigloc

         else

           if(dx < EPS) then

             ldens = .FALSE.

           endif

         endif

         if(nlinesloc < 1) then
           write(6,*)
           write(6,'("  input error in out_band_get_circuit:    ",       &
     &        "number of lines = ",i10)') nlinesloc

           stop

         endif

         if(neigloc < 1) then
           write(6,*)
           write(6,'("  input error in out_band_get_circuit:    ",       &
     &         "number of bands = ",i10)') neigloc

           stop

         endif

       else
         write(6,*)
         write(6,*) '  error in out_band_get_circuit:    ',              &
     &          'problem opening   ',filename 

         stop

       endif
       
       if(nlinesloc /= nlines) then
         write(6,'("   error in out_band_get_circuit:  inconsistent ",   &
     &     "number of lines ",2i5)') nlinesloc,nlines

         stop

       endif

       if(neigloc /= neig) then
         write(6,'("   error in out_band_get_circuit:  inconsistent ",   &
     &     "number of bands ",2i5)') neigloc, neig

         stop

       endif


!      initialize labels

       do n=1,nlines
         labbeg(n) = '      '
         labend(n) = '      '
         lablines(n) = '      '
       enddo

       nrkloc = 0
       nvertloc = 0
       do n=1,nlines

         read(iotape,*,IOSTAT=ioerr) (rkbegin(j,n),j=1,3),               &
     &       (rkend(j,n),j=1,3),nkstep(n),                               &
     &       labbeg(n),lablines(n),labend(n)
          
         if(ioerr /= 0) then
           backspace(iotape)
           read(iotape,*,IOSTAT=ioerr) (rkbegin(j,n),j=1,3),             &
     &                         (rkend(j,n),j=1,3),nkstep(n)
         endif

         rkdist(n) = ZERO
         do i=1,3
         do j=1,3
           rkdist(n) = rkdist(n) + (rkend(i,n)-rkbegin(i,n))*           &
     &                   bdot(i,j)*(rkend(j,n)-rkbegin(j,n))
         enddo
         enddo
         rkdist(n) = sqrt(rkdist(n))

         if(ldens) then
           nkstep(n) = nint(rkdist(n)/dx)
         endif

!        reduces by ninterp

         nkstep(n) = nkstep(n)/ninterp 

         if(nkstep(n) < 1) nkstep(n) = 2

         dist = 1.0
         if(n /= 1) then
           dist = ZERO
           do i=1,3
           do j=1,3
             dist = dist + (rkbegin(i,n)-rkend(i,n-1))*                 &
     &           bdot(i,j)*(rkbegin(j,n)-rkend(j,n-1))
           enddo
           enddo
         endif
         if(dist < 0.0001) then
           ljump(n) = .FALSE.
           nrkloc = nrkloc + nkstep(n)
           nvertloc = nvertloc + 1
         else
           ljump(n) = .TRUE.
           nrkloc = nrkloc + nkstep(n) + 1
           nvertloc = nvertloc + 2
         endif

       enddo
                                                                   
       close(unit=iotape)

       
       if(nvertloc /= nvert) then
         write(6,'("   error in out_band_get_circuit:  inconsistent ",   &
     &     "number of vertical lines ",2i5)') nvertloc, nvert

         stop

       endif

       if(nrkloc /= nrk) then
         write(6,'("   error in out_band_get_circuit:  inconsistent ",   &
     &     "number of k-points ",2i5)') nrkloc, nrk

         stop

       endif

       xjump = rkdist(1)
       do n=1,nlines
         if(rkdist(n) < xjump) xjump = rkdist(n)
       enddo
       xjump = xjump / 5
       if(xjump < 0.1) xjump = 0.1

!      first line

       irk = 0
       n = 1
       xcvert(1) = ZERO
       xklab(1) = ZERO
       label(1) = labbeg(1)

       do i=0,nkstep(n)
         irk = irk+1
         xk(irk) = (i*rkdist(n))/nkstep(n)
         do j=1,3
           rk(j,irk) = (i*rkend(j,n)) + ((nkstep(n)-i)*rkbegin(j,n))
           rk(j,irk) = rk(j,irk) / nkstep(n)
         enddo
       enddo

       xcvert(2) = xk(irk)
       xklab(2) = xk(irk)
       label(2) = labend(1)
       label(nvert+1) = lablines(1)
       xklab(nvert+1) = (xklab(2) + xklab(1))/2
       jrk = 2

       if(nlines > 1) then

         do n=2,nlines

           if(ljump(n)) then

             irk = irk + 1
             xk(irk) = xk(irk-1) + xjump
             jrk = jrk + 1
             xcvert(jrk) = xk(irk)
             xklab(jrk) = xk(irk)
             label(jrk) = labbeg(n)
             do j=1,3
               rk(j,irk) = rkbegin(j,n)
             enddo

           endif

           do i=1,nkstep(n)
             irk = irk+1
             xk(irk) = xk(irk-1) + rkdist(n)/nkstep(n)
             do j=1,3
               rk(j,irk) = (i*rkend(j,n)) + ((nkstep(n)-i)*rkbegin(j,n))
               rk(j,irk) = rk(j,irk) / nkstep(n)
             enddo
           enddo
           jrk = jrk + 1
           xcvert(jrk) = xk(irk)  
           xklab(jrk) = xk(irk)
           label(jrk) = labend(n)
           xklab(nvert+n) = (xklab(jrk) + xklab(jrk-1))/2
           label(nvert+n) = lablines(n)

         enddo

       endif


       deallocate(rkbegin)
       deallocate(rkend)
       deallocate(rkdist)
       deallocate(labbeg)
       deallocate(labend)
       deallocate(lablines)


       return
       end subroutine out_band_get_circuit
