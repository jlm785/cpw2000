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

!>     This subroutine calculates the dimensions required to
!>     store the band structure along a path from
!>     dat on filename (default BAND_LINES.DAT) or invents default.

       subroutine out_band_circuit_size(filename,iotape,ninterp,adot,            &
     &                  neig,nrk,nlines,nvert)

!      Written, April 10, 2014. jlm
!      Modified (density of points) 5 August 2014. JLM
!      Modified, documentation, 4 February 2020. JLM
!      Modified, ninterp, 12 June 2020. JLM
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
      
!      output

       integer, intent(out)               ::  neig                       !<  number of eigenvectors required
       integer, intent(out)               ::  nrk                        !<  number of k-points in the circuit
       integer, intent(out)               ::  nlines                     !<  number of lines in reciprocal space
       integer, intent(out)               ::  nvert                      !<  number of vertical lines

!      local variables

       real(REAL64)      ::  rkbegin(3)               !  begin of line in reciprocal space (lattice coordinates)
       real(REAL64)      ::  rkend(3)                 !  end of line in reciprocal space (lattice coordinates)
       real(REAL64)      ::  rkold(3)                 !  previous end of line in reciprocal space (lattice coordinates)
       integer           ::  nkstep                   !  number of steps in line
       real(REAL64)      ::  rkdist
       real(REAL64)      ::  dist
       real(REAL64)      ::  dx                       !  target average distance between k-points, if it is negative or zero or non-existent use other data

       logical           ::  ldens
       integer           ::  ioerr, ioerr3

       real(REAL64)      ::  vcell, bdot(3,3)

       character(len=120) ::  fline 

!      constants

       real(REAL64), parameter  :: ZERO = 0.0_REAL64
       real(REAL64), parameter  :: EPS = 0.00000001_REAL64

!      counters

       integer    ::  i, j, n


       call adot_to_bdot(adot,vcell,bdot)

!      opens file
     
       open(unit=iotape,file=filename, status='old',iostat=ioerr,        &
     &      form = 'formatted')

       if(ioerr == 0) then

!        reads the first line which may have different formats

         read(iotape,'(a120)') fline

         ldens = .TRUE.

         read(fline,*,iostat=ioerr3) nlines,neig,dx

         if(ioerr3 /= 0) then

           ldens = .FALSE.

           read(fline,*) nlines,neig

         endif

         if(nlines < 1) then
           write(6,*)
           write(6,'("  input error in out_band_circuit_size:    ",      &
     &        "number of lines = ",i10)') nlines

           stop

         else

           if(dx < EPS) then

             ldens = .FALSE.

           endif

         endif

         if(neig < 1) then
           write(6,*)
           write(6,'("  input error in out_band_circuit_size:    ",      &
     &         "number of bands = ",i10)') neig

           stop

         endif

       else
         write(6,*)
         write(6,*) '  error in out_band_circuit_size:    ',             &
     &          'problem opening   ',filename 

         stop

       endif

!      reads the other nlines

       nrk = 0
       nvert = 0
       do n=1,nlines
         read(iotape,*) (rkbegin(j),j=1,3),(rkend(j),j=1,3),nkstep
         
         if(ldens) then
           rkdist = ZERO
           do i=1,3
           do j=1,3
             rkdist = rkdist + (rkend(i)-rkbegin(i))*                    &
     &               bdot(i,j)*(rkend(j)-rkbegin(j))
           enddo
           enddo
           rkdist = sqrt(rkdist)
           nkstep = nint(rkdist/dx)
         endif

!        reduces by ninterp

         nkstep = nkstep/ninterp 

         if(nkstep < 1) nkstep = 2

         dist = 1.0
         if(n /= 1) then
           dist = ZERO
           do i=1,3
           do j=1,3
             dist = dist + (rkbegin(i)-rkold(i))*                        &
     &           bdot(i,j)*(rkbegin(j)-rkold(j))
           enddo
           enddo
         endif
         if(dist < 0.0001) then
           nrk = nrk + nkstep
           nvert = nvert + 1
         else
           nrk = nrk + nkstep + 1
           nvert = nvert + 2
         endif

         do i=1,3
           rkold(i) = rkend(i)
         enddo

       enddo
                                                                   
       close(unit=iotape)


       return
       end subroutine out_band_circuit_size
