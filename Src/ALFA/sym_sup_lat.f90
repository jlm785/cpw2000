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

!>     searches if the structure is a perfect superlattice

       subroutine sym_sup_lat(ipr,tol,isup,jsup,nsup,                    &
     & adot,ntype,natom,rat,                                             &
     & mxdtyp,mxdatm)

!      Written 2 April 2004. JLM
!      Modified 8 June 2014, f90. JLM
!      Modified 30 November 2016.  JLM  unrelated atoms
!      Modified, documentation, December 2019. JLM
!      copyright INESC-MN/Jose Luis Martins

!      version 4.94


       implicit none

       integer, parameter          :: REAL64 = selected_real_kind(12)

!      input

       integer, intent(in)                ::  mxdtyp                     !<  array dimension of types of atoms
       integer, intent(in)                ::  mxdatm                     !<  array dimension of types of atoms

       integer, intent(in)                ::  ipr                        !<  print flag
       real(REAL64), intent(in)           ::  tol                        !<  tolerance for symmetry recognition

       real(REAL64), intent(in)           ::  adot(3,3)                  !<  metric in direct space
       integer, intent(in)                ::  ntype                      !<  number of types of atoms
       integer, intent(in)                ::  natom(mxdtyp)              !<  number of atoms of type i
       real(REAL64), intent(in)           ::  rat(3,mxdatm,mxdtyp)       !<  k-th component (in lattice coordinates) of the position of the n-th atom of type i

!      output

       integer, intent(out)               ::  isup                       !<  isup=0 no superlattice; isup = -1 failure; type of atoms with minimal number
       integer, intent(out)               ::  nsup                       !<  number of atoms not related by superlattice translations (atoms with minimal number)
       integer, intent(out)               ::  jsup(mxdatm)               !<  indicates the nsup atoms not related by superlattice translations (atoms with minimal number)

!      local variables

       integer         ::  jsucc,ntotal,ncount,itymin,natmin
       real(REAL64)    ::  frac(3),diff(3),sdif
       integer         ::  ntr                                           !  number of translations
       real(REAL64)    ::  r0(3),rj(3),dist,distmin

!      allocatable variables

       real(REAL64), allocatable    ::  fsup(:,:)                        !  fractional superlattice translation
       integer, allocatable         ::  jpoint(:,:)                      !  atoms connected by translation
       integer, allocatable         ::  mark(:)                          !  marks accounted atoms

!      parameters

       real(REAL64), parameter :: ZERO = 0.0_REAL64

!      counters

       integer i, j, k, j1, j2, jnear



!      identifies type of atoms with smallest number

       itymin = 1
       natmin = natom(1)
       do i=1,ntype
         if(natom(i) < natmin) then
           itymin = i
           natmin = natom(i)
         endif
       enddo

       isup = 0

       if(natmin > 1) then

         allocate(fsup(3,natmin))
         allocate(jpoint(natmin,natmin+1))

         fsup(1,1) = ZERO
         fsup(2,1) = ZERO
         fsup(3,1) = ZERO

         do j = 1,natmin
           jpoint(j,1) = j
         enddo

         ntr = 1

         do j=2,natmin
           do k=1,3
             frac(k) = rat(k,j,itymin) - rat(k,1,itymin)
             frac(k) = frac(k) - nint(frac(k))
           enddo

!          tests operation

           ntotal = 0
           do i=1,ntype
             ntotal = ntotal + natom(i)
           enddo

           ncount = 0
           do i=1,ntype

             do j1=1,natom(i)
               jsucc = 0
               do j2=1,natom(i)

                 sdif = zero
                 do k=1,3
                   diff(k) = rat(k,j2,i) - rat(k,j1,i) - frac(k)
                   diff(k) = diff(k) - nint(diff(k))
                   sdif = sdif + abs(diff(k))
                 enddo

                 if(sdif < 3*tol) then
                   jsucc = j2

                   exit

                 endif
               enddo

               if(jsucc /= 0) then
                 ncount = ncount+1
                 if(i == itymin) jpoint(j1,ntr+1) = j2
               else

                 exit

               endif

             enddo

           enddo

           if(ncount == ntotal) then
             
             ntr = ntr + 1
             
             fsup(1,ntr) = frac(1)
             fsup(2,ntr) = frac(2)
             fsup(3,ntr) = frac(3)

           endif

         enddo

!        checks consistency

         do i = 1,ntype
           if(mod(natom(i),ntr) /= 0) then

             isup = -1
             write(6,*)
             write(6,*) '   WARNING   inconsistency in sym_sup_lat'
             write(6,*)

             exit

           endif
         enddo

         allocate(mark(natmin))

         if(isup == 0) then

           if(ntr > 1) then

             isup = itymin

             if(ntr == natmin) then

               nsup = 1
               jsup(1) = 1

             else

               do j = 1,natmin
                 mark(j) = 0
               enddo

               nsup = 0
               r0(1) = rat(1,1,itymin)
               r0(2) = rat(2,1,itymin)
               r0(3) = rat(3,1,itymin)

               do i = 1,natmin/ntr
                 do j1 = 1,natmin
                   if(mark(j1) == 0) then
                     nsup = nsup + 1
                     jnear = j1
                     rj(1) = rat(1,jpoint(j1,1),itymin)
                     rj(2) = rat(2,jpoint(j1,1),itymin)
                     rj(3) = rat(3,jpoint(j1,1),itymin)
                     call near_dist(distmin,adot,r0,rj)
                     do j = 1,ntr
                       rj(1) = rat(1,jpoint(j1,j),itymin)
                       rj(2) = rat(2,jpoint(j1,j),itymin)
                       rj(3) = rat(3,jpoint(j1,j),itymin)
                       call near_dist(dist,adot,r0,rj)
                       if(dist < distmin) then
                         jnear = jpoint(j1,j)
                         distmin = dist
                       endif
                     enddo

                     jsup(nsup) = jnear

                     do j = 1,ntr
                       mark(jpoint(j1,j)) = jnear
                     enddo
                   endif
                 enddo
               enddo

!               deallocate(mark)

             endif

           endif

         endif

         if(isup > 0 .and. ipr > 0) then
           write(6,*)
           write(6,'("  WARNING    there is a perfect superlattice")')
           write(6,*)
         endif

         deallocate(fsup)
         deallocate(jpoint)
         deallocate(mark)

       endif

       return

       end subroutine sym_sup_lat

