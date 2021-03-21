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

!>     Sets up the structure factors at the prototype g-vectors.
!>     it also checks the phase factors by comparing them with
!>     the structure factors for all g vectors.

       subroutine structure_factor(sfact,icmplx,                         &
     & ng,kgv,phase,conj,inds,                                           &
     & ntype,natom,rat,                                                  &
     & mxdtyp,mxdatm,mxdgve,mxdnst)

!      adapted from Sverre Froyen plane wave program
!      written May 29 1987. jlm
!      modified July 31 1987. jlm
!      modified (conj) 20 March 1999. jlm
!      Modified documentation, January 2020. JLM
!      copyright INESC-MN/Jose Luis Martins

!      version 4.94


       implicit none

       integer, parameter          :: REAL64 = selected_real_kind(12)

!      input

       integer, intent(in)                ::  mxdtyp                     !<  array dimension of types of atoms
       integer, intent(in)                ::  mxdatm                     !<  array dimension of types of atoms
       integer, intent(in)                ::  mxdgve                     !<  array dimension for g-space vectors
       integer, intent(in)                ::  mxdnst                     !<  array dimension for g-space stars

       integer, intent(in)                ::  ntype                      !<  number of types of atoms
       integer, intent(in)                ::  natom(mxdtyp)              !<  number of atoms of type i
       real(REAL64), intent(in)           ::  rat(3,mxdatm,mxdtyp)       !<  k-th component (in lattice coordinates) of the position of the n-th atom of type i

       integer, intent(in)                ::  ng                         !<  total number of g-vectors with length less than gmax
       integer, intent(in)                ::  kgv(3,mxdgve)              !<  i-th component (reciprocal lattice coordinates) of the n-th g-vector ordered by stars of increasing length
       complex(REAL64), intent(in)        ::  phase(mxdgve)              !<  phase factor of G-vector n
       real(REAL64), intent(in)           ::  conj(mxdgve)               !<  is -1 if one must take the complex conjugate of x*phase
       integer, intent(in)                ::  inds(mxdgve)               !<  star to which g-vector n belongs

!      output

       complex(REAL64), intent(out)       ::  sfact(mxdtyp,mxdnst)       !<  imaginary part of the structure factor
       integer, intent(out)               ::  icmplx                     !<  indicates if the structure factor is complex

!      local variables
       integer       ::  ierror
!       real(REAL64)  ::  errmax, strr, stri, test, fi
       real(REAL64)  ::  errmax, test, fi
       complex(REAL64)  ::  str, sftest

!      parameters
       real(REAL64), parameter :: PI = 3.14159265358979323846_REAL64
       real(REAL64), parameter :: EPS = 0.00000001_REAL64
       real(REAL64), parameter :: ZERO = 0.0_REAL64

!      counters
       integer       ::  i, nt, ja


       icmplx = 0
       ierror = 0
       errmax = ZERO

!      start loop over g-vectors

       do i = 1,ng

!        start loop over atomic types

         do nt = 1,ntype

!          compute structure factors

           str = cmplx(ZERO,ZERO,REAL64)
           do ja = 1,natom(nt)
             fi = kgv(1,i)*rat(1,ja,nt) + kgv(2,i)*rat(2,ja,nt) +        &
     &            kgv(3,i)*rat(3,ja,nt)
             fi = 2*PI*fi
             str = str + cmplx(cos(fi),-sin(fi),REAL64)
           enddo
           if (abs(aimag(str)) > eps) icmplx=1

!          save structure factor if G is prototype

           if (i == 1) then
             sfact(nt,inds(i)) = str
           else
             if (inds(i) /= inds(i-1)) sfact(nt,inds(i)) = str
           endif

           sftest = sfact(nt,inds(i))*conjg(phase(i))
           if(conj(i) < ZERO) sftest = conjg(sftest)
           test = abs(sftest -str)
           if (test > eps) ierror=ierror+1
           if (test > eps) then
             write(6,*) 'i,nt,test', i,nt,test
             write(6,*) sftest
             write(6,*) str
           endif
           if (test > errmax) errmax=test
             
           if(ierror > 5) exit
             
         enddo
             
         if(ierror > 5) exit
             
       enddo

       if(ierror > 0) then

         write(6,*)
         write(6,'("   STOPPED in structure_factor.  The structure ",    &
     &       "and phase factors do not agree.")')
         write(6,'("   Max error = ",e14.6)') errmax

         stop

       endif

       return
       end subroutine structure_factor
