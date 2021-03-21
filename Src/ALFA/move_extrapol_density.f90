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

!>     extrapolation of the charge density for the first
!>     self consistent cycle

       subroutine move_extrapol_density(newcalc, flgcal, nsave, chdsave, &
     & dens, dend, dend1, den,                                           &
     & ng, kgv, phase, conj, ns, mstar,                                  &
     & mxdgve,mxdnst)

!      written 18 September 2002. jlm
!      modified 18 February 2008. jlm
!      modified  8 May 2014. epilbf. jlm
!      modified October 20, 2015, f90. JLM
!      modified documentation, added EPILNG, August 2019. JLM
!      Copyright INESC-MN/Jose Luis Martins

!      version 4.94

       implicit none
       integer, parameter          :: REAL64 = selected_real_kind(12)


!      input

       integer, intent(in)                ::  mxdgve                     !<  array dimension for g-space vectors
       integer, intent(in)                ::  mxdnst                     !<  array dimension for g-space stars

       integer, intent(in)                ::  nsave(3)                   !<  dimensions of chdsave

       logical, intent(in)                ::  newcalc                    !<  indicates that it is a new calculation
       character(len=6), intent(in)       ::  flgcal                     !<  type of calculation

       complex(REAL64), intent(in)        ::  dens(mxdnst)               !<  spherical atomic valence charge density for the prototype G-vector in star j
       complex(REAL64), intent(in)        ::  dend1(mxdnst)              !<  bonding charge density from second previous md step

       integer, intent(in)                ::  ng                         !<  total number of g-vectors with length less than gmax
       integer, intent(in)                ::  kgv(3,mxdgve)              !<  i-th component (reciprocal lattice coordinates) of the n-th g-vector ordered by stars of increasing length
       complex(REAL64), intent(in)        ::  phase(mxdgve)              !<  imaginary part of the phase factor of G-vector n
       real(REAL64), intent(in)           ::  conj(mxdgve)               !<  is -1 if one must take the complex conjugate of x*phase
       integer, intent(in)                ::  ns                         !<  number os stars with length less than gmax
       integer, intent(in)                ::  mstar(mxdnst)              !<  number of g-vectors in the j-th star

       complex(REAL64), intent(in)        ::                             &
     & chdsave(-nsave(1):nsave(1),-nsave(2):nsave(2),-nsave(3):nsave(3)) !<  quantity in reciprocal point i,j,k
      
!      output

       complex(REAL64), intent(out)       ::  den(mxdnst)                !<  density for the prototype G-vector

!      input and output

       complex(REAL64), intent(inout)     ::  dend(mxdnst)               !<  bonding charge density from previous md step 

!      constants

       real(REAL64), parameter  :: ZERO = 0.0_REAL64
       complex(REAL64), parameter  ::  C_ZERO = cmplx(ZERO,ZERO,REAL64)
       
!      counters

       integer       ::  i


       if(newcalc) then

!        first time in this run

         do i=1,ns
           den(i) = dens(i)
           dend(i) = C_ZERO
         enddo

       else

!        previously calculated densities are available

         if(flgcal == 'VCSLNG' .or. flgcal == 'EPILBF' .or.              &
     &      flgcal == 'VCSLBF' .or. flgcal == 'VCSMIC' .or.              &
     &      flgcal == 'EPILNG') then

           call cube_to_star(dend, nsave,chdsave,                        &
     &     ng, kgv, phase, conj, ns, mstar,                              &
     &     mxdgve, mxdnst)

           do i=1,ns
             den(i) = dens(i) + dend(i)
           enddo

         else

           do i=1,ns
             den(i) = dens(i) + 2*dend(i)-dend1(i)
           enddo

         endif
       endif

       return
       end subroutine move_extrapol_density
