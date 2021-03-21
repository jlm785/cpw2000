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

!>     saves density for future extrapolation

       subroutine move_save_density(newcalc, flgcal, nsave, chdsave,     &
     & den, dens, dend, dend1,                                           &
     & ng, kgv, phase, conj, inds, ns, mstar,                            &
     & mxdgve, mxdnst)

!      written 18 september 2002. jlm
!      modified 18 february 2008. jlm
!      modified  8 may 2014. epilbf. jlm
!      modified October 20, 2015, f90. JLM
!      modified January 2019, newcalc intent. JLM
!      Modified, documentation, August 2019, EPILNG. JLM
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
       complex(REAL64), intent(in)        ::  den(mxdnst)                !<  density for the prototype G-vector

       integer, intent(in)                ::  ng                         !<  total number of g-vectors with length less than gmax
       integer, intent(in)                ::  kgv(3,mxdgve)              !<  i-th component (reciprocal lattice coordinates) of the n-th g-vector ordered by stars of increasing length
       complex(REAL64), intent(in)        ::  phase(mxdgve)              !<  imaginary part of the phase factor of G-vector n
       real(REAL64), intent(in)           ::  conj(mxdgve)               !<  is -1 if one must take the complex conjugate of x*phase
       integer, intent(in)                ::  inds(mxdgve)               !<  star to which g-vector n belongs
       integer, intent(in)                ::  ns                         !<  number os stars with length less than gmax
       integer, intent(in)                ::  mstar(mxdnst)              !<  number of g-vectors in the j-th star
       
!      output

       complex(REAL64), intent(out)       ::                             &
     & chdsave(-nsave(1):nsave(1),-nsave(2):nsave(2),-nsave(3):nsave(3)) !<  quantity in reciprocal point i,j,k

       complex(REAL64), intent(out)       ::  dend1(mxdnst)              !<  bonding charge density from second previous md step

!      input and output

       complex(REAL64), intent(inout)     ::  dend(mxdnst)               !<  bonding charge density from previous md step 

!      local variables

       complex(REAL64)  ::  dd
       integer          ::  ic

!      constants

       real(REAL64), parameter  :: ZERO = 0.0_REAL64
       complex(REAL64), parameter  ::  C_ZERO = cmplx(ZERO,ZERO,REAL64)
       
!      counters

       integer          ::  i


!      stores previous bonding charge

       if( .not. newcalc) then
         do i=1,ns
           dend1(i) = dend(i)
         enddo
       endif

!      calculates current bonding charge

       do i=1,ns
         dend(i) = den(i) - dens(i)
       enddo

       if(flgcal == 'VCSLNG' .or. flgcal == 'EPILBF' .or.                &
     &    flgcal == 'VCSLBF' .or. flgcal == 'VCSMIC' .or.                &
     &    flgcal == 'EPILNG') then

         call star_to_cube(dend, nsave, chdsave,                         &
     &   ng, kgv, phase, conj, inds, ns, mstar,                          &
     &   mxdgve,mxdnst)

       else

!        makes sure the update is real

         ic = 1

         dend(1) = C_ZERO
         do i = 2,ns
           if(conj(ic+2) > ZERO) then
             dd = dend(i)*conjg(phase(ic+2))
             dend(i) = (dend(i) + conjg(dd)) / 2
           endif
           ic = ic + mstar(i)
         enddo

!        first time or l-bfgs minimization, next extrapolation
!        should be simpler.

         if(newcalc .or. flgcal == 'LBFSYM') then
           do i=1,ns
             dend1(i) = dend(i)
           enddo
         endif

       endif

       return
       end subroutine move_save_density
