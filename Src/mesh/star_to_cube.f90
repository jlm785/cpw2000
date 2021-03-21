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

!>     distribute information of a star to the whole g-space

       subroutine star_to_cube(den, nsave, chdsave,                      &
     & ng, kgv, phase, conj, inds, ns, mstar,                            &
     & mxdgve, mxdnst)

!      written april 22, 1999. jlm
!      modified july 24, 2002. jlm
!      modified October 21, 2015, f90. JLM
!      Modified, 12 August 2019, dimension of den, documentation. JLM
!      Copyright INESC-MN/Jose Luis Martins

!      version 4.94

       implicit none
       integer, parameter          :: REAL64 = selected_real_kind(12)


!      input

       integer, intent(in)                ::  mxdgve                     !<  array dimension for g-space vectors
       integer, intent(in)                ::  mxdnst                     !<  array dimension for g-space stars

       integer, intent(in)                ::  nsave(3)                   !<  dimensions of chdsave

       complex(REAL64), intent(in)        ::  den(mxdnst)                !<  quantity for the prototype g-vector in star

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

!      local variables

       integer       ::   ngmax
       logical       ::   lout

!      parameters

       real(REAL64), parameter :: ZERO = 0.0_REAL64
       complex(REAL64), parameter  ::  C_ZERO = cmplx(ZERO,ZERO,REAL64)

!      counters

       integer       ::  i, j, k1, k2, k3


!      paranoid check, stops complaint about unused ng...

       if(ng > mxdgve) stop

       do k3 = -nsave(3),nsave(3)
       do k2 = -nsave(2),nsave(2)
       do k1 = -nsave(1),nsave(1)
         chdsave(k1,k2,k3) = C_ZERO
       enddo
       enddo
       enddo

       ngmax = 0
       lout = .FALSE.
       do i=1,ns
         do j=ngmax+1,ngmax+mstar(i)
           if (iabs(kgv(1,j)) > nsave(1)) lout = .TRUE.
           if (iabs(kgv(2,j)) > nsave(2)) lout = .TRUE.
           if (iabs(kgv(3,j)) > nsave(3)) lout = .TRUE.
         enddo

         if(lout) exit

         ngmax = ngmax + mstar(i)
       enddo

       do i=1,ngmax
         k1 = kgv(1,i)
         k2 = kgv(2,i)
         k3 = kgv(3,i)
         if(conj(i) > ZERO) then
           chdsave(k1,k2,k3) = den(inds(i))*conjg(phase(i))
         else
           chdsave(k1,k2,k3) = conjg(den(inds(i)))*phase(i)
         endif
       enddo

       return
       end subroutine star_to_cube
