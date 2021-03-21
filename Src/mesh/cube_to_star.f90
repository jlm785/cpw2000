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

!>     Consolidates information from the whole G-space
!>     into a star with full symetry.

       subroutine cube_to_star(den, nsave, chdsave,                      &
     & ng, kgv, phase, conj, ns, mstar,                                  &
     & mxdgve,mxdnst)

!      Written April 22, 1999. jlm
!      modified July 24, 2002. jlm
!      modified December 19, 2013, denr(1). jlm
!      modified January 8, 2014, f90. jlm
!      Modified, 12 August 2019, dimension of den, documentation. JLM
!      Modified 12 December 2019.  Paranoid check.  JLM
!      Copyright INESC-MN/Jose Luis Martins

!      version 4.94

       implicit none
       integer, parameter          :: REAL64 = selected_real_kind(12)


!      input

       integer, intent(in)                ::  mxdgve                     !<  array dimension for g-space vectors
       integer, intent(in)                ::  mxdnst                     !<  array dimension for g-space stars

       integer, intent(in)                ::  nsave(3)                   !<  dimensions of chdsave

       complex(REAL64), intent(in)        ::                             &
     & chdsave(-nsave(1):nsave(1),-nsave(2):nsave(2),-nsave(3):nsave(3)) !<  quantity in reciprocal point i,j,k

       integer, intent(in)                ::  ng                         !<  total number of g-vectors with length less than gmax
       integer, intent(in)                ::  kgv(3,mxdgve)              !<  i-th component (reciprocal lattice coordinates) of the n-th g-vector ordered by stars of increasing length
       complex(REAL64), intent(in)        ::  phase(mxdgve)              !<  imaginary part of the phase factor of G-vector n
       real(REAL64), intent(in)           ::  conj(mxdgve)               !<  is -1 if one must take the complex conjugate of x*phase
       integer, intent(in)                ::  ns                         !<  number os stars with length less than gmax
       integer, intent(in)                ::  mstar(mxdnst)              !<  number of g-vectors in the j-th star

!      output

       complex(REAL64), intent(out)       ::  den(mxdnst)                !<  quantity for the prototype g-vector in star

!      local variables

       integer       ::   nsmax, ngmax, jmax, jvec
       complex(REAL64)  ::   denu(48)
       logical       ::   lout

!      parameters

       real(REAL64), parameter :: ZERO = 0.0_REAL64
       complex(REAL64), parameter :: C_ZERO = cmplx(ZERO,ZERO,REAL64)

!      counters

       integer       ::  i, j, k1, k2, k3


!      paranoid check, stops complaint about unused ng...

       if(ng > mxdgve) stop

       do i=1,ns
         den(i) = C_ZERO
       enddo

       nsmax = 0
       ngmax = 0
       lout = .FALSE.
       do i=1,ns

!        paranoid check
         if(mstar(i) > 48) then
           write(6,'("   stopped in mesh_cube_to_star,  mstar(",i7,      &
     &             ") = ")') i,mstar(i)

           stop

         endif

         do j=ngmax+1,ngmax+mstar(i)
           if (iabs(kgv(1,j)) > nsave(1)) lout = .TRUE.
           if (iabs(kgv(2,j)) > nsave(2)) lout = .TRUE.
           if (iabs(kgv(3,j)) > nsave(3)) lout = .TRUE.
         enddo

         if(lout) exit

         nsmax = i
         ngmax = ngmax + mstar(i)
       enddo


       den(1) = chdsave(0,0,0)

       if(nsmax == 1) return

       jmax = 1
       do i=2,nsmax
         do j=1,mstar(i)
           k1 = kgv(1,jmax+j)
           k2 = kgv(2,jmax+j)
           k3 = kgv(3,jmax+j)
           denu(j) = chdsave(k1,k2,k3)
         enddo
         do j=1,mstar(i)
           jvec = jmax + j
           if(conj(jvec) > ZERO) then
             den(i) = den(i) + denu(j)*phase(jvec) / mstar(i)
           else
             den(i) = den(i) + conjg(denu(j))*phase(jvec) / mstar(i)
           endif
         enddo
         jmax = jmax + mstar(i)
       enddo

       return
       end subroutine cube_to_star
