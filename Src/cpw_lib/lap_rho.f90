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

!>     calculates the laplacian on the prototype G-vectors
!>     of the density (or another quantity that is real in real space) 
!>     represented in the prototype G-vectors 

       subroutine lap_rho(rho,rholap,adot,                               &
     & ng,kgv,phase,conj,ns,inds,mstar,                                  &
     & mxdgve,mxdnst)

!      Written September 30, 2015 from v_hartree_xc
!      Modified October 7, 2015, output on stars. JLM
!      Modified, documentation, January 2020. JLM
!      copyright INESC-MN/Jose Luis Martins

!      version 4.94


       implicit none

       integer, parameter          :: REAL64 = selected_real_kind(12)

!      input

       integer, intent(in)                ::  mxdgve                     !<  array dimension for g-space vectors
       integer, intent(in)                ::  mxdnst                     !<  array dimension for g-space stars

       complex(REAL64), intent(in)        ::  rho(mxdnst)                !<  density or other quantity in prototype G-vector
       real(REAL64), intent(in)           ::  adot(3,3)                  !<  metric in direct space
      
       integer, intent(in)                ::  ng                         !<  size of g-space
       integer, intent(in)                ::  kgv(3,mxdgve)              !<  G-vectors in reciprocal lattice coordinates 
       complex(REAL64), intent(in)        ::  phase(mxdgve)              !<  phase factor of G-vector n
       real(REAL64), intent(in)           ::  conj(mxdgve)               !<  is -1 if one must take the complex conjugate of x*phase
       integer, intent(in)                ::  ns                         !<  number os stars with length less than gmax
       integer, intent(in)                ::  inds(mxdgve)               !<  star to which g-vector n belongs
       integer, intent(in)                ::  mstar(mxdnst)              !<  number of g-vectors in the j-th star
 
!      output

       complex(REAL64), intent(out)       ::  rholap(mxdnst)             !<  density or other quantity in prototype G-vector

!      local allocatable arrays

       complex(REAL64), allocatable       ::  deng(:)                    !  density (in G)

!      local variables

       real(REAL64)    ::  vcell,bdot(3,3)
       real(REAL64)    ::  xp

!      counters

       integer         ::  i, j ,k

!      parameters

       real(REAL64), parameter :: ZERO = 0.0_REAL64


       call adot_to_bdot(adot,vcell,bdot)


       allocate(deng(mxdgve))

!      unfolds from prototype-G

       call star_of_g_unfold(deng,rho,.FALSE.,                           &
     & ng,phase,conj,inds,                                               &
     & mxdgve,mxdnst)

!      calculates laplacian

       do i = 1,ng
         xp = ZERO
         do j = 1,3
         do k = 1,3
           xp = xp + kgv(j,i)*bdot(j,k)*kgv(k,i)
         enddo
         enddo
         deng(i) = -xp*deng(i)
       enddo

!      symmetrizes by folding into G-prototype

       call star_of_g_fold(rholap,deng,.FALSE.,                          &
     & ng,phase,conj,ns,inds,mstar,                                      &
     & mxdgve,mxdnst)
         

       deallocate(deng)

       return
       end subroutine lap_rho
