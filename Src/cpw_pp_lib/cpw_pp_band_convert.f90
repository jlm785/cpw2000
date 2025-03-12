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

!>  Converts the effective potential or density from one g-space
!>  representation (kgv_in,....) to another (kgv,...).
!>
!>  \author       Jose Luis Martins
!>  \version      5.11
!>  \date         12 March 2025.
!>  \copyright    GNU Public License v2

subroutine cpw_pp_convert(veff, kmax, veff_in,                           &
      ng, kgv, phase, conj, ns, mstar,                                   &
      kgv_in, phase_in, conj_in, ns_in, mstar_in,                        &
      mxdgve, mxdnst, mxdgve_in, mxdnst_in)

! Extracted from cpw_pp_band_prepare. 12 March 2025. JLM

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)


! input

  integer, intent(in)                ::  mxdgve                          !<  array dimension for g-space vectors
  integer, intent(in)                ::  mxdnst                          !<  array dimension for g-space stars

  integer, intent(in)                ::  mxdgve_in                       !<  array dimension for g-space vectors
  integer, intent(in)                ::  mxdnst_in                       !<  array dimension for g-space stars

  integer, intent(in)                ::  kmax(3)                         !<  dimensions of chd

  integer, intent(in)                ::  ng                              !<  total number of g-vectors with length less than gmax
  integer, intent(in)                ::  kgv(3,mxdgve)                   !<  i-th component (reciprocal lattice coordinates) of the n-th g-vector ordered by stars of increasing length
  complex(REAL64), intent(in)        ::  phase(mxdgve)                   !<  imaginary part of the phase factor of G-vector n
  real(REAL64), intent(in)           ::  conj(mxdgve)                    !<  is -1 if one must take the complex conjugate of x*phase
  integer, intent(in)                ::  ns                              !<  number os stars with length less than gmax
  integer, intent(in)                ::  mstar(mxdnst)                   !<  number of g-vectors in the j-th star

  complex(REAL64), intent(in)        ::  veff_in(mxdnst)                 !<  effective potential or charge dinsity for the prototype g-vector in star

  integer, intent(in)                ::  kgv_in(3,mxdgve_in)             !<  i-th component (reciprocal lattice coordinates) of the n-th g-vector ordered by stars of increasing length
  complex(REAL64), intent(in)        ::  phase_in(mxdgve_in)             !<  imaginary part of the phase factor of G-vector n
  real(REAL64), intent(in)           ::  conj_in(mxdgve_in)              !<  is -1 if one must take the complex conjugate of x*phase
  integer, intent(in)                ::  ns_in                           !<  number os stars with length less than gmax
  integer, intent(in)                ::  mstar_in(mxdnst_in)             !<  number of g-vectors in the j-th star

! output

  complex(REAL64), intent(out)       ::  veff(mxdnst)                    !<  effective potential or charge dinsity for the prototype g-vector in star

! allocatable local arrays

  complex(REAL64), allocatable       ::  chd(:,:,:)                      !  effective potential on the grid

! local variables

  integer    ::  ngmax

! counters

  integer    ::  i,j,k

! constants

  real(REAL64), parameter    :: ZERO = 0.0_REAL64
  complex(REAL64), parameter :: C_ZERO = cmplx(ZERO,ZERO,REAL64)

  allocate(chd(-kmax(1):kmax(1),-kmax(2):kmax(2),-kmax(3):kmax(3)))

  do k = -kmax(3),kmax(3)
  do j = -kmax(2),kmax(2)
  do i = -kmax(1),kmax(1)
    chd(i,j,k) = C_ZERO
  enddo
  enddo
  enddo

  veff(:) = C_ZERO

  ngmax = 0
  do i = 1,min(ns,ns_in)
    do j = ngmax+1,ngmax+mstar_in(i)
      if(abs(kgv_in(1,j)) <= kmax(1) .and.                 &
         abs(kgv_in(2,j)) <= kmax(2) .and.                 &
         abs(kgv_in(3,j)) <= kmax(3) ) then

        if(conj_in(j) > ZERO) then
          chd(kgv_in(1,j),kgv_in(2,j),kgv_in(3,j)) =      &
                        veff_in(i)*conjg(phase_in(j))
        else
          chd(kgv_in(1,j),kgv_in(2,j),kgv_in(3,j)) =      &
                        conjg(veff_in(i))*phase_in(j)
        endif

      endif
    enddo
    ngmax = ngmax + mstar_in(i)

  enddo


! collects v_effective in stars

  call cube_to_star(veff, kmax, chd,                                     &
      ng, kgv, phase, conj, ns, mstar,                                   &
      mxdgve, mxdnst)

  deallocate(chd)

  return

end subroutine cpw_pp_convert

