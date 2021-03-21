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

!>  Calculates the energy according to the 
!>  Harris-Weinert-Foulkes functional

subroutine harris_weinert_foulkes(eharrfou, eband, exc, enerew,          &
  ns, mstar, ek,                                                         &
  vhxc, den,                                                             &
  adot,                                                                  &
  mxdnst)

! Adapted from total_ks_energy
! written 2 December 2020. JLM
! copyright INESC-MN/Jose Luis Martins

! Version 4.99

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  mxdnst                          !<  array dimension for g-space stars

  real(REAL64), intent(in)           ::  eband                           !<  integrated band energy.
  real(REAL64), intent(in)           ::  exc                             !<  exchange+correlation energy (previous run)
  real(REAL64), intent(in)           ::  enerew                          !<  Ewald energy

  integer, intent(in)                ::  ns                              !<  number os stars with length less than gmax
  integer, intent(in)                ::  mstar(mxdnst)                   !<  number of g-vectors in the j-th star
  real(REAL64), intent(in)           ::  ek(mxdnst)                      !<  kinetic energy (hartree) of g-vectors in star j

  complex(REAL64), intent(in)        ::  vhxc(mxdnst)                    !<  input hxc potential (hartree) for the prototype G-vector 
  complex(REAL64), intent(in)        ::  den(mxdnst)                     !<  density for the prototype G-vector (old density!)

  real(REAL64), intent(in)           ::  adot(3,3)                       !<  metric in direct space

! output

  real(REAL64), intent(out)          ::  eharrfou                        !<  total energy of the Harris-Weinert-Foulkes functional

! local variables

  real(REAL64)        ::  evhxc, ehart
  real(REAL64)        ::  vcell, bdot(3,3)

! counters

  integer       ::  i

! parameters

  real(REAL64), parameter :: PI = 3.14159265358979323846_REAL64
  real(REAL64), parameter :: ZERO = 0.0_REAL64


  call adot_to_bdot(adot,vcell,bdot)

! initialize hxc-correction, local pseudo, hartree (g>0) energies
! includes correction for average potential

  evhxc = real(vhxc(1),REAL64)*real(den(1),REAL64)
  ehart = ZERO

! start loop over stars

  do i=2,ns

!   Compute the hxc correction from the potential of the
!   valence electrons. This will be subtracted from the sum
!   of the eigenvalues.
!   Then compute the contribution to the total energy from
!   the interaction with the ions and the electrostatic
!   (Hartree) interaction between the electrons.

    evhxc = evhxc + mstar(i)*real(vhxc(i)*conjg(den(i)),REAL64)
    ehart = ehart + mstar(i)*real(den(i)*conjg(den(i)),REAL64) / (2*ek(i))

  enddo

  ehart = 4*PI*ehart/(2*vcell)

! the total energy. modified for average potential

  eharrfou = eband - evhxc + ehart + exc + enerew

  return
end subroutine harris_weinert_foulkes
