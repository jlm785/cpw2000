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

!>  Driver for the calculation of quantum geometric properties
!>  at chosen k-points.
!>
!>  \author       Jose Luis Martins
!>  \version      5.11
!>  \date         4 April 2024.
!>  \copyright    GNU Public License v2

subroutine cpw_pp_qgeom(ioreplay,                                        &
           dims_, flags_, crys_, recip_, pseudo_, atorb_,                &
           pwexp_, strfac_,  vcomp_,                                     &
           efermi, epspsi, icmax)

! adapted from cpw_pp_mass. 4 April 2024.

  use cpw_variables

  implicit none

!  integer, parameter          :: REAL64 = selected_real_kind(12)

  integer, intent(in)                     ::  ioreplay                   !<  tape number for reproducing calculations

  type(dims_t)                            ::  dims_                      !<  array dimensions
  type(flags_t)                           ::  flags_                     !<  computational flags
  type(crys_t)                            ::  crys_                      !<  crystal structure
  type(recip_t)                           ::  recip_                     !<  reciprocal space information
  type(pseudo_t)                          ::  pseudo_                    !<  pseudo-potential (Kleinman-Bylander)
  type(atorb_t)                           ::  atorb_                     !<  atomic orbitals in G-space
  type(pwexp_t)                           ::  pwexp_                     !<  plane-wave expansion choices
  type(strfac_t)                          ::  strfac_                    !<  structure factors
  type(vcomp_t)                           ::  vcomp_                     !<  local potential contributions

 ! other input

  real(real64), intent(in)                ::  epspsi                     !<  accuracy of eigenvalues
  integer, intent(in)                     ::  icmax                      !<  maximum number of iterations for diagonalization
  real(REAL64), intent(in)                ::  efermi                     !<  eigenvalue of highest occupied state (T=0) or fermi energy (T/=0), Hartree

! local variables

  integer            ::  iguess


  iguess = 0

  call out_qgeom(ioreplay,                                               &
        pwexp_%emax, flags_%flgdal, flags_%flgpsd,                       &
        iguess, epspsi, icmax, pseudo_%ztot, efermi,                     &
        crys_%adot, crys_%ntype, crys_%natom, crys_%rat,                 &
        recip_%ng, recip_%kgv, recip_%phase, recip_%conj,                &
        recip_%ns, recip_%inds, recip_%kmax,                             &
        recip_%indv, recip_%ek,                                          &
        strfac_%sfact, strfac_%icmplx,                                   &
        vcomp_%veff,                                                     &
        pseudo_%nq, pseudo_%delq, pseudo_%vkb, pseudo_%nkb,              &
        atorb_%latorb, atorb_%norbat, atorb_%nqwf, atorb_%delqwf,        &
        atorb_%wvfao,atorb_%lorb,                                        &
        dims_%mxdtyp, dims_%mxdatm, dims_%mxdgve, dims_%mxdnst,          &
        dims_%mxdlqp, dims_%mxdcub, dims_%mxdlao)


  return

end subroutine cpw_pp_qgeom

