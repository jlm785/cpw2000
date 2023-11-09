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

!>  Driver for the calculation of the effective masses
!>
!>  \author       Jose Luis Martins
!>  \version      5.08
!>  \date         6 November 2023.
!>  \copyright    GNU Public License v2

subroutine cpw_pp_mass(ioreplay,                                         &
           dims_, flags_, crys_, recip_, pseudo_, atorb_,                &
           pwexp_, strfac_,  vcomp_,                                     &
           epspsi, icmax)

! adapted from cpw_pp_dos. 6 November 2023.

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

  real(real64), intent(in)                :: epspsi                      !<  accuracy of eigenvalues
  integer, intent(in)                     :: icmax                       !<  maximum number of iterations for diagonalization

! local variables

  integer            ::  imethod, ios
  integer            ::  iguess

  write(6,*)
  write(6,*) '  What method do you want to use to calculate effective masses?'
  write(6,*)
  write(6,*) '  0:  exit this section of code '
  write(6,*)
  write(6,*) '  1:  Topological tensor (accurate).'
  write(6,*) '  2:  k.p method (faster).'
  write(6,*) '  3:  Finite differences (just in case),'

  read(5,*,iostat=ios) imethod
  write(ioreplay,*) imethod,'   effective mass method'

  if(ios /= 0) then
    imethod = 0
    write(6,*)
    write(6,*) '  error processing default input '
    write(6,*) '  exiting this program section'
    write(6,*)
  endif

  if(imethod < 0 .or. imethod > 3) then
    imethod = 0
    write(6,*) '  invalid choice '
    write(6,*) '  exiting this program section'
    write(6,*)
  endif

  if(imethod == 1) then

    iguess = 0

    call out_mass_berry(ioreplay,                                        &
          pwexp_%emax, flags_%flgdal, flags_%flgpsd,                     &
          iguess, epspsi, icmax, pseudo_%ztot,                           &
          crys_%adot, crys_%ntype, crys_%natom, crys_%rat,               &
          recip_%ng, recip_%kgv, recip_%phase, recip_%conj,              &
          recip_%ns, recip_%inds, recip_%kmax,                           &
          recip_%indv, recip_%ek,                                        &
          strfac_%sfact, strfac_%icmplx,                                 &
          vcomp_%veff,                                                   &
          pseudo_%nq, pseudo_%delq, pseudo_%vkb, pseudo_%nkb,            &
          atorb_%latorb, atorb_%norbat, atorb_%nqwf, atorb_%delqwf,      &
          atorb_%wvfao,atorb_%lorb,                                      &
          dims_%mxdtyp, dims_%mxdatm, dims_%mxdgve, dims_%mxdnst,        &
          dims_%mxdlqp, dims_%mxdcub, dims_%mxdlao)

  else if(imethod == 2) then

    iguess = 0

    call out_mass_kdotp(ioreplay,                                        &
          pwexp_%emax, flags_%flgdal, flags_%flgpsd,                     &
          iguess, epspsi, icmax, pseudo_%ztot,                           &
          crys_%adot, crys_%ntype, crys_%natom, crys_%rat,               &
          recip_%ng, recip_%kgv, recip_%phase, recip_%conj,              &
          recip_%ns, recip_%inds, recip_%kmax,                           &
          recip_%indv, recip_%ek,                                        &
          strfac_%sfact, strfac_%icmplx,                                 &
          vcomp_%veff,                                                   &
          pseudo_%nq, pseudo_%delq, pseudo_%vkb, pseudo_%nkb,            &
          atorb_%latorb, atorb_%norbat, atorb_%nqwf, atorb_%delqwf,      &
          atorb_%wvfao,atorb_%lorb,                                      &
          dims_%mxdtyp, dims_%mxdatm, dims_%mxdgve, dims_%mxdnst,        &
          dims_%mxdlqp, dims_%mxdcub, dims_%mxdlao)

  else if(imethod == 3) then


    iguess = 0

    call out_mass_fd(ioreplay,                                           &
          pwexp_%emax, flags_%flgdal, flags_%flgpsd,                     &
          iguess, epspsi, icmax, pseudo_%ztot,                           &
          crys_%adot, crys_%ntype, crys_%natom, crys_%rat,               &
          recip_%ng, recip_%kgv, recip_%phase, recip_%conj,              &
          recip_%ns, recip_%inds, recip_%kmax,                           &
          recip_%indv, recip_%ek,                                        &
          strfac_%sfact, strfac_%icmplx,                                 &
          vcomp_%veff,                                                   &
          pseudo_%nq, pseudo_%delq, pseudo_%vkb, pseudo_%nkb,            &
          atorb_%latorb, atorb_%norbat, atorb_%nqwf, atorb_%delqwf,      &
          atorb_%wvfao,atorb_%lorb,                                      &
          dims_%mxdtyp, dims_%mxdatm, dims_%mxdgve, dims_%mxdnst,        &
          dims_%mxdlqp, dims_%mxdcub, dims_%mxdlao)

  endif

  return

end subroutine cpw_pp_mass

