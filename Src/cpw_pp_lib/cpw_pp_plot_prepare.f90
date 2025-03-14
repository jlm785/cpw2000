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

!>  Converts the charge densities for subsequent plotting
!>
!>  \author       Jose Luis Martins
!>  \version      5.11
!>  \date         13 March 2025.
!>  \copyright    GNU Public License v2


subroutine cpw_pp_plot_prepare(dims_, recip_, vcomp_, chdens_,           &
    crys_, strfac_,  pseudo_,                                            &
    dims_in_, recip_in_, vcomp_in_, chdens_in_)

! written 13 March 2025. JLM


  use cpw_variables

  implicit none

!  integer, parameter          :: REAL64 = selected_real_kind(12)

  type(dims_t)                       ::  dims_                           !<  array dimensions

  type(crys_t)                       ::  crys_                           !<  crystal structure
  type(pseudo_t)                     ::  pseudo_                         !<  pseudo-potential (Kleinman-Bylander)
  type(recip_t)                      ::  recip_                          !<  reciprocal space information
!  type(spaceg_t)                     ::  spaceg_                         !<  space group information
!  type(pwexp_t)                      ::  pwexp_                          !<  plane-wave expansion choices
  type(strfac_t)                     ::  strfac_                         !<  structure factors
  type(chdens_t)                     ::  chdens_                         !  charge densities
  type(vcomp_t)                      ::  vcomp_                          !<  local potential contributions

  type(dims_t)                       ::  dims_in_                        !<  input array dimensions

  type(recip_t)                      ::  recip_in_                       !<  input reciprocal space information
  type(chdens_t)                     ::  chdens_in_                      !  charge densities
  type(vcomp_t)                      ::  vcomp_in_                       !<  input local potential contributions

! local variables

  real(REAL64)    ::  vcell, bdot(3,3)
  real(REAL64)    ::  fac
  real(REAL64)    ::  ealpha

! counters

  integer       ::  i, j, n, m

! parameters

  real(REAL64), parameter     ::  PI = 3.14159265358979323846_REAL64
  real(REAL64), parameter     ::  ZERO = 0.0_REAL64
  complex(REAL64), parameter  ::  C_ZERO = cmplx(ZERO,ZERO,REAL64)



  allocate(chdens_%den(dims_%mxdnst))
  allocate(chdens_%dend(dims_%mxdnst))

! converts density

  call cpw_pp_convert(chdens_%den, recip_%kmax, chdens_in_%den,          &
      recip_%ng, recip_%kgv, recip_%phase, recip_%conj, recip_%ns,       &
      recip_%mstar,                                                      &
      recip_in_%kgv, recip_in_%phase, recip_in_%conj, recip_in_%ns,      &
      recip_in_%mstar,                                                   &
      dims_%mxdgve, dims_%mxdnst, dims_in_%mxdgve, dims_in_%mxdnst)

! converts bonding density

  call cpw_pp_convert(chdens_%dend, recip_%kmax, chdens_in_%dend,        &
      recip_%ng, recip_%kgv, recip_%phase, recip_%conj, recip_%ns,       &
      recip_%mstar,                                                      &
      recip_in_%kgv, recip_in_%phase, recip_in_%conj, recip_in_%ns,      &
      recip_in_%mstar,                                                   &
      dims_%mxdgve, dims_%mxdnst, dims_in_%mxdgve, dims_in_%mxdnst)

! atomic terms

  allocate(recip_in_%ek(dims_in_%mxdnst))
  allocate(vcomp_in_%vion(dims_in_%mxdnst))
  allocate(chdens_in_%denc(dims_in_%mxdnst))
  allocate(chdens_in_%dens(dims_in_%mxdnst))
  allocate(pseudo_%vql(dims_%mxdtyp,dims_in_%mxdnst))
  allocate(pseudo_%dvql(dims_in_%mxdnst))
  allocate(pseudo_%dnc(dims_%mxdtyp,dims_in_%mxdnst))
  allocate(pseudo_%ddc(dims_in_%mxdnst))

  call adot_to_bdot(crys_%adot, vcell, bdot)

  recip_in_%ek(:) = ZERO
  m = 1
  do n = 1,recip_in_%ns
    do i = 1,3
    do j = 1,3
      recip_in_%ek(n) = recip_in_%ek(n) + recip_in_%kgv(i,m)*bdot(i,j)*recip_in_%kgv(j,m)
    enddo
    enddo
    recip_in_%ek(n) = recip_in_%ek(n)/2
    m = m + recip_in_%mstar(n)
  enddo

  call v_first(recip_in_%ns, recip_in_%ek, strfac_%sfact, ealpha,        &
      pseudo_%ealraw, pseudo_%nq, pseudo_%delq, pseudo_%vloc,            &
      pseudo_%dcor, pseudo_%dval,                                        &
      crys_%ntype, crys_%adot,                                           &
      vcomp_in_%vion, chdens_in_%denc, chdens_in_%dens,                  &
      pseudo_%vql, pseudo_%dvql, pseudo_%dnc, pseudo_%ddc,               &
      dims_%mxdtyp, dims_%mxdlqp, dims_in_%mxdnst)

! Hartree potential

  allocate(vcomp_%vhar(dims_%mxdnst))

  call adot_to_bdot(crys_%adot, vcell, bdot)

  vcomp_%vhar(:) = C_ZERO
  fac = 2*PI/vcell
  do i = 2,min(recip_in_%ns,recip_%ns)
    vcomp_%vhar(i) = fac*chdens_%den(i)/recip_in_%ek(i)
  enddo

! converts local atomic potential

  allocate(vcomp_%vion(dims_%mxdnst))

  call cpw_pp_convert(vcomp_%vion, recip_%kmax, vcomp_in_%vion,          &
      recip_%ng, recip_%kgv, recip_%phase, recip_%conj, recip_%ns,       &
      recip_%mstar,                                                      &
      recip_in_%kgv, recip_in_%phase, recip_in_%conj, recip_in_%ns,      &
      recip_in_%mstar,                                                   &
      dims_%mxdgve, dims_%mxdnst, dims_in_%mxdgve, dims_in_%mxdnst)

! Exchange correlation potential

  allocate(vcomp_%vxc(dims_%mxdnst))

  do i = 1,recip_%ns
    vcomp_%vxc(i) = vcomp_%veff(i) - vcomp_%vhar(i) - vcomp_%vion(i)
  enddo

! Converts core charge

  allocate(chdens_%denc(dims_%mxdnst))

  call cpw_pp_convert(chdens_%denc, recip_%kmax, chdens_in_%denc,        &
      recip_%ng, recip_%kgv, recip_%phase, recip_%conj, recip_%ns,       &
      recip_%mstar,                                                      &
      recip_in_%kgv, recip_in_%phase, recip_in_%conj, recip_in_%ns,      &
      recip_in_%mstar,                                                   &
      dims_%mxdgve, dims_%mxdnst, dims_in_%mxdgve, dims_in_%mxdnst)

! Converts sum of spherical atomic charges

  allocate(chdens_%dens(dims_%mxdnst))

  call cpw_pp_convert(chdens_%dens, recip_%kmax, chdens_in_%dens,        &
      recip_%ng, recip_%kgv, recip_%phase, recip_%conj, recip_%ns,       &
      recip_%mstar,                                                      &
      recip_in_%kgv, recip_in_%phase, recip_in_%conj, recip_in_%ns,      &
      recip_in_%mstar,                                                   &
      dims_%mxdgve, dims_%mxdnst, dims_in_%mxdgve, dims_in_%mxdnst)

  return

end subroutine cpw_pp_plot_prepare

