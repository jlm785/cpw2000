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

!>  Prepares the data for the self-consistency iterations
!>  Initial potential etc...
!>
!>  \author       José Luís Martins
!>  \version      5.12
!>  \date         around 2020, 24 November 2025.
!>  \copyright    GNU Public License v2

subroutine cpw_scf_prepare(ealpha,iprglob,newcalc,                       &
     nsave, chdsave,                                                     &
     exc,strxc,rhovxc,                                                   &
     dims_,crys_,recip_,strfac_,pseudo_,chdens_,vcomp_,flags_,           &
     ewald_,xc_)

! Adapted from the code without cpw_variables. around 2020. JLM
! New variable for v_Hartree_xc. Indentation. 24 November 2025. JLM

  use cpw_variables

  implicit none

  type(dims_t)                       ::  dims_                           !<  array dimensions
  type(crys_t)                       ::  crys_                           !<  crystal structure
  type(flags_t)                      ::  flags_                          !<  computational flags
  type(recip_t)                      ::  recip_                          !<  reciprocal space information
  type(strfac_t)                     ::  strfac_                         !<  structure factors
  type(chdens_t)                     ::  chdens_                         !<  charge densities
  type(vcomp_t)                      ::  vcomp_                          !<  Componemts of local potential
  type(pseudo_t)                     ::  pseudo_                         !<  pseudo-potential (Kleinman-Bylander)
  type(enfrst_t)                     ::  ewald_                          !<  Ewald energy force stress
  type(xc_t)                         ::  xc_                             !<  exchange and correlation choice

  integer, intent(in)                ::  iprglob                         !<  controls the amount of printing by subroutines
  real(REAL64), intent(out)          ::  ealpha                          !<  G=0 contribution to the total energy (Hartree)
  logical, intent(in)                ::  newcalc                         !<  indicates that it is a new calculation (equivalent iter = ist0+1)

  integer, intent(in)                ::  nsave(3)                        !<  dimensions of chdsave
  complex(REAL64), intent(inout)     ::                                  &
       chdsave(-nsave(1):nsave(1),-nsave(2):nsave(2),-nsave(3):nsave(3)) !<  quantity in reciprocal point i,j,k

  real(REAL64), intent(out)          ::  exc                             !<  Exchange and correlation energy (Hartree)
  real(REAL64), intent(out)          ::  rhovxc                          !<  Integral of rho times vxc (Hartree)
  real(REAL64), intent(out)          ::  strxc(3,3)                      !<  d exc / d adot,  exchange and correlatio contribution to stress tensor (contravariant components)

  integer      ::  ipr

  complex(REAL64), allocatable        ::  rholap(:)                      !  Laplacian of charge density
  complex(REAL64), allocatable        ::  tau(:)                         !  Kinetic energy density (whatever that means)

  integer       ::  i

  real(REAL64), parameter  ::  ZERO = 0.0D0
  complex(REAL64), parameter  ::  C_ZERO = cmplx(ZERO,ZERO,REAL64)

! calculates the structure factors

  call structure_factor(strfac_%sfact, strfac_%icmplx,                   &
      recip_%ng, recip_%kgv, recip_%phase, recip_%conj, recip_%inds,     &
      crys_%ntype, crys_%natom, crys_%rat,                               &
      dims_%mxdtyp, dims_%mxdatm, dims_%mxdgve, dims_%mxdnst)


! Starting quantities on G-vectors


  call v_first(recip_%ns, recip_%ek, strfac_%sfact, ealpha,              &
      pseudo_%ealraw, pseudo_%nq, pseudo_%delq, pseudo_%vloc,            &
      pseudo_%dcor, pseudo_%dval,                                        &
      crys_%ntype, crys_%adot,                                           &
      vcomp_%vion, chdens_%denc, chdens_%dens, pseudo_%vql, pseudo_%dvql,&
      pseudo_%dnc, pseudo_%ddc,                                          &
      dims_%mxdtyp, dims_%mxdlqp, dims_%mxdnst)



! Ewald sums

  call ewald_sum(ewald_%energy, ewald_%force, ewald_%stress,             &
      crys_%adot, crys_%ntype, crys_%natom, crys_%rat, pseudo_%zv,       &
      dims_%mxdtyp, dims_%mxdatm)


  ipr = 0
  if(iprglob > 2) ipr = 1

  call print_energy(ipr, 'Ewald', ewald_%energy, ewald_%force,           &
      ewald_%stress,                                                     &
      crys_%adot, crys_%ntype, crys_%natom, crys_%nameat,                &
      dims_%mxdtyp, dims_%mxdatm)


! density extrapolation

  call move_extrapol_density(newcalc, flags_%flgcal,                     &
      nsave, chdsave,                                                    &
      chdens_%dens, chdens_%dend, chdens_%dend1, chdens_%den,            &
      recip_%ng, recip_%kgv, recip_%phase, recip_%conj,                  &
      recip_%ns, recip_%mstar,                                           &
      dims_%mxdgve, dims_%mxdnst)


  allocate(rholap(dims_%mxdnst))
  allocate(tau(dims_%mxdnst))

  do i = 1,recip_%ns
    tau(i) = C_ZERO
    rholap(i) = C_ZERO
  enddo

  ipr = 0
  if(iprglob > 1) ipr = 1

  call v_hartree_xc(ipr, xc_%author, xc_%tblaha, .FALSE., crys_%adot,    &
      exc, strxc, rhovxc,                                                &
      vcomp_%vhar, vcomp_%vxc, chdens_%den, chdens_%denc, rholap, tau,   &
      recip_%ng, recip_%kgv, recip_%phase, recip_%conj, recip_%ns,       &
      recip_%inds, recip_%kmax, recip_%mstar, recip_%ek,                 &
      dims_%mxdgve, dims_%mxdnst)

  deallocate(rholap)
  deallocate(tau)

  return

end subroutine cpw_scf_prepare
