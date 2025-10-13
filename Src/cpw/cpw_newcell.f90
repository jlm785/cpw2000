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

!>  Performs required steps when a new primitive cell
!>  appears in a minimization/molecular dynamics simulation
!>
!>  \author       Jose Luis Martins
!>  \version      5.12
!>  \date         Before February 2020, 12 October 2025
!>  \copyright    GNU Public License v2


subroutine cpw_newcell(lkpg, symkip, iprglob, symtol, kmscr, lnewnrk,    &
     dims_, recip_, crys_, spaceg_, pwexp_, strfac_, pseudo_,            &
     chdens_, vcomp_, flags_, kpoint_, atorb_,                           &
     hamallk_, psiallk_, filename_)

! Written from previous code before February 2020 (November 2019?). JLM
! Modified, option to read/write wave-functions from/to disk. 12 October 2025. JLM


  use cpw_variables

  implicit none

  type(dims_t)                       ::  dims_                           !<  array dimensions
  type(recip_t)                      ::  recip_                          !<  reciprocal space information
  type(crys_t)                       ::  crys_                           !<  crystal structure
  type(spaceg_t)                     ::  spaceg_                         !<  space group information
  type(pwexp_t)                      ::  pwexp_                          !<  plane-wave expansion choices
  type(strfac_t)                     ::  strfac_                         !<  structure factors
  type(chdens_t)                     ::  chdens_                         !<  charge densities
  type(vcomp_t)                      ::  vcomp_                          !<  Componemts of local potential
  type(pseudo_t)                     ::  pseudo_                         !<  pseudo-potential (Kleinman-Bylander)
  type(kpoint_t)                     ::  kpoint_                         !<  k-point data
  type(atorb_t)                      ::  atorb_                          !<  atomic orbitals in G-space
  type(hamallk_t)                    ::  hamallk_                        !<  hamiltonian size and indexation for all k-points
  type(psiallk_t)                    ::  psiallk_                        !<  psi for all k-points
  type(flags_t)                      ::  flags_                          !<  computational flags
  type(filename_t)                   ::  filename_                       !<  filenames

  logical, intent(inout)             ::  lkpg                            !<  If true use the previous G-vectors (same mtxd and isort)
  integer, intent(inout)             ::  kmscr(7)                        !<  max value of kgv(i,n) used for the potential FFT mesh (DUAL APPROXIMATION TYPE)
  logical, intent(inout)             ::  lnewnrk                         !<  new number of k-points


  integer, intent(in)                ::  iprglob                         !<  controls the amount of printing by subroutines
  logical, intent(in)                ::  symkip                          !<  calculate symmetry
  real(REAL64), intent(in)           ::  symtol                          !<  tolerance for symmetry recognition

! local variables

  integer              ::  ipr
  real(REAL64)         ::  tol
  integer              ::  isym

  integer              ::  istatus                         !  istatus = 0, successful; 1 not closed; 2 no inverse; 3 inconsistent with atomic positions

  real(REAL64)         ::  adotnew(3,3)


! do not adjust adot for symmetry

  adotnew(:,:) = crys_%adot(:,:)

  if(lkpg) then

    ipr = 0
    tol = symtol

    call sym_test(ipr, tol, istatus,                                &
       spaceg_%ntrans, spaceg_%mtrx, spaceg_%tnp,                   &
       crys_%ntype, crys_%natom, crys_%rat, adotnew,                &
       dims_%mxdtyp, dims_%mxdatm)

    if(istatus /= 0 .and. istatus /= 3) then
      lkpg = .FALSE.
      write(6,*)
      write(6,*)   '  Switched off the fixed k+G mode'
      write(6,*)
    endif

  endif

  if(lkpg) then

    call gsp_update_ek(recip_%ns, recip_%mstar, recip_%ek,          &
    recip_%kgv, crys_%adot,                                         &
    dims_%mxdgve, dims_%mxdnst)


  else

!   calculates the symmetry operations

    isym = 0
    if(symkip) isym = 1
    ipr = 0
    tol = symtol
    if(iprglob > 0) ipr = 1

    call sym_identify(isym, ipr, tol,                               &
    spaceg_%ntrans, spaceg_%mtrx, spaceg_%tnp,                      &
    adotnew, crys_%ntype, crys_%natom, crys_%rat,                   &
    dims_%mxdtyp, dims_%mxdatm)


!   reconstructs reciprocal space and allocates/deallocates recip_, strfac_, pseudo_, chdens_, vcomp_


    call cpw_gspace(iprglob, kmscr,                                 &
        dims_, crys_, spaceg_, pwexp_, recip_, strfac_, pseudo_,    &
        chdens_, vcomp_, flags_)


!   calculates integration k-points

    call cpw_bzint(iprglob, lnewnrk,                                &
              dims_, kpoint_, crys_, spaceg_,pwexp_)



!   hamiltonian matrix dimensions

!    maxdim = 0

    call cpw_size_alloc_hampsi(lnewnrk, dims_,                      &
      crys_, kpoint_, pwexp_, recip_, atorb_,                       &
      hamallk_, psiallk_, filename_)


  endif

  return

end subroutine cpw_newcell
