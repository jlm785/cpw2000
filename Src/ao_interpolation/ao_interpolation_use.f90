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

!>  Reads the modified tight-binding (MTB) data in hsr2.dat and
!>  calculates band structures, densities of states, files for the
!>  dielectric function, Wannier90 and WannierTools.  It does not
!>  need PW_RHO_V.DAT.
!>
!>  \author       claude
!>  \version      5.13
!>  \date         23 October 2026.
!>  \copyright    GNU Public License v2

subroutine ao_interpolation_use(ioreplay)

! Written 23 October 2026 from ao_interpolation_sub (Carlos Loia Reis,
! Jose Luis Martins and others). claude

  use NonOrthoInterp

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  ioreplay                        !<  tape number for reproducing calculations

! dimensions

  integer                            ::  mxdtyp                          !  array dimension of types of atoms
  integer                            ::  mxdatm                          !  array dimension of number of atoms of a given type
  integer                            ::  mxdlao                          !  array dimension of orbital per atom type
  integer                            ::  mxdlqp                          !  array dimension of the radial functions
  integer                            ::  mxdgve                          !  array dimension for g-space vectors
  integer                            ::  mxdnst                          !  array dimension for g-space stars
  integer                            ::  mxdcub                          !  array dimension for 3-index g-space

! MTB data

  type(noiData_t)                    ::  noiData                         !  see NonOrthoInterp

! header of hsr2.dat

  character(len=250)                 ::  meta_cpw2000                    !  metadata from cpw2000 (superlattice information for unfolding)
  character(len=50)                  ::  title                           !  title for plots
  character(len=140)                 ::  subtitle                        !  subtitle for plots
  real(REAL64)                       ::  adot(3,3)                       !  metric in direct space
  real(REAL64)                       ::  ztot                            !  total charge density (electrons/cell)
  real(REAL64)                       ::  emax                            !  kinetic energy cutoff of plane wave expansion (hartree)
  integer                            ::  ntype                           !  number of types of atoms
  integer, allocatable               ::  natom(:)                        !  number of atoms of type i
  character(len=2), allocatable      ::  nameat(:)                       !  chemical symbol for the type i
  real(REAL64), allocatable          ::  rat(:,:,:)                      !  lattice coordinates of atom j of type i
  integer                            ::  ntrans                          !  number of symmetry operations in the factor group
  integer                            ::  mtrx(3,3,48)                    !  rotation matrix (in reciprocal lattice coordinates) of the symmetry operations
  real(REAL64)                       ::  tnp(3,48)                       !  2*pi* fractional translation of the symmetry operations
  integer, allocatable               ::  n_bsets(:)                      !  number of atomic basis sets (1)
  integer, allocatable               ::  jbset(:)                        !  basis set chosen (1)
  integer, allocatable               ::  norbat(:,:)                     !  number of atomic orbitals of atom type nt
  integer, allocatable               ::  lorb(:,:,:)                     !  angular momentum of orbital n of atom type nt
  integer, allocatable               ::  nqwf(:)                         !  number of points for wavefunction interpolation
  real(REAL64), allocatable          ::  delqwf(:)                       !  step used in the wavefunction interpolation
  real(REAL64), allocatable          ::  wvfao(:,:,:,:)                  !  (1/q**l) * wavefunction for atom type nt, orbital n
  integer                            ::  lso, loptical                   !  spin-orbit, optical matrix elements (1 if present)
  integer                            ::  imtb                            !  MTB method, 1 projection, 2 matching
  integer                            ::  mp_grid(3)                      !  grid of k-points of the MTB calculation
  integer                            ::  nband                           !  number of basis functions

! G-space (only for unfolding and Wannier90)

  logical                            ::  lgspace                         !  G-space already calculated
  integer                            ::  ng                              !  total number of g-vectors with length less than gmax
  integer, allocatable               ::  kgv(:,:)                        !  i-th component (reciprocal lattice coordinates) of the n-th g-vector ordered by stars of increasing length
  complex(REAL64), allocatable       ::  phase(:)                        !  phase factor of G-vector n
  real(REAL64), allocatable          ::  conj(:)                         !  is -1 if one must take the complex conjugate of x*phase
  integer                            ::  ns                              !  number os stars with length less than gmax
  integer, allocatable               ::  inds(:)                         !  star to which g-vector n belongs
  integer                            ::  kmax(3)                         !  max value of kgv(i,n)
  integer, allocatable               ::  indv(:)                         !  kgv(i,indv(jadd)) is the g-vector associated with jadd
  integer, allocatable               ::  mstar(:)                        !  number of g-vectors in the j-th star
  real(REAL64), allocatable          ::  ek(:)                           !  kinetic energy (Hartree) of g-vectors in star j
  integer, allocatable               ::  izstar(:)                       !  is 0 if the phase=0

! other variables

  logical            ::  lok
  integer            ::  ipr
  character(len=1)   ::  yesno

! parameters

  integer, parameter :: IOHSR = 777                                      !  tape number of hsr2.dat


! reads hsr2.dat

  open(unit=IOHSR, file='hsr2.dat', form='unformatted', status='old')

  call ao_int_hsr2_size(IOHSR, lok, ntype, mxdatm, mxdlao, mxdlqp)

  if(.not. lok) then
    write(6,*)
    write(6,*) '   STOPPED in ao_interpolation_use:  hsr2.dat has another version'
    write(6,*)

    stop

  endif

  mxdtyp = ntype

  allocate(natom(mxdtyp), nameat(mxdtyp), rat(3,mxdatm,mxdtyp))
  allocate(n_bsets(mxdtyp), jbset(mxdtyp))
  allocate(norbat(1,mxdtyp), lorb(mxdlao,1,mxdtyp))
  allocate(nqwf(mxdtyp), delqwf(mxdtyp))
  allocate(wvfao(-2:mxdlqp,mxdlao,1,mxdtyp))

  call ao_int_hsr2_read(IOHSR, meta_cpw2000, title, subtitle,            &
      adot, ztot, emax, ntype, natom, nameat, rat,                       &
      ntrans, mtrx, tnp,                                                 &
      norbat, lorb, nqwf, delqwf, wvfao,                                 &
      lso, loptical, imtb, mp_grid, nband,                               &
      mxdtyp, mxdatm, mxdlqp, mxdlao)

  call NonOrthoInterpReadFromFile(noiData, IOHSR)

  close(unit=IOHSR)

! the basis set in hsr2.dat is the only one

  n_bsets(:) = 1
  jbset(:) = 1

  lgspace = .FALSE.

!-------------------------------------------------------------------

  write(6,*)
  write(6,*) '  do you want a band structure with unfolding'
  write(6,*) '  using ao interpolation? (y/n)'
  write(6,*)
  read(5,*) yesno
  write(ioreplay,*) yesno,'   band unfold'

  if(yesno == 'y' .or. yesno == 'Y') then

    write(6,*)
    write(6,*) '     band structure with unfolding  '
    write(6,*) '         (ao interpolation)         '

    call ao_int_use_gspace()

    call ao_interpolation_out_band_fold_full(ioreplay,                   &
        meta_cpw2000, title, subtitle, noiData,                          &
        emax, ztot,                                                      &
        adot, ntype, natom, nameat,                                      &
        rat,                                                             &
        ng, kgv,                                                         &
        n_bsets, jbset,                                                  &
        norbat, nqwf, delqwf,                                            &
        wvfao, lorb,                                                     &
        mxdtyp, mxdatm, mxdgve, mxdlqp,                                  &
        mxdlao, 1)

    write(6,*) '     '
    write(6,*) '     completed'
    write(6,*)

  endif


  write(6,*)
  write(6,*) '  do you want a band structure '
  write(6,*) '  using ao interpolation? (y/n)'
  write(6,*)
  read(5,*) yesno
  write(ioreplay,*) yesno,'   band'

  if(yesno == 'y' .or. yesno == 'Y') then

    write(6,*) '     band structure (ao interpolation)'

    call ao_interpolation_out_band(title, subtitle,                      &
        noiData, ztot, adot)
    write(6,*) '     '
    write(6,*) '     completed'
    write(6,*)

  endif


  write(6,*)
  write(6,*) '  do you want prepare file for density of states '
  write(6,*) '  using ao interpolation? (y/n)'
  write(6,*)
  read(5,*) yesno
  write(ioreplay,*) yesno,'   dos prepare'

  if(yesno == 'y' .or. yesno == 'Y') then

    write(6,*) '     files for density of states (ao interpolation)'

    call ao_interpolation_out_dos(noiData, ztot, adot, ntrans, mtrx)

    write(6,*) '     '
    write(6,*) '     completed'
    write(6,*)

  endif


  if(noiData % loptical == 1) then

    write(6,*)
    write(6,*) '  do you want prepare file for the dielectric function '
    write(6,*) '  using ao interpolation? (y/n)'
    write(6,*)
    read(5,*) yesno
    write(ioreplay,*) yesno,'   ie_prepare'

    if(yesno == 'y' .or. yesno == 'Y') then
      write(6,*) '     files for dielectric function calculation (ao interpolation)'

      call ao_interpolation_out_ie(noiData, ztot, adot, ntrans, mtrx)

      write(6,*) '     '
      write(6,*) '     completed'
      write(6,*)

    endif

  endif


  write(6,*)
  write(6,*) '  do you want the files for Wannier90 run? '
  write(6,*)
  read(5,*) yesno
  write(ioreplay,*) yesno,'   w90_interface'

  if(yesno == 'y' .or. yesno == 'Y') then

    call ao_int_use_gspace()

    call ao_interpolation_w90(noiData,                                   &
        emax, ztot,                                                      &
        adot,                                                            &
        ntype, natom, nameat, rat,                                       &
        ng, kgv,                                                         &
        n_bsets, jbset,                                                  &
        norbat, nqwf, delqwf,                                            &
        wvfao, lorb,                                                     &
        mxdtyp, mxdatm, mxdgve, mxdlqp,                                  &
        mxdlao, 1)

  endif


  write(6,*)
  write(6,*) '  do you want prepare file compatible with WannierTools? '
  write(6,*)
  read(5,*) yesno
  write(ioreplay,*) yesno,'   write hr'

  if(yesno == 'y' .or. yesno == 'Y') then

     call ao_interpolation_write_hr(noiData, adot, ioreplay)

     call ao_interpolation_write_wt_in(noiData, ztot,                    &
          adot, ntype, natom,                                            &
          nameat, rat, n_bsets, jbset,                                   &
          norbat, lorb,                                                  &
          mxdtyp, mxdatm, mxdlao, 1)

     write(6,*) '     '
     write(6,*) '     completed'

  endif


  call NonOrthoInterpFinish(noiData)

  deallocate(natom, nameat, rat)
  deallocate(n_bsets, jbset)
  deallocate(norbat, lorb)
  deallocate(nqwf, delqwf)
  deallocate(wvfao)

  if(lgspace) then
    deallocate(kgv, phase, conj, inds, indv, mstar, ek, izstar)
  endif

  write(6,*)
  write(6,*) '  All Done with ao interpolation!'
  write(6,*)

  return

contains

!>  G-space of the plane-wave cutoff, calculated once

  subroutine ao_int_use_gspace()

    if(lgspace) return

    call size_g_space(emax, adot, ntrans, mtrx,                          &
        mxdgve, mxdnst, mxdcub)

    allocate(kgv(3,mxdgve))
    allocate(phase(mxdgve))
    allocate(conj(mxdgve))
    allocate(inds(mxdgve))
    allocate(indv(mxdcub))
    allocate(mstar(mxdnst))
    allocate(ek(mxdnst))
    allocate(izstar(mxdnst))

    ipr = 1
    call g_space(ipr, emax,                                              &
        adot, ntrans, mtrx, tnp,                                         &
        ng, kgv, phase, conj,                                            &
        inds, kmax, indv, ns, mstar, ek, izstar,                         &
        mxdgve, mxdnst, mxdcub)

    lgspace = .TRUE.

    return

  end subroutine ao_int_use_gspace

end subroutine ao_interpolation_use
