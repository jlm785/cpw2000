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

!>  Driver for the atomic orbital interpolation, MTB (Modified Tight
!>  Binding) method.  If the file hsr2.dat exists it may be used,
!>  otherwise the MTB hamiltonian is calculated from PW_RHO_V.DAT and
!>  written to hsr2.dat (ao_interpolation_create).  The results are
!>  then calculated from hsr2.dat (ao_interpolation_use).
!>
!>  \author       Carlos Loia Reis, Jose Luis Martins
!>  \version      5.13
!>  \date         December 18, 2013, 23 October 2026.
!>  \copyright    GNU Public License v2


subroutine ao_interpolation_sub(ioreplay)

! written December 18, 2013. jlm
! modified March 31, 2014.
! modified May 22, 2014.
! modified, mxdlao, December 1, 2015.  JLM
! Modified, ioreplay, 8 January 2017. JLM
! Modified, cpw_variables. 7 January 2020. JLM
! Modified, epspsi, icmax, 14 June 2020. JLM
! Modified, efermi, 29 November 2021. JLM
! Modified, size of author, 13 January 2024.
! Modified, indentation, duplicate ao_interpolation_prepare. 2 October 2024. JLM
! Modified, cpw_pp_band_dos_init/prepare. 12 March 2025. JLM
! Modified, inconsistent cpw_pp_band_prepare, 27 July 2026. Lukas Bauer
! Generalized Kohn-Sham meta-GGA, lgks from author. 7 October 2026. JLM+claude
! Logical assignments written as if statements. 12 October 2026. JLM+claude
! pseudo_%nq, pseudo_%delq, pseudo_%vkb renamed nqnl, delqnl, vkbraw. 14 October 2026. JLM+claude
! atorb_%n_bsets, atorb_%jbset, dims_%mxdset (several atomic basis sets). 16 October 2026. JLM+claude
! Choice of atomic basis set (cpw_pp_atorb_questions), consistency with hsr.dat. 20 October 2026. JLM+claude
! tnp passed to ao_interpolation_prepare (projection method). 20 October 2026. JLM+claude
! icmplx no longer passed to ao_interpolation_out_band_fold_full and ao_interpolation_w90. 23 October 2026. JLM+claude
! Split in ao_interpolation_create and ao_interpolation_use, file hsr2.dat. 23 October 2026. JLM+claude

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  ioreplay                        !<  tape number for reproducing calculations

! header of hsr2.dat

  integer                            ::  mxdtyp                          !  array dimension of types of atoms
  integer                            ::  mxdatm                          !  array dimension of number of atoms of a given type
  integer                            ::  mxdlao                          !  array dimension of orbital per atom type
  integer                            ::  mxdlqp                          !  array dimension of the radial functions

  character(len=250)                 ::  meta_cpw2000                    !  metadata from cpw2000
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
  integer                            ::  mtrx(3,3,48)                    !  rotation matrix of the symmetry operations
  real(REAL64)                       ::  tnp(3,48)                       !  2*pi* fractional translation of the symmetry operations
  integer, allocatable               ::  norbat(:,:)                     !  number of atomic orbitals of atom type nt
  integer, allocatable               ::  lorb(:,:,:)                     !  angular momentum of orbital n of atom type nt
  integer, allocatable               ::  nqwf(:)                         !  number of points for wavefunction interpolation
  real(REAL64), allocatable          ::  delqwf(:)                       !  step used in the wavefunction interpolation
  real(REAL64), allocatable          ::  wvfao(:,:,:,:)                  !  (1/q**l) * wavefunction for atom type nt, orbital n
  integer                            ::  lso, loptical                   !  spin-orbit, optical matrix elements (1 if present)
  integer                            ::  imtb                            !  MTB method, 1 projection, 2 matching
  integer                            ::  mp_grid(3)                      !  grid of k-points of the MTB calculation
  integer                            ::  nband                           !  number of basis functions

! other variables

  logical            ::  lexist, lok, lcreate
  character(len=1)   ::  yesno

! parameters

  integer, parameter :: IOHSR = 777                                      !  tape number of hsr2.dat


  write(6,*)
  write(6,'("  Analysis of the electronic structure ")')
  write(6,'("  with the atomic orbital interpolation (MTB)")')
  write(6,*)

  lcreate = .TRUE.

  inquire(file='hsr2.dat', exist=lexist)

  if(lexist) then

    open(unit=IOHSR, file='hsr2.dat', form='unformatted', status='old')

    call ao_int_hsr2_size(IOHSR, lok, ntype, mxdatm, mxdlao, mxdlqp)

    if(lok) then

      mxdtyp = ntype

      allocate(natom(mxdtyp), nameat(mxdtyp), rat(3,mxdatm,mxdtyp))
      allocate(norbat(1,mxdtyp), lorb(mxdlao,1,mxdtyp))
      allocate(nqwf(mxdtyp), delqwf(mxdtyp))
      allocate(wvfao(-2:mxdlqp,mxdlao,1,mxdtyp))

      call ao_int_hsr2_read(IOHSR, meta_cpw2000, title, subtitle,        &
          adot, ztot, emax, ntype, natom, nameat, rat,                   &
          ntrans, mtrx, tnp,                                             &
          norbat, lorb, nqwf, delqwf, wvfao,                             &
          lso, loptical, imtb, mp_grid, nband,                           &
          mxdtyp, mxdatm, mxdlqp, mxdlao)

      call ao_int_hsr2_print(title, subtitle, ztot, emax,                &
          ntype, natom, nameat, norbat, lorb,                            &
          lso, loptical, imtb, mp_grid, nband,                           &
          mxdtyp, mxdlao)

      deallocate(natom, nameat, rat)
      deallocate(norbat, lorb)
      deallocate(nqwf, delqwf)
      deallocate(wvfao)

      write(6,*) '  Do you want to use it (y/n)?'
      write(6,*)
      read(5,*) yesno
      write(ioreplay,*) yesno,'   use hsr2.dat'

      if(yesno == 'y' .or. yesno == 'Y') then
        lcreate = .FALSE.
        write(6,*)
        write(6,*) '  using the existing hsr2.dat'
        write(6,*)
      else
        write(6,*)
        write(6,*) '  recalculates the MTB data'
        write(6,*)
      endif

    else

      write(6,*)
      write(6,*) '  hsr2.dat was written by another version, recalculates'
      write(6,*)

    endif

    close(unit=IOHSR)

  endif

  if(lcreate) call ao_interpolation_create(ioreplay)

  call ao_interpolation_use(ioreplay)

  return

end subroutine ao_interpolation_sub
