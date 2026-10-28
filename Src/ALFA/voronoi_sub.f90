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

!>  Performs the analysis of the crystal geometry
!>  from the PW_RHO_V.DAT file
!>
!>  \author       Jose Luis Martins
!>  \version      5.13
!>  \date         22 April 2021. 28 October 2026.
!>  \copyright    GNU Public License v2

subroutine voronoi_sub(ioreplay)

! Added call to sym_space_group_name. 23 January 2022. JLM
! Added filenme to call write_cif. 9 June 2025. JLM
! mxdset in pw_rho_v_in_size. 7 October 2026. JLM+claude
! Analysis moved to voronoi_analysis. 28 October 2026. JLM+claude

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  ioreplay                        !<  tape number for reproducing calculations

! dimensions

  integer                            ::  mxdtyp                          !  array dimension of types of atoms
  integer                            ::  mxdatm                          !  array dimension of number of atoms of a given type

! crystal variables

  real(REAL64)                       ::  adot(3,3)                       !  metric in direct space in atomic units (Bohr radius)

  integer                            ::  ntype                           !  number of types of atoms
  integer, allocatable               ::  natom(:)                        !  number of atoms of type i
  real(REAL64), allocatable          ::  rat(:,:,:)                      !  lattice coordinates of atom j of type i
  character(len=2), allocatable      ::  nameat(:)                       !  chemical symbol for the type i

! other variables

  character(len=60)       ::  filename                                   !
  integer                 ::  io

! unused variables

  integer        ::  mxdgvein, mxdnstin, mxdlqp, mxdlao, mxdset

! counters

  integer        ::  i, j, k

  filename = 'PW_RHO_V.DAT'
  io = 21

  call pw_rho_v_in_size(filename, io,                                    &
     mxdtyp, mxdatm, mxdgvein, mxdnstin, mxdlqp, mxdlao, mxdset)

  allocate(natom(mxdtyp))
  allocate(rat(3,mxdatm,mxdtyp))
  allocate(nameat(mxdtyp))

! copyed/adapted from pw_rho_v_in_crystal_calc

  open(unit=io,file=trim(filename),status='old',form='UNFORMATTED')

  read(io) ntype
  if(ntype > mxdtyp) then
    write(6,*)
    write(6,'("  STOPPED in voronoi_sub")')
    write(6,'("  ntype = ",i8," is greater than mxdtyp =",i8)')          &
         ntype,mxdtyp

    stop

  endif

  read(io) (natom(i),i=1,ntype)
  do i=1,ntype
    if(natom(i) > mxdatm) then
      write(6,*)
      write(6,'("  STOPPED in voronoi_sub")')
      write(6,'("  for atom type ",i4," natom = ",i8," is ",             &
        & "greater than mxdatm =",i8)') i,natom(i),mxdatm

      stop

    endif
  enddo

  read(io)                                                               !  bdate
  read(io)                                                               !  author
  read(io)                                                               !  emax
  read(io)                                                               !  meta

  read(io) ((adot(i,j),i=1,3),j=1,3)

  read(io) (nameat(i),i=1,ntype)

  do i=1,ntype
    read(io) ((rat(j,k,i),j=1,3),k=1,natom(i))
  enddo

  close(unit=io)

  call voronoi_analysis(ioreplay,                                        &
      adot, ntype, natom, nameat, rat,                                   &
      mxdtyp, mxdatm)

  deallocate(natom)
  deallocate(rat)
  deallocate(nameat)

  return
end subroutine voronoi_sub
