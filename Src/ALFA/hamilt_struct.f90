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

!>  Calculates the indexes, kinetic energies, length of k+G
!>
!>  \author       Jose Luis Martins
!>  \version      4.94
!>  \date         8 February 1990,
!>  \copyright    GNU Public License v2

subroutine hamilt_struct(emax, rkpt, mtxd, isort, qmod, ekpg, lkplusg,   &
      ng, kgv, adot,                                                     &
      mxdgve, mxddim)

! written February 8 1990. jlm
! modified 24 March 1999. jlm
! modified June 3, 1999. jlm
! modified February 21, 2000. jlm
! modified 19 November 2013 (veffr1). jlm
! modified (f90) 14 January 2014.  jlm
! Split matrix_diag_kb into ham_struct and ham_diag_kb, 7 February 2014. JLM
! Modified lkplusg, October 2018. JLM
! Modified documentation, August 2019. JLM

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  mxdgve                          !<  array dimension for g-space vectors
  integer, intent(in)                ::  mxddim                          !<  array dimension of plane-waves

  real(REAL64), intent(in)           ::  emax                            !<  largest kinetic energy included in hamiltonian diagonal. (hartree).
  real(REAL64), intent(in)           ::  rkpt(3)                         !<  component in lattice coordinates of the k-point
  integer, intent(in)                ::  ng                              !<  total number of g-vectors with length less than gmax
  integer, intent(in)                ::  kgv(3,mxdgve)                   !<  i-th component (reciprocal lattice coordinates) of the n-th g-vector ordered by stars of increasing length
  real(REAL64), intent(in)           ::  adot(3,3)                       !<  metric in direct space

  logical, intent(in)                ::  lkplusg                         !<  If true use the previous G-vectors (same mtxd and isort)

! output

  real(REAL64), intent(out)          ::  qmod(mxddim)                    !<  length of k+g-vector of row/column i
  real(REAL64), intent(out)          ::  ekpg(mxddim)                    !<  kinetic energy (hartree) of k+g-vector of row/column i

! input and output

  integer, intent(inout)             ::  mtxd                            !<  dimension of the hamiltonian
  integer, intent(inout)             ::  isort(mxddim)                   !<  g-vector associated with row/column i of hamiltonian

! allocatable work arrays

  real(REAL64), allocatable   ::   ekin(:)
  integer, allocatable        ::   irow(:)

! local varaibles

  integer         ::  mtry
  real(real64)    ::  qk(3),vcell,bdot(3,3)
  real(real64)    ::  xsum

! parameters
  real(REAL64), parameter ::  EPS = 0.1E-20

! counters

  integer   ::  i

  allocate(ekin(ng))
  allocate(irow(ng))

  call adot_to_bdot(adot, vcell, bdot)

! calculate the kinetic energies

  do i = 1,ng
    qk(1) = rkpt(1) + kgv(1,i)
    qk(2) = rkpt(2) + kgv(2,i)
    qk(3) = rkpt(3) + kgv(3,i)
    ekin(i) = (qk(1)*bdot(1,1) + qk(2)*bdot(2,1) + qk(3)*bdot(3,1))*qk(1) +   &
              (qk(1)*bdot(1,2) + qk(2)*bdot(2,2) + qk(3)*bdot(3,2))*qk(2) +   &
              (qk(1)*bdot(1,3) + qk(2)*bdot(2,3) + qk(3)*bdot(3,3))*qk(3)
    ekin(i) = ekin(i) / 2
  enddo

  if(lkplusg) then

    do i = 1,mtxd
      irow(i) = isort(i)
    enddo

  else

! sort ekin

    xsum = (rkpt(1)*bdot(1,1) + rkpt(2)*bdot(2,1) + rkpt(3)*bdot(3,1))*rkpt(1) + &
           (rkpt(1)*bdot(1,2) + rkpt(2)*bdot(2,2) + rkpt(3)*bdot(3,2))*rkpt(2) + &
           (rkpt(1)*bdot(1,3) + rkpt(2)*bdot(2,3) + rkpt(3)*bdot(3,3))*rkpt(3)

!   sorts if it is not the gamma point

    if(xsum < eps) then
      do i = 1,ng
        irow(i) = i
      enddo
    else
      call sort(ng, ekin, irow)
    endif

!   find mtxd
!   mtxd is the size of the matrix

    mtxd = 0
    do i = 1,ng
      if(ekin(irow(i)) > emax) exit
      mtxd = i
    enddo

  endif

  if(mtxd > mxddim) then
    mtry = nint(mxddim * (emax / ekin(irow(mxddim)))**1.5)
    write(6,'("   STOPPED in hamilt_struct:    Try to increase",    &
       " mxddim to about",i8)') mtry

    stop

  endif
  if(mtxd == ng) then
    write(6,'("   STOPPED in hamilt_struct:    Something is",       &
        " wrong, the size of g-space",i8," is the same as the",     &
        " hamiltonian size")') ng

    stop

  endif

  do i = 1,mtxd
    ekpg(i) = ekin(irow(i))
    qmod(i) = sqrt(2*ekin(irow(i)))
    isort(i) = irow(i)
  enddo


  deallocate(ekin)
  deallocate(irow)


  return

end subroutine hamilt_struct
