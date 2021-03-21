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

!>  checks if a lattice vector translation is needed
!>  for a k-point that is too far away from the 1st Brillouin zone

subroutine h_kb_dia_check_k(emax, rkpt,                                  &
    lkshift, qualpsi, mtxd, isort,                                       &
    ng, kgv, adot, kmscr,                                                &
    mxdgve, mxddim)

! written 30 January 2021, based previous code.  JLM
! copyright  Jose Luis Martins/Carlos Loia Reis/INESC-MN

! version 4.99

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

  integer, intent(in)                ::  kmscr(7)                        !<  max value of kgv(i,n) used for the potential fft mesh

! output

  logical, intent(out)               ::  lkshift                         !<  a shift to the 1st Brillouin zone is necessary
  real(REAL64), intent(out)          ::  qualpsi                         !<  estimate of the quality of psi after applying the shift.  1 good, 0 terrible.
  integer, intent(out)               ::  mtxd                            !<  dimension of the hamiltonian
  integer, intent(out)               ::  isort(mxddim)                   !<  g-vector associated with row/column i of hamiltonian

!      allocatable work arrays

  real(REAL64), allocatable   ::   ekin(:)
  integer, allocatable        ::   irow(:)

! other variables

  real(real64)    ::  qk(3),vcell,bdot(3,3)

  integer         ::  n1, n2, n3
  integer         ::  kd1, kd2, kd3

  real(real64)    ::  qual
  
! parameters
  
  real(REAL64), parameter    ::  UM = 1.0_REAL64, ZERO = 0.0_REAL64

! counters

  integer        ::  i


  allocate(ekin(ng))
  allocate(irow(ng))

  call adot_to_bdot(adot,vcell,bdot)

! calculate the kinetic energies

  kd1 = 0
  kd2 = 0
  kd3 = 0

  do i=1,ng
    qk(1) = rkpt(1) + kgv(1,i)
    qk(2) = rkpt(2) + kgv(2,i)
    qk(3) = rkpt(3) + kgv(3,i)
    ekin(i) = (qk(1)*bdot(1,1) + qk(2)*bdot(2,1) + qk(3)*bdot(3,1))*qk(1)   &
            + (qk(1)*bdot(1,2) + qk(2)*bdot(2,2) + qk(3)*bdot(3,2))*qk(2)   &
            + (qk(1)*bdot(1,3) + qk(2)*bdot(2,3) + qk(3)*bdot(3,3))*qk(3)
    ekin(i) = ekin(i) / 2
    if(emax > ekin(i)) then
      if (iabs(kgv(1,i)) > kd1) kd1 = iabs(kgv(1,i))
      if (iabs(kgv(2,i)) > kd2) kd2 = iabs(kgv(2,i))
      if (iabs(kgv(3,i)) > kd3) kd3 = iabs(kgv(3,i))
    endif
  enddo

! find mtxd, the size of the matrix

  call sort(ng,ekin,irow)

  mtxd = 0
  do i=1,ng
    if(ekin(irow(i)) > emax) exit 
    mtxd = i
  enddo

  if(mtxd > mxddim) mtxd = mxddim

  do i=1,mtxd
    isort(i) = irow(i)
  enddo
  

  n1 = kmscr(4)
  n2 = kmscr(5)
  n3 = kmscr(6)

  if(kd1 > (n1-1)/2 .or. kd2 > (n2-1)/2 .or. kd3 > (n3-1)/2) then
    lkshift = .TRUE.
    qualpsi = ZERO
    if(kd1 > (n1-1)/2) then
      qual = UM*(kd1 - (n1-1)/2) / (UM*((n1-1)/2))
      if(qual > UM) qual = UM
      qualpsi = UM - qual*qual
    endif
    if(kd2 > (n2-1)/2) then
      qual = UM*(kd2 - (n2-1)/2) / (UM*((n2-1)/2))
      if(qual > UM) qual = UM
      qualpsi = max(UM - qual*qual, qualpsi)
    endif
    if(kd3 > (n3-1)/2) then
      qual = UM*(kd3 - (n3-1)/2) / (UM*((n3-1)/2))
      if(qual > UM) qual = UM
      qualpsi = max(UM - qual*qual, qualpsi)
    endif
  else
    lkshift = .FALSE.
    qualpsi = UM
  endif

end subroutine h_kb_dia_check_k
