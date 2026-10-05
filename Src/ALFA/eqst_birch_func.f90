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

! Birch (third order Birch-Murnaghan) equation of state:
!   eqst_eofv   E(V)
!   eqst_pofv   p(V)
!   eqst_bofv   B(V)
!   eqst_vofp   V(p)
!   eqst_hofp   H(p)
! Merged from separate files into eqst_birch_func.f90. 3 October 2026. JLM+claude


!>  Birch equation of state E(V)
!>
!>  \author       Jose Luis Martins
!>  \version      5.13
!>  \date         September 1989, 3 October 2026.
!>  \copyright    GNU Public License v2

subroutine eqst_eofv(e, v, ezero, vzero, bzero, bprim)

! Prefix eqst_ for all subroutines and functions, split from eqst.f90. 3 October 2026. JLM+claude
! One argument per declaration, intent and description of all arguments. 3 October 2026. JLM+claude
! Integer literal in csi. 3 October 2026. JLM+claude

  implicit none
  integer, parameter  :: REAL64 = selected_real_kind(12)

! input

  real(REAL64), intent(in)    ::  v                                      !<  volume (bohr^3)
  real(REAL64), intent(in)    ::  ezero                                  !<  energy at equilibrium (Rydberg)
  real(REAL64), intent(in)    ::  vzero                                  !<  equilibrium volume (bohr^3)
  real(REAL64), intent(in)    ::  bzero                                  !<  bulk modulus at vzero (Rydberg/bohr^3)
  real(REAL64), intent(in)    ::  bprim                                  !<  pressure derivative of the bulk modulus

! output

  real(REAL64), intent(out)   ::  e                                      !<  energy at volume v (Rydberg)

! other variables

  real(REAL64)                ::  x, csi, x23

! constants

  real(REAL64), parameter    ::  UM = 1.0_REAL64

  x = vzero/v
  csi = 3 - 3*bprim/(4*UM)
  x23 = x**((2*UM)/(3*UM))
  e = ezero + ((9*UM)/(8*UM))*bzero*vzero*                             &
      (1+2*csi/(3*UM)+x23*(-2-2*csi+x23*(1+2*csi-2*csi*x23/(3*UM))))

  return

end subroutine eqst_eofv


!>  Birch equation of state p(V)
!>
!>  \author       Jose Luis Martins
!>  \version      5.13
!>  \date         September 1989, 3 October 2026.
!>  \copyright    GNU Public License v2

subroutine eqst_pofv(p, v, vzero, bzero, bprim)

! Prefix eqst_ for all subroutines and functions, split from eqst.f90. 3 October 2026. JLM+claude
! One argument per declaration, intent and description of all arguments. 3 October 2026. JLM+claude
! Integer literal in csi. 3 October 2026. JLM+claude

  implicit none
  integer, parameter  :: REAL64 = selected_real_kind(12)

! input

  real(REAL64), intent(in)    ::  v                                      !<  volume (bohr^3)
  real(REAL64), intent(in)    ::  vzero                                  !<  equilibrium volume (bohr^3)
  real(REAL64), intent(in)    ::  bzero                                  !<  bulk modulus at vzero (Rydberg/bohr^3)
  real(REAL64), intent(in)    ::  bprim                                  !<  pressure derivative of the bulk modulus

! output

  real(REAL64), intent(out)   ::  p                                      !<  pressure at volume v (Rydberg/bohr^3)

! other variables

  real(REAL64)                ::  x, csi, f

! constants

  real(REAL64), parameter    ::  UM = 1.0_REAL64

  x = vzero/v
  csi = 3 - 3*bprim/(4*UM)
  f = (x**((2*UM)/(3*UM))-UM) / 2
  p = 3*bzero*f*((UM+2*f)**((5*UM)/(2*UM)))*(UM-2*csi*f)

  return

end subroutine eqst_pofv


!>  Birch equation of state B(V)
!>
!>  \author       Jose Luis Martins
!>  \version      5.13
!>  \date         September 1989, 3 October 2026.
!>  \copyright    GNU Public License v2

subroutine eqst_bofv(b, v, vzero, bzero, bprim)

! Prefix eqst_ for all subroutines and functions, split from eqst.f90. 3 October 2026. JLM+claude
! One argument per declaration, intent and description of all arguments. 3 October 2026. JLM+claude
! Integer literal in csi. 3 October 2026. JLM+claude

  implicit none
  integer, parameter  :: REAL64 = selected_real_kind(12)

! input

  real(REAL64), intent(in)    ::  v                                      !<  volume (bohr^3)
  real(REAL64), intent(in)    ::  vzero                                  !<  equilibrium volume (bohr^3)
  real(REAL64), intent(in)    ::  bzero                                  !<  bulk modulus at vzero (Rydberg/bohr^3)
  real(REAL64), intent(in)    ::  bprim                                  !<  pressure derivative of the bulk modulus

! output

  real(REAL64), intent(out)   ::  b                                      !<  bulk modulus at volume v (Rydberg/bohr^3)

! other variables

  real(REAL64)                ::  x, csi, x23

! constants

  real(REAL64), parameter    ::  UM = 1.0_REAL64

  x = vzero/v
  csi = 3 - 3*bprim/(4*UM)
  x23 = x**((2*UM)/(3*UM))
  b = -bzero*(x**((5*UM)/(3*UM)))*                                      &
               (5-7*x23+csi*(5+x23*(-14+x23*9))) / (2*UM)

  return

end subroutine eqst_bofv


!>  Birch equation of state V(p)
!>
!>  The bulk modulus B(V) of the Birch equation of state is zero when
!>    9 csi y^2 - (7 + 14 csi) y + 5 (1 + csi) = 0,  y = (vzero/V)^(2/3)
!>  with csi = 3 - 3 bprim / 4.  At y = 1 the left side is -2 (B = bzero > 0).
!>  The mechanically stable branch is the interval ylo < y < yhi around 1
!>  limited by the nearest roots.  For bprim < 4 the pressure has a maximum
!>  at yhi, and for any bprim it has a minimum (negative) at ylo.
!>  In that interval p(V) is monotonic, so a bracketed root is unique.
!>
!>  \author       Jose Luis Martins+claude
!>  \version      5.13
!>  \date         September 1989, 3 October 2026.
!>  \copyright    GNU Public License v2

subroutine eqst_vofp(v, p, vzero, bzero, bprim, ierr)

! Prefix eqst_ for all subroutines and functions, split from eqst.f90. 3 October 2026. JLM+claude
! Rewritten for robustness: stays on the stable branch (B > 0),
! bracketing and num_zbrent instead of Newton, ierr instead of stop. 3 October 2026. JLM+claude

  implicit none
  integer, parameter  :: REAL64 = selected_real_kind(12)

! input

  real(REAL64), intent(in)    ::  p                                      !<  pressure
  real(REAL64), intent(in)    ::  vzero                                  !<  equilibrium volume
  real(REAL64), intent(in)    ::  bzero                                  !<  bulk modulus at vzero
  real(REAL64), intent(in)    ::  bprim                                  !<  pressure derivative of the bulk modulus

! output

  real(REAL64), intent(out)   ::  v                                      !<  volume at pressure p
  integer, intent(out)        ::  ierr                                   !<  0: OK; 1: p outside the stable range, v is the limiting volume; 2: vzero or bzero not positive, v = vzero

! local variables

  real(REAL64)                ::  csi
  real(REAL64)                ::  a2, a1, a0                             !  coefficients of the quadratic in y
  real(REAL64)                ::  disc, q
  real(REAL64)                ::  root(2)
  integer                     ::  nroot
  real(REAL64)                ::  ylo, yhi                               !  limits of the stable branch
  real(REAL64)                ::  vlo, vhi                               !  bracket for the volume
  real(REAL64)                ::  plim                                   !  pressure at the end of the stable branch
  real(REAL64)                ::  pv

  real(REAL64)                ::  xval, fval                             !  for num_zbrent
  integer                     ::  iflag

! constants

  real(REAL64), parameter     ::  ZERO = 0.0_REAL64, UM = 1.0_REAL64
  real(REAL64), parameter     ::  TOL = 1.0E-12_REAL64                   !  relative tolerance in volume
  integer, parameter          ::  MXHALF = 200                           !  maximum number of halvings of the volume

! counters

  integer                     ::  i


  ierr = 0

  if(vzero <= ZERO .or. bzero <= ZERO) then
    v = vzero
    ierr = 2

    return

  endif

  if(p == ZERO) then
    v = vzero

    return

  endif

! roots of the quadratic, numerically stable form (the discriminant is always positive)

  csi = 3 - 3*bprim/(4*UM)
  a2 = 9*csi
  a1 = -(7 + 14*csi)
  a0 = 5*(1 + csi)

  disc = a1*a1 - 4*a2*a0
  q = -(a1 + sign(sqrt(disc),a1)) / 2

  nroot = 0
  if(q /= ZERO) then
    nroot = nroot + 1
    root(nroot) = a0 / q
  endif
  if(a2 /= ZERO) then
    nroot = nroot + 1
    root(nroot) = q / a2
  endif

! nearest roots on each side of y = 1.  ylo = 0 or yhi = 0 mean no limit.

  ylo = ZERO
  yhi = ZERO
  do i = 1,nroot
    if(root(i) > ZERO .and. root(i) < UM) then
      if(root(i) > ylo) ylo = root(i)
    elseif(root(i) > UM) then
      if(yhi == ZERO .or. root(i) < yhi) yhi = root(i)
    endif
  enddo

  if(p > ZERO) then

!   compression, vlo < v < vzero

    vhi = vzero

    if(yhi > ZERO) then
      vlo = vzero / yhi**((3*UM)/(2*UM))
      call eqst_pofv(plim, vlo, vzero, bzero, bprim)
      if(p > plim) then
        v = vlo
        ierr = 1

        return

      endif
    else
      vlo = vzero
      do i = 1,MXHALF
        vlo = vlo / 2
        call eqst_pofv(pv, vlo, vzero, bzero, bprim)
        if(pv >= p) exit
      enddo
      if(pv < p) then
        v = vlo
        ierr = 1

        return

      endif
    endif

  else

!   expansion, vzero < v < vhi

    vlo = vzero

    if(ylo > ZERO) then
      vhi = vzero / ylo**((3*UM)/(2*UM))
      call eqst_pofv(plim, vhi, vzero, bzero, bprim)
      if(p < plim) then
        v = vhi
        ierr = 1

        return

      endif
    else

!     cannot happen for the Birch equation (p -> 0 for large v), for safety

      v = vzero
      ierr = 1

      return

    endif

  endif

! solves p(v) = p in the bracket

  iflag = 0
  fval = ZERO
  do
    call num_zbrent(xval, fval, vlo, vhi, TOL*vzero, iflag)
    if(iflag == 0) exit
    call eqst_pofv(pv, xval, vzero, bzero, bprim)
    fval = pv - p
  enddo

  v = xval

  return

end subroutine eqst_vofp


!>  Birch equation of state H(p)
!>
!>  \author       Jose Luis Martins
!>  \version      5.13
!>  \date         September 1989, 3 October 2026.
!>  \copyright    GNU Public License v2

subroutine eqst_hofp(h, p, ezero, vzero, bzero, bprim, ierr)

! Prefix eqst_ for all subroutines and functions, split from eqst.f90. 3 October 2026. JLM+claude
! ierr from eqst_vofp. 3 October 2026. JLM+claude
! One argument per declaration, intent and description of all arguments. 3 October 2026. JLM+claude

  implicit none

  integer, parameter  :: REAL64 = selected_real_kind(12)

! input

  real(REAL64), intent(in)    ::  p                                      !<  pressure (Rydberg/bohr^3)
  real(REAL64), intent(in)    ::  ezero                                  !<  energy at equilibrium (Rydberg)
  real(REAL64), intent(in)    ::  vzero                                  !<  equilibrium volume (bohr^3)
  real(REAL64), intent(in)    ::  bzero                                  !<  bulk modulus at vzero (Rydberg/bohr^3)
  real(REAL64), intent(in)    ::  bprim                                  !<  pressure derivative of the bulk modulus

! output

  real(REAL64), intent(out)   ::  h                                      !<  enthalpy at pressure p (Rydberg)
  integer, intent(out)        ::  ierr                                   !<  0: OK; 1: p outside the stable range of the equation of state; 2: invalid parameters

! other variables

  real(REAL64)                ::  v, e

  call eqst_vofp(v, p, vzero, bzero, bprim, ierr)
  call eqst_eofv(e, v, ezero, vzero, bzero, bprim)
  h = e + p*v

  return

end subroutine eqst_hofp
