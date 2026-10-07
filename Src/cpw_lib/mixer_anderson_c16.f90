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

!>  Computes a new input vector (prototype G-vector of each star)
!>  in an iterative scheme using Anderson's extrapolation scheme,
!>  eqs 4.1-4.9, 4.15-4.18 of
!>  D.G.Anderson J.Assoc.Computing Machinery, 12, 547 (1965).
!>  The scalar product is the real space one, sum_j mstar(j) Re(x_j conj(y_j)).
!>  The history is kept by the calling program in vmem.
!>
!>  icy = 1, vin is replaced by vin + vout (as in atom_atm_dmixp)
!>  icy = 2, linear mixing, starts (or restarts) the history.
!>  icy = 3, Anderson with one previous iteration.
!>  icy > 3, Anderson with two previous iterations if id = 3.
!>
!>  \author       Jose Luis Martins
!>  \version      5.13
!>  \date         1980s, 7 October 2026.
!>  \copyright    GNU Public License v2

subroutine mixer_anderson_c16(icy, id, beta, vout, vin, vmem,            &
    ns, mstar, mxdnst)

! Adapted from atom_atm_dmixp of the pseudopotential code (adapted from K.C.Pandey).
! Complex version for stars with real space metric. 7 October 2026. JLM+claude

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  mxdnst                          !<  array dimension for g-space stars

  integer, intent(in)                ::  icy                             !<  cycle number, icy = 2 (re)starts the history with linear mixing
  integer, intent(in)                ::  id                              !<  method: 1 linear mixing, 2 Anderson with one, 3 with two previous iterations
  real(REAL64), intent(in)           ::  beta                            !<  mixing coefficient

  integer, intent(in)                ::  ns                              !<  number os stars with length less than gmax
  integer, intent(in)                ::  mstar(mxdnst)                   !<  number of g-vectors in the j-th star

  complex(REAL64), intent(in)        ::  vout(mxdnst)                    !<  output vector for the prototype G-vector

! input and output

  complex(REAL64), intent(inout)     ::  vin(mxdnst)                     !<  input: old input vector, output: new input vector
  complex(REAL64), intent(inout)     ::  vmem(mxdnst,4)                  !<  history (residual-1, residual-2, input-1, input-2). Do not change outside!

! local allocatable arrays

  complex(REAL64), allocatable       ::  res(:)                          !  residual vout - vin
  complex(REAL64), allocatable       ::  c(:)                            !  difference of residuals with iteration-1
  complex(REAL64), allocatable       ::  d(:)                            !  difference of residuals with iteration-2

! local variables

  integer                ::  in
  real(REAL64)           ::  d11, d12, d22
  real(REAL64)           ::  rd1m, rd2m
  real(REAL64)           ::  t1, t2
  real(REAL64)           ::  x
  real(REAL64)           ::  a2, det, dett

! parameters

  real(REAL64), parameter            ::  ZERO = 0.0_REAL64
  real(REAL64), parameter            ::  UM = 1.0_REAL64
  real(REAL64), parameter            ::  DETOL = 1.0E-9_REAL64

! counters

  integer       ::  i


  in = icy - 1

  if(in == 0) then

    do i = 1,ns
      vin(i) = vin(i) + vout(i)
    enddo

    return

  endif

  allocate(res(ns))

  do i = 1,ns
    res(i) = vout(i) - vin(i)
  enddo

  if(id == 1) then

    do i = 1,ns
      vin(i) = vin(i) + beta*res(i)
    enddo

  elseif(in == 1) then

!   linear mixing, starts history

    do i = 1,ns
      vmem(i,1) = res(i)
      vmem(i,3) = vin(i)
      vin(i) = vin(i) + beta*res(i)
    enddo

  else

    allocate(c(ns), d(ns))

    do i = 1,ns
      c(i) = vmem(i,1) - res(i)
    enddo
    if(id == 3 .and. in > 2) then
      do i = 1,ns
        d(i) = vmem(i,2) - res(i)
      enddo
    endif
    if(id > 2) then
      do i = 1,ns
        vmem(i,2) = vmem(i,1)
      enddo
    endif
    do i = 1,ns
      vmem(i,1) = res(i)
    enddo

    d11 = ZERO
    rd1m = ZERO
    do i = 1,ns
      d11 = d11 + mstar(i)*real(c(i)*conjg(c(i)),REAL64)
      rd1m = rd1m + mstar(i)*real(res(i)*conjg(c(i)),REAL64)
    enddo

    if(d11 <= ZERO) then

!     identical residuals, linear mixing

      t1 = ZERO
      t2 = ZERO

    elseif(in <= 2 .or. id <= 2) then

      t1 = -rd1m / d11
      t2 = ZERO

    else

      d22 = ZERO
      d12 = ZERO
      rd2m = ZERO
      do i = 1,ns
        d22 = d22 + mstar(i)*real(d(i)*conjg(d(i)),REAL64)
        d12 = d12 + mstar(i)*real(c(i)*conjg(d(i)),REAL64)
        rd2m = rd2m + mstar(i)*real(res(i)*conjg(d(i)),REAL64)
      enddo
      a2 = d11*d22
      det = a2 - d12*d12
      dett = ZERO
      if(a2 > ZERO) dett = det / a2

      if(abs(dett) > DETOL) then
        t1 = (-rd1m*d22 + rd2m*d12) / det
        t2 = ( rd1m*d12 - rd2m*d11) / det
      else
        t1 = -rd1m / d11
        t2 = ZERO
      endif

    endif

!   new vin = x (vin + beta res) + t1 (vin_1 + beta res_1) + t2 (vin_2 + beta res_2)

    x = UM - t1 - t2
    do i = 1,ns
      res(i) = beta*(res(i) + t1*c(i)) + t1*vmem(i,3)
    enddo
    if(t2 /= ZERO) then
      do i = 1,ns
        res(i) = res(i) + t2*(beta*d(i) + vmem(i,4))
      enddo
    endif
    do i = 1,ns
      vmem(i,4) = vmem(i,3)
      vmem(i,3) = vin(i)
      vin(i) = x*vin(i) + res(i)
    enddo

    deallocate(c, d)

  endif

  deallocate(res)

  return

end subroutine mixer_anderson_c16
