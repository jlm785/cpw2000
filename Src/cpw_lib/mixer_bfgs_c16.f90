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

!>     Computes the new exchange correlation potential
!>     from the input and the output potential.
!>
!>  \author       Jose Luis Martins
!>  \version      5.01
!>  \date         May 20,1999. Modified
!>  \copyright    GNU Public License v2


subroutine mixer_bfgs_c16(itmix, adot, ztot, bandwid, penngap, energy,   &
  ng, phase, conj, ns, mstar, ek,                                        &
  vhxc, vhxcout, delvhxc,                                                &
  mxdgve, mxdnst, mxdupd, mxdscf)

! Uses the Broyden-Fletcher-Goldfarb-Shanno method
! Numerical Recipes (2nd ed), Eq. 10.7.8.
! A good guess of the dielectric matrix is the key to its success.
! Stores only the rank-1 updates for mxdupd stars,
! for the other stars it uses a dielectric model.

! Adapted from Sverre Froyen plane wave program.
! Written May 20,1999.
! Modified 17 april 2006. JLM
! Modified, complex and f90, 22 September 2015. JLM
! Modified documentation, January 2020. JLM
! Modified mxdupd,mxdscf, February 2020. JLM
! Modified, dgrfi, vout inconsistent if not set to zero. 20 June 2021. JLM
! Modernized. 20 June 2021. JLM
! copyright INESC-MN/Jose Luis Martins

! Version 5.01

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  mxdgve                          !<  array dimension for g-space vectors
  integer, intent(in)                ::  mxdnst                          !<  array dimension for g-space stars

  integer, intent(in)                ::  mxdupd                          !<  array dimension for Jacobian update of scf potential
  integer, intent(in)                ::  mxdscf                          !<  array dimension of old scf-iteration information

  integer, intent(in)                ::  itmix                           !<  iteration number of mixer. Initial efault is iter in cpw_scf. itmix = 1 restarts.
  real(REAL64), intent(in)           ::  adot(3,3)                       !<  metric in direct space
  real(REAL64), intent(in)           ::  ztot                            !<  total charge density (electrons/cell)
  real(REAL64), intent(in)           ::  bandwid                         !<  occupied band width estimate (Hartree)
  real(REAL64), intent(in)           ::  penngap                         !<  Penn gap estimate (Hartree)
  real(REAL64), intent(in)           ::  energy                          !<  Kohn-Sham energy (hartree)

  integer, intent(in)                ::  ng                              !<  total number of g-vectors with length less than gmax
  complex(REAL64), intent(in)        ::  phase(mxdgve)                   !<  phase factor of G-vector n
  real(REAL64), intent(in)           ::  conj(mxdgve)                    !<  is -1 if one must take the complex conjugate of x*phase
  integer, intent(in)                ::  ns                              !<  number os stars with length less than gmax
  integer, intent(in)                ::  mstar(mxdnst)                   !<  number of g-vectors in the j-th star
  real(REAL64), intent(in)           ::  ek(mxdnst)                      !<  kinetic energy (hartree) of g-vectors in star j

  complex(REAL64), intent(in)        ::  vhxc(mxdnst)                    !<  input hxc potential (hartree) for the prototype G-vector
  complex(REAL64), intent(in)        ::  vhxcout(mxdnst)                 !<  output hxc potential (hartree) for the prototype G-vector

! output

  complex(REAL64), intent(out)       ::  delvhxc(mxdnst)                 !<  correction to input hxc potential (hartree) for the prototype G-vector

! local allocatable arrays

  real(REAL64), save, allocatable    ::  hessd(:)                        !  diagonal hessian (metric)
  complex(REAL64), save, allocatable ::  vold(:)                         !  old values of vhxc
  complex(REAL64), save, allocatable ::  voldout(:)                      !  old values of vhxcout

  complex(REAL64), save, allocatable ::  hy(:,:)                         !  H_i x (grad f_i+1 - gr!ad f_i)
  complex(REAL64), save, allocatable ::  pvec(:,:)                       !  x_i+1 - x_i (normalized or not normalized)
  complex(REAL64), save, allocatable ::  uvec(:,:)                       !  vector u (Eq. 10.7.10  Num Rec - 2nd edition)
  complex(REAL64), allocatable       ::  dgrfi(:)                        !  cf_i+1 - grad f_i

! local variables

  integer             ::  icount, ic
  integer             ::  nj
  real(REAL64)        ::  vcell,bdot(3,3)
  real(REAL64)        ::  xmix, xmix0, xmixmin
  real(REAL64), save  ::  xmixmax = 1.0_REAL64
  real(REAL64), save  ::  oldener
  real(REAL64)        ::  xsum
  real(REAL64)        ::  vnorm1,vnorm2
  complex(REAL64)     ::  vnew

! counters

  integer       ::  i, j, k

! parameters

  real(REAL64), parameter :: EPS = 0.00001_REAL64
  real(REAL64), parameter :: EPSVN = 0.00000001_REAL64
  real(REAL64), parameter :: ZERO = 0.0_REAL64, UM = 1.0_REAL64
  complex(REAL64), parameter  ::  C_ZERO = cmplx(ZERO,ZERO,REAL64)

! hard coded parameters, could be input parameters

!  mxdupd = max(min(mxdnst/10,1000),min(mxdnst,100))
!  mxdupd = max(mxdnst/2,min(mxdnst,100))
!  mxdupd = mxdnst
!  mxdscf = 15


! paranoid check, stops complaint about unused ng...

  if(ng > mxdgve) stop

  icount = itmix - 1

  if(icount == 0) oldener = energy

! allocations on first iteration

  if(icount == 0) then
    if(allocated(hessd)) deallocate(hessd)
    allocate(hessd(mxdupd))
    if(allocated(vold)) deallocate(vold)
    allocate(vold(mxdupd))
    if(allocated(voldout)) deallocate(voldout)
    allocate(voldout(mxdupd))
    if(allocated(hy)) deallocate(hy)
    allocate(hy(mxdupd,mxdscf))
    if(allocated(pvec)) deallocate(pvec)
    allocate(pvec(mxdupd,mxdscf))
    if(allocated(uvec)) deallocate(uvec)
    allocate(uvec(mxdupd,mxdscf))
  endif

  if(icount == 0) xmixmax = UM

  call adot_to_bdot(adot, vcell, bdot)

! Thomas-Fermi wave-vector squared, Slater alpha, xc-screening, Fermi wave-vector
! are initialized for later screening

  call mixer_screening(xmix, .TRUE., ek(1), ztot, vcell)

! estimate of xmix limits with some "empirical" parameters

  xmixmin = penngap*penngap / (penngap*penngap + bandwid*bandwid)
  if(xmixmin > 0.3) xmixmin = 0.3
  xmix0 = 0.3

! reduces xmixmax

  if(icount > 0) then
    if(energy > oldener) then
      xmixmax = 0.8*xmixmax
      if(xmixmax < 2.0*xmixmin) xmixmax = 2.0*xmixmin
      if(xmixmax < 0.1) xmixmax = 0.1
    endif
  endif

  nj = ns

  if(nj > mxdupd) then
    do i = 1,mxdupd
      nj = mxdupd + 1 - i

      if(ek(nj)+eps < ek(mxdupd+1)) exit

    enddo

    if(nj == 1) then
      write(6,*)
      write(6,*) '   WARNING    mixer is using simple mixing '
      write(6,'("increase mxdupd from ",i8)') mxdupd
      write(6,*)
    endif
  endif

! not enough memory, forget older vectors

  if(icount > mxdscf) then
    icount = mxdscf

    if(icount > 1) then
      do k=2,icount
        do j=2,nj
          hy(j,k-1) = hy(j,k)
          pvec(j,k-1) = pvec(j,k)
          uvec(j,k-1) = uvec(j,k)
        enddo
      enddo
    endif

  endif

  if (nj > 1) then

!   initialize diagonal of inverse jacobian matrix hess

    if (icount == 0) then
      hessd(1) = ZERO
      do i=2,nj
        call mixer_screening(xmix, .FALSE., ek(i), ztot, vcell)
        if(xmix < xmixmin) xmix = xmixmin
        if(xmix > xmix0) xmix = xmix0
        hessd(i) = - xmix
      enddo
    endif


!   update inverse jacobian matrix h - p x hy

    if (icount > 0) then

      allocate(dgrfi(nj))

      do i = 2,nj
        dgrfi(i) = (vhxcout(i)-vhxc(i)) - (voldout(i)-vold(i))
      enddo

      do i=2,nj
          hy(i,icount) = hessd(i) * dgrfi(i)
          pvec(i,icount) = (vhxc(i)-vold(i))
      enddo

      if(icount > 1) then

        do k = 1,icount-1
          xsum = ZERO
          do j = 2,nj
            xsum = xsum + real(hy(j,k) * conjg(dgrfi(j)),REAL64)
          enddo
          do i = 2,nj
            hy(i,icount) = hy(i,icount) + hy(i,k) * xsum
          enddo
        enddo

        do k = 1,icount-1
          xsum = ZERO
          do j = 2,nj
            xsum = xsum + real(pvec(j,k) * conjg(dgrfi(j)),REAL64)
          enddo
          do i = 2,nj
            hy(i,icount) = hy(i,icount) - xsum * pvec(i,k)
          enddo
        enddo

        do k = 1,icount-1
          xsum = ZERO
          do j = 2,nj
            xsum = xsum + real(uvec(j,k) * conjg(dgrfi(j)),REAL64)
          enddo
          do i = 2,nj
            hy(i,icount) = hy(i,icount) - xsum * uvec(i,k)
          enddo
        enddo

      endif

      vnorm1 = ZERO
      do i = 2,nj
        vnorm1 = vnorm1 + real((vhxc(i)-vold(i)) * conjg(dgrfi(i)),REAL64)
      enddo

      vnorm2 = ZERO
      do i = 2,nj
        vnorm2 = vnorm2 + real(hy(i,icount) * conjg(dgrfi(i)),REAL64)
      enddo

      deallocate(dgrfi)

      if(-vnorm1 <= EPSVN .or. -vnorm2 <= EPSVN) then

!       we are in the noise

        do i = 2,nj
          uvec(i,icount) = C_ZERO
          pvec(i,icount) = C_ZERO
          hy(i,icount) = C_ZERO
        enddo

      else

!       renormalizes stuff

        do i = 2,nj
          uvec(i,icount) = pvec(i,icount) / vnorm1 - hy(i,icount) / vnorm2
        enddo

        do i = 2,nj
          pvec(i,icount) = pvec(i,icount) / sqrt(-vnorm1)
        enddo

        do i = 2,nj
          hy(i,icount) = hy(i,icount) / sqrt(-vnorm2)
        enddo

        do i = 2,nj
          uvec(i,icount) = uvec(i,icount) * sqrt(-vnorm2)
        enddo

      endif

    endif                                                                !  icount > 0


!   save vhxc and vhxcout for next iteration

    do i = 2,nj
      vold(i) = vhxc(i)
      voldout(i) = vhxcout(i)
    enddo

!   compute correction to vhxc from inverse jacobian

    do i = 2,nj
      delvhxc(i) = -hessd(i) * (vhxcout(i)-vhxc(i))
    enddo


    if(icount > 0) then
      do k = 1,icount

        xsum = ZERO
        do j=2,nj
          xsum = xsum + real(hy(j,k) * conjg(vhxcout(j)-vhxc(j)),REAL64)
        enddo

        do i = 2,nj
          delvhxc(i) = delvhxc(i) - xsum * hy(i,k)
        enddo

        xsum = ZERO
        do j = 2,nj
          xsum = xsum + real(pvec(j,k) * conjg(vhxcout(j)-vhxc(j)),REAL64)
        enddo

        do i = 2,nj
          delvhxc(i) = delvhxc(i) + xsum * pvec(i,k)
        enddo

        xsum = ZERO
        do j = 2,nj
          xsum = xsum + real(uvec(j,k) * conjg(vhxcout(j)-vhxc(j)),REAL64)
        enddo

        do i = 2,nj
          delvhxc(i) = delvhxc(i) + xsum * uvec(i,k)
        enddo

      enddo

    endif

  endif


! loop over remaining stars

  delvhxc(1) = C_ZERO
  if (ns > nj) then
    do i=nj+1,ns
      call mixer_screening(xmix, .FALSE., ek(i), ztot, vcell)
      if(xmix < xmixmin) xmix = xmixmin
      if(xmix > xmixmax) xmix = xmixmax
      delvhxc(i) = xmix*(vhxcout(i)-vhxc(i))
    enddo
  endif

! makes sure the update is real

  ic = 1
  do i = 2,ns
    if(conj(ic+2) > 0.0) then
      vnew = delvhxc(i)*conjg(phase(ic+2))
      delvhxc(i) = (delvhxc(i) + conjg(vnew)) / 2
    endif
    ic = ic + mstar(i)
  enddo

  oldener = energy

  return
end subroutine mixer_bfgs_c16
