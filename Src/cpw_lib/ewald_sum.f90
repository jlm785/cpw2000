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

!>  Computes the ewald contribution to the
!>  total energy, forces and stresses.
!>
!>  \author       Jose Luis Martins,Sverre Froyen
!>  \version      4.94
!>  \date         ecember 16 1987.  6 June 2024.
!>  \copyright    GNU Public License v2

subroutine ewald_sum(enerew, forcew, strew,                              &
    adot, ntype, natom, rat, zv,                                         &
    mxdtyp, mxdatm)

! Adapted from Sverre Froyen plane wave program
! written December 16 1987. jlm
! modified December 9 1993. jlm
! modified March 1999. jlm
! modified for f90, 22 October 2015. JLM
! Modified documentation, January 2020. JLM
! Indentation, 5 June 2024. JLM


  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  mxdtyp                          !<  array dimension of types of atoms
  integer, intent(in)                ::  mxdatm                          !<  array dimension of number of atoms of a given type

  real(REAL64), intent(in)           ::  adot(3,3)                       !<  metric in direct space
  integer, intent(in)                ::  ntype                           !<  number of types of atoms
  integer, intent(in)                ::  natom(mxdtyp)                   !<  number of atoms of type i
  real(REAL64), intent(in)           ::  rat(3,mxdatm,mxdtyp)            !<  lattice coordinates of atom j of type i

  real(REAL64), intent(in)           ::  zv(mxdtyp)                      !<  valence of atom with type i

! output

  real(REAL64), intent(out)          ::  enerew                          !<  Ewald energy (Hartree)
  real(REAL64), intent(out)          ::  forcew(3,mxdatm,mxdtyp)         !<  d enerew / d rat,  Ewald contribution to force (contravariant components)
  real(REAL64), intent(out)          ::  strew(3,3)                      !<  d enerew / d adot,  Ewald contribution to stress tensor (contravariant components)

! local allocatable arrays

  real(REAL64), allocatable          ::  zz(:)
  real(REAL64), allocatable          ::  rc(:,:)
  real(REAL64), allocatable          ::  fsumr(:,:)
  real(REAL64), allocatable          ::  fsumg(:,:)
  real(REAL64), allocatable          ::  zcgdr(:)
  real(REAL64), allocatable          ::  zsgdr(:)

! local variables

  real(REAL64)             ::  bdot(3,3),vcell

  real(REAL64)             ::  expgi,dexpgi
  real(REAL64)             ::  eki
  real(REAL64)             ::  rp(3),xir(3)
  real(REAL64)             ::  ssumg(3,3),ssumr(3,3)
  real(REAL64)             ::  fsub(3),ssub(3,3),ssum0(3,3)
  real(REAL64)             ::  esub,esumg,esumr
  real(REAL64)             ::  rmod,arg,exp1,exp2,esum0
  real(REAL64)             ::  tol,eps,seps,sepi,rmax
  real(REAL64)             ::  zz1,zz2,sumc,sums,fac,gdr

  real(REAL64)             ::  gmax, xkg(3)

  integer                  ::  imx,jmx,kmx,natot,im
  integer                  ::  img, jmg, kmg

! parameters

  real(REAL64), parameter  :: PI = 3.14159265358979323846_REAL64
  real(REAL64), parameter  :: SMALL = 0.000000000001_REAL64
  real(REAL64), parameter  :: ZERO = 0.0_REAL64, UM = 1.0_REAL64

! counters

  integer       ::   i, j, k, l, m, n, ii, jj, kk

! functions

  real(REAL64),external    ::  erfc


  allocate(zz(mxdatm*mxdtyp))
  allocate(rc(3,mxdatm*mxdtyp))
  allocate(fsumr(3,mxdatm*mxdtyp))
  allocate(fsumg(3,mxdatm*mxdtyp))
  allocate(zcgdr(mxdtyp*mxdatm))
  allocate(zsgdr(mxdtyp*mxdatm))

! calculates cell volume

  bdot(1,1) = adot(2,2)*adot(3,3) - adot(3,2)*adot(2,3)
  bdot(2,1) = adot(3,2)*adot(1,3) - adot(1,2)*adot(3,3)
  bdot(3,1) = adot(1,2)*adot(2,3) - adot(2,2)*adot(1,3)
  bdot(1,2) = adot(2,3)*adot(3,1) - adot(3,3)*adot(2,1)
  bdot(2,2) = adot(3,3)*adot(1,1) - adot(1,3)*adot(3,1)
  bdot(3,2) = adot(1,3)*adot(2,1) - adot(2,3)*adot(1,1)
  bdot(1,3) = adot(2,1)*adot(3,2) - adot(3,1)*adot(2,2)
  bdot(2,3) = adot(3,1)*adot(1,2) - adot(1,1)*adot(3,2)
  bdot(3,3) = adot(1,1)*adot(2,2) - adot(2,1)*adot(1,2)

! cell volume (squared)

  vcell = bdot(1,1)*adot(1,1) + bdot(2,1)*adot(2,1) + bdot(3,1)*adot(3,1)
  if(vcell < SMALL) then
    write(6,'("    Stopped in ewald_sum   cell volume squared= ", e12.4)') vcell

    stop

  endif

  do j = 1,3
  do i = 1,3
    bdot(i,j) = 4*PI*PI*bdot(i,j)/vcell
  enddo
  enddo
  vcell = sqrt(vcell)

! compute various parameters

  tol = -log(SMALL)
  gmax = sqrt(max(bdot(1,1),bdot(2,2),bdot(3,3)))*tol/PI
  eps = gmax*gmax/(4*tol)
  img = int(gmax*sqrt(adot(1,1))/(2*PI)) + 1
  jmg = int(gmax*sqrt(adot(2,2))/(2*PI)) + 1
  kmg = int(gmax*sqrt(adot(3,3))/(2*PI)) + 1

  seps = sqrt(eps)
  sepi = 2*seps/sqrt(PI)
  rmax = sqrt(tol/eps)
  imx = int(rmax*sqrt(bdot(1,1))/(2*PI)) + 1
  jmx = int(rmax*sqrt(bdot(2,2))/(2*PI)) + 1
  kmx = int(rmax*sqrt(bdot(3,3))/(2*PI)) + 1

! store data for atoms into new arrays

  natot = 0
  do i = 1,ntype
    do j = 1,natom(i)
      natot = natot + 1
      zz(natot) = zv(i)
      rc(1,natot) = rat(1,j,i)
      rc(2,natot) = rat(2,j,i)
      rc(3,natot) = rat(3,j,i)
    enddo
  enddo

! translates to unit cell at origin

  do i = 1,natot
    rc(1,i) = rc(1,i) - real(int(rc(1,i)),REAL64)
    if(rc(1,i) < 0) rc(1,i) = rc(1,i) + UM
    rc(2,i) = rc(2,i) - real(int(rc(2,i)),REAL64)
    if(rc(2,i) < 0) rc(2,i) = rc(2,i) + UM
    rc(3,i) = rc(3,i) - real(int(rc(3,i)),REAL64)
    if(rc(3,i) < 0) rc(3,i) = rc(3,i) + UM
  enddo

  zz1 = ZERO
  zz2 = ZERO
  do i = 1,natot
    zz1 = zz1 + zz(i)
  enddo
  do i = 1,natot
    zz2 = zz2 + zz(i)*zz1
  enddo

! initialize sums : esum(g,r)   -   energy
!                   fsum(g,r)   -   force
!                   ssum(g,r)   -   stress

! note that the sums are in units of basis vectors
! (g) in units of bi's, and (r) in units of ai's

  esumg = ZERO
  esumr = ZERO
  do i = 1,3
  do j = 1,3
    ssumg(j,i)   = ZERO
    ssumr(j,i)   = ZERO
  enddo
  enddo
  do j = 1,natot
    fsumg(1,j) = ZERO
    fsumr(1,j) = ZERO
    fsumg(2,j) = ZERO
    fsumr(2,j) = ZERO
    fsumg(3,j) = ZERO
    fsumr(3,j) = ZERO
  enddo

! start sum in g space

! loop over g vectors

  do ii = -img,img
  xkg(1) = real(ii,REAL64)
  do jj = -jmg,jmg
  xkg(2) = real(jj,REAL64)
  do kk = -kmg,kmg
  xkg(3) = real(kk,REAL64)

    eki = ZERO
    do i = 1,3
    do j = 1,3
      eki = eki + xkg(j)*bdot(j,i)*xkg(i)
    enddo
    enddo

   if(ii == 0 .and. jj == 0 .and. kk == 0) then
      expgi = ZERO
      dexpgi = ZERO
    else
      expgi = exp(-eki/(4*eps))/eki
      dexpgi = - (UM/(2*eps) + 2*UM/eki)*expgi
    endif

    sumc = ZERO
    sums = ZERO

!   start loop over atoms in cell

    do i = 1,natot
      gdr = 2*PI*(xkg(1)*rc(1,i) + xkg(2)*rc(2,i) + xkg(3)*rc(3,i))
      zcgdr(i) = zz(i)*cos(gdr)
      zsgdr(i) = zz(i)*sin(gdr)
      sumc = sumc + zcgdr(i)
      sums = sums + zsgdr(i)
    enddo

    esumg = esumg + expgi*(sumc*sumc+sums*sums)

    do i = 1,natot
      fac = 2*expgi*(sumc*zsgdr(i) - sums*zcgdr(i))
      fsumg(1,i) = fsumg(1,i) + xkg(1)*fac
      fsumg(2,i) = fsumg(2,i) + xkg(2)*fac
      fsumg(3,i) = fsumg(3,i) + xkg(3)*fac
    enddo

    fac = dexpgi*(sumc*sumc + sums*sums)
    do i = 1,3
    do j = 1,3
      ssumg(j,i) = ssumg(j,i) + fac * xkg(j) * xkg(i)
    enddo
    enddo

  enddo
  enddo
  enddo

! end g sum

  esumg = 2*PI*esumg/vcell
  do i = 1,3
  do j = 1,3
    ssumg(j,i) = 2*PI*ssumg(j,i)/vcell
  enddo
  enddo
  do j = 1,natot
    fsumg(1,j) = 2*PI*fsumg(1,j)/vcell
    fsumg(2,j) = 2*PI*fsumg(2,j)/vcell
    fsumg(3,j) = 2*PI*fsumg(3,j)/vcell
  enddo


! start sum in r space

  esum0 = ZERO
  do i = 1,3
  do j = 1,3
    ssum0(j,i) = ZERO
  enddo
  enddo

  do ii = -imx,imx
  xir(1) = real(ii,REAL64)
  do jj = -jmx,jmx
  xir(2) = real(jj,REAL64)
  do kk = -kmx,kmx
  xir(3) = real(kk,REAL64)
    rmod = ZERO
    do i = 1,3
    do j = 1,3
      rmod = rmod + xir(j)*adot(j,i)*xir(i)
    enddo
    enddo
    if (ii /= 0 .or. jj /= 0 .or. kk /= 0) then
      rmod = sqrt(rmod)
      arg = seps*rmod
      if (arg < 25.0) then
        exp1 = erfc(arg) / rmod
        exp2 = (exp1 + sepi*exp(-arg*arg))/(rmod*rmod)
        esum0 = esum0 + exp1
        do i = 1,3
        do j = 1,3
          ssum0(j,i) = ssum0(j,i) + exp2 * xir(i)*xir(j)
        enddo
        enddo
      endif
    endif
  enddo
  enddo
  enddo
  esum0 = esum0 - pi/(eps*vcell) - sepi

! start loop over atoms in cell

  do i = 1,natot

!   term with a=b

    esumr = esumr + zz(i)*zz(i)*esum0/(2*UM)
    do j = 1,3
    do k = 1,3
      ssumr(k,j) = ssumr(k,j) + zz(i)*zz(i) * ssum0(k,j) / (2*UM)
    enddo
    enddo
    im = i-1
    if (im /= 0) then

!     terms with a#b

      do j = 1,im

!       loop over lattice points

        esub = ZERO
        do k = 1,3
          fsub(k) = ZERO
          do l = 1,3
            ssub(l,k) = ZERO
          enddo
        enddo

        do ii = -imx,imx
        xir(1) = real(ii,REAL64)
        do jj = -jmx,jmx
        xir(2) = real(jj,REAL64)
        do kk = -kmx,kmx
        xir(3) = real(kk,REAL64)
          rp(1) = xir(1) + rc(1,i) - rc(1,j)
          rp(2) = xir(2) + rc(2,i) - rc(2,j)
          rp(3) = xir(3) + rc(3,i) - rc(3,j)
          rmod = ZERO
          do l = 1,3
          do m = 1,3
            rmod = rmod + rp(l)*adot(l,m)*rp(m)
          enddo
          enddo
          rmod = sqrt(rmod)
          arg = seps*rmod
          if (arg < 25.0) then
            exp1 = erfc(arg) / rmod
            exp2 = (exp1 + sepi*exp(-arg*arg))/(rmod*rmod)
            esub = esub + exp1
            fsub(1) = fsub(1) + rp(1) * exp2
            fsub(2) = fsub(2) + rp(2) * exp2
            fsub(3) = fsub(3) + rp(3) * exp2
            do l = 1,3
            do m = 1,3
              ssub(m,l) = ssub(m,l) + rp(m) * exp2 * rp(l)
            enddo
            enddo
          endif
        enddo
        enddo
        enddo

        esub = esub - PI/(eps*vcell)
        esumr = esumr + zz(i)*zz(j)*esub
        do l = 1,3
        do m = 1,3
          ssumr(m,l) = ssumr(m,l) + zz(i)*zz(j)*ssub(m,l)
        enddo
        enddo
        do k = 1,3
          fsumr(k,i) = fsumr(k,i) + zz(i)*zz(j)*fsub(k)
          fsumr(k,j) = fsumr(k,j) - zz(i)*zz(j)*fsub(k)
        enddo
      enddo
    endif
  enddo

! end r sum

  enerew = esumg + esumr

! force

  n = 0
  do i = 1,ntype
    do j = 1,natom(i)
      n = n + 1
      do k = 1,3
        forcew(k,j,i) = fsumr(k,n) + (bdot(k,1)*fsumg(1,n) +             &
            bdot(k,2)*fsumg(2,n)+ bdot(k,3)*fsumg(3,n))/(2*PI)
      enddo
    enddo
  enddo

! stress

  do k = 1,3
  do j = 1,3
    strew(j,k) = ssumr(j,k) +                                            &
         ( bdot(j,1)*ssumg(1,1)*bdot(k,1)                                &
         + bdot(j,2)*ssumg(2,2)*bdot(k,2)                                &
         + bdot(j,3)*ssumg(3,3)*bdot(k,3)                                &
         + bdot(j,1)*ssumg(1,2)*bdot(k,2)                                &
         + bdot(j,2)*ssumg(2,3)*bdot(k,3)                                &
         + bdot(j,3)*ssumg(3,1)*bdot(k,1)                                &
         + bdot(j,2)*ssumg(2,1)*bdot(k,1)                                &
         + bdot(j,3)*ssumg(3,2)*bdot(k,2)                                &
         + bdot(j,1)*ssumg(1,3)*bdot(k,3) ) / (4*PI*PI)                  &
            + (esumg - zz2*pi/(2*eps*vcell))*bdot(j,k) / (4*PI*PI)
  enddo
  enddo

  deallocate(zz)
  deallocate(rc)
  deallocate(fsumr)
  deallocate(fsumg)
  deallocate(zcgdr)
  deallocate(zsgdr)

  return

end subroutine ewald_sum
