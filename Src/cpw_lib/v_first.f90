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

!>  Calculates the starting potential and the values of the
!>  pseudopotential in the G-vectors. Scales by vcell
!>
!>  \author       Jose Luis Martins
!>  \version      5.11
!>  \date         October 1994, March 5 2025.
!>  \copyright    GNU Public License v2

  subroutine v_first(ns, ek, sfact, ealpha, ealraw,                      &
      nq, delq, vloc, dcor, dval,                                        &
      ntype, adot,                                                       &
      vion, denc, dens, vql, dvql, dnc, ddc,                             &
      mxdtyp, mxdlqp, mxdnst)

! written october 94. jlm
! modified 22 march 1999. jlm
! modified 21 february 2000. jlm
! modified 29 january 2008 (wvfraw). jlm
! modified 28 september 2013, vionr(1). jlm
! Modified 25 October 2015. f90. removed vkb.  JLM
! Modified documentation, January 2020. JLM
! Modified, indentation, dens(1)=0. JLM

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  mxdtyp                          !<  array dimension of types of atoms
  integer, intent(in)                ::  mxdlqp                          !<  array dimension for local potential
  integer, intent(in)                ::  mxdnst                          !<  array dimension for g-space stars

  integer, intent(in)                ::  ns                              !<  number os stars with length less than gmax
  real(REAL64), intent(in)           ::  ek(mxdnst)                      !<  kinetic energy (hartree) of g-vectors in star j
  complex(REAL64), intent(in)        ::  sfact(mxdtyp,mxdnst)            !<  real part of the structure factor

  real(REAL64), intent(in)           ::  ealraw                          !<  G=0 contrib. to the total energy. (non norm. to vcell, Hartree)

  integer, intent(in)                ::  nq(mxdtyp)                      !<  number of points for the local pseudopotential interpolation
  real(REAL64), intent(in)           ::  delq(mxdtyp)                    !<  step used in the interpolation
  real(REAL64), intent(in)           ::  vloc(-1:mxdlqp,mxdtyp)          !<  local pseudopotential for atom k (hartree)
  real(REAL64), intent(in)           ::  dcor(-1:mxdlqp,mxdtyp)          !<  core charge density for atom k
  real(REAL64), intent(in)           ::  dval(-1:mxdlqp,mxdtyp)          !<  valence charge density for atom k

  integer, intent(in)                ::  ntype                           !<  number of types of atoms
  real(REAL64), intent(in)           ::  adot(3,3)                       !<  metric in real space

! output

  real(REAL64), intent(out)          ::  ealpha                          !<  G=0 contribution to the total energy (Hartree)

  complex(REAL64), intent(out)       ::  vion(mxdnst)                    !<  ionic potential for the prototype G-vector in star j
  complex(REAL64), intent(out)       ::  denc(mxdnst)                    !<  core charge density for the prototype G-vector
  complex(REAL64), intent(out)       ::  dens(mxdnst)                    !<  spherical atomic valence charge density for the prototype G-vector in star j

  real(REAL64), intent(out)          ::  vql(mxdtyp,mxdnst)              !<  local pseudopotential for atom type i and prototype g-vector in star j
  real(REAL64), intent(out)          ::  dnc(mxdtyp,mxdnst)              !<  core charge for atom type i and prototype g-vector in star j
  complex(REAL64), intent(out)       ::  dvql(mxdnst)                    !<  derivative of the local pseudopotential for the prototype g-vector in star j
  complex(REAL64), intent(out)       ::  ddc(mxdnst)                     !<  derivative of the core charge for the prototype g-vector in star j

! local variables

  real(REAL64)           ::  vionr1
  integer                ::  n, nql
  real(REAL64)           ::  vcell, bdot(3,3)
  real(REAL64)           ::  gmax,glmax,delql
  real(REAL64)           ::  qj,xn,q2vn,q2vp,q2vm
  real(REAL64)           ::  vqj,dvqj,dcj,ddcj,dvj

! constants

  real(REAL64), parameter     ::  ZERO = 0.0_REAL64, UM = 1.0_REAL64
  complex(REAL64), parameter  ::  C_ZERO = cmplx(ZERO,ZERO,REAL64)
  real(REAL64), parameter     ::  SMALL = 1.0E-10_REAL64

! counters

  integer       ::  i, j, nt


! scaling

  call adot_to_bdot(adot,vcell,bdot)

  ealpha = ealraw/vcell
!  fac = UM / sqrt(vcell)

! initialize arrays

  do i = 1,ns
    vion(i) = C_ZERO
    denc(i) = C_ZERO
    dens(i) = C_ZERO
    dvql(i) = C_ZERO
    ddc(i) = C_ZERO
  enddo

  do nt = 1,ntype

!   check argument limits and compare with gmax

    delql = delq(nt)
    nql = nq(nt)
    gmax = sqrt(2*ek(ns))
    glmax = delql * (nql-1)
    if (glmax < gmax) then
      write(6,'("  *** warning  in v_first tape contains ",              &
         &  "information about local potential up to gmax = ",e10.3)')   &
            glmax
    endif

!   compute ionic local potential

    do j = 2,ns
      vql(nt,j) = ZERO

!     quadratic interpolation of q**2 * potential

      qj = sqrt(2*ek(j))
      xn = qj / delql
      n  =  int(xn + UM/2)
      if (n < nql) then
        if (n < 0) n = 0
        xn = xn - n*UM
        q2vn = (n*n)*vloc(n,nt)
        q2vp = ((n+1)*(n+1))*vloc(n+1,nt)
        q2vm = ((n-1)*(n-1))*vloc(n-1,nt)
        vqj  = q2vn * (UM-xn) * (UM+xn)                                  &
             + (UM/2) * (q2vp*(UM+xn) - q2vm*(UM-xn)) * xn
        vqj  = delql*delql * vqj / (2* vcell * ek(j))
        dvqj = -2 * q2vn * xn                                            &
             + q2vp * ((UM/2)+xn) - q2vm * ((UM/2)-xn)
        dvqj = (delql*dvqj/(2*vcell*qj) - vqj) / (2*ek(j))
        vion(j) = vion(j) + vqj*sfact(nt,j)
        dvql(j) = dvql(j) + dvqj*sfact(nt,j)
        vql(nt,j) = vqj
      ENDIF
    ENDDO

!   compute core charge

    do j = 1,ns
      dnc(nt,j) = ZERO

!     interpolate vda

      qj = sqrt(2*ek(j))
      xn = qj/delql
      n = int(xn + UM/2)
      if (n < nql) then
        if(n < 0) n = 0
        xn = xn - n*UM
        dcj = dcor(n,nt) * (UM+xn) * (UM-xn)                             &
            + (UM/2) * (dcor(n+1,nt)*(UM+xn)                             &
                       -dcor(n-1,nt)*(UM-xn)) * xn
        dnc(nt,j) = dcj
        denc(j) = denc(j) + dcj*sfact(nt,j)
        if (j /= 1) then
          ddcj = -2*dcor(n,nt) * xn                                      &
               + dcor(n+1,nt)*((UM/2)+xn)                                &
               - dcor(n-1,nt)*((UM/2)-xn)
          ddcj = ddcj/(2*delql*qj)
          ddc(j) = ddc(j) + ddcj*sfact(nt,j)
        endif
      endif
    enddo
    ddcj = (dcor(1,nt)-dcor(0,nt))  / (delql*delql)
    ddcj = ddcj * real(sfact(nt,1),REAL64)
    ddc(1) = ddc(1) + cmplx(ddcj,ZERO,REAL64)

!   compute valence charge

    do j = 1,ns

!     interpolate vda

      xn = sqrt(2*ek(j))/delql
      n = int(xn + UM/2)
      if (n < nql) then

        if(n < 0) n = 0
        xn = xn - n*UM
        dvj = dval(n,nt) * (UM+xn) * (UM-xn)                             &
            + (UM/2) * (dval(n+1,nt)*(UM+xn)                             &
            - dval(n-1,nt)*(UM-xn)) * xn

!       sum up the charge density

        dens(j) = dens(j) + dvj*sfact(nt,j)

      endif
    enddo

!   end loop over atomic types

  enddo

! average value of electrostatic potential
! beware of dens(1)=0 (to allow empty lattice calculations)

  if(abs(dens(1)) < SMALL) then
    vionr1 = ZERO
    if(abs(ealraw) > SMALL) then
      write(6,*)
      write(6,*)  '    WARNING   WARNING   WATNING   '
      write(6,*)  '    There are no electrons but raw ealpha is not zero'
      write(6,*)
    endif
  else
    vionr1 = ealraw / (vcell*real(dens(1),REAL64))
  endif
  vion(1) = cmplx(vionr1,ZERO,REAL64)

  return

end subroutine v_first
