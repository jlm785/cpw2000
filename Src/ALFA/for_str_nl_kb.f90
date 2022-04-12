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

!>  computes and adds contributions to the
!>  Hellman-Feynman forces and stress from the eigenvectors
!>  non-local Kleinman and Bylander type of pseudopotential
!>
!>  \author       Jose Luis Martins
!>  \version      5.0.5
!>  \date         May 12 1990, 12 April 2022.
!>  \copyright    GNU Public License v2

subroutine for_str_nl_kb(fnlkb, strnlkb,                                 &
    mtxd, rkpt, neig, occp,                                              &
    isort, psi,                                                          &
    kgv,                                                                 &
    nqnl, delqnl, vkb, nkb,                                              &
    ntype, natom, rat ,adot,                                             &
    mxdtyp,mxdatm,mxdgve,mxdlqp,mxddim,mxdbnd)

! Written May 12 1990. jlm
! Modified March 7 1999. jlm
! Modified January 5, 2017, f90. JLM
! Modified, documentation, January 2020. JLM
! Added the f non-local contributions. 2 December 2021. JLM
! Corrected lmx bug. 12 April 2022. JLM


  implicit none

  integer, parameter :: REAL64 = selected_real_kind(12)

  integer, parameter :: LMAX = 3                                         !  hard coded max l in pseudo

! input

  integer, intent(in)                ::  mxdtyp                          !<  array dimension of types of atoms
  integer, intent(in)                ::  mxdatm                          !<  array dimension of number of atoms of a given type
  integer, intent(in)                ::  mxdgve                          !<  array dimension for g-space vectors
  integer, intent(in)                ::  mxdlqp                          !<  array dimension for local potential
  integer, intent(in)                ::  mxddim                          !<  array dimension of plane-waves
  integer, intent(in)                ::  mxdbnd                          !<  array dimension for number of bands

  integer, intent(in)                ::  mtxd                            !<  dimension of the hamiltonian
  real(REAL64), intent(in)           ::  rkpt(3)                         !<  component in lattice coordinates of the k-point
  integer, intent(in)                ::  neig                            !<  number of eigenvectors (requested on input, modified by degeneracies on output)
  real(REAL64), intent(in)           ::  occp(mxdbnd)                    !<  fractional ocupation of level j

  integer, intent(in)                ::  isort(mxddim)                   !<  g-vector associated with row/column i of hamiltonian
  complex(REAL64), intent(in)        ::  psi(mxddim,mxdbnd)              !<  component j of vector i

  integer, intent(in)                ::  kgv(3,mxdgve)                   !<  i-th component (reciprocal lattice coordinates) of the n-th g-vector ordered by stars of increasing length

  integer, intent(in)                ::  nqnl(mxdtyp)                    !<  number of points for the non-local pseudopotential interpolation
  real(REAL64), intent(in)           ::  delqnl(mxdtyp)                  !<  step used in the interpolation
  real(REAL64), intent(in)        ::  vkb(-2:mxdlqp,0:LMAX,-1:1,mxdtyp)  !<  (1/q**l) * KB nonlocal pseudo. for atom k, ang. mom. l. normalized to vcell, hartree
  integer, intent(in)                ::  nkb(0:LMAX,-1:1,mxdtyp)         !<  KB pseudo.  normalization for atom k, ang. mom. l

  integer, intent(in)                ::  ntype                           !<  number of types of atoms
  integer, intent(in)                ::  natom(mxdtyp)                   !<  number of atoms of type i
  real(REAL64), intent(in)           ::  rat(3,mxdatm,mxdtyp)            !<  lattice coordinates of atom j of type i
  real(REAL64), intent(in)           ::  adot(3,3)                       !<  metric in real space

! output

  real(REAL64), intent(out)          ::  fnlkb(3,mxdatm,mxdtyp)          !<  unsymmetrized  non-local pseudopotential contribution to the force (Hartree/Bohr)
  real(REAL64), intent(out)          ::  strnlkb(3,3)                    !<  unsymmetrized non-local pseudopotential contribution to the  covariant stress tensor (Hartree,Bohr)

! allocatable arrays

  complex(REAL64), allocatable       ::  anl(:,:,:)
  real(REAL64), allocatable          ::  xnlkb(:)
  complex(REAL64), allocatable       ::  dhd(:,:,:)

! local variables

  real(REAL64)         ::  vcell,bdot(3,3),avec(3,3),bvec(3,3)
  integer              ::  nanl
  real(REAL64)         ::  qk(3),qcar(3),qi,qinv
  real(REAL64)         ::  dqi(6),dqcar(6,3)
  real(REAL64)         ::  fi,xni,xi
  complex(REAL64)      ::  st, cfac
  real(REAL64)         ::  fac, vfac
  real(REAL64)         ::  vq,dvq
  real(REAL64)         ::  enl,enonlo
  integer              ::  ni,lm,lmmin,lmmax,ind
  real(REAL64)         ::  vint, dvintdxi
  integer              ::  ifail
  integer              ::  lmx, lsize

! variables associated with spherical harmonics

  real(REAL64)         ::  flm((LMAX+1)*(LMAX+1))
  real(REAL64)         ::  dflm(6,(LMAX+1)*(LMAX+1))
  real(REAL64)         ::  qitol(-1:LMAX)

  real(REAL64)         ::  rylm(0:LMAX,-LMAX:LMAX)                  !  r**l sqrt(4 pi / 2l+1) Y_lm(theta,phi), REAL VERSION
  real(REAL64)         ::  drylmdrc(3,0:LMAX,-LMAX:LMAX)            !  d rylm / d rc
  real(REAL64)         ::  d2rylmdrc2(3,3,0:LMAX,-LMAX:LMAX)        !  d^2 rylm / d rc^2

  complex(REAL64)      ::  zylm(0:LMAX,-LMAX:LMAX)                  !  r**l sqrt(4 pi / 2l+1) Y_lm(theta,phi)
  complex(REAL64)      ::  dzylmdrc(3,0:LMAX,-LMAX:LMAX)            !  d zylm / d rc
  complex(REAL64)      ::  d2zylmdrc2(3,3,0:LMAX,-LMAX:LMAX)        !  d^2 zylm / d rc^2

! constants

  real(REAL64), parameter   :: ZERO = 0.0_REAL64, UM = 1.0_REAL64
  real(REAL64), parameter   :: PI = 3.14159265358979323846_REAL64
  complex(REAL64), parameter  ::  C_UM = cmplx(UM,ZERO,REAL64)
  complex(REAL64), parameter  ::  C_ZERO = cmplx(ZERO,ZERO,REAL64)
  complex(REAL64), parameter  ::  C_I = cmplx(ZERO,UM,REAL64)
  real(REAL64), parameter  :: EPS = 1.0D-12


! counters

  integer       ::  i, j, k, l, m, kk, kl, n


! initialize fnlkb and strnlkb

  enonlo = ZERO

  do i = 1,3
  do j = 1,3
    strnlkb(j,i) = ZERO
  enddo
  enddo

  do i = 1,ntype
    do j = 1,natom(i)
    do k = 1,3
      fnlkb(k,j,i) = ZERO
    enddo
    enddo
  enddo

  if (neig < 1) return

! finds real l maximum and allocates arrays

  lmx = 0
  do k = 1,ntype
    do l = 0,LMAX
      if(nkb(l,0,k) /= 0 .and. l > lmx) lmx = l
    enddo
  enddo
  lsize = (lmx+1)*(lmx+1)

  allocate(anl(10,lsize,mtxd))
  allocate(dhd(10,lsize,neig))
  allocate(xnlkb(lsize))

  call adot_to_bdot(adot,vcell,bdot)
  call adot_to_avec(adot,avec,bvec)
  vfac = UM / sqrt(vcell)


  do k = 1,ntype

    nanl = 0
    do l = 0,lmx
      if(nkb(l,0,k) /= 0) nanl = nanl + 2*l+1
    enddo

    do kk = 1,natom(k)

!     calculates the matrix elements

      do i = 1,mtxd

        qk(1) = rkpt(1) + UM*kgv(1,isort(i))
        qk(2) = rkpt(2) + UM*kgv(2,isort(i))
        qk(3) = rkpt(3) + UM*kgv(3,isort(i))
        qi = (qk(1)*bdot(1,1) + qk(2)*bdot(2,1) + qk(3)*bdot(3,1))*qk(1) +    &
             (qk(1)*bdot(1,2) + qk(2)*bdot(2,2) + qk(3)*bdot(3,2))*qk(2) +    &
             (qk(1)*bdot(1,3) + qk(2)*bdot(2,3) + qk(3)*bdot(3,3))*qk(3)

!       derivatives and angular functions

        if(qi > EPS) then

          qi = sqrt(qi)
          qinv = UM / qi

!         derivatives of qi with respect to strain
!         in lattice coordinates

          dqi(1) = - qk(1)*qk(1)*qinv
          dqi(2) = - qk(2)*qk(2)*qinv
          dqi(3) = - qk(3)*qk(3)*qinv
          dqi(4) = - qk(1)*qk(2)*qinv
          dqi(5) = - qk(2)*qk(3)*qinv
          dqi(6) = - qk(3)*qk(1)*qinv

!         qvector in cartesian coordinates and normalized

          qcar(1) = bvec(1,1)*qk(1)+bvec(1,2)*qk(2)+bvec(1,3)*qk(3)
          qcar(2) = bvec(2,1)*qk(1)+bvec(2,2)*qk(2)+bvec(2,3)*qk(3)
          qcar(3) = bvec(3,1)*qk(1)+bvec(3,2)*qk(2)+bvec(3,3)*qk(3)
          qcar(1) = qcar(1)*qinv
          qcar(2) = qcar(2)*qinv
          qcar(3) = qcar(3)*qinv

!         derivatives of qcar with respect to strain
!         in lattice coordinates

          dqcar(1,1) = - avec(1,1) * qk(1) / (2*PI)
          dqcar(2,1) = - avec(1,2) * qk(2) / (2*PI)
          dqcar(3,1) = - avec(1,3) * qk(3) / (2*PI)
          dqcar(4,1) = - avec(1,1) * qk(2) / (2*PI)
          dqcar(5,1) = - avec(1,2) * qk(3) / (2*PI)
          dqcar(6,1) = - avec(1,3) * qk(1) / (2*PI)
          dqcar(1,2) = - avec(2,1) * qk(1) / (2*PI)
          dqcar(2,2) = - avec(2,2) * qk(2) / (2*PI)
          dqcar(3,2) = - avec(2,3) * qk(3) / (2*PI)
          dqcar(4,2) = - avec(2,1) * qk(2) / (2*PI)
          dqcar(5,2) = - avec(2,2) * qk(3) / (2*PI)
          dqcar(6,2) = - avec(2,3) * qk(1) / (2*PI)
          dqcar(1,3) = - avec(3,1) * qk(1) / (2*PI)
          dqcar(2,3) = - avec(3,2) * qk(2) / (2*PI)
          dqcar(3,3) = - avec(3,3) * qk(3) / (2*PI)
          dqcar(4,3) = - avec(3,1) * qk(2) / (2*PI)
          dqcar(5,3) = - avec(3,2) * qk(3) / (2*PI)
          dqcar(6,3) = - avec(3,3) * qk(1) / (2*PI)

          do kl = 1,6
            dqcar(kl,1) = (dqcar(kl,1) - dqi(kl)*qcar(1)) * qinv
            dqcar(kl,2) = (dqcar(kl,2) - dqi(kl)*qcar(2)) * qinv
            dqcar(kl,3) = (dqcar(kl,3) - dqi(kl)*qcar(3)) * qinv
          enddo

!         angular functions

          call ylm_rc(qcar,rylm,drylmdrc,d2rylmdrc2,.FALSE.,             &
                   zylm,dzylmdrc,d2zylmdrc2,LMAX,1,ifail)

          do l = 0,lmx
          do m = -l,l
            ind = 1+l*l+l-m
            flm(ind) = rylm(l,m)
            do kl = 1,6
              dflm(kl,ind) = drylmdrc(1,l,m)*dqcar(kl,1) +               &
                             drylmdrc(2,l,m)*dqcar(kl,2) +               &
                             drylmdrc(3,l,m)*dqcar(kl,3)
            enddo
          enddo
          enddo

        else

          qi = ZERO
          do kl = 1,6
            dqi(kl) = ZERO
          enddo
          do lm = 1,lsize
            flm(lm) = ZERO
            do kl = 1,6
              dflm(kl,lm) = ZERO
            enddo
          enddo
          flm(1) = UM

        endif

!       "powers" of qi qitol(-1) = 0

        qitol(-1) = ZERO
        qitol(0) = UM
        do l = 1,lmx
          qitol(l) = qi*qitol(l-1)
        enddo

!       structure factor

        fi = kgv(1,isort(i))*rat(1,kk,k) +                               &
             kgv(2,isort(i))*rat(2,kk,k) +                               &
             kgv(3,isort(i))*rat(3,kk,k)
        st = exp(2*PI*C_I * fi)

        ind = 0
        do l = 0,lmx

          if(nkb(l,0,k) /= 0) then
            xni = qi/delqnl(k)
            ni  = nint(xni)
            if(ni < 0) ni = 0
            vint = ZERO
            dvintdxi = ZERO
            if(ni < nqnl(k)-1) then
              xi  = xni - ni*UM
              vint = vkb(ni,l,0,k) * (UM+xi) * (UM-xi)                   &
                 + (UM/2) * (vkb(ni+1,l,0,k)*(UM+xi)                     &
                           -vkb(ni-1,l,0,k)*(UM-xi)) * xi
              dvintdxi = -2 * xi *vkb(ni,l,0,k) +                        &
                       (xi + UM/2) * vkb(ni+1,l,0,k) +                   &
                       (xi - UM/2) * vkb(ni-1,l,0,k)

                vq = vfac * vint * qitol(l)
                dvq = vfac * dvintdxi * qitol(l) / delqnl(k) +           &
                      vfac * vint * l * qitol(l-1)

            endif

            lmmin = l*l + 1
            lmmax = (l+1)*(l+1)

            do lm = lmmin,lmmax
              ind = ind + 1
              anl(1,ind,i) = st*flm(lm)*vq
              xnlkb(ind) = nkb(l,0,k)*UM

!             forces

              cfac = -C_I*st*flm(lm)*vq
              anl(2,ind,i) =  cfac*kgv(1,isort(i))
              anl(3,ind,i) =  cfac*kgv(2,isort(i))
              anl(4,ind,i) =  cfac*kgv(3,isort(i))

!             stresses

              do kl = 1,6
                anl(4+kl,ind,i) = st * (dflm(kl,lm)*vq + flm(lm)*dvq*dqi(kl))
              enddo

            enddo

          endif
!
        enddo

      enddo

!     calculates the matrix vector product

      call zgemm('n','n',10*nanl,neig,mtxd,C_UM,anl,10*lsize,            &
                 psi,mxddim,C_ZERO,dhd,10*lsize)

      do n = 1,neig
        do i = 1,nanl


          enl = occp(n) * xnlkb(i) * real(dhd(1,i,n)*conjg(dhd(1,i,n)))

          enonlo = enonlo + enl

          fnlkb(1,kk,k) = fnlkb(1,kk,k) + 2 * occp(n) * xnlkb(i) *       &
                    real(dhd(2,i,n)*conjg(dhd(1,i,n)),REAL64)
          fnlkb(2,kk,k) = fnlkb(2,kk,k) + 2 * occp(n) * xnlkb(i) *       &
                    real(dhd(3,i,n)*conjg(dhd(1,i,n)),REAL64)
          fnlkb(3,kk,k) = fnlkb(3,kk,k) + 2 * occp(n) * xnlkb(i) *       &
                    real(dhd(4,i,n)*conjg(dhd(1,i,n)),REAL64)

          fac = 8*PI*PI * occp(n) * xnlkb(i)
          strnlkb(1,1) = strnlkb(1,1) + enl * adot(1,1) -                &
                   fac* real(dhd(5,i,n)*conjg(dhd(1,i,n)),REAL64)
          strnlkb(2,2) = strnlkb(2,2) + enl * adot(2,2) -                &
                   fac* real(dhd(6,i,n)*conjg(dhd(1,i,n)),REAL64)
          strnlkb(3,3) = strnlkb(3,3) + enl * adot(3,3) -                &
                   fac* real(dhd(7,i,n)*conjg(dhd(1,i,n)),REAL64)
          strnlkb(1,2) = strnlkb(1,2) + enl * adot(1,2) -                &
                   fac * real(dhd(8,i,n)*conjg(dhd(1,i,n)),REAL64)
          strnlkb(2,3) = strnlkb(2,3) + enl * adot(2,3) -                &
                   fac * real(dhd(9,i,n)*conjg(dhd(1,i,n)),REAL64)
          strnlkb(3,1) = strnlkb(3,1) + enl * adot(3,1) -                &
                   fac * real(dhd(10,i,n)*conjg(dhd(1,i,n)))

        enddo
      enddo

    enddo
  enddo

  deallocate(anl)
  deallocate(dhd)
  deallocate(xnlkb)

  strnlkb(2,1) = strnlkb(1,2)
  strnlkb(3,2) = strnlkb(2,3)
  strnlkb(1,3) = strnlkb(3,1)

  do i = 1,ntype
    do j = 1,natom(i)
    do k = 1,3
      fnlkb(k,j,i) = 2*PI * fnlkb(k,j,i)
    enddo
    enddo
  enddo

  return
  end subroutine for_str_nl_kb
