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

!>  Sets up the special representation of the
!>  non local part of the hamiltonian for a given k-point
!>  Kleinman and Bylander pseudo-potential
!>  and its first and second derivatives with respect to coordinates.

subroutine proj_nl_kb_der_c16(rkpt, mtxd, isort, nanl, nder,             &
    ng, kgv,                                                             &
    nqnl, delqnl, vkb, nkb,                                              &
    ntype, natom, rat, adot,                                             &
    anlga, xnlkb, danlgadrk, d2anlgadrk2,                                &
    mxdtyp, mxdatm, mxdlqp, mxddim, mxdanl, mxdgve)

! written f90/complex  19 June 2012
! Modified 7 January 2014, style. JLM
! modified, dimensions vkb, real rylm, March 31, 2014. jlm
! Adapted Summer 2015. vkb normalization. JLM
! Modified documentation, January 2020. JLM
! openmp version 19 January 2020. JLM
! copyright INESC-MN/Jose Luis Martins

! version 4.99


  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  mxdtyp                          !<  array dimension of types of atoms
  integer, intent(in)                ::  mxdatm                          !<  array dimension of number of atoms of a given type
  integer, intent(in)                ::  mxddim                          !<  array dimension of plane-waves
  integer, intent(in)                ::  mxdlqp                          !<  array dimension for local potential
  integer, intent(in)                ::  mxdanl                          !<  array dimension of number of projectors
  integer, intent(in)                ::  mxdgve                          !<  array dimension of G-space vectors

  real(REAL64), intent(in)           ::  rkpt(3)                         !<  k-point reciprocal lattice coordinates
  integer, intent(in)                ::  mtxd                            !<  wavefunction dimension
  integer, intent(in)                ::  isort(mxddim)                   !<  G-vector corresponding to coefficient i of wavefunction

  integer, intent(in)                ::  nder                            !<  Indicates order of derivatives to be calculated, nder = 0,1,2

  integer, intent(in)                ::  ng                              !<  size of g-space
  integer, intent(in)                ::  kgv(3,mxdgve)                   !<  G-vectors in reciprocal lattice coordinates
  integer, intent(in)                ::  ntype                           !<  number of types of atoms
  integer, intent(in)                ::  natom(mxdtyp)                   !<  number of atoms of type i
  real(REAL64), intent(in)           ::  rat(3,mxdatm,mxdtyp)            !<  lattice coordinates of atom j of type i
  real(REAL64), intent(in)           ::  adot(3,3)                       !<  metric in real space
  integer, intent(in)                ::  nqnl(mxdtyp)                    !<  number of points for pseudo interpolation for atom k
  real(REAL64), intent(in)           ::  delqnl(mxdtyp)                  !<  step used in the pseudo interpolation for atom k
  real(REAL64), intent(in)           ::  vkb(-2:mxdlqp,0:3,-1:1,mxdtyp)  !<  (1/q**l) * KB nonlocal pseudo. for atom k, ang. mom. l. NOT normalized to vcell, hartree
  integer, intent(in)                ::  nkb(0:3,-1:1,mxdtyp)            !<  KB pseudo.  normalization for atom k, ang. mom. l

! output

  integer, intent(out)               ::  nanl                            !<  half of number of projectors without spin
  complex(REAL64), intent(out)       ::  anlga(mxddim,mxdanl)            !<  KB projectors without spin-orbit
  real(REAL64), intent(out)          ::  xnlkb(mxdanl)                   !<  KB normalization without spin-orbit
  complex(REAL64), intent(out)       ::  danlgadrk(mxddim,mxdanl,3)      !<  d anlga / d rkpt
  complex(REAL64), intent(out)       ::  d2anlgadrk2(mxddim,mxdanl,3,3)  !<  d^2 anlga / d rkpt^2

! allocatable local variables

  complex(REAL64), allocatable       ::  st(:,:)
  complex(REAL64), allocatable       ::  vqil(:,:)
  complex(REAL64), allocatable       ::  dvqildq(:,:)
  complex(REAL64), allocatable       ::  d2vqildq2(:,:)

  integer, allocatable               ::  k_ind(:)
  integer, allocatable               ::  kk_ind(:)
  integer, allocatable               ::  l_ind(:)
  integer, allocatable               ::  m_ind(:)


! local variables

  real(REAL64)    ::  avec(3,3), bvec(3,3)
  real(REAL64)    ::  vcell, bdot(3,3),fac
  real(REAL64)    ::  qi,qk(3), qcar(3), fi, xni, xi, vq0

  real(REAL64)    ::  rylm(0:3,-3:3)                 !   r**l sqrt(4 pi / 2l+1) "Y_lm(l,m)+-Y_lm(l,-m)"
  real(REAL64)    ::  drylmdqc(3,0:3,-3:3)           !   d rylm / d r
  real(REAL64)    ::  d2rylmdqc2(3,3,0:3,-3:3)       !   d^2 rylm / d rc^2

  complex(REAL64) ::  zylm(0:3,-3:3)                 !   unused
  complex(REAL64) ::  dzylmdqc(3,0:3,-3:3)           !   unused
  complex(REAL64) ::  d2zylmdqc2(3,3,0:3,-3:3)       !   unused

  real(REAL64)    ::  dqidrk(3)
  real(REAL64)    ::  d2qidrk2(3,3)
  real(REAL64)    ::  drylmdrk(3,0:3,-3:3)
  real(REAL64)    ::  d2rylmdrk2(3,3,0:3,-3:3)
  real(REAL64)    ::  dvq0dqi, d2vq0dqi2

  integer         ::  lmax, ifail
  logical         ::  lqi

! constants

  real(REAL64), parameter    :: ZERO = 0.0_REAL64, UM = 1.0_REAL64
  complex(REAL64), parameter :: C_I = cmplx(ZERO,UM,REAL64)
  real(REAL64), parameter    :: PI = 3.14159265358979323846_REAL64
  real(REAL64), parameter    :: EPS = 1.0E-8_REAL64

! counters

  integer         ::  ind
  integer         ::  j, k, kk, ik, l, n, m
  integer         ::  i, ni


! paranoid check, stops complaint about unused ng...

  if(ng > mxdgve) stop

  if(nder < 0 .OR. nder > 2) then
    write(6,'("  proj_nl_kb_der_c16:    cannot calculate",               &
           &  " derivatives up to ",i4)')  nder

    stop

  endif


  allocate(st(mxdatm,mxdtyp))
  allocate(vqil(0:3,mxdtyp))
  allocate(dvqildq(0:3,mxdtyp))
  allocate(d2vqildq2(0:3,mxdtyp))

  allocate(k_ind(nanl), kk_ind(nanl), l_ind(nanl), m_ind(nanl))

! cell stuff

  call adot_to_avec(adot, avec, bvec)
  call adot_to_bdot(adot, vcell, bdot)
  fac = UM / sqrt(vcell)

! checks dimensions and maximum angular momentum

  ind = 0
  lmax = 0
  do k=1,ntype
    do l = 0,3
      if(nkb(l,0,k) /= 0) then
        ind = ind + (2*l+1)*natom(k)
        lmax = max(lmax,l)
      endif
    enddo
  enddo

  nanl = ind
  if(nanl > mxdanl) then
    write(6,'("  proj_nl_kb:    increase mxdanl from ",i8,               &
           &  " to ",2i8)')  mxdanl, nanl, ng

    stop

  endif

  if(lmax > 3) then
    write(6,'("  proj_nl_kb:    lmax = ",i8,                            &
           &  " > 3 (max allowed in code)")')  lmax

    stop

  endif

! fills indexation array

  ind = 0
  do k = 1,ntype
    do l = 0,3
      if(nkb(l,0,k) /= 0) then

        do m = -l,l
          do kk = 1,natom(k)
            ind = ind + 1
            k_ind(ind) = k
            kk_ind(ind) = kk
            l_ind(ind) = l
            m_ind(ind) = m
          enddo
        enddo

      endif
    enddo    !  l=0,3
  enddo      !  k=1,ntype


!      starts loop over g-vectors

!$omp  parallel do default(private)                                      &
!$omp& shared(kgv, bvec, bdot, ntype, natom, mtxd, rkpt, isort, rat)     &
!$omp& shared(nder, lmax, delqnl, nqnl, vkb, fac)                        &
!$omp& shared(nanl, anlga, danlgadrk, d2anlgadrk2)                       &
!$omp& shared(k_ind, kk_ind, l_ind, m_ind)
  do i = 1,mtxd

    qk(1) = rkpt(1) + UM*(kgv(1,isort(i)))
    qk(2) = rkpt(2) + UM*(kgv(2,isort(i)))
    qk(3) = rkpt(3) + UM*(kgv(3,isort(i)))
    qcar(1) = bvec(1,1)*qk(1)+bvec(1,2)*qk(2)+bvec(1,3)*qk(3)
    qcar(2) = bvec(2,1)*qk(1)+bvec(2,2)*qk(2)+bvec(2,3)*qk(3)
    qcar(3) = bvec(3,1)*qk(1)+bvec(3,2)*qk(2)+bvec(3,3)*qk(3)
    qi = qcar(1)*qcar(1) + qcar(2)*qcar(2) + qcar(3)*qcar(3)
    qi = sqrt(qi)

    if(nder > 0) then
      if(qi > EPS) then

        lqi = .TRUE.

        do n=1,3
          dqidrk(n) = ( bdot(n,1)*qk(1) + bdot(n,2)*qk(2) +              &
                        bdot(n,3)*qk(3) ) / qi
        enddo

        do n=1,3
        do m=1,3
          d2qidrk2(n,m) = bdot(n,m) / qi - dqidrk(n)*dqidrk(m) / qi
        enddo
        enddo

      else

        lqi = .FALSE.

!       these quantities are ill defined for q -> 0

        dqidrk(1) = ZERO
        dqidrk(2) = ZERO
        dqidrk(3) = ZERO
        do n=1,3
        do m=1,3
          d2qidrk2(n,m) = ZERO
        enddo
        enddo

      endif
    endif

!   angular functions

    call ylm_rc(qcar, rylm, drylmdqc, d2rylmdqc2, .FALSE.,               &
        zylm, dzylmdqc, d2zylmdqc2, lmax, nder, ifail)

!   convert to lattice coordinates

    if(nder > 0) then
      do l = 0,lmax
      do m =-l,l
        do n = 1,3
          drylmdrk(n,l,m) = ZERO
          do k = 1,3
            drylmdrk(n,l,m) = drylmdrk(n,l,m) + drylmdqc(k,l,m) * bvec(k,n)
          enddo
        enddo
      enddo
      enddo
    endif

    if(nder > 1) then
      do l = 0,lmax
      do m =-l,l
        do n = 1,3
        do j = 1,3
          d2rylmdrk2(j,n,l,m) = ZERO
          do k=1,3
          do ik=1,3
            d2rylmdrk2(j,n,l,m) = d2rylmdrk2(j,n,l,m) +                   &
                  bvec(ik,j)*d2rylmdqc2(ik,k,l,m) * bvec(k,n)
          enddo
          enddo
        enddo
        enddo
      enddo
      enddo
    endif

!   structure phase factor

    do k = 1,ntype
      do kk = 1,natom(k)
        fi = kgv(1,isort(i))*rat(1,kk,k) +                               &
             kgv(2,isort(i))*rat(2,kk,k) +                               &
             kgv(3,isort(i))*rat(3,kk,k)
        fi = -2*PI*fi
        st(kk,k) = cmplx(cos(fi),sin(fi),REAL64)
      enddo
    enddo

!   i**l * v(q)

    do k = 1,ntype
      do l = 0,3

        xni = qi/delqnl(k)
        ni  = nint(xni)
        if(ni < 0) ni = 0

        vq0 = ZERO
        dvq0dqi = ZERO
        d2vq0dqi2 = ZERO

        if(ni < nqnl(k)-1) then

          xi  = xni - UM*ni
          vq0 = vkb(ni,l,0,k) * (UM+xi) * (UM-xi)                         &
               + (UM/2) * (vkb(ni+1,l,0,k)*(UM+xi)                        &
                          -vkb(ni-1,l,0,k)*(UM-xi)) * xi
          vq0 = fac * vq0

          if(nder > 0) then
            dvq0dqi = (-2*xi*vkb(ni,l,0,k)                               &
               + (UM/2) * (vkb(ni+1,l,0,k)*(UM+2*xi)                     &
                          -vkb(ni-1,l,0,k)*(UM-2*xi))) / delqnl(k)
            dvq0dqi = fac * dvq0dqi
          endif

          if(nder > 1) then
            d2vq0dqi2 = (-2*vkb(ni,l,0,k) +  vkb(ni+1,l,0,k)             &
                        +vkb(ni-1,l,0,k)) / (delqnl(k)*delqnl(k))
            d2vq0dqi2 = fac * d2vq0dqi2
          endif

        endif

        vqil(l,k) = vq0*C_I**l
        dvqildq(l,k) = dvq0dqi*C_I**l
        d2vqildq2(l,k) = d2vq0dqi2*C_I**l

      enddo

    enddo

!   starts loop over second index

    do ind = 1,nanl

      kk = kk_ind(ind)
      k = k_ind(ind)
      l = l_ind(ind)
      m = m_ind(ind)

      anlga(i,ind) = st(kk,k)*rylm(l,m)*vqil(l,k)

      if(nder > 0) then
        do n = 1,3
          danlgadrk(i,ind,n) = st(kk,k)*                                 &
             ( drylmdrk(n,l,m)*vqil(l,k)                                 &
               + rylm(l,m)*dvqildq(l,k)*dqidrk(n) )
        enddo
      endif

      if(nder > 1) then

!       second derivative of the projector. d2qidrk2 -> 1/q,
!       for small q, d2qidrk2 -> 1/q, dqidrk undefined, rylm -> q**l

        if(lqi) then                                                    !      if(qi > EPS). avoids float imprecision
          do n=1,3
          do j=1,3
            d2anlgadrk2(i,ind,j,n) = st(kk,k)*                           &
              (  d2rylmdrk2(j,n,l,m)*vqil(l,k)                           &
               + drylmdrk(n,l,m)*dvqildq(l,k)*dqidrk(j)                  &
               + drylmdrk(j,l,m)*dvqildq(l,k)*dqidrk(n)                  &
               + rylm(l,m)*d2vqildq2(l,k)*dqidrk(n)*dqidrk(j)            &
               + rylm(l,m)*dvqildq(l,k)*d2qidrk2(j,n)   )
          enddo
          enddo
        else
          do n=1,3
          do j=1,3
            d2anlgadrk2(i,ind,j,n) = st(kk,k)*                           &
                ( d2rylmdrk2(j,n,l,m)*vqil(l,k) +                        &
                  rylm(l,m)*d2vqildq2(l,k)*bdot(j,n) )
          enddo
          enddo
        endif

      endif
    enddo

  enddo
!$omp end parallel do

! sign of projector normalization

  do ind = 1,nanl
    xnlkb(ind) = UM*nkb(l_ind(ind),0,k_ind(ind))
  enddo

  deallocate(st)
  deallocate(vqil, dvqildq, d2vqildq2)

  deallocate(k_ind, kk_ind, l_ind, m_ind)

  return
end subroutine proj_nl_kb_der_c16
