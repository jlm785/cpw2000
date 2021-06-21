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
!>  and its derivatives

!>  spin-orbit perturbation only

subroutine proj_nl_kb_so_pert_der_c16(rkpt, mtxd, isort,                 &
    nanl, nanlso, nder,                                                  &
    ng, kgv,                                                             &
    nqnl, delqnl, vkb, nkb,                                              &
    ntype, natom, rat, adot,                                             &
    anlspin, anlsop, anlsom, xnlkbspin, xnlkbso,                         &
    danlspindrk, danlsopdrk, danlsomdrk,                                 &
    d2anlspindrk2, d2anlsopdrk2, d2anlsomdrk2,                           &
    mxdtyp, mxdatm, mxdlqp, mxddim, mxdanl, mxdaso, mxdgve)


! written 28 January 2020 from the non-perturbation version. JLM
! openmp version 19 January 2020. JLM

! copyright INESC-MN/Jose Luis Martins

! version 4.95


  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  mxdtyp                          !<  array dimension of types of atoms
  integer, intent(in)                ::  mxdatm                          !<  array dimension of number of atoms of a given type
  integer, intent(in)                ::  mxddim                          !<  array dimension of plane-waves
  integer, intent(in)                ::  mxdlqp                          !<  array dimension for local potential
  integer, intent(in)                ::  mxdanl                          !<  array dimension of number of projectors
  integer, intent(in)                ::  mxdaso                          !<  array dimension of number of projectors with spin-orbit
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
  integer, intent(out)               ::  nanlso                          !<  number of projectors witn spin

  complex(REAL64), intent(out)       ::  anlspin(mxddim,mxdanl)          !<  KB projectors without spin-orbit
  complex(REAL64), intent(out)       ::  anlsop(mxddim,mxdaso)           !<  KB projectors with spin-orbit m_s = + 1/2
  complex(REAL64), intent(out)       ::  anlsom(mxddim,mxdaso)           !<  KB projectors with spin-orbit m_s = - 1/2

  real(REAL64), intent(out)          ::  xnlkbspin(mxdanl)               !<  KB normalization without spin-orbit
  real(REAL64), intent(out)          ::  xnlkbso(mxdaso)                 !<  KB normalization with spin-orbit

  complex(REAL64), intent(out)     ::  danlspindrk(mxddim,mxdanl,3)      !<  d anlspin / d rkpt
  complex(REAL64), intent(out)     ::  d2anlspindrk2(mxddim,mxdanl,3,3)  !<  d^2 anlspin / d rkpt^2
  complex(REAL64), intent(out)     ::  danlsopdrk(mxddim,mxdaso,3)       !<  d anlsop / d rkpt
  complex(REAL64), intent(out)     ::  d2anlsopdrk2(mxddim,mxdaso,3,3)   !<  d^2 anlsop / d rkpt^2
  complex(REAL64), intent(out)     ::  danlsomdrk(mxddim,mxdaso,3)       !<  d anlsom / d rkpt
  complex(REAL64), intent(out)     ::  d2anlsomdrk2(mxddim,mxdaso,3,3)   !<  d^2 anlsom / d rkpt^2


! allocatable local variables

  complex(REAL64), allocatable  ::  st(:,:)
  complex(REAL64), allocatable  ::  vqil(:,:,:)
  complex(REAL64), allocatable  ::  dvqildq(:,:,:)
  complex(REAL64), allocatable  ::  d2vqildq2(:,:,:)

  integer, allocatable          ::  k_ind(:)                             !  index of atom type
  integer, allocatable          ::  kk_ind(:)                            !  index of atom of type k
  integer, allocatable          ::  l_ind(:)                             !  quantum number l
  integer, allocatable          ::  m_ind(:)                             !  quantum number m

  integer, allocatable          ::  k_ind_so(:)                          !  index of atom type
  integer, allocatable          ::  kk_ind_so(:)                         !  index of atom of type k
  integer, allocatable          ::  l_ind_so(:)                          !  quantum number l
  integer, allocatable          ::  j_ind_so(:)                          !  2* quantum number j
  integer, allocatable          ::  mj_ind_so(:)                         !  2* quantum number m_j

! local variables

  real(REAL64)    ::  avec(3,3), bvec(3,3)
  real(REAL64)    ::  vcell, bdot(3,3),fac
  real(REAL64)    ::  qi,qk(3), qcar(3), fi, xni, xi

  real(REAL64)    ::  rylm(0:3,-3:3)                                     !   r**l sqrt(4 pi / 2l+1) "Y_lm(l,m)+-Y_lm(l,-m)"
  real(REAL64)    ::  drylmdqc(3,0:3,-3:3)                               !   d rylm / d r
  real(REAL64)    ::  d2rylmdqc2(3,3,0:3,-3:3)                           !   d^2 rylm / d rc^2

  complex(REAL64) ::  zylm(0:3,-3:3)                                     !  r**l sqrt(4 pi / 2l+1) Y_lm(theta,phi)
  complex(REAL64) ::  dzylmdqc(3,0:3,-3:3)                               !  d zylm / d rc
  complex(REAL64) ::  d2zylmdqc2(3,3,0:3,-3:3)                           !  d^2 zylm / d rc^2

  real(REAL64)    ::  dqidrk(3)
  real(REAL64)    ::  d2qidrk2(3,3)

  real(REAL64)    ::  vq(-1:1)
  real(REAL64)    ::  dvqdqi(-1:1), d2vqdqi2(-1:1)

!  real(REAL64)    ::  drylmdrk(3,0:3,-3:3)
!  real(REAL64)    ::  d2rylmdrk2(3,3,0:3,-3:3)
  complex(REAL64) ::  dzylmdrk(3,0:3,-3:3)
  complex(REAL64) ::  d2zylmdrk2(3,3,0:3,-3:3)

  real(REAL64)    ::  cg(0:3,-7:7)                                       !  Clebsch-Gordan coefficient (0:lmax,-2*lmax-1:2*lmax+1)

  integer         ::  lmax, ifail
  logical         ::  lqi

! constants

  real(REAL64), parameter    :: ZERO = 0.0_REAL64, UM = 1.0_REAL64
  complex(REAL64), parameter :: C_I = cmplx(ZERO,UM,REAL64)
  complex(REAL64), parameter :: C_ZERO = cmplx(ZERO,ZERO,REAL64)
  real(REAL64), parameter    :: PI = 3.14159265358979323846_REAL64
  real(REAL64), parameter    :: EPS = 1.0E-8_REAL64

! counters

  integer         ::  ind
  integer         ::  k, kk, l, m, j, mj, ms
  integer         ::  i, ni, n, ik, n1,n2


! paranoid check, stops complaint about unused ng...

  if(ng > mxdgve) stop

  if(nder < 0 .OR. nder > 2) then
    write(6,'("  proj_nl_kb_der_c16:    cannot calculate",               &
           &  " derivatives up to ",i4)')  nder

    stop

  endif


  allocate(st(mxdatm,mxdtyp))
  allocate(vqil(0:3,-3:3,mxdtyp))
  allocate(dvqildq(0:3,-3:3,mxdtyp))
  allocate(d2vqildq2(0:3,-3:3,mxdtyp))

  allocate(k_ind(nanl), kk_ind(nanl), l_ind(nanl), m_ind(nanl))
  allocate(k_ind_so(nanlso), kk_ind_so(nanlso), l_ind_so(nanlso))
  allocate(j_ind_so(nanlso), mj_ind_so(nanlso))

! cell stuff

  call adot_to_avec(adot, avec, bvec)
  call adot_to_bdot(adot, vcell, bdot)
  fac = UM / sqrt(vcell)

! checks dimensions and maximum angular momentum

  ind = 0
  lmax = 0
  do k = 1,ntype
    do l = 0,3
      if(nkb(l,0,k) /= 0) then
        ind = ind + (2*l+1)*natom(k)
        lmax = max(lmax,l)
      endif
    enddo
  enddo
  nanl = ind
  if(nanl > mxdanl) then
    write(6,'("  proj_nl_kb_so_der:    increase mxdanl from ",i8,        &
           &  " to ",2i8)')  mxdanl, nanl, ng

    stop

  endif

  if(lmax > 3) then
    write(6,'("  proj_nl_kb_so_der:    lmax = ",i8,                      &
           &  " > 3 (max allowed in code)")')  lmax

    stop

  endif

  ind = 0
  do k=1,ntype
    do l = 1,3
      if(nkb(l,-1,k) /= 0) then
        ind = ind + (2*l)*natom(k)
        lmax = max(lmax,l)
      endif
    enddo
    do l = 0,3
      if(nkb(l,1,k) /= 0) then
        ind = ind + (2*l+2)*natom(k)
        lmax = max(lmax,l)
      endif
    enddo
  enddo

  nanlso = ind
  if(nanlso > mxdaso) then
    write(6,'("  proj_nl_kb_so_der:    increase mxdaso from ",i8,        &
           &  " to ",i8)')  mxdaso,nanlso

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

  ind = 0
  do k = 1,ntype
    do l = 0,3
      if(nkb(l,1,k) /= 0) then
        do mj = -2*l-1,2*l+1,2
          do kk = 1,natom(k)
            ind = ind + 1
            k_ind_so(ind) = k
            kk_ind_so(ind) = kk
            l_ind_so(ind) = l
            j_ind_so(ind) = 2*l+1
            mj_ind_so(ind) = mj
          enddo
        enddo
      endif
      if(nkb(l,-1,k) /= 0) then
        if(l > 0) then
          do mj = -2*l+1,2*l-1,2
            do kk = 1,natom(k)
              ind = ind + 1
              k_ind_so(ind) = k
              kk_ind_so(ind) = kk
              l_ind_so(ind) = l
              j_ind_so(ind) = 2*l-1
              mj_ind_so(ind) = mj
            enddo
          enddo
        endif
      endif
    enddo
  enddo

! Clebsch-Gordan coefficients

  do l = 0,lmax
  do j = -2*lmax-1,2*lmax+1
    cg(l,j) = ZERO
  enddo
  enddo

  do l = 0,lmax
  do j = -2*l-1,2*l+1
    cg(l,j) = sqrt( ((2*l+1 + j)*UM) / ((4*l+2)*UM) )
  enddo
  enddo

! starts loop over g-vectors

!$omp  parallel do default(private)                                      &
!$omp& shared(kgv, bvec, bdot, ntype, natom, mtxd, rkpt, isort, rat)     &
!$omp& shared(nder, lmax, delqnl, nqnl, vkb, fac)                        &
!$omp& shared(nanl, anlspin, danlspindrk, d2anlspindrk2)                 &
!$omp& shared(nanlso, anlsop, danlsopdrk, d2anlsopdrk2)                  &
!$omp& shared(anlsom, danlsomdrk, d2anlsomdrk2)                          &
!$omp& shared(k_ind, kk_ind, l_ind, m_ind)                               &
!$omp& shared(k_ind_so, kk_ind_so, l_ind_so, j_ind_so, mj_ind_so)        &
!$omp& shared(cg)
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

    call ylm_rc(qcar, rylm, drylmdqc, d2rylmdqc2, .TRUE.,                &
        zylm, dzylmdqc, d2zylmdqc2, lmax, nder, ifail)

!   convert to lattice coordinates

    if(nder > 0) then
      do l = 0,lmax
      do m =-l,l
        do n = 1,3
          dzylmdrk(n,l,m) = C_ZERO
          do k = 1,3
            dzylmdrk(n,l,m) = dzylmdrk(n,l,m) + dzylmdqc(k,l,m) * bvec(k,n)
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
          d2zylmdrk2(j,n,l,m) = C_ZERO
          do k=1,3
          do ik=1,3
            d2zylmdrk2(j,n,l,m) = d2zylmdrk2(j,n,l,m) +                  &
                  bvec(ik,j)*d2zylmdqc2(ik,k,l,m) * bvec(k,n)
          enddo
          enddo
        enddo
        enddo
      enddo
      enddo
    endif


!   structure phase factor

    do k=1,ntype
      do kk=1,natom(k)
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

        do m=-1,1
          vq(m) = ZERO
          if(ni < nqnl(k)-1) then
            xi  = xni - UM*ni
            vq(m) = vkb(ni,l,m,k) * (UM+xi) * (UM-xi)                    &
               + (UM/2) * (vkb(ni+1,l,m,k)*(UM+xi)                       &
                          -vkb(ni-1,l,m,k)*(UM-xi)) * xi
            vq(m) = fac * vq(m)

            if(nder > 0) then
              dvqdqi(m) = (-2*xi*vkb(ni,l,m,k)                           &
               + (UM/2) * (vkb(ni+1,l,m,k)*(UM+2*xi)                     &
                          -vkb(ni-1,l,m,k)*(UM-2*xi))) / delqnl(k)
              dvqdqi(m) = fac * dvqdqi(m)
            endif

            if(nder > 1) then
              d2vqdqi2(m) = (-2*vkb(ni,l,m,k) +  vkb(ni+1,l,m,k)         &
                        +vkb(ni-1,l,m,k)) / (delqnl(k)*delqnl(k))
              d2vqdqi2(m) = fac * d2vqdqi2(m)
            endif

          endif

          vqil(l,m,k) = vq(m)*C_I**l
          dvqildq(l,m,k) = dvqdqi(m)*C_I**l
          d2vqildq2(l,m,k) = d2vqdqi2(m)*C_I**l

        enddo

      enddo

    enddo

!   loop over second index, spin-orbit
!   mj, ms are twice hte half-integers

    do ind = 1,nanlso

      k = k_ind_so(ind)
      kk = kk_ind_so(ind)
      l = l_ind_so(ind)
      j = j_ind_so(ind)
      mj = mj_ind_so(ind)
      ms = j - 2*l

      if( mj == -2*l-1 ) then

        anlsop(i,ind) = C_ZERO
        anlsom(i,ind) = st(kk,k)*zylm(l,-l)*vqil(l,ms,k)

      elseif( mj == 2*l+1 ) then

        anlsop(i,ind) = st(kk,k)*zylm(l, l)*vqil(l,ms,k)
        anlsom(i,ind) = ZERO


      elseif( l > 0 ) then

        anlsop(i,ind) =    st(kk,k)*cg(l, ms*mj)*zylm(l,(mj-1)/2)*vqil(l,ms,k)
        anlsom(i,ind) = ms*st(kk,k)*cg(l,-ms*mj)*zylm(l,(mj+1)/2)*vqil(l,ms,k)

      endif

!     Calculates first derivatives

      if(nder > 0) then

        if( mj == -2*l-1 ) then

          do n = 1,3

            danlsopdrk(i,ind,n) = C_ZERO

            danlsomdrk(i,ind,n) = st(kk,k)*                              &
                  (  dzylmdrk(n,l,-l)*vqil(l,ms,k) +                     &
                     zylm(l,-l)*dvqildq(l,ms,k)*dqidrk(n) )

          enddo

        elseif( mj == 2*l+1 ) then

          do n = 1,3

            danlsopdrk(i,ind,n) = st(kk,k)*                              &
                  (  dzylmdrk(n,l, l)*vqil(l,ms,k) +                     &
                     zylm(l, l)*dvqildq(l,ms,k)*dqidrk(n) )

            danlsomdrk(i,ind,n) = C_ZERO

          enddo

        elseif( l > 0 ) then

          do n = 1,3
            danlsopdrk(i,ind,n) =    st(kk,k)*cg(l, ms*mj)*              &
                  (  dzylmdrk(n,l,(mj-1)/2)*vqil(l,ms,k) +               &
                     zylm(l,(mj-1)/2)*dvqildq(l,ms,k)*dqidrk(n) )

            danlsomdrk(i,ind,n) = ms*st(kk,k)*cg(l,-ms*mj)*              &
                  (  dzylmdrk(n,l,(mj+1)/2)*vqil(l,ms,k) +               &
                     zylm(l,(mj+1)/2)*dvqildq(l,ms,k)*dqidrk(n) )
          enddo
        endif

      endif

!     Calculates second derivatives

      if(nder > 1) then

!       second derivative of the projector. d2qidrk2 -> 1/q,
!       for small q, d2qidrk2 -> 1/q, dqidrk undefined, rylm -> q**l
!       so (qi < EPS) is treated separately.

        if(lqi) then                                                     !  qi > EPS

          if( mj == -2*l-1 ) then

            do n1 = 1,3
            do n2 = 1,3

              d2anlsopdrk2(i,ind,n2,n1) = C_ZERO

              d2anlsomdrk2(i,ind,n2,n1) = st(kk,k)*                      &
                (  d2zylmdrk2(n2,n1,l,-l)*vqil(l,ms,k)                   &
                 + dzylmdrk(n1,l,-l)*dvqildq(l,ms,k)*dqidrk(n2)          &
                 + dzylmdrk(n2,l,-l)*dvqildq(l,ms,k)*dqidrk(n1)          &
                 + zylm(l,-l)*d2vqildq2(l,ms,k)*dqidrk(n1)*dqidrk(n2)    &
                 + zylm(l,-l)*dvqildq(l,ms,k)*d2qidrk2(n2,n1)   )

            enddo
            enddo

          elseif( mj == 2*l+1 ) then

            do n1 = 1,3
            do n2 = 1,3

              d2anlsopdrk2(i,ind,n2,n1) = st(kk,k)*                      &
                (  d2zylmdrk2(n2,n1,l,l)*vqil(l,ms,k)                    &
                 + dzylmdrk(n1,l,l)*dvqildq(l,ms,k)*dqidrk(n2)           &
                 + dzylmdrk(n2,l,l)*dvqildq(l,ms,k)*dqidrk(n1)           &
                 + zylm(l,l)*d2vqildq2(l,ms,k)*dqidrk(n1)*dqidrk(n2)     &
                 + zylm(l,l)*dvqildq(l,ms,k)*d2qidrk2(n2,n1)   )

              d2anlsomdrk2(i,ind,n2,n1) = C_ZERO

            enddo
            enddo

          elseif( l > 0 ) then

            do n1 = 1,3
            do n2 = 1,3

              d2anlsopdrk2(i,ind,n2,n1) =    st(kk,k)*cg(l, ms*mj)*      &
                (  d2zylmdrk2(n2,n1,l,(mj-1)/2)*vqil(l,ms,k)             &
                 + dzylmdrk(n1,l,(mj-1)/2)*dvqildq(l,ms,k)*dqidrk(n2)    &
                 + dzylmdrk(n2,l,(mj-1)/2)*dvqildq(l,ms,k)*dqidrk(n1)    &
                 + zylm(l,(mj-1)/2)*d2vqildq2(l,ms,k)*                   &
                             dqidrk(n1)*dqidrk(n2)                       &
                 + zylm(l,(mj-1)/2)*dvqildq(l,ms,k)*d2qidrk2(n2,n1)   )

              d2anlsomdrk2(i,ind,n2,n1) = ms*st(kk,k)*cg(l,-ms*mj)*      &
                (  d2zylmdrk2(n2,n1,l,(mj+1)/2)*vqil(l,ms,k)             &
                 + dzylmdrk(n1,l,(mj+1)/2)*dvqildq(l,ms,k)*dqidrk(n2)    &
                 + dzylmdrk(n2,l,(mj+1)/2)*dvqildq(l,ms,k)*dqidrk(n1)    &
                 + zylm(l,(mj+1)/2)*d2vqildq2(l,ms,k)*                   &
                             dqidrk(n1)*dqidrk(n2)                       &
                 + zylm(l,(mj+1)/2)*dvqildq(l,ms,k)*d2qidrk2(n2,n1)   )

            enddo
            enddo

          endif

        else                                                             !  qi < EPS

          if( mj == -2*l-1 ) then

            do n1 = 1,3
            do n2 = 1,3

              d2anlsopdrk2(i,ind,n2,n1) = C_ZERO

              d2anlsomdrk2(i,ind,n2,n1) = st(kk,k)*                      &
                  ( d2zylmdrk2(n2,n1,l,-l)*vqil(l,ms,k) +                &
                    zylm(l,-l)*d2vqildq2(l,ms,k)*bdot(n2,n1) )

            enddo
            enddo

          elseif( mj == 2*l+1 ) then

            do n1 = 1,3
            do n2 = 1,3

              d2anlsopdrk2(i,ind,n2,n1) = st(kk,k)*                      &
                  ( d2zylmdrk2(n2,n1,l, l)*vqil(l,ms,k) +                &
                    zylm(l, l)*d2vqildq2(l,ms,k)*bdot(n2,n1) )

              d2anlsomdrk2(i,ind,n2,n1) = C_ZERO

            enddo
            enddo

          elseif( l > 0 ) then

            do n1 = 1,3
            do n2 = 1,3

              d2anlsopdrk2(i,ind,n2,n1) =    st(kk,k)*cg(l, ms*mj)*      &
                ( d2zylmdrk2(n2,n1,l,(mj-1)/2)*vqil(l,ms,k) +            &
                  zylm(l,(mj-1)/2)*d2vqildq2(l,ms,k)*bdot(n2,n1) )

              d2anlsomdrk2(i,ind,n2,n1) = ms*st(kk,k)*cg(l,-ms*mj)*      &
                ( d2zylmdrk2(n2,n1,l,(mj+1)/2)*vqil(l,ms,k) +            &
                  zylm(l,(mj+1)/2)*d2vqildq2(l,ms,k)*bdot(n2,n1) )

            enddo
            enddo
          endif

        endif

      endif

    enddo

!   loop over second index, without spin-orbit

    do ind = 1,nanl

      k = k_ind(ind)
      kk = kk_ind(ind)
      l = l_ind(ind)
      m = m_ind(ind)

      anlspin(i,ind) = st(kk,k)*zylm(l,m)*vqil(l,0,k)

      if(nder > 0) then
        do n = 1,3
          danlspindrk(i,ind,n) = st(kk,k)*                               &
             ( dzylmdrk(n,l,m)*vqil(l,0,k)                               &
               + zylm(l,m)*dvqildq(l,0,k)*dqidrk(n) )
        enddo
      endif

      if(nder > 1) then

!       second derivative of the projector. d2qidrk2 -> 1/q,
!       for small q, d2qidrk2 -> 1/q, dqidrk undefined, rylm -> q**l
!       if(qi > EPS). avoids float imprecision

        if(lqi) then
          do n1 = 1,3
          do n2 = 1,3
            d2anlspindrk2(i,ind,n2,n1) = st(kk,k)*                       &
              (  d2zylmdrk2(n2,n1,l,m)*vqil(l,0,k)                       &
               + dzylmdrk(n1,l,m)*dvqildq(l,0,k)*dqidrk(n2)              &
               + dzylmdrk(n2,l,m)*dvqildq(l,0,k)*dqidrk(n1)              &
               + zylm(l,m)*d2vqildq2(l,0,k)*dqidrk(n1)*dqidrk(n2)        &
               + zylm(l,m)*dvqildq(l,0,k)*d2qidrk2(n2,n1)   )
          enddo
          enddo
        else
          do n1 = 1,3
          do n2 = 1,3
            d2anlspindrk2(i,ind,n2,n1) = st(kk,k)*                       &
                ( d2zylmdrk2(n2,n1,l,m)*vqil(l,0,k) +                    &
                  zylm(l,m)*d2vqildq2(l,0,k)*bdot(n2,n1) )
          enddo
          enddo
        endif

      endif

    enddo


  enddo        !  i=1,mtxd
!$omp end parallel do



  do ind = 1,nanl
    xnlkbspin(ind) = UM*nkb(l_ind(ind),0,k_ind(ind))
  enddo

  do ind = 1,nanlso
    l = l_ind_so(ind)
    j = j_ind_so(ind)
    ms = j - 2*l
    k = k_ind_so(ind)
    xnlkbso(ind) = UM*nkb(l,ms,k)
  enddo

  deallocate(st)
  deallocate(vqil)

  deallocate(dvqildq)
  deallocate(d2vqildq2)

  deallocate(k_ind, kk_ind, l_ind, m_ind)
  deallocate(k_ind_so, kk_ind_so, l_ind_so)
  deallocate(j_ind_so, mj_ind_so)

  return
end subroutine proj_nl_kb_so_pert_der_c16
