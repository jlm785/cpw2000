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

!>  Sets up the special representation of the projectors of the
!>  non local part of the hamiltonian for a given k-point
!>  Kleinman and Bylander pseudo-potential.
!>  It includes both non-relativistic and relativistic
!>  pseudopotentials with the same representation.
!>
!>  For spin-orbit perturbation use proj_nl_kb_so_pert_c16
!>  that generates the projectors in both relativistic and non-relativistic
!>
!>
!>  \author       Jose Luis Martins
!>  \version      5.0.9
!>  \date         3 December 2023
!>  \copyright    GNU Public License v2

subroutine proj_nl_kb_so_c16(rkpt, mtxd, isort,                          &
    nanlsp,                                                              &
    ng, kgv,                                                             &
    nqnl, delqnl, vkb, nkb,                                              &
    ntype, natom, rat, adot,                                             &
    anlsp, xnlkbsp,                                                      &
    mxdtyp, mxdatm, mxdlqp, mxddim, mxdasp, mxdgve)


! Written 3 December 2023 from previous code. JLM
! It is based on the perturbative subroutine and the old-subroutine.


  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)
  integer, parameter          :: LMAX = 3                                !  hard coded max ang. mom.

! input

  integer, intent(in)                ::  mxdtyp                          !<  array dimension of types of atoms
  integer, intent(in)                ::  mxdatm                          !<  array dimension of number of atoms of a given type
  integer, intent(in)                ::  mxddim                          !<  array dimension of plane-waves
  integer, intent(in)                ::  mxdlqp                          !<  array dimension for local potential
  integer, intent(in)                ::  mxdasp                          !<  array dimension of number of projectors with or without spin-orbit
  integer, intent(in)                ::  mxdgve                          !<  array dimension of G-space vectors

  real(REAL64), intent(in)           ::  rkpt(3)                         !<  k-point reciprocal lattice coordinates
  integer, intent(in)                ::  mtxd                            !<  wavefunction dimension
  integer, intent(in)                ::  isort(mxddim)                   !<  G-vector corresponding to coefficient i of wavefunction

  integer, intent(in)                ::  ng                              !<  size of g-space
  integer, intent(in)                ::  kgv(3,mxdgve)                   !<  G-vectors in reciprocal lattice coordinates
  integer, intent(in)                ::  ntype                           !<  number of types of atoms
  integer, intent(in)                ::  natom(mxdtyp)                   !<  number of atoms of type i
  real(REAL64), intent(in)           ::  rat(3,mxdatm,mxdtyp)            !<  lattice coordinates of atom j of type i
  real(REAL64), intent(in)           ::  adot(3,3)                       !<  metric in real space
  integer, intent(in)                ::  nqnl(mxdtyp)                    !<  number of points for pseudo interpolation for atom k
  real(REAL64), intent(in)           ::  delqnl(mxdtyp)                  !<  step used in the pseudo interpolation for atom k
  real(REAL64), intent(in)        ::  vkb(-2:mxdlqp,0:LMAX,-1:1,mxdtyp)  !<  (1/q**l) * KB nonlocal pseudo. for atom k, ang. mom. l. NOT normalized to vcell, hartree
  integer, intent(in)                ::  nkb(0:LMAX,-1:1,mxdtyp)         !<  KB pseudo.  normalization for atom k, ang. mom. l

! output

  integer, intent(out)               ::  nanlsp                          !<  number of Kleinman-Bylander projectors with or without spin-orbit depending on pseudopotential

  complex(REAL64), intent(out)       ::  anlsp(2*mxddim,mxdasp)          !<  KB projectors in spin representation

  real(REAL64), intent(out)          ::  xnlkbsp(mxdasp)                 !<  KB normalization without spin-orbit

! allocatable local variables

  complex(REAL64), allocatable  ::  st(:,:)
  complex(REAL64), allocatable  ::  vqil(:,:,:)

  integer, allocatable          ::  k_ind(:)                             !  index of atom type
  integer, allocatable          ::  kk_ind(:)                            !  index of atom of type k
  integer, allocatable          ::  l_ind(:)                             !  quantum number l
  integer, allocatable          ::  m_ind(:)                             !  quantum number m
  integer, allocatable          ::  ms_ind(:)                            !  2*quantum number ms
  integer, allocatable          ::  j_ind(:)                             !  2* quantum number j
  integer, allocatable          ::  mj_ind(:)                            !  2* quantum number m_j

  logical, allocatable          ::  lso(:)                               !  indicates if the pseudopotential has relativistic components

! local variables

  real(REAL64)    ::  avec(3,3),bvec(3,3)
  real(REAL64)    ::  vcell, bdot(3,3),fac
  real(REAL64)    ::  qi,qk(3),qcar(3),fi,xni,xi

  real(REAL64)    ::  rylm(0:LMAX,-LMAX:LMAX)                            !   r**l sqrt(4 pi / 2l+1) "Y_lm(l,m)+-Y_lm(l,-m)"
  real(REAL64)    ::  drylmdqc(3,0:LMAX,-LMAX:LMAX)                      !   unused
  real(REAL64)    ::  d2rylmdqc2(3,3,0:LMAX,-LMAX:LMAX)                  !   unused

  complex(REAL64) ::  zylm(0:LMAX,-LMAX:LMAX)                            !  r**l sqrt(4 pi / 2l+1) Y_lm(theta,phi)
  complex(REAL64) ::  dzylmdqc(3,0:LMAX,-LMAX:LMAX)                      !  unused
  complex(REAL64) ::  d2zylmdqc2(3,3,0:LMAX,-LMAX:LMAX)                  !  unused


  real(REAL64)    ::  vq(-1:1)

  real(REAL64)    ::  cg(0:LMAX,-(2*LMAX+1):2*LMAX+1)                    !  Clebsch-Gordan coefficient (0:lmax,-2*lmax-1:2*lmax+1)

  integer         ::  lmx, ifail

! constants

  real(REAL64), parameter    :: ZERO = 0.0_REAL64, UM = 1.0_REAL64
  complex(REAL64), parameter :: C_I = cmplx(ZERO,UM,REAL64)
  complex(REAL64), parameter :: C_ZERO = cmplx(ZERO,ZERO,REAL64)
  real(REAL64), parameter    :: PI = 3.14159265358979323846_REAL64

! counters

  integer         ::  ind
  integer         ::  k, kk, l, m, j, mj, ms
  integer         ::  i, ni, ic


! paranoid check, stops complaint about unused ng...

  if(ng > mxdgve) stop

  allocate(st(mxdatm,mxdtyp))
  allocate(vqil(0:LMAX,-LMAX:LMAX,mxdtyp))

  allocate(lso(mxdtyp))


! cell stuff

  call adot_to_avec(adot, avec, bvec)
  call adot_to_bdot(adot, vcell, bdot)
  fac = UM / sqrt(vcell)

! checks dimensions and maximum angular momentum

  nanlsp = 0
  lmx = 0
  do k = 1,ntype

    ic = 0
    if(nkb(0, 1,k) /= 0) ic = ic+1
    do l = 1,LMAX
      if(nkb(l,-1,k) /= 0) ic = ic+1
      if(nkb(l, 1,k) /= 0) ic = ic+1
    enddo

    if(ic == 0) then

      lso(k) = .FALSE.

      do l = 0,LMAX
        if(nkb(l,0,k) /= 0) then
          nanlsp = nanlsp + 2*(2*l+1)*natom(k)
          lmx = max(lmx,l)
        endif
      enddo

    else

      lso(k) = .TRUE.
      if(nkb(0, 1,k) /= 0) nanlsp = nanlsp + 2*natom(k)

      do l = 1,LMAX
        if(nkb(l,-1,k) /= 0) then
          nanlsp = nanlsp + (2*l)*natom(k)
          lmx = max(lmx,l)
        endif
        if(nkb(l, 1,k) /= 0) then
          nanlsp = nanlsp + (2*l+2)*natom(k)
          lmx = max(lmx,l)
        endif
      enddo

    endif
  enddo

  if(nanlsp > mxdasp) then
    write(6,'("  proj_nl_kb_so_c16(1):    increase mxdasp from ",i8,     &
           &  " to ",2i8)')  mxdasp, nanlsp, ng

    stop

  endif

  if(lmx > LMAX) then
    write(6,'("  proj_nl_kb_so_c16:    lmx = ",i8,                       &
           &  " > ",i3," (max allowed in code)")')  lmx, LMAX

    stop

  endif

! fills indexation array

  allocate(k_ind(nanlsp), kk_ind(nanlsp), l_ind(nanlsp))
  allocate(m_ind(nanlsp), ms_ind(nanlsp))
  allocate(j_ind(nanlsp), mj_ind(nanlsp))

  k_ind(:) = 0
  kk_ind(:) = 0
  l_ind(:) = 0
  m_ind(:) = 0
  j_ind(:) = 0
  mj_ind(:) = 0

  ind = 0
  do k = 1,ntype
    do l = 0,lmx
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
    enddo    !  l=0,lmx
  enddo      !  k=1,ntype


  ind = 0
  do k = 1,ntype

    if(lso(k)) then

      do l = 0,lmx

        if(nkb(l,1,k) /= 0) then
          do mj = -2*l-1,2*l+1,2
            do kk = 1,natom(k)
              ind = ind + 1
              k_ind(ind) = k
              kk_ind(ind) = kk
              l_ind(ind) = l
              j_ind(ind) = 2*l+1
              mj_ind(ind) = mj
            enddo
          enddo
        endif

        if(nkb(l,-1,k) /= 0) then
          if(l > 0) then
            do mj = -2*l+1,2*l-1,2
              do kk = 1,natom(k)
                ind = ind + 1
                k_ind(ind) = k
                kk_ind(ind) = kk
                l_ind(ind) = l
                j_ind(ind) = 2*l-1
                mj_ind(ind) = mj
              enddo
            enddo
          endif
        endif

      enddo    !  l=0,lmx

    else

      do l = 0,lmx
        if(nkb(l,0,k) /= 0) then

          do m = -l,l
            do kk = 1,natom(k)
              do ms = -1,1,2
                ind = ind + 1
                k_ind(ind) = k
                kk_ind(ind) = kk
                l_ind(ind) = l
                m_ind(ind) = m
                ms_ind(ind) = ms
              enddo
            enddo
          enddo

        endif
      enddo    !  l=0,lmx

    endif

  enddo        !  k=1,ntype

! Clebsch-Gordan coefficients

  do l = 0,lmx
  do j = -2*lmx-1,2*lmx+1
    cg(l,j) = ZERO
  enddo
  enddo
  do l = 0,lmx
  do j = -2*l-1,2*l+1
    cg(l,j) = sqrt( ((2*l+1 + j)*UM) / ((4*l+2)*UM) )
  enddo
  enddo

! starts loop over g-vectors

!$omp  parallel do default(private)                                      &
!$omp& shared(kgv, bvec, bdot, ntype, natom, mtxd, rkpt, isort, rat)     &
!$omp& shared(lmx, delqnl, nqnl, vkb, fac, lso)                          &
!$omp& shared(anlsp, nanlsp)                                             &
!$omp& shared(k_ind, kk_ind, l_ind, m_ind, ms_ind)                       &
!$omp& shared(j_ind, mj_ind)                                             &
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

!   angular functions

    call ylm_rc(qcar, rylm, drylmdqc, d2rylmdqc2, .TRUE.,                &
        zylm, dzylmdqc, d2zylmdqc2, lmx, 0, ifail)

!   structure phase factor

    do k = 1,ntype
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

      do l = 0,lmx

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

          endif

          vqil(l,m,k) = vq(m)*C_I**l

        enddo

      enddo

    enddo

!   loop over second index

    do ind = 1,nanlsp

      k = k_ind(ind)
      kk = kk_ind(ind)
      l = l_ind(ind)

      if(lso(k)) then

!       pseudopotential is relativistic

        j = j_ind(ind)
        mj = mj_ind(ind)
        ms = j - 2*l

        if( mj == -2*l-1 ) then

          anlsp(2*i-1,ind) = C_ZERO
          anlsp(2*i  ,ind) = st(kk,k)*zylm(l,-l)*vqil(l,ms,k)

        elseif( mj == 2*l+1 ) then

          anlsp(2*i-1,ind) = st(kk,k)*zylm(l, l)*vqil(l,ms,k)
          anlsp(2*i  ,ind) = C_ZERO


        elseif( l > 0 ) then

          anlsp(2*i-1,ind) =    st(kk,k)*zylm(l,(mj-1)/2)*vqil(l,ms,k)*cg(l, ms*mj)
          anlsp(2*i  ,ind) = ms*st(kk,k)*vqil(l,ms,k)*zylm(l,(mj+1)/2)*cg(l,-ms*mj)

        endif

      else

!       pseudopotential is non-relativistic

        l = l_ind(ind)
        m = m_ind(ind)
        ms = ms_ind(ind)

        if(ms == 1) then
          anlsp(2*i-1,ind) = st(kk,k)*zylm(l,m)*vqil(l,0,k)
          anlsp(2*i  ,ind) = C_ZERO
        else
          anlsp(2*i-1,ind) = C_ZERO
          anlsp(2*i  ,ind) = st(kk,k)*zylm(l,m)*vqil(l,0,k)
        endif

      endif

    enddo


  enddo        !  i=1,mtxd
!$omp end parallel do

  do ind = 1,nanlsp
    k = k_ind(ind)
    l = l_ind(ind)
    if(lso(k)) then
      j = j_ind(ind)
      ms = j - 2*l
      xnlkbsp(ind) = UM*nkb(l,ms,k)
    else
      xnlkbsp(ind) = UM*nkb(l,0,k)
    endif
  enddo

  deallocate(st)
  deallocate(vqil)

  deallocate(lso)

  deallocate(k_ind, kk_ind, l_ind)
  deallocate(m_ind, ms_ind)
  deallocate(j_ind, mj_ind)

  return

end subroutine proj_nl_kb_so_c16
