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
!>
!>  \author       Jose Luis Martins
!>  \version      5.0.3
!>  \date         May 12 1990, 7 December 2021.
!>  \copyright    GNU Public License v2

subroutine proj_nl_kb_c16(rkpt, mtxd, isort, nanl,                       &
    ng, kgv,                                                             &
    nqnl, delqnl, vkb, nkb,                                              &
    ntype, natom, rat, adot,                                             &
    anlga, xnlkb,                                                        &
    mxdtyp, mxdatm, mxdlqp, mxddim, mxdanl, mxdgve)

! written f90/complex  19 June 2012
! Modified 7 January 2014, style. JLM
! modified, dimensions vkb, real rylm, March 31, 2014. jlm
! Adapted Summer 2015. vkb normalization. JLM
! Modified documentation, January 2020. JLM
! openmp version 19 January 2020. JLM
! LMAX,lmx, 7 December 2021. JLM
! copyright INESC-MN/Jose Luis Martins


  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)
  integer, parameter          :: LMAX = 3                                !  hard coded max ang. mom.

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

  integer, intent(in)                ::  nanl                            !<  number of projectors without spin
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

  complex(REAL64), intent(out)       ::  anlga(mxddim,mxdanl)            !<  KB projectors without spin-orbit
  real(REAL64), intent(out)          ::  xnlkb(mxdanl)                   !<  KB normalization without spin-orbit

! allocatable local variables

  complex(REAL64), allocatable  ::  st(:,:)
  complex(REAL64), allocatable  ::  vqil(:,:)

  integer, allocatable          ::  k_ind(:)
  integer, allocatable          ::  kk_ind(:)
  integer, allocatable          ::  l_ind(:)
  integer, allocatable          ::  m_ind(:)


! local variables

  real(REAL64)    ::  avec(3,3),bvec(3,3)
  real(REAL64)    ::  vcell, bdot(3,3),fac
  real(REAL64)    ::  qi,qk(3),qcar(3),fi,xni,xi,vq0

  real(REAL64)    ::  rylm(0:LMAX,-LMAX:LMAX)                 !   r**l sqrt(4 pi / 2l+1) "Y_lm(l,m)+-Y_lm(l,-m)"
  real(REAL64)    ::  drylmdqc(3,0:LMAX,-LMAX:LMAX)           !   unused
  real(REAL64)    ::  d2rylmdqc2(3,3,0:LMAX,-LMAX:LMAX)       !   unused

  complex(REAL64) ::  zylm(0:LMAX,-LMAX:LMAX)                 !   unused
  complex(REAL64) ::  dzylmdqc(3,0:LMAX,-LMAX:LMAX)           !   unused
  complex(REAL64) ::  d2zylmdqc2(3,3,0:LMAX,-LMAX:LMAX)       !   unused

  integer         ::  lmx, ifail

! constants

  real(REAL64), parameter    :: ZERO = 0.0_REAL64, UM = 1.0_REAL64
  complex(REAL64), parameter :: C_I = cmplx(ZERO,UM,REAL64)
  real(REAL64), parameter    :: PI = 3.14159265358979323846_REAL64

! counters

  integer         ::  ind
  integer         ::  k, kk, l, m
  integer         ::  i, ni


! paranoid check, stops complaint about unused ng...

  if(ng > mxdgve) stop

  allocate(st(mxdatm,mxdtyp))
  allocate(vqil(0:LMAX,mxdtyp))

  allocate(k_ind(nanl), kk_ind(nanl), l_ind(nanl), m_ind(nanl))

! cell stuff

  call adot_to_avec(adot, avec, bvec)
  call adot_to_bdot(adot, vcell, bdot)
  fac = UM / sqrt(vcell)

! checks dimensions and maximum angular momentum

  ind = 0
  lmx = 0
  do k=1,ntype
    do l = 0,LMAX
      if(nkb(l,0,k) /= 0) then
        ind = ind + (2*l+1)*natom(k)
        lmx = max(lmx,l)
      endif
    enddo
  enddo

  if(ind /= nanl) then
    write(6,'("  proj_nl_kb:    inconsistent value of nanl ",i8,         &
           &  " should be ",2i8)')  nanl, ind, ng

    stop

  endif

  if(nanl > mxdanl) then
    write(6,'("  proj_nl_kb:    increase mxdanl from ",i8," to ",i8)')   &
               mxdanl,nanl

    stop

  endif

  if(lmx > LMAX) then
    write(6,'("  proj_nl_kb:    lmx = ",i8,                               &
           &  " > ",i3," (max allowed in code)")')  lmx, LMAX

    stop

  endif

! fills indexation array

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

! starts loop over g-vectors

!$omp  parallel do default(private)                                      &
!$omp& shared(kgv, bvec, ntype, natom, mtxd, rkpt, isort, rat)           &
!$omp& shared(lmx, delqnl, nqnl, vkb, fac)                               &
!$omp& shared(nanl, anlga, k_ind, kk_ind, l_ind, m_ind)
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

    call ylm_rc(qcar, rylm, drylmdqc, d2rylmdqc2, .FALSE.,               &
        zylm, dzylmdqc, d2zylmdqc2, lmx, 0, ifail)

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
      do l = 0,lmx

        xni = qi/delqnl(k)
        ni  = nint(xni)
        if(ni < 0) ni = 0

        vq0 = ZERO
        if(ni < nqnl(k)-1) then
          xi  = xni - UM*ni
          vq0 = vkb(ni,l,0,k) * (UM+xi) * (UM-xi)                         &
               + (UM/2) * (vkb(ni+1,l,0,k)*(UM+xi)                        &
                          -vkb(ni-1,l,0,k)*(UM-xi)) * xi
          vq0 = fac * vq0
        endif

        vqil(l,k) = vq0*C_I**l

      enddo

    enddo

!   starts loop over second index

    do ind = 1,nanl

      kk = kk_ind(ind)
      k = k_ind(ind)
      l = l_ind(ind)
      m = m_ind(ind)

      anlga(i,ind) = st(kk,k)*rylm(l,m)*vqil(l,k)

    enddo

  enddo
!$omp end parallel do


! sign of projector normalization

  do ind = 1,nanl
    xnlkb(ind) = UM*nkb(l_ind(ind),0,k_ind(ind))
  enddo

  deallocate(st)
  deallocate(vqil)

  deallocate(k_ind, kk_ind, l_ind, m_ind)

  return
end subroutine proj_nl_kb_c16
