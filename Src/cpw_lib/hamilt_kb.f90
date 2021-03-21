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

!>     sets up the hamiltonian for a given k-point
!>     for a Kleinman-Bylander type pseudopotential

       subroutine hamilt_kb(mtxd,hdiag,isort,qmod,                       &
     & ng,kgv,phase,conj,inds,kmax,indv,ek,                              &
     & sfact,veff,nqnl,delqnl,vkb,nkb,                                   &
     & ntype,adot,hamk,                                                  &
     & mxdtyp,mxdgve,mxdnst,mxdcub,mxdlqp,mxddim)

!      written february 8 1990. jlm
!      modified (conj) 20 march 1999. jlm
!      modified (4f) 21 february 2000. jlm
!      modified (f90) 15 January 2014.  jlm
!      modified, vkb,nkb 16 April 2014.  JLM
!      modified, f90, hamk storage, September 2015. JLM
!      modified, vkb normalization, 26 September 2015. JLM
!      Modified, documentation, January 2020. JLM
!      copyright INESC-MN/Jose Luis Martins

!      version 4.94

       implicit none

       integer, parameter          :: REAL64 = selected_real_kind(12)

!      input

       integer, intent(in)                ::  mxdtyp                     !<  array dimension of types of atoms
       integer, intent(in)                ::  mxdgve                     !<  array dimension for g-space vectors
       integer, intent(in)                ::  mxdnst                     !<  array dimension for g-space stars
       integer , intent(in)               ::  mxdcub                     !<  array dimension for 3-index g-space
       integer, intent(in)                ::  mxdlqp                     !<  array dimension for local potential
       integer, intent(in)                ::  mxddim                     !<  array dimension of plane-waves

       integer, intent(in)                ::  mtxd                       !<  dimension of the hamiltonian
       real(REAL64), intent(in)           ::  hdiag(mxddim)              !<  hamiltonian diagonal
       integer, intent(in)                ::  isort(mxddim)              !<  g-vector associated with row/column i of hamiltonian
       real(REAL64), intent(in)           ::  qmod(mxddim)               !<  length of k+g-vector of row/column i

       integer, intent(in)                ::  ng                         !<  total number of g-vectors with length less than gmax
       integer, intent(in)                ::  kgv(3,mxdgve)              !<  i-th component (reciprocal lattice coordinates) of the n-th g-vector ordered by stars of increasing length
       complex(REAL64), intent(in)        ::  phase(mxdgve)              !<  real part of the phase factor of G-vector n
       real(REAL64), intent(in)           ::  conj(mxdgve)               !<  is -1 if one must take the complex conjugate of x*phase
       integer, intent(in)                ::  inds(mxdgve)               !<  star to which g-vector n belongs
       integer, intent(in)                ::  kmax(3)                    !<  max value of kgv(i,n)
       integer, intent(in)                ::  indv(mxdcub)               !<  kgv(i,indv(jadd)) is the g-vector associated with jadd. jadd is defined by the g-vector components and kmax
       real(REAL64), intent(in)           ::  ek(mxdnst)                 !<  kinetic energy (hartree) of g-vectors in star j

       complex(REAL64), intent(in)        ::  veff(mxdnst)               !<  real part of the ionic potential (hartree) for the prototype g-vector in star j
       complex(REAL64), intent(in)        ::  sfact(mxdtyp,mxdnst)       !<  real part of the structure factor

       integer, intent(in)                ::  nqnl(mxdtyp)               !<  number of points for the non-local pseudopotential interpolation
       real(REAL64), intent(in)           ::  delqnl(mxdtyp)             !<  step used in the interpolation
       real(REAL64), intent(in)    ::    vkb(-2:mxdlqp,0:3,-1:1,mxdtyp)  !<  (1/q**l) * KB nonlocal pseudo. for atom k, ang. mom. l. normalized to vcell, hartree
       integer, intent(in)                ::  nkb(0:3,-1:1,mxdtyp)       !<  KB pseudo.  normalization for atom k, ang. mom. l

       integer, intent(in)                ::  ntype                      !<  number of types of atoms
       real(REAL64), intent(in)           ::  adot(3,3)                  !<  metric in real space

!      output

       complex(REAL64), intent(out)       ::  hamk(mxddim,mxddim)        !<  small hamiltonian

!      local varaibles

       real(REAL64)   :: qi,qj
       real(REAL64)   :: vqnl
       real(REAL64)   :: gamma
       real(REAL64)   :: plg(0:3)
       real(REAL64)   :: vqi(0:3)
       real(REAL64)   :: vqj(0:3)
       real(REAL64)   :: xi,xni,xj,xnj
       complex(REAL64)  :: hh
       real(REAL64)   ::  vcell, bdot(3,3),fac

!      parameters

       real(REAL64), parameter ::  EPS = 0.1E-8
       real(REAL64), parameter ::  ZERO = 0.0_REAL64, UM = 1.0_REAL64

!      counters

       integer   ::  i, j
       integer   ::  l, nt
       integer   ::  im, irowi, iadd, jadd, kadd
       integer   ::  ni, nj
       integer   ::  imx,jmx,kmx
       integer   ::  kxd,kyd,kzd


!      paranoid check, stops complaint about unused ng...

       if(ng > mxdgve) stop

       call adot_to_bdot(adot,vcell,bdot)
       fac = UM / sqrt(vcell)

!      set up matrix

!      it is stored in a one dimensional array in the packed form:

!      a11,a21,a22,a31,a32,a33,a41...   a(i,j)=a(i*(i-1)/2 +j)


       imx = 2*kmax(1) + 1
       jmx = 2*kmax(2) + 1
       kmx = 2*kmax(3) + 1

!      start loop over columns

       do i=1,mtxd
         irowi = (i*i-i)/2
         iadd = isort(i)
         qi = qmod(i)

!        diagonal of hamiltonian

         hamk(i,i) = cmplx(hdiag(i),ZERO,REAL64)

!        start loop over rows (off diagonal elements)

         if (i > 1) then
           im = i-1
           do j=1,im
             jadd = isort(j)
             hamk(i,j) = cmplx(ZERO,ZERO,REAL64)
             qj = qmod(j)
             kxd = kgv(1,iadd)-kgv(1,jadd)+kmax(1)+1
             if (kxd >= 1 .and. kxd <= imx) then
             kyd = kgv(2,iadd)-kgv(2,jadd)+kmax(2)+1
             if (kyd >= 1 .and. kyd <= jmx) then
             kzd = kgv(3,iadd)-kgv(3,jadd)+kmax(3)+1
             if (kzd >= 1 .and. kzd <= kmx) then

!              calculate address for arrays vloc(r)(i) and s(r)(i)

               kadd = ((kxd-1)*jmx+kyd-1)*kmx+kzd
               kadd = indv(kadd)
               if (kadd /= 0) then

!                local part of potential

                 hamk(i,j) = veff(inds(kadd))*conjg(phase(kadd))
                 if(conj(kadd) < ZERO) hamk(i,j) = conjg(hamk(i,j))

!                nonlocal part of potential

!                compute p (gamma)
!                         l

                 gamma = UM
                 if (qi > eps .and. qj > eps) then
                   gamma = (qi*qi + qj*qj - 2*ek(inds(kadd))) /          &
     &                     (qi*qj*2)
                 endif
                 plg(0) = UM
                 plg(1) = gamma
                 plg(2) = (3*gamma*gamma - UM)/ 2
                 plg(3) = gamma*(5*gamma*gamma - 3*UM)/ 2

!                add nonlocal part

                 do nt = 1,ntype
                   xni = qi/delqnl(nt)
                   ni  = nint(xni)
                   if(ni < 0) ni = 0
                   do l = 0,3
                     vqi(l) = ZERO
                     if (ni < nqnl(nt)-1) then
                       xi  = xni - UM*ni
                       if(nkb(l,0,nt) /= 0) then
                         vqi(l) = vkb(ni,l,0,nt) * (UM+xi) * (UM-xi)     &
     &                    +  (vkb(ni+1,l,0,nt)*(UM+xi)                   &
     &                       -vkb(ni-1,l,0,nt)*(UM-xi)) * xi / 2
                         vqi(l) = fac * vqi(l) * qi**l
                       endif
                     endif
                   enddo

                   xnj = qj/delqnl(nt)
                   nj  = nint(xnj)
                   if(nj < 0) nj = 0
                   do l = 0,3
                     vqj(l) = ZERO
                     if (nj < nqnl(nt)-1) then
                       xj  = xnj - UM*nj
                       if(nkb(l,0,nt) /= 0) then
                         vqj(l) = vkb(nj,l,0,nt) * (UM+xj) * (UM-xj)     &
     &                    +  (vkb(nj+1,l,0,nt)*(UM+xj)                   &
     &                       -vkb(nj-1,l,0,nt)*(UM-xj)) * xj / 2
                         vqj(l) = fac * vqj(l) * qj**l
                       endif
                     endif
                   enddo

                   vqnl = ZERO
                   do l = 0,3
                     vqnl = vqnl + vqi(l)*vqj(l)*nkb(l,0,nt)*plg(l)
                   enddo

                   hh = vqnl*sfact(nt,inds(kadd))*conjg(phase(kadd))
                   if(conj(kadd) < ZERO) hh = conjg(hh)
                   hamk(i,j) = hamk(i,j) + hh

                 enddo
               endif
             endif
             endif
             endif

           enddo
         endif
       enddo

!      completes the matrix

       do i = 1,mtxd
         if (i > 1) then
           im = i-1
           do j=1,im
             hamk(j,i) = conjg(hamk(i,j))
           enddo
         endif
       enddo

       return
       end subroutine hamilt_kb
