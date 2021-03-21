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

       subroutine proj_nl_kb_so_c16(rkpt,mtxd,isort,nanl,nanlso,        &
     & ng,kgv,                                                          &
     & nqnl,delqnl,vkb,nkb,                                             &
     & ntype,natom,rat,adot,                                            &
     & anlspin,anlso,xnlkbspin,xnlkbso,                                 &
     & mxdtyp,mxdatm,mxdlqp,mxddim,mxdanl,mxdaso,mxdgve)

!      subroutine sets up the special representation of the
!      non local part of the hamiltonian for a given k-point
!      Kleinman and Bylander pseudo-potential

!      SPIN-ORBIT VERSION

!      written f90/complex  19 June 2012
!      Modified 7 January 2014, style. JLM
!      modified, dimensions vkb, March 31, 2014. jlm
!      modified lmax bug, 16 July 2014. jlm

!      copyright INESC-MN/Jose Luis Martins

!      version 4.70


       implicit none
       integer, parameter          :: REAL64 = selected_real_kind(12)

!      input

       integer, intent(in)                ::  mxdtyp                     !  array dimension of types of atoms
       integer, intent(in)                ::  mxdatm                     !  array dimension of number of atoms of a given type
       integer, intent(in)                ::  mxddim                     !  array dimension of plane-waves
       integer, intent(in)                ::  mxdlqp                     !  array dimension for local potential
       integer, intent(in)                ::  mxdanl                     !  array dimension of number of projectors
       integer, intent(in)                ::  mxdaso                     !  array dimension of number of projectors with spin-orbit
       integer, intent(in)                ::  mxdgve                     !  array dimension of G-space vectors

       real(REAL64), intent(in)           ::  rkpt(3)                    !  k-point reciprocal lattice coordinates
       integer, intent(in)                ::  mtxd                       !  wavefunction dimension
       integer, intent(in)                ::  isort(mxddim)              !  G-vector corresponding to coefficient i of wavefunction
       integer, intent(in)                ::  ng                         !  size of g-space
       integer, intent(in)                ::  kgv(3,mxdgve)              !  G-vectors in reciprocal lattice coordinates
       integer, intent(in)                ::  ntype                      !  number of types of atoms
       integer, intent(in)                ::  natom(mxdtyp)              !  number of atoms of type i
       real(REAL64), intent(in)           ::  rat(3,mxdatm,mxdtyp)       !  lattice coordinates of atom j of type i
       real(REAL64), intent(in)           ::  adot(3,3)                  !  metric in real space
       integer, intent(in)                ::  nqnl(mxdtyp)               !  number of points for pseudo interpolation for atom k
       real(REAL64), intent(in)           ::  delqnl(mxdtyp)             !  step used in the pseudo interpolation for atom k
       real(REAL64), intent(in)    ::    vkb(-2:mxdlqp,0:3,-1:1,mxdtyp)  !  (1/q**l) * KB nonlocal pseudo. for atom k, ang. mom. l. NOT normalized to vcell, hartree
       integer, intent(in)                ::  nkb(0:3,-1:1,mxdtyp)       !  KB pseudo.  normalization for atom k, ang. mom. l

!      output

       integer, intent(out)               ::  nanl                       !  half of number of projectors without spin
       integer, intent(out)               ::  nanlso                     !  number of projectors witn spin
       complex(REAL64), intent(out)       ::  anlspin(2*mxddim,2*mxdanl) !  KB projectors without spin-orbit
       complex(REAL64), intent(out)       ::  anlso(2*mxddim,mxdaso)     !  KB projectors with spin-orbit
       real(REAL64), intent(out)          ::  xnlkbspin(2*mxdanl)        !  KB normalization without spin-orbit
       real(REAL64), intent(out)          ::  xnlkbso(mxdaso)            !  KB normalization with spin-orbit

!      local variables

       complex(REAL64), allocatable  ::  st(:)

       real(REAL64)    ::  avec(3,3),bvec(3,3)
       real(REAL64)    ::  vcell, bdot(3,3),fac
       integer         ::  n2mj
       real(REAL64)    ::  qi,qk(3),qcar(3),fi,xni,xi,vq(-1:1)

       real(REAL64)    ::  rylm(0:3,-3:3)                 !  r**l sqrt(4 pi / 2l+1) "Y_lm(l,m)+-Y_lm(l,-m)"
       real(REAL64)    ::  drylmdqc(3,0:3,-3:3)           !  unused
       real(REAL64)    ::  d2rylmdqc2(3,3,0:3,-3:3)       !  unused

       complex(REAL64) ::  zylm(0:3,-3:3)                 !  r**l sqrt(4 pi / 2l+1) Y_lm(theta,phi)
       complex(REAL64) ::  dzylmdqc(3,0:3,-3:3)           !  unused
       complex(REAL64) ::  d2zylmdqc2(3,3,0:3,-3:3)       !  unused

       complex(REAL64) ::  pitol                            !   i**l

       integer         ::  lmax, ifail

!      constants

       real(REAL64), parameter :: ZERO = 0.0_REAL64, UM = 1.0_REAL64
       real(REAL64), parameter :: PI = 3.14159265358979323846_REAL64

!      counters

       integer         ::  ind,i,k,l,m,kk,ni,ind_nl


!      paranoid check, stops complaint about unused ng...

       if(ng > mxdgve) stop

       allocate(st(mxdatm))

!      cell stuff

       call adot_to_avec(adot,avec,bvec)
       call adot_to_bdot(adot,vcell,bdot)
       fac = UM / sqrt(vcell)

!      checks dimensions and maximum angular momentum

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
         write(6,'("  proj_nl_kb_so(1):    increase mxdanl from ",i8,    &
     &            " to ",i8)')  mxdanl,nanl

         stop

       endif

       if(lmax > 3) then
         write(6,'("  proj_nl_kb:    lmax = ",i8,                        &
     &            " > 3 (max allowed in code)")')  lmax

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
         write(6,'("  proj_nl_kb_so:    increase mxdaso from ",i8,       &
     &            " to ",i8)')  mxdaso,nanlso

         stop

       endif

!      starts loop over g-vectors

       do i=1,mtxd
         qk(1) = rkpt(1) + UM*(kgv(1,isort(i)))
         qk(2) = rkpt(2) + UM*(kgv(2,isort(i)))
         qk(3) = rkpt(3) + UM*(kgv(3,isort(i)))
         qcar(1) = bvec(1,1)*qk(1)+bvec(1,2)*qk(2)+bvec(1,3)*qk(3)
         qcar(2) = bvec(2,1)*qk(1)+bvec(2,2)*qk(2)+bvec(2,3)*qk(3)
         qcar(3) = bvec(3,1)*qk(1)+bvec(3,2)*qk(2)+bvec(3,3)*qk(3)
         qi = qcar(1)*qcar(1) + qcar(2)*qcar(2) + qcar(3)*qcar(3)
         qi = sqrt(qi)

!        angular functions

         call ylm_rc(qcar,rylm,drylmdqc,d2rylmdqc2,.TRUE.,               &
     &            zylm,dzylmdqc,d2zylmdqc2,lmax,0,ifail)

!        starts loop over second index

         ind = 0
         ind_nl = 0

         do k=1,ntype

!          structure phase factor

           do kk=1,natom(k)
             fi = kgv(1,isort(i))*rat(1,kk,k) +                          &
     &            kgv(2,isort(i))*rat(2,kk,k) +                          &
     &            kgv(3,isort(i))*rat(3,kk,k)
             fi = -2*PI*fi
             st(kk) = cmplx(cos(fi),sin(fi),REAL64)
           enddo

           do l=0,3

             xni = qi/delqnl(k)
             ni  = nint(xni)
             if(ni < 0) ni = 0

             do m=-1,1
               vq(m) = ZERO
               if(ni < nqnl(k)-1) then
                 xi  = xni - UM*ni
                 vq(m) = vkb(ni,l,m,k) * (UM+xi) * (UM-xi)               &
     &              + (UM/2) * (vkb(ni+1,l,m,k)*(UM+xi)                  &
     &                         -vkb(ni-1,l,m,k)*(UM-xi)) * xi
                 vq(m) = fac * vq(m)
               endif
             enddo

             pitol = cmplx(ZERO,UM,REAL64)**l

             if(nkb(l,0,k) /= 0) then
               do m=-l,l
                 do kk=1,natom(k)
                   ind = ind + 1
                   anlspin(2*i-1,2*ind-1) = pitol*st(kk)*zylm(l,m)*vq(0)
                   anlspin(2*i  ,2*ind-1) = ZERO
                   anlspin(2*i-1,2*ind  ) = ZERO
                   anlspin(2*i  ,2*ind  ) = pitol*st(kk)*zylm(l,m)*vq(0)
                 enddo
               enddo
             endif

             if(nkb(l,1,k) /= 0) then
!              j = l+1/2  m_j = -l-1/2
               do kk=1,natom(k)
                 ind_nl = ind_nl + 1
                 anlso(2*i-1,ind_nl) = ZERO
                 anlso(2*i  ,ind_nl) = pitol*st(kk)*zylm(l,-l)*vq(1)
               enddo
               if(l > 0) then
                 do n2mj = -(2*l-1),2*l-1,2
!                  j = l+1/2   m_j = n2mj /2
                   do kk=1,natom(k)
                     ind_nl = ind_nl + 1
                     anlso(2*i-1,ind_nl) = pitol*st(kk)*vq(1)            &
     &                     *sqrt( ((2*l+1+n2mj)*UM) / (4*l+2) )          &
     &                               *zylm(l,(n2mj-1)/2)
                     anlso(2*i  ,ind_nl) = pitol*st(kk)*vq(1)            &
     &                     *sqrt( ((2*l+1-n2mj)*UM) / (4*l+2) )          &
     &                               *zylm(l,(n2mj+1)/2)
                   enddo
                 enddo
               endif
!              j = l+1/2  m_j = l+1/2
               do kk=1,natom(k)
                 ind_nl = ind_nl + 1
                 anlso(2*i-1,ind_nl) = pitol*st(kk)*zylm(l,l)*vq(1)
                 anlso(2*i  ,ind_nl) = ZERO
               enddo
             endif                     ! nkb(l,1,k) /= 0

             if(nkb(l,-1,k) /= 0) then
               if(l > 0) then
                 do n2mj = -(2*l-1),2*l-1,2
!                  j = l-1/2   m_j = n2mj /2
                   do kk=1,natom(k)
                     ind_nl = ind_nl + 1
                     anlso(2*i-1,ind_nl) = pitol*st(kk)*vq(-1)           &
     &                     *sqrt( ((2*l+1-n2mj)*UM) / (4*l+2) )          &
     &                               *zylm(l,(n2mj-1)/2)
                     anlso(2*i  ,ind_nl) =-pitol*st(kk)*vq(-1)           &
     &                     *sqrt( ((2*l+1+n2mj)*UM) / (4*l+2) )          &
     &                               *zylm(l,(n2mj+1)/2)
                   enddo
                 enddo
               endif
             endif                     ! nkb(l,-1,k) /= 0

           enddo    !  l=0,3

         enddo      !  k=1,ntype

       enddo        !  i=1,mtxd

!      sign of projector normalization

       ind = 0
       do k=1,ntype
         do l=0,3
           if(nkb(l,0,k) /= 0) then
             do m=-l,l
               do kk=1,natom(k)
                 ind = ind + 1
                 xnlkbspin(2*ind-1) = UM*nkb(l,0,k)
                 xnlkbspin(2*ind  ) = UM*nkb(l,0,k)
               enddo
             enddo
           endif
         enddo
       enddo

       ind = 0
       do k=1,ntype
         do l=0,3
           if(nkb(l,1,k) /= 0) then
             do m=-2*l-1,2*l+1,2
               do kk=1,natom(k)
                 ind = ind + 1
                 xnlkbso(ind) = UM*nkb(l,1,k)
               enddo
             enddo
           endif
           if(nkb(l,-1,k) /= 0) then
             if(l > 0) then
               do m=-2*l+1,2*l-1,2
                 do kk=1,natom(k)
                   ind = ind + 1
                   xnlkbso(ind) = UM*nkb(l,-1,k)
                 enddo
               enddo
             endif
           endif
         enddo
       enddo

       deallocate(st)

       return

       end subroutine proj_nl_kb_so_c16
