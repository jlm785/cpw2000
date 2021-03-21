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

!>     This subroutine calculates the band structure along a path
!>     using atomic orbital interpolation

       subroutine ao_interpolation_out_band_fold_full(ioreplay,         &
     & pwline,title,subtitle,noiData,                                   &
     & emax,ztot,                                                       &
     & adot,ntype,natom,rat,                                            &
     & ng,kgv,                                                          &
     & icmplx,                                                          &
     & norbat,nqwf,delqwf,wvfao,lorb,                                   &
     & mxdtyp,mxdatm,mxdgve,mxdlqp,mxdlao)

!      Written by Carlos Loia Reis, May 2020, based on previous code,
!      ao_interpolation_out_band_fold.
!      Modified, documentation, May 2020. JLM

!      copyright  Carlos Loia Reis /INESC-MN.

!      version 4.98

       use NonOrthoInterp

       implicit none

       integer, parameter                 :: REAL64 = selected_real_kind(12)

!      input
       integer, intent(in)                ::  mxdtyp                     !<  array dimension of types of atoms
       integer, intent(in)                ::  mxdatm                     !<  array dimension of number of atoms of a given type
       integer, intent(in)                ::  mxdgve                     !<  array dimension for g-space vectors
       integer, intent(in)                ::  mxdlqp                     !<  array dimension for local potential
       integer, intent(in)                ::  mxdlao                     !<  array dimension of orbital per atom type

       integer, intent(in)                :: ioreplay                    !<  tape number for reproducing calculations

       character(len=50), intent(in)      ::  title                      !<  title for plots
       character(len=140), intent(in)     ::  subtitle                   !<  subtitle for plots
       character(len=60), intent(in)      ::  pwline                     !<  line with unfolding data

       type(noiData_t)                    ::  noiData                    !<  see NonOrthoInterp

       real(REAL64), intent(in)           ::  emax                       !<  largest kinetic energy included in hamiltonian diagonal. (hartree).
       real(REAL64), intent(in)           ::  ztot                       !<  total charge density (electrons/cell)

       real(REAL64), intent(in)           ::  adot(3,3)                  !<  metric in direct space
       integer, intent(in)                ::  ntype                      !<  number of types of atoms
       integer, intent(in)                ::  natom(mxdtyp)              !<  number of atoms of type i
       real(REAL64), intent(in)           ::  rat(3,mxdatm,mxdtyp)       !<  k-th component (in lattice coordinates) of the position of the n-th atom of type i

       integer, intent(in)                ::  ng                         !<  total number of g-vectors with length less than gmax
       integer, intent(in)                ::  kgv(3,mxdgve)              !<  i-th component (reciprocal lattice coordinates) of the n-th g-vector ordered by stars of increasing length

       integer, intent(in)                ::  icmplx                     !<  indicates if the structure factor is complex

       integer, intent(in)                ::  norbat(mxdtyp)             !<  number of atomic orbitals for atom k
       integer, intent(in)                ::  nqwf(mxdtyp)               !<  number of points for wavefunction interpolation for atom k
       real(REAL64), intent(in)           ::  delqwf(mxdtyp)             !<  step used in the wavefunction interpolation for atom k
       integer, intent(in)                ::  lorb(mxdlao,mxdtyp)        !<  angular momentum of orbital n of atom k
       real(REAL64), intent(in)      ::  wvfao(-2:mxdlqp,mxdlao,mxdtyp)  !<  (1/q**l) * wavefunction for atom k, ang. mom. l (unnormalized to vcell)

!      allocatable arrays

       complex(REAL64), allocatable       ::  psi_ao(:,:)
       complex(REAL64), allocatable       ::  ao_basis(:,:)
       complex(REAL64), allocatable       ::  ao_basis_so(:,:)

       real(REAL64), allocatable          ::  hdiag(:)                   !  hamiltonian diagonal
       integer, allocatable               ::  isort(:)                   !  g-vector associated with row/column i of hamiltonian
       real(REAL64), allocatable          ::  qmod(:)                    !  length of k+g-vector of row/column i
       real(REAL64), allocatable          ::  ekpg(:)                    !  kinetic energy (hartree) of k+g-vector of row/column i

       integer, allocatable               ::  infolcao(:,:)              !  information about the original atomic orbital.  (type of atom, atom of that type, n,l,m)


!      allocatable arrays for Brillouin zone path

       integer                            ::  nlines                     !  number of lines in reciprocal space
       integer, allocatable               ::  nkstep(:)                  !  number of steps in line
       logical, allocatable               ::  ljump(:)                   !  indicates if the new line contains a jump from the preceeding
       integer                            ::  nvert                      !  number of vertical lines in plot
       real(REAL64), allocatable          ::  xcvert(:)                  !  x coordinate of vertical line
       real(REAL64), allocatable          ::  xk(:)                      !  x coordinate of k-point in plot
       real(REAL64), allocatable          ::  rk(:,:)                    !  x coordinate of k-point in plot
       real(REAL64), allocatable          ::  rk_fld(:,:)                !  x coordinate of k-point in plot
       real(REAL64), allocatable          ::  e_of_k(:,:)                !  band energies of k-point in plot
       real(REAL64), allocatable          ::  e_of_k_so(:,:)             !  spin-orbit band energies of k-point in plot
       character(len=6), allocatable      ::  label(:)                   !  label of symmetry k-points
       real(REAL64), allocatable          ::  xklab(:)                   !  x coordinate of label

       real(REAL64),allocatable           ::  pkn(:,:)                   !  weight of unfolded band
       real(REAL64),allocatable           ::  pkn_so(:,:)                !  weight of unfolded band (spin-orbit)

!      local variables

       integer                            ::  mxddim                     !  array dimension for the hamiltonian
       real(REAL64)                       ::  bdot(3,3),vcell
       integer                            ::  mtxd                       !  dimension of the hamiltonian
       integer                            ::  nd, mxdorb, lmax, k
       logical                            ::  lnewanl

       logical                            ::  true_pkn

       integer           ::  neig                                        !  number of eigenvectors required (maybe modified on output)
       real(REAL64)      ::  rkpt(3)                                     !  j-th component in lattice coordinates of the k-point

       real(REAL64)      ::  eref                                        !  reference energy for plot
       integer           ::  nocc                                        !  number of occupied states (different color)
       integer           ::  nstyle                                      !  choice of plot style

       integer           ::  irk
       integer           ::  iotape
       integer           ::  nrk2

       integer           ::  iMinv(3,3)
       integer           ::  idet
       real(REAL64)      ::  avec(3,3),bvec(3,3), adot_pc(3,3)

       integer            ::  idum(3,3)
       character(len=6)   ::  fdum
       integer            ::  ioerr
       character(len=60)  ::  pwlinloc                                   !  local version of pwline

       real(REAL64), allocatable :: ev_interp(:)

       logical       ::  lkplusg                                         !  If true use the previous G-vectors (same mtxd and isort)

       integer neig_try
       character(len=1)            ::  yesno

!      counters

       integer    ::   i, j, n

!      constants

       real(REAL64), parameter     :: ZERO = 0.0_REAL64, UM = 1.0_REAL64
       complex(REAL64), parameter  ::  C_ZERO = cmplx(ZERO,ZERO,REAL64)
       complex(REAL64), parameter  ::  C_UM = cmplx(UM,ZERO,REAL64)


!-----------------------------------------------------------------------
!!!!      pwline from old PW_RHO_V.DAT does not work due to a bug in out_rho_v
!!!!      this is a workaround

!      checks it is the debugged version

       true_pkn = .false.

       read(pwline,'(3(3i4,2x),2x,a6)',IOSTAT=ioerr)                     &
     &           ((idum(i,j),i=1,3),j=1,3),fdum

       if(ioerr ==0 .and. adjustl(trim(fdum)) == 'fcc SL') then

         pwlinloc = pwline

       else

         write(6,*)
         write(6,*) '  Trying to use PW.DAT for the file identifier'
         write(6,*)

         open(unit=11,file='PW.DAT',status='OLD',form='FORMATTED',IOSTAT=ioerr)

         if(ioerr == 0) then

           read(11,'(a60)',IOSTAT=ioerr) pwlinloc
           close(unit=11)

         else

           write(6,*)
           write(6,*)  '  Not using unfolding'
           write(6,*)

           do i = 1,60
             pwlinloc(i:i) = ' '
           enddo

           write(pwlinloc,'(3(3i4,2x),2x,a6)') 1,0,0,  0,1,0,  0,0,1,  'fcc SL'

         endif

       endif
!!!!
!!!!        PrepFold needs avec
       call adot_to_avec_sym(adot,avec,bvec)
!!!!        ploting routines need metric of primitive cell !
       call Fold_Get_adot_pc(pwlinloc, avec, adot_pc)
!-----------------------------------------------------------------------

       iotape = 13
       call out_band_circuit_size('BAND_LINES.DAT', iotape, 1, adot_pc,  &    ! note call with adot_pc
     &                  neig, nrk2, nlines, nvert)

       allocate(xk(nrk2))
       allocate(rk(3,nrk2))
       allocate(rk_fld(3,nrk2))
       allocate(xcvert(nvert))
       allocate(ljump(nlines))
       allocate(nkstep(nlines))
       allocate(label(nvert+nlines))
       allocate(xklab(nvert+nlines))

       call out_band_get_circuit('BAND_LINES.DAT',iotape,1,adot_pc,     &     ! note call with adot_pc
     &                  xk,rk,xcvert,ljump,nkstep,label,xklab,          &
     &                  neig,nrk2,nlines,nvert)


      neig_try = nint(1.5*(ztot/2))
      if(neig_try <(ztot/2) + 8) then
        neig_try = (ztot/2) + 8
      endif

      write(6,*) ' Interpolated  band structure with unfolding '
      write(6,*) ' version 4.931 updated '
      write(6,*)
      write(6,*) ' I recomend the following value '
      write(6,*) ' for the number of bands:', neig_try
      write(6,*) ' do you want to accept? (y/n)'

      read(5,*) yesno
      write(ioreplay,*) yesno,'   accept suggested number of bands'

      if(yesno == 'y' .or. yesno == 'Y') then
        neig=neig_try
      else
        write(6,*) '  Using',neig,' bands'
      endif

!---
      write(6,*)
      write(6,*) 'do you want to compute the spectral weights? (y/n)'
      write(6,*) 'this requires additional computing resources'

      read(5,*) yesno
      write(ioreplay,*) yesno,'   true_pkn'

      if(yesno == 'y' .or. yesno == 'Y') then
        true_pkn = .true.
      endif



       allocate(e_of_k(neig,nrk2))
       allocate(e_of_k_so(2*neig,nrk2))

!-----------------------------------------------------------------------
       call Fold_Prep(pwlinloc, avec,rk, iMinv,idet, rk_fld, nrk2)
!-----------------------------------------------------------------------

!      finds mxddim, mxdbnd

!       mxdbnd = neig
!-----------------------------------------------------------------------
       allocate(pkn(nrk2,neig))
       allocate(pkn_so(nrk2,2*neig))
       allocate(ev_interp(noiData%nband))
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!      Spectral weights are not calculated for now
       pkn(:,:) = 1.0D0
       pkn_so(:,:) = 1.0D0
!-----------------------------------------------------------------------

!
!--------------- atomic orbital preparation section---------------------


      if (true_pkn) then


       call adot_to_bdot(adot,vcell,bdot)

!      finds mxdddim
       mxddim = 1
       do irk=1,nrk2
         do j=1,3
           rkpt(j) = rk_fld(j,irk)                                       ! computed in folded k-point
         enddo
         call size_mtxd(emax,rkpt,adot,ng,kgv,nd)
         if(nd > mxddim) mxddim = nd
       enddo

!      finds mxdorb
       mxdorb = 0
       lmax = 0
       do k=1,ntype
         do j=1,norbat(k)
           mxdorb = mxdorb + (2*lorb(j,k)+1)*natom(k)
           lmax = max(lorb(j,k),lmax)
         enddo
       enddo

       allocate(hdiag(mxddim))
       allocate(isort(mxddim))
       allocate(qmod(mxddim))
       allocate(ekpg(mxddim))


!      allocate wavefunctions needded for pkn calculation
       allocate(ao_basis(mxddim,mxdorb))

       if(noiData%lso==1) then
        allocate(psi_ao(2*mxddim,2*mxdorb))
        allocate(ao_basis_so(2*mxddim,2*mxdorb))
       else
        allocate(psi_ao(mxddim,mxdorb))
       endif

       endif

       do irk=1,nrk2

!        loop over k-points

         do j=1,3
           rkpt(j) = rk_fld(j,irk)                                       ! computed in folded k-point
         enddo

         write(*,'(i5,3f8.5)') irk, rkpt(1),rkpt(2),rkpt(3)

         call NonOrthoInterpRun(noiData,rkpt,ev_interp)


         if(noiData%lso==1) then
          do j=1,2*neig
            e_of_k_so(j,irk) = ev_interp(j)
          enddo
         else
          do j=1,neig
            e_of_k(j,irk) = ev_interp(j)
          enddo
         endif


         if (true_pkn) then
         ! -pkn calclulation section---------

         lkplusg = .false.

         call hamilt_struct(emax,rkpt,mtxd,isort,qmod,ekpg, lkplusg,     &
     &   ng,kgv,adot,                                                    &
     &   mxdgve,mxddim)


         lnewanl = .TRUE.

         allocate(infolcao(5,mxdorb))

         call atomic_orbital_c16(rkpt,mtxd,isort,icmplx,                 &
     &   mxdorb,ao_basis,infolcao,                                       &
     &   ng,kgv,                                                         &
     &   norbat,nqwf,delqwf,wvfao,lorb,                                  &
     &   ntype,natom,rat,adot,                                           &
     &   mxdtyp,mxdatm,mxdlqp,mxddim,mxdorb,mxdgve,mxdlao)

         deallocate(infolcao)


        if(noiData%lso==1) then
          do i=1,mxdorb
          do j=1,mtxd
            ao_basis_so(2*j-1,2*i         ) = ao_basis(j,i)
            ao_basis_so(2*j  ,2*i         ) = C_ZERO
            ao_basis_so(2*j-1,2*i-1       ) = C_ZERO
            ao_basis_so(2*j  ,2*i-1       ) = ao_basis(j,i)
          enddo
          enddo


          ! noiData%wrk has the eigenvectors corresponding to ev_interp

          call zgemm('N','N',2*mtxd,2*mxdorb,2*mxdorb,C_UM,ao_basis_so,2*mxddim,        &
      &        noiData%wrk,2*mxdorb,C_ZERO,psi_ao,2*mxddim)

          call Fold_GetPknSO(pkn_so,iMinv,idet,irk,kgv,isort,psi_ao, nrk2, 2*neig, mtxd, ng, mxdgve,2*mxddim,2*mxdorb)

        else

          call zgemm('N','N',mtxd,mxdorb,mxdorb,C_UM,ao_basis,mxddim,                   &
     &         noiData%wrk,mxdorb,C_ZERO,psi_ao,mxddim)

          call Fold_GetPkn(pkn,iMinv,idet,irk,kgv,isort,psi_ao, nrk2, neig, mtxd, ng, mxdgve,mxddim,mxdorb)

        endif

        ! end pkn calclulation section---------
         endif

       enddo

!      writes the output files for xmgrace

       iotape = 15
       nstyle = 2

      if(noiData%lso==1) then
!-----------------------------------------------------------------------
       n = min(nint(ztot + 0.01),2*neig)
       eref = e_of_k_so(n,1)
       do irk = 1,nrk2
       do j=1,n
         if(e_of_k_so(j,irk) > eref) eref = e_of_k_so(j,irk)
       enddo
       enddo

       noiData%eref  =eref

       nocc = n
!       do irk = 1,nrk2
!         if(e_of_k(n+1,irk) < eref) then
!           nocc = 0
!
!           exit
!
!         endif
!       enddo

!-----------------------------------------------------------------------
       call out_band_fold_xmgrace('band_so.agr',iotape,                  &
     &        title,subtitle,nstyle,                                     &
     &        pkn_so,2*neig,nrk2,xk,e_of_k_so,eref,nocc,                 &
     &        nvert,xcvert,nlines,ljump,nkstep,label,xklab)
!-----------------------------------------------------------------------
      else
!-----------------------------------------------------------------------
       n = min(nint(0.5*ztot + 0.01),neig)
       eref = e_of_k(n,1)
       do irk = 1,nrk2
       do j=1,n
         if(e_of_k(j,irk) > eref) eref = e_of_k(j,irk)
       enddo
       enddo

       noiData%eref  =eref

!       write(*,*) 'EREF IS' , noiData%eref

       nocc = n
!       do irk = 1,nrk2
!         if(e_of_k(n+1,irk) < eref) then
!           nocc = 0
!
!           exit
!
!         endif
!       enddo

!-----------------------------------------------------------------------
       call out_band_fold_xmgrace('band.agr',iotape,                     &
     &        title,subtitle,nstyle,                                     &
     &        pkn,neig,nrk2,xk,e_of_k,eref,nocc,                         &
     &        nvert,xcvert,nlines,ljump,nkstep,label,xklab)
!-----------------------------------------------------------------------
      endif

!-----------------------------------------------------------------------
       deallocate(rk_fld)
       deallocate(pkn)
       deallocate(pkn_so)
!-----------------------------------------------------------------------

       deallocate(nkstep)
       deallocate(ljump)

       deallocate(xcvert)
       deallocate(xk)
       deallocate(rk)
       deallocate(e_of_k)
       deallocate(e_of_k_so)
       deallocate(label)
       deallocate(xklab)
       deallocate(ev_interp)



       if (true_pkn) then
        deallocate(hdiag)
        deallocate(isort)
        deallocate(qmod)
        deallocate(ekpg)

        deallocate(ao_basis)

        if(noiData%lso==1) then
          deallocate(psi_ao)
          deallocate(ao_basis_so)
        else
          deallocate(psi_ao)
        endif
       endif


       write(*,*) 'Band Structure Calcultion Done.'

       return
       end subroutine ao_interpolation_out_band_fold_full
