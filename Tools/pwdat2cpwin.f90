       program pwdat2cpwin

!      converts the old PW.DAT file format to the new cpw.in format

!      copyright Jose Luis Martins/INESC-MN

!      version 4.93   October 2018

       implicit none
       
       integer, parameter          :: REAL64 = selected_real_kind(12)

!      dimensions

       integer                            ::  mxdtyp                     !  array dimension of types of atoms
       integer                            ::  mxdatm                     !  array dimension of number of atoms of a given type

!      crystal structure

       real(REAL64)                       ::  alatt                      !  lattice constant

       real(REAL64)                       ::  adot(3,3)                  !  metric in direct space in atomic units (Bohr radius)

       integer                            ::  ntype                      !  number of types of atoms
       integer, allocatable               ::  natom(:)                   !  number of atoms of type i
       real(REAL64), allocatable          ::  rat(:,:,:)                 !  lattice coordinates of atom j of type i
       real(REAL64), allocatable          ::  atmass(:)                  !  atomic mass of atoms of type i
       character(len=2), allocatable      ::  nameat(:)                  !  chemical symbol for the type i

!      reciprocal space

       real(REAL64)                       ::  emax                       !  kinetic energy cutoff of plane wave expansion (Hartree).

!      k-point data

       integer                            ::  nx, ny, nz                 !  size of the integration mesh in k-space (nx*ny*nz)
       real(REAL64)                       ::  sx, sy, sz                 !  offset of the integration mesh (usually 0.5)

       integer                            ::  nbandin                    !  target for number of bands      
                 
       character(len=250)                 ::  meta_pwdat                 !  metadata from cpw_in or PW.DAT
       integer             ::  ioerr                                     !  input/output error

       call size_mxdtyp_mxdatm('PW.DAT',21,mxdtyp,mxdatm)

       open(UNIT=5,FILE='PW.DAT',STATUS='OLD',FORM='FORMATTED')

       read(5,'(a250)',iostat=ioerr) meta_pwdat
       if(ioerr == 0) then
         write(6,*) meta_pwdat(1:60)
         write(6,*) meta_pwdat(61:110)
         write(6,*) meta_pwdat(111:250)
       else
         write(6,*) '  cpw:   Unable to read titles'
       endif
       write(6,*)
       write(6,*)

       allocate(natom(mxdtyp))
       allocate(rat(3,mxdatm,mxdtyp))
       allocate(atmass(mxdtyp))
       allocate(nameat(mxdtyp))

!      gets crystal data

       call read_data(adot,ntype,natom,nameat,rat,atmass,alatt,        &
     & mxdtyp,mxdatm)

       read(5,'(1x,f9.4)') emax

       read(5,'(i5,15x,3i5,5x,3f5.2)') nbandin,nx,ny,nz,sx,sy,sz

       close(UNIT=5)
     
       call write_cpwout(meta_pwdat,                                     &
     & adot,ntype,natom,nameat,rat,atmass,alatt,                         &
     & emax,nbandin,nx,ny,nz,sx,sy,sz,                                   &
     & .FALSE.,.FALSE.,                                                  &
     & mxdtyp,mxdatm)

       stop
       end program pwdat2cpwin
