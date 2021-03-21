!>  generates a formatted PW_try.DAT file to use in the old cpw
!>  program. The hermetic and formatted PW.DAT file was
!>  supposed to restrict the use of cpw to the initiated...
!>  For the new code cpw2000 it generates a cpw.in input file
!>  which has a modern format (SIESTA style) that can be parsed.

program gen_PW

! written by J.L.Martins and J.M.Pacheco
! July-October 2001
! Modernized and cpw.in, 15 February 2021. JLM

! version 4.99

  implicit none

  integer, parameter          ::  REAL64 = selected_real_kind(12)

  integer                       ::  mxdtyp
  integer                       ::  mxdatm

  integer                       ::  ibravais
  real(REAL64)                  ::  avec(3,3)
  real(REAL64)                  ::  conv(3,3)
  character(len = 20)           ::  title

  integer                       ::  nband
  integer                       ::  nx,ny,nz
  real(REAL64)                  ::  emax
  real(REAL64)                  ::  sx,sy,sz

  integer                       ::  ntype

  character(len = 200)          ::  source

! allocatable arrays

  integer, allocatable          ::  natom(:)
  character(len=2), allocatable ::  nameat(:)
  real(REAL64), allocatable     ::  rat(:,:,:)

! other variables

  real(REAL64)                  ::  vec(3,3)
  real(REAL64)                  ::  alatt
  character(len = 1)            ::  iop(3,3)
  integer                       ::  ncount

  integer                       ::  iopw
  character(len = 10)           ::  filepw

  integer                       ::  ioreplay
  character(len = 16)           ::  filereplay

  integer                       ::  iocpw
  character(len = 6)            ::  filecpw

! counter

  integer               ::  i, j


  iopw = 10
  filepw = 'PW_try.DAT'
  open(unit = iopw, file = filepw, form = 'FORMATTED', status='UNKNOWN')

  ioreplay = 35
  filereplay = 'replay_genpw.dat'
  open(unit = ioreplay, file = filereplay, form = 'FORMATTED', status='UNKNOWN')

  iocpw = 11
  filecpw = 'cpw.in'
  open(unit = iocpw, file = filecpw, status='UNKNOWN', form='FORMATTED')

  write(6,*) '  type in a title (20 characters) for the calculation'
  write(6,*) '  that title is arbitrary...'

  read(5,*) title
  write(ioreplay,'(a20,10x,"title")' ) title
  write(iopw,'(a20)' ) title

! now do the lattice

  call gen_latt(ioreplay, alatt, iop, vec, avec, conv)

  write(iopw,'(f12.6)') alatt
  do i=1,3
    write(iopw,'(3(2x,a1,f12.6))') (iop(j,i),vec(j,i),j=1,3)
  enddo


  call gen_basis_ntype(ioreplay, ntype)

  mxdtyp = ntype

  allocate(natom(mxdtyp))
  allocate(nameat(mxdtyp))

  call gen_basis_comp(ioreplay, ntype, natom, nameat, ncount, mxdtyp)

  mxdatm = ncount

  allocate(rat(3,mxdatm,mxdtyp))

  call gen_basis(ioreplay, avec, conv, ntype, natom, nameat, rat,        &
      mxdatm, mxdtyp)

  write(6,*)' enter cut-off energy in Hartree'
  read(5,*) emax
  write(ioreplay,'(1x,f9.4,5x,"emax")') emax

  write(iopw,'(1x,f9.4)') emax

  write(6,*)' enter number of one-electron states'
  read(5,*) nband
  write(ioreplay,'(1x,i10,5x,"nband")') nband

  write(6,*)' enter size of mesh (nx,ny,nz) of integration k-points'
  read(5,*) nx,ny,nz
  write(ioreplay,'(5x,3i5,5x,"nkxyz")') nx,ny,nz

  write(6,*)' enter shifts (x,y,z) of k-point mesh'
  read(5,*) sx,sy,sz
  write(ioreplay,'(5x,3f14.8,5x,"shiftxyz")') sx,sy,sz

  write(iopw,'(i5,15x,3i5,5x,3f5.2)') nband, nx,ny,nz, sx,sy,sz

  write(6,*)' enter (in quotes if URL) the source of your structure (max 200 char)'
  read(5,*) source
  write(ioreplay,*) '"',trim(adjustl(source)),'"'

  call write_cpwin(iocpw, title, trim(adjustl(source)), alatt,           &
      avec, ntype, natom, nameat, rat,                                   &
      emax, nband, nx,ny,nz, sx,sy,sz,                                   &
      mxdtyp, mxdatm)


  close(unit = iopw)
  close(unit = ioreplay)
  close(unit = iocpw)

  stop
end program gen_PW



subroutine gen_basis_ntype(ioreplay, ntype)

! gets the number of types of atoms

  implicit none

  integer, parameter          ::  REAL64 = selected_real_kind(12)

! input

  integer, intent(in)           ::  ioreplay

! output

  integer, intent(out)          ::  ntype



  write(6,*)' enter number of diff elements'
  read(5,*) ntype
  write(10,'(i5)') ntype
  write(ioreplay,'(1x,i5,5x,"ntype")') ntype


  return
end subroutine gen_basis_ntype


subroutine gen_basis_comp(ioreplay, ntype, natom, nameat, ncount, mxdtyp)

! generate the basis (atomic positions)

  implicit none

  integer, parameter          ::  REAL64 = selected_real_kind(12)

! input

  integer, intent(in)           ::  ioreplay
  integer, intent(in)           ::  mxdtyp
  integer, intent(in)           ::  ntype

! output

  integer, intent(out)          ::  natom(mxdtyp)
  character(len=2),intent(out)  ::  nameat(ntype)

  integer, intent(out)          ::  ncount

! counters

  integer                       ::  nt


! enter composition

  ncount = 0
  do nt = 1,ntype

    nameat(nt) = '  '
    write(6,*)' enter name of element type (without blanks)',nt
    read(5,*) nameat(nt)
    write(ioreplay,'(1x,a2,5x,"nameat for ",i5)') nameat(nt),nt

    write(6,'(" enter number of atoms of element ",a2)') nameat(nt)
    read(5,*) natom(nt)
    write(ioreplay,'(1x,i5,5x,"natom for ",a2)') natom(nt), nameat(nt)

    ncount = max(ncount,natom(nt))

  enddo

  return
end subroutine gen_basis_comp




subroutine gen_basis(ioreplay, avec, conv, ntype, natom, nameat, rat,    &
         mxdatm, mxdtyp)

! generate the basis (atomic positions)

  implicit none

  integer, parameter          ::  REAL64 = selected_real_kind(12)

! input

  integer, intent(in)           ::  ioreplay
  integer, intent(in)           ::  mxdatm
  integer, intent(in)           ::  mxdtyp

  real(REAL64), intent(in)      ::  avec(3,3)
  real(REAL64), intent(in)      ::  conv(3,3)

  integer, intent(in)           ::  ntype
  integer, intent(in)           ::  natom(mxdtyp)
  character(len=2), intent(in)  ::  nameat(mxdtyp)

! output

  real(REAL64), intent(out)     ::  rat(3,mxdatm,mxdtyp)

! other variables

  real(REAL64)                  ::  bvec(3,3), vcell
  integer                       ::  iflag

  character(len=1)              ::  jop(3) ,yesno
  character(len=2)              ::  elem
  character(len=40)             ::  filnam
  real(REAL64)                  ::  coor(3),cart(3)
  real(REAL64)                  ::  xx
  integer                       ::  ixx
  integer                       ::  ncount

! constants

  real(REAL64), parameter  ::  ZERO = 0.0_REAL64, UM = 1.0_REAL64
  real(REAL64), parameter  ::  EPS = 1.0E-10_REAL64
  real(REAL64), parameter  ::  EPS2 = 1.0E-04_REAL64

! counters

  integer                       ::  nt, i, j, k


! reciprocal lattice

  bvec(1,1) = avec(2,2)*avec(3,3) - avec(3,2)*avec(2,3)
  bvec(2,1) = avec(3,2)*avec(1,3) - avec(1,2)*avec(3,3)
  bvec(3,1) = avec(1,2)*avec(2,3) - avec(2,2)*avec(1,3)
  bvec(1,2) = avec(2,3)*avec(3,1) - avec(3,3)*avec(2,1)
  bvec(2,2) = avec(3,3)*avec(1,1) - avec(1,3)*avec(3,1)
  bvec(3,2) = avec(1,3)*avec(2,1) - avec(2,3)*avec(1,1)
  bvec(1,3) = avec(2,1)*avec(3,2) - avec(3,1)*avec(2,2)
  bvec(2,3) = avec(3,1)*avec(1,2) - avec(1,1)*avec(3,2)
  bvec(3,3) = avec(1,1)*avec(2,2) - avec(2,1)*avec(1,2)

  vcell = bvec(1,1)*avec(1,1) + bvec(2,1)*avec(2,1) + bvec(3,1)*avec(3,1)

  if(abs(vcell) < EPS) then
    write(6,'("  zero cell volume",e12.3)') vcell
    stop
  endif

  do j = 1,3
    bvec(1,j) = bvec(1,j)/vcell
    bvec(2,j) = bvec(2,j)/vcell
    bvec(3,j) = bvec(3,j)/vcell
  enddo

! type coordinates

  write(6,*) ' What coordinates will you be entering? '
  write(6,*) ' PRIMITIVE LATTICE (1)'
  write(6,*) ' CARTESIAN (2)'
  write(6,*) ' CONVENTIONAL LATTICE (3)'
  read(5,*) iflag
  write(ioreplay,'(1x,i5,5x,"coordinate type")') iflag

  if(iflag == 1) then
    write(6,*) '  Remember your primitive vectors are '
    write(6,*)
    write(6,'(a5,3x,3f10.4)') ' a1: ',(avec(j,1),j=1,3)
    write(6,'(a5,3x,3f10.4)') ' a2: ',(avec(j,2),j=1,3)
    write(6,'(a5,3x,3f10.4)') ' a3: ',(avec(j,3),j=1,3)
    write(6,*)
  elseif(iflag == 3) then
    write(6,*) '  Remember your conventional vectors are '
    write(6,*) '  (enter only atoms in primitive cell!)'
    write(6,*)
    write(6,'(a5,3x,3f10.4)') ' c1: ',(conv(j,1),j=1,3)
    write(6,'(a5,3x,3f10.4)') ' c2: ',(conv(j,2),j=1,3)
    write(6,'(a5,3x,3f10.4)') ' c3: ',(conv(j,3),j=1,3)
    write(6,*)
  else
    iflag = 2
    write(6,*) '  Expecting cartesian coordinates '
    write(6,*)
  endif

! enter coordinates

  do nt = 1,ntype
    write(6,*) ' for atoms of ', nameat(nt)

    write(10,'(i5,3x,a2)') natom(nt), nameat(nt)

    do j = 1,natom(nt)

      write(6,*)' enter coordinates of atom number',j
      read(5,*) (coor(k),k=1,3)
      write(ioreplay,'(1x,3f14.6,5x,"coordinates for ",a2,2x,i5)')       &
               (coor(k),k=1,3), nameat(nt), j

      if(iflag == 1) then
        do k=1,3
          rat(k,j,nt) = coor(k)
        enddo
      elseif(iflag == 2) then
        do k=1,3
          rat(k,j,nt) = ZERO
          do i=1,3
            rat(k,j,nt) = rat(k,j,nt) + bvec(i,k)*coor(i)
          enddo
        enddo
      else
        do k=1,3
          cart(k) = ZERO
          do i=1,3
            cart(k) = cart(k) + conv(k,i)*coor(i)
          enddo
        enddo
        do k=1,3
          rat(k,j,nt) = ZERO
          do i=1,3
            rat(k,j,nt) = rat(k,j,nt) + bvec(i,k)*cart(i)
          enddo
        enddo
      endif

!     brings all atoms close to the origin

      do k=1,3
        rat(k,j,nt) = rat(k,j,nt) - nint(rat(k,j,nt))*UM
      enddo

      do k=1,3
        jop(k) = ' '
        coor(k) = rat(k,j,nt)
        xx = 3*rat(k,j,nt)
        ixx = nint(xx)
        if(abs(xx - ixx*UM) < EPS2 .AND. mod(ixx,3) /= 0) then
          coor(k) = ixx*UM
          jop(k) = 'T'
        endif
      enddo

      write(10,'(3(2x,a1,f12.6))')(jop(k),coor(k),k=1,3)

    enddo
  enddo

  return
end subroutine gen_basis



subroutine gen_latt(ioreplay, alatt, iop, vec, avec, conv)

! writes all the stuff that is relevant to the lattice

  implicit none

! version 4.99

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)           ::  ioreplay

! output

  real(REAL64), intent(out)     ::  alatt
  character(len=1), intent(out) ::  iop(3,3)
  real(REAL64), intent(out)     ::  vec(3,3)

  real(REAL64), intent(out)     ::  avec(3,3)
  real(REAL64), intent(out)     ::  conv(3,3)

! other variables

  real(REAL64)                  ::  a, b, c
  real(REAL64)                  ::  alfa, beta, gama
  real(REAL64)                  ::  a1, c1, c2, c3

  integer                       ::  ibravais

  character(len = 1)            ::  yesno
  real(REAL64)                  ::  sgn

! constants

  real(REAL64), parameter  ::  ZERO = 0.0_REAL64, UM = 1.0_REAL64
  real(REAL64), parameter  ::  PI = 3.14159265358979323846_REAL64

! counter

  integer               ::  loop, i, j

! initialization  you have 20 tries


  do loop = 1,20

!   initialization

    do i = 1,3
    do j = 1,3
      vec(j,i) = ZERO
      conv(j,i) = ZERO
      avec(j,i) = ZERO
      iop(j,i) = ' '
    enddo
    enddo

    write(6,*) 'what is the bravais lattice?'
    write(6,*) ' 1 = simple cubic'
    write(6,*) ' 2 = body centered cubic'
    write(6,*) ' 3 = face centerd cubic'
    write(6,*) ' 4 = simple tetragonal'
    write(6,*) ' 5 = centered tetragonal '
    write(6,*) ' 6 = orthorhombic'
    write(6,*) ' 7 = base centered orthorhombic'
    write(6,*) ' 8 = body centered orthorhombic'
    write(6,*) ' 9 = face centered orthorhombic'
    write(6,*) ' 10 = monoclinic'
    write(6,*) ' 11 = centered monoclinic'
    write(6,*) ' 12 = triclinic'
    write(6,*) ' 13 = trigonal'
    write(6,*) ' 14 = hexagonal'

    read(5,*) ibravais
    write(ioreplay,'(i5,5x,"ibravais")') ibravais

    if(ibravais < 1 .OR. ibravais > 14) then
      ibravais= 12
      write(6,*) '  bad choice, assume triclinic:-)'
    endif

    if(ibravais == 1) then

      write(6,*) '  simple cubic '

      write(6,*) 'enter the lattice constant a'
      read(5,*) a
      write(ioreplay,'(1x,f9.4,5x,"latt")') a

      alatt = a
      do i=1,3
        vec(i,i) = UM
        conv(i,i) = a
      enddo

    elseif(ibravais == 2) then

      write(6,*) ' body centered cubic '

      write(6,*) ' enter the lattice constant a'
      write(ioreplay,'(1x,f9.4,5x,"latt")') a
      read(5,*) a

      alatt = a
      do i=1,3
      do j=1,3
        vec(j,i) = UM / 2
      enddo
      enddo

      do i=1,3
        vec(i,i) = -UM / 2
        conv(i,i) = a
      enddo

    elseif(ibravais == 3) then

      write(6,*) ' face centered cubic '

      write(6,*) ' enter the lattice constant a'
      read(5,*) a
      write(ioreplay,'(1x,f9.4,5x,"latt")') a

      alatt = a
      do i=1,3
      do j=1,3
        vec(j,i) = UM / 2
      enddo
      enddo

      do i=1,3
        vec(i,i) = ZERO
        conv(i,i) = a
      enddo

    elseif(ibravais == 4) then

      write(6,*) '  simple tetragonal  '

      write(6,*) ' enter the lattice constants a,c'
      write(6,*) ' c is the tetragonal axis '
      read(5,*) a, c
      write(ioreplay,'(1x,2f9.4,5x,"latt")') a, c

      alatt = a

      vec(1,1) = UM
      vec(2,2) = UM
      vec(3,3) = c/a

      conv(1,1) = a
      conv(2,2) = a
      conv(3,3) = c

    elseif(ibravais == 5) then

      write(6,*) '   body centered tetragonal  '

      write(6,*) ' enter the lattice constants a,c'
      write(6,*) ' c is the tetragonal axis '
      read(5,*) a, c
      write(ioreplay,'(1x,2f9.4,5x,"latt")') a, c

      alatt = a

      vec(1,1) = UM
      vec(2,2) = UM
      vec(1,3) = UM / 2
      vec(2,3) = UM / 2
      vec(3,3) = (UM/2)*c/a

      conv(1,1) = a
      conv(2,2) = a
      conv(3,3) = c

    elseif(ibravais == 6) then

      write(6,*) '   simple orthorhombic  '

      write(6,*) ' enter the lattice constants a,b,c'
      read(5,*) a, b, c
      write(ioreplay,'(1x,3f9.4,5x,"latt")') a, b, c

      alatt = a

      vec(1,1) = UM
      vec(2,2) = b/a
      vec(3,3) = c/a

      conv(1,1) = a
      conv(2,2) = b
      conv(3,3) = c

    elseif(ibravais == 7) then

      write(6,*) '  base centered orthorhombic'

      write(6,*) ' enter the lattice constants a,b,c'
      write(6,*) ' center is in the a-b plane, normal to c'
      read(5,*) a, b, c
      write(ioreplay,'(1x,3f9.4,5x,"latt")') a, b, c

      alatt = a

      vec(1,1) = UM / 2
      vec(2,1) = -(UM/2)*b/a
      vec(1,2) = UM / 2
      vec(2,2) = (UM/2)*b/a
      vec(3,3) = c/a

      conv(1,1) = a
      conv(2,2) = b
      conv(3,3) = c

    elseif(ibravais == 8) then

      write(6,*) '  body centered orthorhombic'

      write(6,*) ' enter the lattice constants a,b,c'
      read(5,*) a, b, c
      write(ioreplay,'(1x,3f9.4,5x,"latt")') a, b, c

      alatt = a

      vec(1,1) = UM
      vec(2,2) = b/a
      vec(1,3) = UM / 2
      vec(2,3) = (UM/2)*b/a
      vec(3,3) = (UM/2)*c/a

      conv(1,1) = a
      conv(2,2) = b
      conv(3,3) = c

    elseif(ibravais == 9) then

      write(6,*) '  face centered orthorhombic'

      write(6,*) ' enter the lattice constants a,b,c'
      read(5,*) a, b, c
      write(ioreplay,'(1x,3f9.4,5x,"latt")') a, b, c

      alatt = a

      vec(2,1) = (UM/2)*b/a
      vec(3,1) = (UM/2)*c/a
      vec(1,2) = UM / 2
      vec(3,2) = (UM/2)*c/a
      vec(1,3) = UM / 2
      vec(2,3) = (UM/2)*b/a

      conv(1,1) = a
      conv(2,2) = b
      conv(3,3) = c

    elseif(ibravais == 10) then

      write(6,*) '  monoclinic'
      write(6,*) ' enter the lattice constants a,b,c'
      read(5,*) a, b, c
      write(ioreplay,'(1x,3f9.4,5x,"latt")') a, b, c

      write(6,*) ' enter angle beta in degrees'
      write(6,*) '   beta is between a and c'
      write(6,*) '   (b is perpendicular to a and c)'
      read(5,*) beta
      write(ioreplay,'(1x,f9.4,5x,"beta")') beta
      beta = PI*beta / 180

      alatt = UM

      vec(1,1) = a
      vec(2,2) = b
      vec(1,3) = cos(beta)*c
      vec(3,3) = sin(beta)*c

      conv(1,1) = a
      conv(2,2) = b
      conv(1,3) = cos(beta)*c
      conv(3,3) = sin(beta)*c

    elseif(ibravais == 11) then

      write(6,*) ' 11 = centered monoclinic'

      write(6,*) ' enter the lattice constants a,b,c'
      read(5,*) a, b, c
      write(ioreplay,'(1x,3f9.4,5x,"latt")') a, b, c

      write(6,*) ' enter angle beta in degrees'
      write(6,*) '   beta is between a and c'
      write(6,*) '   (b is perpendicular to a and c)'
      write(6,*) '   center is in the a-b plane'
      read(5,*) beta
      write(ioreplay,'(1x,f9.4,5x,"beta")') beta
      beta = PI*beta / 180

      alatt = UM

      vec(1,1) = (UM/2)*a
      vec(2,1) = -(UM/2)*b
      vec(1,2) = (UM/2)*a
      vec(2,2) = (UM/2)*b
      vec(1,3) = cos(beta)*c
      vec(3,3) = sin(beta)*c

      conv(1,1) = a
      conv(2,2) = b
      conv(1,3) = cos(beta)*c
      conv(3,3) = sin(beta)*c

    elseif(ibravais == 12) then

      write(6,*) '  triclinic'
      write(6,*) ' enter the lattice constants a,b,c'
      read(5,*) a,b,c
      write(ioreplay,'(1x,3f9.4,5x,"latt")') a, b, c

      write(6,*) ' enter angles alpha,beta,gamma in degrees'
      write(6,*) ' alpha is between b and c, etc...'
      read(5,*) alfa,beta,gama
      write(ioreplay,'(1x,3f9.4,5x,"angles")') alfa, beta, gama
      alfa = PI*alfa / 180
      beta = PI*beta / 180
      gama = PI*gama / 180
      c1 = cos(beta)
      c2 = ((cos(alfa)-cos(gama)*cos(beta))/sin(gama))
      c3 = sqrt(UM-c1*c1-c2*c2)

      alatt = UM

      vec(1,1) = a
      vec(1,2) = cos(gama)*b
      vec(2,2) = sin(gama)*b
      vec(1,3) = c1*c
      vec(2,3) = c2*c
      vec(3,3) = c3*c

      conv(1,1) = a
      conv(1,2) = cos(gama)*b
      conv(2,2) = sin(gama)*b
      conv(1,3) = c1*c
      conv(2,3) = c2*c
      conv(3,3) = c3*c

    elseif(ibravais == 13) then

      write(6,*) '   trigonal'

      write(6,*) ' do you want to use pseudo-hexagonal cell? (y/n)'
      read(5,*) yesno
      write(ioreplay,'(1x,a1,5x,"pseudo-hexagonal")') yesno
      if(yesno == 'y') then
        write(6,*) ' enter the lattice constants a,c'
        read(5,*) a,c
        write(ioreplay,'(1x,2f9.4,5x,"latt")') a, c

        alatt = a

        vec(1,1) = UM
        vec(1,2) = UM/2
        vec(2,2) = (3*UM)/(4*UM)
        vec(3,2) = ZERO
        vec(1,3) = UM/2
        vec(2,3) = UM/2
        vec(3,3) = c/((3*UM)*a)

        iop(2,2) = 'S'
        iop(2,3) = 'H'

        conv(1,1) = a/2
        conv(2,1) =-sqrt((3*UM)/(4*UM))*a
        conv(3,1) = ZERO
        conv(1,2) = a/2
        conv(2,2) = sqrt((3*UM)/(4*UM))*a
        conv(3,2) = ZERO
        conv(1,3) = ZERO
        conv(2,3) = ZERO
        conv(3,3) = c

      else

        write(6,*) ' enter the lattice constant a'
        read(5,*) a
        write(ioreplay,'(1x,f9.4,5x,"latt")') a

        write(6,*) ' enter angle alpha in degrees'
        read(5,*) alfa
        write(ioreplay,'(1x,f9.4,5x,"angle")') alfa
        alfa = PI*alfa / 180

        a1 = sqrt((2*UM)*(UM-cos(alfa))/(3*UM))
        c1 = sqrt((UM+(2*UM)*cos(alfa))/(3*UM))

        alatt = a*a1

        vec(1,1) = UM
        vec(2,1) = ZERO
        vec(3,1) = c1/a1
        vec(1,2) =-UM / 2
        vec(2,2) = (3*UM)/(4*UM)
        vec(3,2) = c1/a1
        vec(1,3) =-UM / 2
        vec(2,3) =-(3*UM)/(4*UM)
        vec(3,3) = c1/a1

        iop(2,2) = 'S'
        iop(2,3) = 'S'

        conv(1,1) = a1*a
        conv(2,1) = ZERO
        conv(3,1) = c1*a
        conv(1,2) =-(UM/2)*a1*a
        conv(2,2) = dsqrt(((3*UM)/(4*UM))*a1*a1)*a
        conv(3,2) = c1*a
        conv(1,3) =-(UM/2)*a1*a
        conv(2,3) =-dsqrt(((3*UM)/(4*UM))*a1*a1)*a
        conv(3,3) = c1*a

      endif

    elseif(ibravais == 14) then

      write(6,*) '    hexagonal'

      write(6,*) ' enter the lattice constants a,c'
      read(5,*) a,c
      write(ioreplay,'(1x,2f9.4,5x,"latt")') a, c

      alatt = a

      vec(1,1) = UM/2
      vec(2,1) =-(3*UM)/(4*UM)
      vec(1,2) = UM/2
      vec(2,2) = (3*UM)/(4*UM)
      vec(3,3) = c/a

      iop(2,1) = 'S'
      iop(2,2) = 'S'

      conv(1,1) = a/2
      conv(2,1) =-sqrt((3*UM)/(4*UM))*a
      conv(1,2) = a/2
      conv(2,2) = sqrt((3*UM)/(4*UM))*a
      conv(3,3) = c
    endif

!   now the real lattice vectors

    do j = 1,3
    do i = 1,3
      sgn = UM
      if(vec(i,j) < ZERO) sgn = -sgn
      avec(i,j) = abs(vec(i,j))
      if(iop(i,j) == 'S') avec(i,j) = sqrt(avec(i,j))
      if(iop(i,j) == 'H') avec(i,j) = avec(i,j)/sqrt(3*UM)
      avec(i,j) = sgn*alatt*avec(i,j)
    enddo
    enddo

    write(6,*)
    write(6,*) '  The following lines will be added to PW.DAT'
    write(6,*)

    write(6,'(f12.6)') alatt
    do i=1,3
      write(6,'(3(2x,a1,f12.6))') (iop(j,i),vec(j,i),j=1,3)
    enddo

    write(6,*)
    write(6,*) '  The primitive lattice vectors are'
    write(6,*)

    do i=1,3
      write(6,'(3(3x,f12.6))') (avec(j,i),j=1,3)
    enddo

    write(6,*)
    write(6,*) '  The conventional lattice vectors are'
    write(6,*)

    do i=1,3
      write(6,'(3(3x,f12.6))') (conv(j,i),j=1,3)
    enddo

    write(6,*)
    write(6,*) ' are you happy? (y/n)'
    read(5,*) yesno
    write(ioreplay,'(1x,a1,5x,"zen question")') yesno

    if(yesno == 'y' .OR. yesno == 'Y') exit

  enddo

  return
end subroutine gen_latt




!  writes the cpw.in file with final geometry

subroutine write_cpwin(io, title, source, alatt,                         &
    avec, ntype, natom, nameat, rat,                                     &
    emax, nband, nx,ny,nz, sx,sy,sz,                                     &
    mxdtyp, mxdatm)

! Adapted June 2017. JLM
! Bug squashed (metadata not from rede) September 2017.
! Adapted from the subroutine write_cpwout
! Copyright  J.L.Martins, INESC-MN.

! version 4.94

  implicit none

  integer, parameter  :: REAL64 = selected_real_kind(12)

! input:

  integer, intent(in)                ::  mxdtyp                     !<  array dimension of types of atoms
  integer, intent(in)                ::  mxdatm                     !<  array dimension of types of atoms

  integer, intent(in)                ::  io                         !  tape number

  character(len=20), intent(in)      ::  title                      !<  title cpw_in or PW.DAT
  character(len=*), intent(in)       ::  source

  real(REAL64), intent(in)           ::  alatt                      !<  lattice constant

  real(REAL64), intent(in)           ::  avec(3,3)                  !<  lattice vectors
  integer, intent(in)                ::  ntype                      !<  number of types of atoms
  integer, intent(in)                ::  natom(mxdtyp)              !<  number of atoms of type i
  character(len=2), intent(in)       ::  nameat(mxdtyp)             !<  chemical symbol for the type i
  real(REAL64), intent(in)           ::  rat(3,mxdatm,mxdtyp)       !<  k-th component (in lattice coordinates) of the position of the n-th atom of type i

  real(REAL64), intent(in)           ::  emax                       !<  kinetic energy cutoff of plane wave expansion (Hartree).
  integer, intent(in)                ::  nband                      !<  target for number of bands

  integer, intent(in)                ::  nx, ny, nz                 !<  size of the integration mesh in k-space (nx*ny*nz)
  real(REAL64), intent(in)           ::  sx, sy, sz                 !<  offset of the integration mesh (usually 0.5)

! local:

  real(REAL64)       ::  bvec(3,3)

  integer            ::  nlatpl                                     ! number of lattice planes
  character(len=3)   ::  vers                                       ! program version
  character(len=140) ::  metadata                                   ! metadata of the calculation

  integer  ::  isl1(3),isl2(3),isl3(3)
  integer  ::  nat
  integer  ::  istart, iend, ioerr


! counters

  integer       ::  nt, n, i, j


! open file


  write(io,'(72a1)') ("#",j=1,72)
  write(io,'("#",70x,"#")')
  write(io,'("#",6x,"cpw.in input file generated by gen_PW",27x,"#")')
  write(io,'("#",70x,"#")')
  write(io,'(72a1)') ("#",j=1,72)


  write(io,*)
  write(io,'("SystemLabel",19x,a20)') title
  write(io,*)

  write(io,*)
  write(io,'("#------------------------------------------------")')
  write(io,'("# Crystal structure")')
  write(io,'("#------------------------------------------------")')
  write(io,*)
  write(io,'("LatticeConstant",10x,f16.8,5x,"bohr")') alatt
  write(io,*)

  write(io,'("%block LatticeVectors")')
  write(io,'(3(3x,f16.8))') avec(1,1)/alatt,avec(2,1)/alatt,avec(3,1)/alatt
  write(io,'(3(3x,f16.8))') avec(1,2)/alatt,avec(2,2)/alatt,avec(3,2)/alatt
  write(io,'(3(3x,f16.8))') avec(1,3)/alatt,avec(2,3)/alatt,avec(3,3)/alatt
  write(io,'("%endblock LatticeVectors")')
  write(io,*)

! Finds total number of atoms

  nat = 0
  do i=1,ntype
    nat = nat + natom(i)
  enddo

  write(io,'("NumberOfSpecies",9x,i8)') ntype
  write(io,*)

  write(io,'("NumberOfAtoms",9x,i8)') nat
  write(io,*)

  write(io,'("%block Chemical_Species_Label")')
  do i = 1,ntype
    call p_tbl_charge(nameat(i),n)
    write(io,'(i6,3x,i5,3x,a2)') i,n,nameat(i)
  enddo
  write(io,'("%endblock Chemical_Species_Label")')
  write(io,*)

  write(io,'("AtomicCoordinatesFormat",5x,"Fractional")')
  write(io,*)


  write(io,'("%block AtomicCoordinatesAndAtomicSpecies")')

  do nt = 1,ntype
  do i = 1,natom(nt)
     write(io,'(3x,3f16.8,4x,i5,5x,"#  ",a2,3x,i5)')                 &
           (rat(j,i,nt),j=1,3),nt,nameat(nt),i
  enddo
  enddo
  write(io,'("%endblock AtomicCoordinatesAndAtomicSpecies")')
  write(io,*)

  write(io,'("StructureSource",15x,200a1)') (source(i:i),i=1,len_trim(source))
  write(io,*)

  write(io,'("#------------------------------------------------")')
  write(io,'("# Energy cutoff, bands,  and Brillouin mesh")')
  write(io,'("#------------------------------------------------")')
  write(io,*)

  write(io,'("PWEnergyCutoff",13x,f12.4,6x,"hartree")')  emax
  write(io,*)

  write(io,'("NumberOfEigenStates",12x,i8)') nband
  write(io,*)

  write(io,'("%block kgrid_Monkhorst_Pack")')
  write(io,'(8x,3(2x,i6),3x,f12.6)')  nx,0,0,sx
  write(io,'(8x,3(2x,i6),3x,f12.6)')  0,ny,0,sy
  write(io,'(8x,3(2x,i6),3x,f12.6)')  0,0,nz,sz
  write(io,'("%endblock kgrid_Monkhorst_Pack")')
  write(io,*)


  write(io,*)
  write(io,'("#------------------------------------------------")')
  write(io,'("# Active options")')
  write(io,'("#------------------------------------------------")')
  write(io,*)

  write(io,'("MD.TypeOfRun                  ONE           ",             &
       "# ONE,EPILBF,MICRO,LANG,LBFSYM,VCSLNG,VCSLBF,RSTRT")')
  write(io,*)

  write(io,'("UseSymmetry                   .true.        ",             &
       "# .true. , .false. ")')
  write(io,*)

  write(io,'("MD.UseKeatingCorrections      .false.       ",             &
       "# .true. , .false. ")')
  write(io,*)

  write(io,'("MD.UseFixedkplusG             .true.        ",             &
       "# .true. , .false. ")')
  write(io,*)


  write(io,'("TypeOfScfDiag                 PW            ",             &
       "# PW,AO,AOJC,AOJCPW")')
  write(io,*)

  write(io,'("DualApproximation             .true.        ",             &
       "#  .true. , .false.")')
  write(io,*)

  write(io,'("XC.Authors                    CA            ",             &
       "# CA, PBE, TBL")')
  write(io,*)

  write(io,'("Xc.TBL.C                      1.04          ",             &
       "# sets Tran-Blaha constant (if negative use calculated)")')
  write(io,*)

  write(io,'("PrintingLevel                 2              ",            &
       "# 1, 2, 3")')
  write(io,*)

  write(io,*)
  write(io,'("#------------------------------------------------")')
  write(io,'("# MD Inactive options")')
  write(io,'("#------------------------------------------------")')
  write(io,*)

  write(io,'("#MD.InitialTemperature        300 K           #")')
  write(io,'("#MD.TargetTemperature         300 K           #")')
  write(io,'("#MD.TargetPressure            0 GPa           #")')
  write(io,*)
  write(io,'("#MD.NumberOfSteps             10              #")')
  write(io,'("#MD.LengthTimeStep            2.4 fs          #")')
  write(io,'("#MD.FrictionFracInvTimeStep   20.0            #")')
  write(io,*)
  write(io,'("#MD.CG.Tolerance         0.0001 ''har/bohr''    #")')
  write(io,'("#MD.CG.StepMax                0.01 bohr       #")')
  write(io,'("#MD.CG.FixedkplusGTol    0.01 ''har/bohr''      #")')
  write(io,*)
  write(io,'("#%block MD.TargetStress                       #")')
  write(io,'("   0.0 0.0 0.0                                #")')
  write(io,'("   0.0 0.0 0.0                                #")')
  write(io,'("   0.0 0.0 0.0                                #")')
  write(io,'("#%endblock   MD.TargetStress                  #")')
  write(io,*)
  write(io,'("#MD.CellMass                  10.0            #")')
  write(io,'("#MD.Seed                      76978           #")')

  write(io,*)
  write(io,'("#------------------------------------------------")')
  write(io,'("# Electronic Structure Inactive options")')
  write(io,'("#------------------------------------------------")')
  write(io,*)

  write(io,'("#MaxSCFIterations             20              #")')
  write(io,*)
  write(io,'("#MaxSCFIterations             20              #")')
  write(io,'("#TypeOfPseudoMixing           BROYD1          #",     &
       " BROYD1, BFGS#")')
  write(io,*)
  write(io,'("#ElectronicTemperature        1000 K          #")')
  write(io,'("#TypeOfPseudopotential        PSEUKB          #",     &
       " PSEUKB")')
  write(io,*)
  write(io,'("#ScfTolerance                 0.00005         #")')
  write(io,'("#DiagTolerance                0.0001          #")')
  write(io,'("#SymmTolerance                1.0E-5          #")')

  return
end subroutine write_cpwin



!>Function determines the nuclear charge of an element.
!>All elements from H to Lr are included.

subroutine p_tbl_charge(name,iz)

! Version dated May 1, 1991 written by  njtj
! Modified documentation, August 2019. JLM
! Modified for integer. July 2013
! copyright inesc-mn/Jose Luis Martins/Ivo Souza

! Version 4.94 of cpw

  implicit none

! input:

  character(len=2),intent(in)        :: name                        !<  element symbol

! output:

  integer,intent(out)                :: iz                          !<  nuclear charge (atomic number)

  if (name == 'H ' .or. name == ' H') then
    iz = 1
  elseif (name == 'He') then
    iz = 2
  elseif (name == 'Li') then
    iz = 3
  elseif (name == 'Be') then
    iz = 4
  elseif (name == 'B ' .or. name == ' B') then
    iz = 5
  elseif (name == 'C ' .or. name == ' C') then
    iz = 6
  elseif (name == 'N ' .or. name == ' N') then
    iz = 7
  elseif (name == 'O ' .or. name == ' O') then
    iz = 8
  elseif (name == 'F ' .or. name == ' F') then
    iz = 9
  elseif (name == 'Ne') then
    iz = 10
  elseif (name == 'Na') then
    iz = 11
  elseif (name == 'Mg') then
    iz = 12
  elseif (name == 'Al') then
    iz = 13
  elseif (name == 'Si') then
    iz = 14
  elseif (name == 'P ' .or. name == ' P') then
    iz = 15
  elseif (name == 'S ' .or. name == ' S') then
    iz = 16
  elseif (name == 'Cl') then
    iz = 17
  elseif (name == 'Ar') then
    iz = 18
  elseif (name == 'K ' .or. name == ' K') then
    iz = 19
  elseif (name == 'Ca') then
    iz = 20
  elseif (name == 'Sc') then
    iz = 21
  elseif (name == 'Ti') then
    iz = 22
  elseif (name == 'V ' .or. name == ' V') then
    iz = 23
  elseif (name == 'Cr') then
    iz = 24
  elseif (name == 'Mn') then
    iz = 25
  elseif (name == 'Fe') then
    iz = 26
  elseif (name == 'Co') then
    iz = 27
  elseif (name == 'Ni') then
    iz = 28
  elseif (name == 'Cu') then
    iz = 29
  elseif (name == 'Zn') then
    iz = 30
  elseif (name == 'Ga') then
    iz = 31
  elseif (name == 'Ge') then
    iz = 32
  elseif (name == 'As') then
    iz = 33
  elseif (name == 'Se') then
    iz = 34
  elseif (name == 'Br') then
    iz = 35
  elseif (name == 'Kr') then
    iz = 36
  elseif (name == 'Rb') then
    iz = 37
  elseif (name == 'Sr') then
    iz = 38
  elseif (name == 'Y ' .or. name == ' Y') then
    iz = 39
  elseif (name == 'Zr') then
    iz = 40
  elseif (name == 'Nb') then
    iz = 41
  elseif (name == 'Mo') then
    iz = 42
  elseif (name == 'Tc') then
    iz = 43
  elseif (name == 'Ru') then
    iz = 44
  elseif (name == 'Rh') then
    iz = 45
  elseif (name == 'Pd') then
    iz = 46
  elseif (name == 'Ag') then
    iz = 47
  elseif (name == 'Cd') then
    iz = 48
  elseif (name == 'In') then
    iz = 49
  elseif (name == 'Sn') then
    iz = 50
  elseif (name == 'Sb') then
    iz = 51
  elseif (name == 'Te') then
    iz = 52
  elseif (name == 'I ' .or. name == ' I') then
    iz = 53
  elseif (name == 'Xe') then
    iz = 54
  elseif (name == 'Cs') then
    iz = 55
  elseif (name == 'Ba') then
    iz = 56
  elseif (name == 'La') then
    iz = 57
  elseif (name == 'Ce') then
    iz = 58
  elseif (name == 'Pr') then
    iz = 59
  elseif (name == 'Nd') then
    iz = 60
  elseif (name == 'Pm') then
    iz = 61
  elseif (name == 'Sm') then
    iz = 62
  elseif (name == 'Eu') then
    iz = 63
  elseif (name == 'Gd') then
    iz = 64
  elseif (name == 'Tb') then
    iz = 65
  elseif (name == 'Dy') then
    iz = 66
  elseif (name == 'Ho') then
    iz = 67
  elseif (name == 'Er') then
    iz = 68
  elseif (name == 'Tm') then
    iz = 69
  elseif (name == 'Yb') then
    iz = 70
  elseif (name == 'Lu') then
    iz = 71
  elseif (name == 'Hf') then
    iz = 72
  elseif (name == 'Ta') then
    iz = 73
  elseif (name == 'W ' .or. name == ' W') then
    iz = 74
  elseif (name == 'Re') then
    iz = 75
  elseif (name == 'Os') then
    iz = 76
  elseif (name == 'Ir') then
    iz = 77
  elseif (name == 'Pt') then
    iz = 78
  elseif (name == 'Au') then
    iz = 79
  elseif (name == 'Hg') then
    iz = 80
  elseif (name == 'Tl') then
    iz = 81
  elseif (name == 'Pb') then
    iz = 82
  elseif (name == 'Bi') then
    iz = 83
  elseif (name == 'Po') then
    iz = 84
  elseif (name == 'At') then
    iz = 85
  elseif (name == 'Rn') then
    iz = 86
  elseif (name == 'Fr') then
    iz = 87
  elseif (name == 'Ra') then
    iz = 88
  elseif (name == 'Ac') then
    iz = 89
  elseif (name == 'Th') then
    iz = 90
  elseif (name == 'Pa') then
    iz = 91
  elseif (name == ' U' .or. name == 'U ') then
    iz = 92
  elseif (name == 'Np') then
    iz = 93
  elseif (name == 'Pu') then
    iz = 94
  elseif (name == 'Am') then
    iz = 95
  elseif (name == 'Cm') then
    iz = 96
  elseif (name == 'Bk') then
    iz = 97
  elseif (name == 'Cf') then
    iz = 98
  elseif (name == 'Es') then
    iz = 99
  elseif (name == 'Fm') then
    iz = 100
  elseif (name == 'Md') then
    iz = 101
  elseif (name == 'No') then
    iz = 102
  elseif (name == 'Lr') then
    iz = 103
  elseif (name == 'Rf') then
    iz = 104
  elseif (name == 'Db') then
    iz = 105
  elseif (name == 'Sg') then
    iz = 106
  elseif (name == 'Bh') then
    iz = 107
  elseif (name == 'Hs') then
    iz = 108
  elseif (name == 'Mt') then
    iz = 109
  elseif (name == 'Ds') then
    iz = 110
  elseif (name == 'Rg') then
    iz = 111
  elseif (name == 'Cn') then
    iz = 112
  elseif (name == 'Nh') then
    iz = 113
  elseif (name == 'Fl') then
    iz = 114
  elseif (name == 'Mc') then
    iz = 115
  elseif (name == 'Lv') then
    iz = 116
  elseif (name == 'Ts') then
    iz = 117
  elseif (name == 'Og') then
    iz = 118
  else
    write(6,*)
    write(6,*)
    write(6,'("  Element ",a2," unknown")') name
    write(6,'("  Using charge 200")')
    write(6,*)
    iz = 200
  endif

  return
end subroutine p_tbl_charge
