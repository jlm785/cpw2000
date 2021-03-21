!  
!  E l e c t r o n i c   S t r u c t u r e   D a t a   F o r m a t
!  ---------------------------------------------------------------
!  
!                             E S D F
!                             =======
!  
!  Author: Chris J. Pickard (c)
!  Email : cp@min.uni-kiel.de
!  Place : Kiel, Germany
!  Date  : 5/6th August 1999
!  
!  Summary
!  -------
!  
!  This module is designed to simplify and enhance the input of data into
!  electronic structure codes (for example, CASTEP). It works from a
!  highly flexible input file. Data is input in a "label <value>"
!  fashion, and input is independent of the ordering of the input
!  file. An important feature is the requirement that most inputs require
!  default settings to be supplied within the main program calling
!  ESDF. This means that rarely used variables will not clutter everyday
!  input files, and, even more usefully, "intelligence" may be built into
!  the main code as the defaults may be dependent of other set
!  variables. Block data may also be read in. Another important feature
!  is the ability to define "physical" values. This means that the input
!  files need not depend on the internal physical units used by the main
!  program.
!  
!  
!  History
!  -------
!  
!  ESDF has been written from scratch in F90, but is heavily based
!  (especially for the concept) on the FDF package developed by Alberto
!  Garcia and Jose Soler. It is not as "flexible" as FDF - there is no
!  provision for branching to further input files. This simplifies the
!  code, and I hope that it is still useful without this feature. Also,
!  the input and defaults are not dumped to a output file currently. I've
!  not found this a hindrance as of now.
!  
!  
!  Future
!  ------ 
!  
!  My intention is to make this release available to Alberto Garcia and
!  Jose Soler for their comments. It might be a good idea to use this as
!  a base for fully converting the FDF package to F90. Or it may remain
!  as a cut down version of FDF. I certainly hope that a package of the
!  FDF sort becomes widely used in the electronic structure community. My
!  experience has been very positive.
!  
!  
!  Usage
!  -----
!  
!  First, "Use esdf" wherever you wish to make use of its features. In
!  the main program call the initialisation routine: Call
!  esdf_init('input.esdf'). "input.esdf" is the name of the input file -
!  it could be anything. This routine opens the input file, and reads
!  into a dynamically allocated storage array. The comments and blank
!  lines are stripped out. You are now ready to use the
!  esdf_functions. For example, if you want to read in the number of
!  atoms in your calculation, you would use: natom =
!  esdf_integer('NumberOfAtoms',1), where 'NumberOfAtoms' is the label to
!  search for in the input file, and '1' is the default. Call esdf_close to
!  deallocate the data arrays. You may then open another input file using
!  esdf_init. It is not currently possible to open more that on input 
!  file at one time.
!  
!  
!  Syntax
!  ------
!  
!  The input file can contain comments. These are defined as anything to
!  the right of, and including, '#', ';', or '!'. It is straightforward
!  to modify the routine to accept further characters. Blank lines are
!  ignored -- use comments and blank lines to make you input file
!  readable.
!  
!  The "labels" are case insensitive (e.g. unitCell is equivalent to
!  UnItceLL) and punctuation insensitive (unit.cell is equivalent to
!  unit_cell is equivalent to unitcell). Punctuation characters are '.',
!  '_', and '-' at the moment. Again - use this feature to improve
!  readability.
!  
!  The following are equivalent ways of defining a physical quantity:
!  
!  "AgeOfUniverse = 24.d0 s" or "AgeOfUniverse : 24.d0 S" or
!  "AgeOfUniverse 24.d0 S"
!   
!  It would be read in by the main program in the following way:
!  
!  aou = esdf_physical('ageofuniverse',77.d0,ns)
!  
!  "aou" is the double precision variable, 77.d0 is the default number of
!  "ns" or nanoseconds. 24s will be converted automatically to its
!  equivalent number of nanoseconds.
!  
!  Block data should be placed in the input file as follows:
!  
!  %block cellvectors 
!  1.0 1.0 0.0 
!  0.0 1.0 1.0 
!  1.0 0.0 1.0 
!  %endblock cellvectors
!  
!  And it may be read:
!  
!    if(esdf_block('CellVectors',nlines))
!      if(nlines.ne.3) then (... break out here if the incorrect number
!  of lines)
!      do i=1,nlines
!        read(block_data(i),*) x,y,z
!      end do
!    endif
!  
!  
!  List of functions
!  -----------------
!  
!  Self explanatory:
!  
!  esdf_string(label,default)
!  esdf_integer(label,default)
!  esdf_single(label,default)
!  esdf_double(label,default)
!  esdf_physical(label,default,unit)
!  
!  A little more explanation:
!  
!  esdf_defined(label) is true if "label" found, false otherwise
!  
!  esdf_boolean(label,default) is true if "label yes/true/t (case/punct.insens)
!                              is false if"label no/false/f (case/punct.insens)
!  
!  The Help feature
!  ----------------
!  
!  The routine "esdf_help(helpword,searchword)" can be used to access the
!  information contained within the "esdf_key_mod" module.
!  
!  If "helpword" is "search" (case insensitive), then all variables whose
!  description contains "searchword" will be output.
!  
!  If "helpword" is "basic", "inter", "expert" or "dummy" the varibles of
!  that type will be displayed.
!  
!  If "helpword" is one of the valid labels, then a description of this
!  label will be output.
!  
!  
!  Finishing off
!  -------------
!  
!  Three routines, "esdf_dump", "esdf_warnout" and "esdf_close", can be
!  used to finish the use of ESDF. "esdf_dump" outputs a file "ESDF.esdf"
!  which could be used as an input file for further runs. "esdf_warnout"
!  outputs ESDF warnings to screen, and "esdf_close" deallocates the
!  allocated ESDF arrays.
!  
!  Contact the Author
!  ------------------
!  
!  This code is under development, and the author would be very happy to
!  receive comments by email. Any use in a commercial software package is 
!  forbidden without prior arrangement with the author (Chris J. Pickard).

Module esdf

  Use esdf_key

  Implicit None

  ! Kind parameters 

  Integer, Private, Parameter :: I4B = Selected_int_kind(9)
  Integer, Private, Parameter :: DP  = Kind(1.d0)
  Integer, Private, Parameter :: SP  = Kind(1.0)

  ! Set the length of the lines

  Integer(I4B), Public, Parameter :: llength=80
  Integer(I4B), Private, Parameter ::  nphys = 56, ndump = 2000
  Integer(I4B), Private :: nrecords,nwarns,ndmp
  Character(llength), Private, Dimension(:), Allocatable :: llist,warns,dump 
  Character(llength), Private, Dimension(:,:), Allocatable :: tlist 

  ! The public block data array

  Character(llength), Public, Dimension(:), Allocatable :: block_data 

  ! Set the physical units database

  Type phys_unit
     Character(10) :: d,n ! d - dimension n - name
     Real(DP)      :: u   ! u - unit
  End Type phys_unit

  Type(phys_unit), Private, Dimension(nphys) :: phy

  !
  !     We allow case variations in the units. This could be dangerous
  !     (meV --> MeV !!) in real life, but not in this restricted 
  !     field.
  !     
  ! m - mass l - length t - time e - energy f - force p - pressure c- charge
  ! d - dipole mom - mom inert ef - efield
  !

  Data phy(1)%d /'m'/;Data phy(1)%n /'kg'/;Data phy(1)%u /1.d0/
  Data phy(2)%d /'m'/;Data phy(2)%n /'g'/;Data phy(2)%u /1.d-3/
  Data phy(3)%d /'m'/;Data phy(3)%n /'amu'/;Data phy(3)%u /1.66054d-27/
  Data phy(4)%d /'l'/;Data phy(4)%n /'m'/;Data phy(4)%u /1.d0/
  Data phy(5)%d /'l'/;Data phy(5)%n /'nm'/;Data phy(5)%u /1.d-9/
  Data phy(6)%d /'l'/;Data phy(6)%n /'ang'/;Data phy(6)%u /1.d-10/
  Data phy(7)%d /'l'/;Data phy(7)%n /'bohr'/;Data phy(7)%u /0.529177d-10/
  Data phy(8)%d /'t'/;Data phy(8)%n /'s'/;Data phy(8)%u /1.d0/
  Data phy(9)%d /'t'/;Data phy(9)%n /'ns'/;Data phy(9)%u /1.d-9/
  Data phy(10)%d /'t'/;Data phy(10)%n /'ps'/;Data phy(10)%u /1.d-12/
  Data phy(11)%d /'t'/;Data phy(11)%n /'fs'/;Data phy(11)%u /1.d-15/
  Data phy(12)%d /'e'/;Data phy(12)%n /'j'/;Data phy(12)%u /1.d0/
  Data phy(13)%d /'e'/;Data phy(13)%n /'erg'/;Data phy(13)%u /1.d-7/
  Data phy(14)%d /'e'/;Data phy(14)%n /'ev'/;Data phy(14)%u /1.60219d-19/
  Data phy(15)%d /'e'/;Data phy(15)%n /'mev'/;Data phy(15)%u /1.60219d-22/
  Data phy(16)%d /'e'/;Data phy(16)%n /'ry'/;Data phy(16)%u /2.17991d-18/
  Data phy(17)%d /'e'/;Data phy(17)%n /'mry'/;Data phy(17)%u /2.17991d-21/
  Data phy(18)%d /'e'/;Data phy(18)%n /'hartree'/;Data phy(18)%u /4.35982d-18/
  Data phy(19)%d /'e'/;Data phy(19)%n /'kcal/mol'/;Data phy(19)%u /6.94780d-21/
  Data phy(20)%d /'e'/;Data phy(20)%n /'mhartree'/;Data phy(20)%u /4.35982d-21/
  Data phy(21)%d /'e'/;Data phy(21)%n /'kj/mol'/;Data phy(21)%u /1.6606d-21/
  Data phy(22)%d /'e'/;Data phy(22)%n /'hz'/;Data phy(22)%u /6.6262d-34/
  Data phy(23)%d /'e'/;Data phy(23)%n /'thz'/;Data phy(23)%u /6.6262d-22/
  Data phy(24)%d /'e'/;Data phy(24)%n /'cm-1'/;Data phy(24)%u /1.986d-23/
  Data phy(25)%d /'e'/;Data phy(25)%n /'cm^-1'/;Data phy(25)%u /1.986d-23/
  Data phy(26)%d /'e'/;Data phy(26)%n /'cm**-1'/;Data phy(26)%u /1.986d-23/
  Data phy(27)%d /'f'/;Data phy(27)%n /'N'/;Data phy(27)%u /1.d0/
  Data phy(28)%d /'f'/;Data phy(28)%n /'ev/ang'/;Data phy(28)%u /1.60219d-9/
  Data phy(29)%d /'f'/;Data phy(29)%n /'ry/bohr'/;Data phy(29)%u /4.11943d-8/
  Data phy(30)%d /'l'/;Data phy(30)%n /'cm'/;Data phy(30)%u /1.d-2/
  Data phy(31)%d /'p'/;Data phy(31)%n /'pa'/;Data phy(31)%u /1.d0/
  Data phy(32)%d /'p'/;Data phy(32)%n /'mpa'/;Data phy(32)%u /1.d6/
  Data phy(33)%d /'p'/;Data phy(33)%n /'gpa'/;Data phy(33)%u /1.d9/
  Data phy(34)%d /'p'/;Data phy(34)%n /'atm'/;Data phy(34)%u /1.01325d5/
  Data phy(35)%d /'p'/;Data phy(35)%n /'bar'/;Data phy(35)%u /1.d5/
  Data phy(36)%d /'p'/;Data phy(36)%n /'mbar'/;Data phy(36)%u /1.d11/
  Data phy(37)%d /'p'/;Data phy(37)%n /'ry/bohr**3'/;Data phy(37)%u /1.47108d13/
  Data phy(38)%d /'p'/;Data phy(38)%n /'ev/ang**3'/;Data phy(38)%u /1.60219d11/
  Data phy(39)%d /'c'/;Data phy(39)%n /'c'/;Data phy(39)%u /1.d0/
  Data phy(40)%d /'c'/;Data phy(40)%n /'e'/;Data phy(40)%u /1.602177d-19/
  Data phy(41)%d /'d'/;Data phy(41)%n /'C*m'/;Data phy(41)%u /1.d0/
  Data phy(42)%d /'d'/;Data phy(42)%n /'D'/;Data phy(42)%u /3.33564d-30/
  Data phy(43)%d /'d'/;Data phy(43)%n /'debye'/;Data phy(43)%u /3.33564d-30/
  Data phy(44)%d /'d'/;Data phy(44)%n /'e*bohr'/;Data phy(44)%u /8.47835d-30/
  Data phy(45)%d /'d'/;Data phy(45)%n /'e*ang'/;Data phy(45)%u /1.602177d-29/
  Data phy(46)%d /'mom'/;Data phy(46)%n /'kg*m**2'/;Data phy(46)%u /1.d0/
  Data phy(47)%d /'mom'/;Data phy(47)%n /'ry*fs**2'/;Data phy(47)%u /2.1799d-48/
  Data phy(48)%d /'ef'/;Data phy(48)%n /'v/m'/;Data phy(48)%u /1.d0/
  Data phy(49)%d /'ef'/;Data phy(49)%n /'v/nm'/;Data phy(49)%u /1.d9/
  Data phy(50)%d /'ef'/;Data phy(50)%n /'v/ang'/;Data phy(50)%u /1.d10/
  Data phy(51)%d /'ef'/;Data phy(51)%n /'v/bohr'/;Data phy(51)%u /1.8897268d10/
  Data phy(52)%d /'ef'/;Data phy(52)%n /'ry/bohr/e'/;Data phy(52)%u /2.5711273d11/
  Data phy(53)%d /'ef'/;Data phy(53)%n /'har/bohr/e'/;Data phy(53)%u /5.1422546d11/
  Data phy(54)%d /'e'/;Data phy(54)%n /'k'/;Data phy(54)%u /1.38066d-23/
  Data phy(55)%d /'f'/;Data phy(55)%n /'har/bohr'/;Data phy(55)%u /8.23886d-8/
  Data phy(56)%d /'t'/;Data phy(56)%n /'autime'/;Data phy(56)%u /2.418884d-17/


Contains

  Subroutine esdf_init(filename)

    Character(*), Intent(in) :: filename

    ! Local

    Integer(I4B), Parameter :: ncomm=3,ndiv=3
    Integer(I4B) :: unit,ierr,i,j,ic,nt,ndef,nread
    Character(llength) :: cjunk,ctemp
    Character(1) :: comment(ncomm),divide(ndiv)
    Logical :: inblock

    ! Define comment characters

    Data comment /'#',';','!'/
    Data divide /' ','=',':'/

    ! "reduce" the keyword list for comparison

    Do i = 1,numkw
       ctemp = kw(i)%label
       kw(i)%label = esdf_reduce(ctemp)
    End Do

    ! Open the esdf file


    Call esdf_file(unit,filename,ierr)
    cjunk = 'Unable to open main input file "'//Trim(filename)//'"'

    If(ierr.Eq.1) Then
       Write(*,*) 'ESDF WARNING: '//Trim(cjunk)//' - using defaults'
       nread = 0
    Else
       nread = Huge(1)
    Endif

    ! Count the number of records (excluding blank lines and commented lines)

    nrecords = 0

    Do i=1,nread
       Read(unit,'(a)',End=100) cjunk
       Do j=1,ncomm
          ic=Index(cjunk,comment(j))
          If(ic.Gt.0) cjunk(ic:) = ' '
       End Do
       If(Len_trim(cjunk).Gt.0) Then
          nrecords = nrecords + 1
       Endif
    End Do
100 Rewind(unit)


    ! Allocate the array to hold the records and tokens

    Allocate(llist(nrecords),block_data(nrecords),tlist(llength,nrecords),&
         warns(nrecords),dump(ndump))

    ! Set the number of warnings to zero

    nwarns = 0 ; warns = ' ' ; ndmp = 0 ; dump = ' '

    ! Read in the records

    nrecords = 0
    Do i=1,nread
       Read(unit,'(a)',End=101) cjunk
       Do j=1,ncomm
          ic=Index(cjunk,comment(j))
          If(ic.Gt.0) cjunk(ic:) = ' '
       End Do

       If(Len_trim(cjunk).Gt.0) Then
          nrecords=nrecords+1
          llist(nrecords) = Adjustl(cjunk)
       Endif
    End Do
101 Close(unit)

    ! Now read in the tokens from llist

    tlist = ' '

    Do i=1,nrecords
       ctemp = llist(i)
       nt=0
       Do While(Len_trim(ctemp).Gt.0)
          ic = Minval(Index(ctemp,divide),mask=Index(ctemp,divide)>0)
          If(ic.Gt.1) Then
             nt=nt+1
             tlist(nt,i) = Adjustl(ctemp(:ic-1))
          Endif
          ctemp = Adjustl(ctemp(ic+1:))
       End Do
    End Do

    ! Check if any of the "labels" in the input file are unrecognised

    inblock=.False.
    Do i=1,nrecords
       ! Check if we are in a block
       If(esdf_reduce(tlist(1,i)).Eq.'%block') Then
          inblock = .True.
          ! Check if block label is recognised
          If((Count(esdf_reduce(tlist(2,i)).Eq.kw%label).Eq.0)) Then
             ctemp='Label "'//Trim(esdf_reduce(tlist(2,i)))//&
                  &'" not in keyword list'
             If(Count(ctemp.Eq.warns).Eq.0) Call esdf_warn(ctemp) 
          Endif
          ! Check if "label" is multiply defined in the input file
          ndef=0
          Do j=1,nrecords
             If(esdf_reduce(tlist(2,i)).Eq.esdf_reduce(tlist(2,j))) ndef=ndef+1
          End Do
          ctemp='Label "'//Trim(esdf_reduce(tlist(2,i)))//&
               &'" is multiply defined in the input file. '
          If((ndef.Gt.2).And.(Count(ctemp.Eq.warns).Eq.0))&
               Call esdf_warn(ctemp)
       Endif
       ! Check it is in the list of keywords
       If((Count(esdf_reduce(tlist(1,i)).Eq.kw%label).Eq.0)&
            .And.(.Not.inblock)) Then
          ctemp='Label "'//Trim(esdf_reduce(tlist(1,i)))//&
               &'" not in keyword list'
          If(Count(ctemp.Eq.warns).Eq.0) Call esdf_warn(ctemp)
       Endif
       If(.Not.inblock) Then
          ! Check if "label" is multiply defined in the input file
          ndef=0
          Do j=1,nrecords
             If(esdf_reduce(tlist(1,i)).Eq.esdf_reduce(tlist(1,j))) ndef=ndef+1
          End Do
          ctemp='Label "'//Trim(esdf_reduce(tlist(1,i)))//&
               &'" is multiply defined in the input file. '
          If((ndef.Gt.1).And.(Count(ctemp.Eq.warns).Eq.0)) &
               Call esdf_warn(ctemp)
       Endif
       ! Check if we have left a block
       If(esdf_reduce(tlist(1,i)).Eq.'%endblock') inblock= .False.

    End Do

  End Subroutine esdf_init

  !  
  ! Return the string attached to the "label"
  !

  Function esdf_string(label,default)

    Character(*), Intent(in) :: label,default
    Character(llength) :: esdf_string

    ! Local

    Integer(I4B) :: i
    Character(llength) :: ctemp
!jlm begin
    Integer(I4B) :: ioerr

    ! Check "label" is defined

    Call esdf_lblchk(label,'T')

    ! Set to default

    esdf_string = default
    ioerr = 0

    Do i=1,nrecords

       ! Search in the first token for "label"
       ! The first instance is returned

       If(esdf_reduce(tlist(1,i)).Eq.esdf_reduce(label)) Then
          esdf_string = llist(i)(Index(llist(i),Trim(tlist(2,i))):)
          Exit
       Endif

    End Do

   ! Dump the string used 

    ndmp=ndmp+1
    Write(dump(ndmp),*,iostat = ioerr)  Trim(esdf_reduce(label)),':',Trim(esdf_string)
    If(Count(dump(ndmp).Eq.dump(1:ndmp-1)).Gt.0) ndmp=ndmp-1

    If(ioerr /= 0) then
       ctemp = 'Unable to dump "'//Trim(esdf_reduce(label))//'" in esdf_string'
       Call esdf_warn(ctemp)
   EndIf

    Return

  End Function esdf_string

  !  
  ! Return the integer attached to the "label"
  !

  Function esdf_integer(label,default)

    Integer(I4B), Intent(in) :: default
    Character(*), Intent(in) :: label
    Integer(I4B) :: esdf_integer

    ! Local

    Integer(I4B) :: i
    Character(llength) :: ctemp
!jlm begin
    Integer(I4B) :: ioerr
    
    ! Check "label" is defined

    Call esdf_lblchk(label,'I')

    ! Set to default

    esdf_integer = default
    ioerr = 0

    Do i=1,nrecords

       ! Search in the first token for "label"
       ! The first instance is returned

       If(esdf_reduce(tlist(1,i)).Eq.esdf_reduce(label)) Then
          Read(tlist(2,i),*,iostat = ioerr) esdf_integer
          Exit
       Endif

    End Do
    
    If(ioerr /= 0) then
       ctemp = 'Unable to parse "'//Trim(esdf_reduce(label))//'" in esdf_integer'
       Call esdf_die(ctemp)
    EndIf

    ! Dump the value used 

    ndmp=ndmp+1
    Write(dump(ndmp),*,iostat = ioerr)  Trim(esdf_reduce(label)),':',esdf_integer
    If(Count(dump(ndmp).Eq.dump(1:ndmp-1)).Gt.0) ndmp=ndmp-1
    
    If(ioerr /= 0) then
       ctemp = 'Unable to dump "'//Trim(esdf_reduce(label))//'" in esdf_integer'
       Call esdf_warn(ctemp)
    EndIf
!jlm end

    Return

  End Function esdf_integer

  !  
  ! Return the single precisioned value attached to the "label"
  !

  Function esdf_single(label,default)

    Real(SP), Intent(in) :: default
    Character(*), Intent(in) :: label
    Real(SP) :: esdf_single

    ! Local

    Integer(I4B) :: i
    Character(llength) :: ctemp
!jlm begin
    Integer(I4B) :: ioerr

    ! Check "label" is defined

    Call esdf_lblchk(label,'S')

    ! Set to default

    esdf_single = default
    ioerr = 0

    Do i=1,nrecords

       ! Search in the first token for "label"
       ! The first instance is returned

       If(esdf_reduce(tlist(1,i)).Eq.esdf_reduce(label)) Then
          Read(tlist(2,i),*,iostat = ioerr) esdf_single
          Exit
       Endif

    End Do
    
    If(ioerr /= 0) then
       ctemp = 'Unable to parse "'//Trim(esdf_reduce(label))//'" in esdf_single'
       Call esdf_die(ctemp)
    EndIf

    ! Dump the value used 

    ndmp=ndmp+1
    Write(dump(ndmp),*,iostat = ioerr)  Trim(esdf_reduce(label)),':',esdf_single
    If(Count(dump(ndmp).Eq.dump(1:ndmp-1)).Gt.0) ndmp=ndmp-1

    If(ioerr /= 0) then
       ctemp = 'Unable to dump "'//Trim(esdf_reduce(label))//'" in esdf_single'
       Call esdf_warn(ctemp)
!jlm end
    EndIf

    Return

  End Function esdf_single

  !  
  ! Return the double precisioned value attached to the "label"
  !

  Function esdf_double(label,default)

    Real(DP), Intent(in) :: default
    Character(*), Intent(in) :: label
    Real(DP) :: esdf_double

    ! Local

    Integer(I4B) :: i
    Character(llength) :: ctemp
!jlm begin
    Integer(I4B) :: ioerr

    ! Check "label" is defined

    Call esdf_lblchk(label,'D')
    
    ! Set to default

    esdf_double = default
    ioerr = 0

    Do i=1,nrecords

       ! Search in the first token for "label"
       ! The first instance is returned

       If(esdf_reduce(tlist(1,i)).Eq.esdf_reduce(label)) Then
          Read(tlist(2,i),*,iostat = ioerr) esdf_double
          Exit
       Endif

    End Do  
    
    If(ioerr /= 0) then
       esdf_double = default
       ctemp = 'Unable to parse "'//Trim(esdf_reduce(label))//'" in esdf_double'
       Call esdf_die(ctemp)
    EndIf

  ! Dump the value used 

    ndmp=ndmp+1
    Write(dump(ndmp),*,iostat = ioerr)  Trim(esdf_reduce(label)),':',esdf_double
    If(Count(dump(ndmp).Eq.dump(1:ndmp-1)).Gt.0) ndmp=ndmp-1
    
    If(ioerr /= 0) then
       ctemp = 'Unable to dump "'//Trim(esdf_reduce(label))//'" in esdf_double'
       Call esdf_warn(ctemp)
    EndIf
!jlm end

    Return

  End Function esdf_double

  !  
  ! Return the double precisioned physical value attached to the "label"
  ! Units converted to "dunit"
  !

  Function esdf_physical(label,default,dunit)

    Real(DP), Intent(in) :: default
    Character(*), Intent(in) :: label,dunit
    Real(DP) :: esdf_physical

    ! Local

    Integer(I4B) :: i
    Character(llength) :: ctemp,iunit
!jlm begin
    Integer(I4B) :: ioerr1,ioerr2

    ! Check "label" is defined

    Call esdf_lblchk(label,'P')

    ! Set to default

    esdf_physical = default
    ioerr1 = 0
    ioerr2 = 0

    Do i=1,nrecords

       ! Search in the first token for "label"
       ! The first instance is returned

       If(esdf_reduce(tlist(1,i)).Eq.esdf_reduce(label)) Then
          Read(tlist(2,i),*,iostat = ioerr1) esdf_physical
          Read(tlist(3,i),*,iostat = ioerr2) iunit
          esdf_physical = esdf_convfac(iunit,dunit)*esdf_physical
          Exit
       Endif

    End Do
    
    If(ioerr1 /= 0 .OR. ioerr2 /= 0) then
       esdf_physical = default
       ctemp = 'Unable to parse "'//Trim(esdf_reduce(label))//'" in esdf_physical'
       Call esdf_die(ctemp)
    EndIf

     ! Dump the value used 

    ndmp=ndmp+1
    Write(dump(ndmp),*,iostat = ioerr1)  Trim(esdf_reduce(label)),':',esdf_physical,' ',Trim(dunit)
    If(Count(dump(ndmp).Eq.dump(1:ndmp-1)).Gt.0) ndmp=ndmp-1
    
    If(ioerr1 /= 0) then
       ctemp = 'Unable to dump "'//Trim(esdf_reduce(label))//'" in esdf_physical'
       Call esdf_warn(ctemp)
   EndIf
!jlm end

    Return


  End Function esdf_physical

  !
  ! Is the "label" defined in the input file
  !

  Function esdf_defined(label)

    Character(*), Intent(in) :: label
    Logical :: esdf_defined

    ! Local

    Integer(I4B) :: i
    Character(llength) :: ctemp
!jlm begin
    Integer(I4B) :: ioerr

    ! Check "label" is defined

    Call esdf_lblchk(label,'E')

    ! Set to default

    esdf_defined = .False.
    ioerr = 0

    Do i=1,nrecords

       ! Search in the first token for "label"
       ! The first instance is returned

       If(esdf_reduce(tlist(1,i)).Eq.esdf_reduce(label)) Then
          esdf_defined = .True.
          Exit
       Endif

    End Do

    ! Dump the value used 

    If(esdf_defined) Then

       ndmp=ndmp+1
       Write(dump(ndmp),*,iostat = ioerr)  Trim(esdf_reduce(label)),':'
       If(Count(dump(ndmp).Eq.dump(1:ndmp-1)).Gt.0) ndmp=ndmp-1

    Endif

    If(ioerr /= 0) then
       ctemp = 'Unable to dump "'//Trim(esdf_reduce(label))//'" in esdf_defined'
       Call esdf_warn(ctemp)
    EndIf
!jlm end

    Return


  End Function esdf_defined
  !
  ! Is the "label" defined in the input file
  !

  Function esdf_boolean(label,default)

    Character(*), Intent(in) :: label
    Logical, Intent(in) :: default
    Logical :: esdf_boolean

    ! Local

    Integer(I4B) :: i
    Character(llength) :: ctemp,positive(3),negative(3)

    Data positive /'yes','true','t'/
    Data negative /'no','false','f'/
!jlm begin
    Integer(I4B) :: ioerr

    ! Check "label" is defined

    Call esdf_lblchk(label,'L')

    ! Set to default

    esdf_boolean = default
    ioerr = 0

    Do i=1,nrecords

       ! Search in the first token for "label"
       ! The first instance is returned

       If(esdf_reduce(tlist(1,i)).Eq.esdf_reduce(label)) Then
          If(Len_trim(tlist(2,i)).Eq.0) Then
             esdf_boolean = .True.
             Exit
          Endif
          If(Any(Index(positive,esdf_reduce(tlist(2,i))).Gt.0)) Then
             esdf_boolean = .True.
             Exit
          Endif
          If(Any(Index(negative,esdf_reduce(tlist(2,i))).Gt.0)) Then
             esdf_boolean = .False.
             Exit
          Endif
          Call esdf_die('Unable to parse boolean value')

       Endif

    End Do

    ! Dump the value used 

    ndmp=ndmp+1
    Write(dump(ndmp),*,iostat = ioerr)  Trim(esdf_reduce(label)),': ',esdf_boolean
    If(Count(dump(ndmp).Eq.dump(1:ndmp-1)).Gt.0) ndmp=ndmp-1
    
    If(ioerr /= 0) then
       ctemp = 'Unable to dump "'//Trim(esdf_reduce(label))//'" in esdf_boolean'
       Call esdf_warn(ctemp)
    EndIf
!jlm end

    Return


  End Function esdf_boolean

  Function esdf_block(label,nlines)

    Character(*), Intent(in) :: label
    Integer(I4B), Intent(out) :: nlines
    Logical :: esdf_block

    ! Local

    Integer(I4B) :: i,j
    Character(llength) :: ctemp
!jlm begin
    Integer(I4B) :: ioerr

    ! Check "label" is defined

    Call esdf_lblchk(label,'B')

    ctemp ='Block "'//Trim(esdf_reduce(label))//'" not closed correctly '
    
    esdf_block=.False.
    ioerr = 0

    nlines = 0

    Do i=1,nrecords
       If((esdf_reduce(tlist(1,i)).Eq.esdf_reduce('%block'))&
            .And.(esdf_reduce(tlist(2,i)).Eq.esdf_reduce(label))) Then
          esdf_block = .True.
          Do While(esdf_reduce(tlist(1,i+nlines+1))&
               .Ne.esdf_reduce('%endblock'))
             nlines=nlines+1
             If(nlines+i.Gt.nrecords) Call esdf_die(ctemp)
             block_data(nlines)=llist(i+nlines)
          End Do

          If(esdf_reduce(tlist(2,i+nlines+1)).Ne.esdf_reduce(label))&
               Call esdf_die(ctemp)
          Exit
       Endif
    End Do
    If(.Not.esdf_block) Return

    ! Dump the block 

    ndmp=ndmp+1
    Write(dump(ndmp),*,iostat = ioerr) '%block ',Trim(esdf_reduce(label)),': '
    
    If(ioerr /= 0) then
       ctemp = 'Unable to dump "'//Trim(esdf_reduce(label))//'" in esdf_block'
       Call esdf_warn(ctemp)
    EndIf

    If(Count(dump(ndmp).Eq.dump(1:ndmp-1)).Gt.0) Then
       ndmp=ndmp-1
       Return
    Endif
    Do j=1,nlines
       ndmp=ndmp+1
       dump(ndmp)=block_data(j)
    End Do
    ndmp=ndmp+1
    Write(dump(ndmp),*,iostat = ioerr) '%endblock ',Trim(esdf_reduce(label)),': '
    
    If(ioerr /= 0) then
       ctemp = 'Unable to dump "'//Trim(esdf_reduce(label))//'" in esdf_block'
       Call esdf_warn(ctemp)
    EndIf

    Return

    
  End Function esdf_block

  !
  ! Reduce the string to lower case and remove punctuation
  !

  Function esdf_reduce(string)

    Character(*), Intent(in) :: string
    Character(llength) :: esdf_reduce

    ! Local

    Integer(I4B), Parameter :: npunct=3
    Integer(I4B) :: iA,iZ,ishift,ic,i,ln
    Character(llength) :: ctemp
    Character(1) :: punct(npunct)
    Logical :: keep_going

    ! Define the punctuation to be removed

    Data punct /'.','_','-'/

    ! Initialise system dependant bounds in collating sequence

    iA = Ichar('A');iZ = Ichar('Z') 
    ishift = Ichar('a')-iA 

    ! Initialise output

    ln = len(string)

    esdf_reduce(1:ln) = string(1:ln) ; esdf_reduce(ln+1:)=' '

    ! Drop all upper case characters to lower case

    Do i=1,llength
       ic = Ichar(esdf_reduce(i:i))
       If((ic.Ge.iA).And.(ic.Le.iZ)) esdf_reduce(i:i) = Char(ishift+ic) 
    Enddo

    ! Now remove punctuation

    Do i=1,npunct

       keep_going=.True.

       Do While(keep_going)
          ic = Index(esdf_reduce,punct(i))
          If(ic.Gt.0) Then
             ctemp = esdf_reduce
             esdf_reduce(ic:)=ctemp(ic+1:)
          Else
             keep_going=.False.
          End If
       End Do

    End Do

    esdf_reduce = Trim(Adjustl(esdf_reduce))

  End Function esdf_reduce

  !
  ! Find the conversion factor between physical units
  !

  Function esdf_convfac(from,to)

    Character(*), Intent(in) :: from,to
    Real(DP) :: esdf_convfac

    ! Local

    Integer(I4B) :: i,ifrom,ito
    Character(llength) :: ctemp

    ! Find the index numbers of the from and to units

    ifrom = 0 ; ito = 0
    Do i=1,nphys
       If(esdf_reduce(from).Eq.phy(i)%n) ifrom = i
       If(esdf_reduce(to).Eq.phy(i)%n) ito = i
    End Do

    ! Check that the units were recognised

    If(ifrom.Eq.0) Then
       ctemp = 'Units not recognised in input file : '//Trim(esdf_reduce(from))
       Call esdf_die(ctemp)
    Endif

    If(ito.Eq.0) Then
       ctemp = 'Units not recognised in Program : '//Trim(esdf_reduce(to))
       Call esdf_die(ctemp)
    Endif

    ! Check that from and to are of the same dimensions

    If(phy(ifrom)%d.Ne.phy(ito)%d) Then
       ctemp = 'Dimensions Do not match : '//Trim(esdf_reduce(from))&
            & //' vs '//Trim(esdf_reduce(to))
       Call esdf_die(ctemp)
    Endif

    ! Set the conversion factor

    esdf_convfac = phy(ifrom)%u/phy(ito)%u

  End Function esdf_convfac

  ! 
  ! Find an unused i/o unit
  !

  Function esdf_unit(ierr)
    Integer(I4B), Intent(out) :: ierr
    Integer(I4B) :: esdf_unit
    ! Local 
    Logical :: op
    ierr=0
    Do esdf_unit=10,99
       Inquire(unit=esdf_unit,opened=op,err=100)
       If(.Not.op) Return
    End Do
    Call esdf_warn('Unable to find a free i/o unit using esdf_unit')
    ierr = 1
    Return
100 Call esdf_die('Error opening files by esdf_unit')
  End Function esdf_unit

  !
  ! Open an old file
  !

  Subroutine esdf_file(unit,filename,ierr)
    Character(*), Intent(in) :: filename
    Integer(I4B), Intent(out) :: unit,ierr
    Logical :: ex
    unit = esdf_unit(ierr)
    If(ierr.Gt.0) Return
    Inquire(file=Trim(filename),exist=ex,err=100)
    If(.Not.ex) Goto 100
    Open(unit=unit,file=Trim(filename),form='formatted',status='old',err=100)
    Return
100 ierr=1
    Return
  End Subroutine esdf_file

  ! Open a new file

  Subroutine esdf_newfile(unit,filename,ierr)
    Character(*), Intent(in) :: filename
    Integer(I4B), Intent(out) :: unit,ierr
    unit = esdf_unit(ierr)
    If(ierr.Gt.0) Return
    Open(unit=unit,file=Trim(filename),form='formatted',status='replace',err=100)
    Return
100 ierr=1
    Return
  End Subroutine esdf_newfile

  !
  ! Check that the label is known, and used correctly
  !

  Subroutine esdf_lblchk(string,typ)
    Character(*), Intent(in) :: string
    Character(1), Intent(in) :: typ
    ! Local
    Character(llength) :: ctemp
    Character(1) :: tp
    Integer(I4B) :: i
    ! Check if label is recognised

    i=Count(esdf_reduce(string).Eq.kw%label)
    ctemp = 'Label "'//Trim(esdf_reduce(string))//'" not recognised in&
         & keyword list'
    If(i.Eq.0) Call esdf_die(ctemp)
    ctemp = 'Label "'//Trim(esdf_reduce(string))//'" is multiply defined'
    If(i.Gt.1) Call esdf_die(ctemp)
    ctemp = 'Label "'//Trim(esdf_reduce(string))//'" has been used with the wrong type'
    tp = ' '
    i=0
    Do While(tp.Eq.' ')
       i=i+1
       If(esdf_reduce(string).Eq.kw(i)%label) tp=kw(i)%typ
    End Do

!   modification by JLM

    If(typ.Ne.tp .AND. typ /= 'E') Call esdf_die(ctemp)
  End Subroutine esdf_lblchk


  Subroutine esdf_help(helpword,searchword)

    Implicit None

    Character(*) :: helpword,searchword

    ! Local

    Integer(I4B)  :: i,indx,indx2,ln
    Character(20) :: ctyp,clev
    Character(60) :: title,fm
    Character(80) :: ctemp
    Character(1)  :: cl

    helpword = esdf_reduce(helpword)
    searchword = esdf_reduce(searchword)

    If(esdf_reduce(helpword).Eq.'search') Then

       If(Len_trim(searchword).Lt.1) Call esdf_die('help: "searchword" is empty')

       ! Search for useful keywords

       Do i=1,numkw
          If((Index(kw(i)%label,Trim(searchword)).Gt.0).Or.&
               (Index(kw(i)%dscrpt,Trim(searchword)).Gt.0)) Then 
             indx=Index(kw(i)%dscrpt,'!*')-1
             If(indx.Eq.-1) &
                  Call esdf_die('help: keyword description incorrectly formatted')
             title = kw(i)%dscrpt(1:indx)
             ln=Len(Trim(title))
             If(ln.Gt.80) Call esdf_die('help: keyword title too long')

             Write (*,*) kw(i)%label,Trim(title)

          End If
       End Do
       Stop


    Endif

    ! All keywords, short description

    If('all'.Eq.helpword) Then
       Do i=1,numkw
          If(Len_trim(kw(i)%label).Gt.0) Then 
             indx=Index(kw(i)%dscrpt,'!*')-1
             If(indx.Eq.-1) &
                  Call esdf_die('help: keyword description incorrectly formatted')
             title = kw(i)%dscrpt(1:indx)
             ln=Len(Trim(title))
             If(ln.Gt.80) Call esdf_die('help: keyword title too long')

             Write (*,*) kw(i)%label,Trim(title)

          End If
       End Do
       Stop
    End If

    ! All specific levels of keywords

    If(Any((/'basic ','inter ','expert','dummy '/).Eq.helpword)) Then

       Select Case(helpword)
       Case('basic')  ; cl = 'B'
       Case('inter')  ; cl = 'I'
       Case('expert') ; cl = 'E'
       Case('dummy')  ; cl = 'D'
       End Select

       Do i=1,numkw
          If(kw(i)%typ(3:3).Eq.cl) Then 
             indx=Index(kw(i)%dscrpt,'!*')-1
             If(indx.Eq.-1) &
                  Call esdf_die('help: keyword description incorrectly formatted')
             title = kw(i)%dscrpt(1:indx)
             ln=Len(Trim(title))
             If(ln.Gt.80) Call esdf_die('help: keyword title too long')

             Write (*,*) kw(i)%label,Trim(title)

          End If
       End Do
       Stop
    End If

    ! More information about a specific keyword

    If(.Not.Any(kw%label.Eq.helpword)) &
         Call esdf_die('help: keyword not recognised')
    If(Count(kw%label.Eq.helpword).Gt.1) &
         Call esdf_die('help: keyword entry duplicated')
    Do i=1,numkw
       If(kw(i)%label.Eq.helpword) Then 
          indx=Index(kw(i)%dscrpt,'!*')+1
          If(indx.Eq.1) &
               Call esdf_die('help: keyword description incorrectly formatted')
          title = kw(i)%dscrpt(1:indx)
          ln=Len(Trim(title))
          If(ln.Gt.80) &
               Call esdf_die('help: keyword title too long')
          If(ln.Le.9) Write(fm,'("(",i2,"x,a",i1,")")') 40-ln/2,ln
          If(ln.Gt.9) Write(fm,'("(",i2,"x,a",i2,")")') 40-ln/2,ln
          Write (*,fm) Trim(title)
          Write (*,*)
          Select Case(kw(i)%typ(1:1))
          Case('I') ; ctyp ='Integer'
          Case('S') ; ctyp ='Single Precision'
          Case('D') ; ctyp ='Double Precision'
          Case('P') ; ctyp ='Physical'
          Case('T') ; ctyp ='String'
          Case('E') ; ctyp ='Defined'
          Case('B') ; ctyp ='Block'
          Case('L') ; ctyp ='Boolean'   
          End Select
          Select Case(kw(i)%typ(3:3))
          Case('B') ; clev ='Basic'
          Case('I') ; clev ='Intermediate'
          Case('E') ; clev ='Expert'
          Case('D') ; clev ='Dummy'
          End Select
          Write (fm,'(a,i2,a)') '("Type: ",a,',&
               78-(7+Len_trim(clev))-(6+Len_trim(ctyp)),'x," Level: ",a)'
          Write (ctemp,fm) Trim(ctyp),Trim(clev)
          Write (*,'(a)') Trim(ctemp)
          Write (*,*)
          indx=indx+1
          ln = Len(Trim(kw(i)%dscrpt))
          Do While (indx.Lt.ln)
             ctemp = kw(i)%dscrpt(indx:Min(indx+80,ln))
             indx2=Index(ctemp,' ',back=.True.)
             Write (*,'(a)') Trim(Adjustl(ctemp(:indx2)))
             indx=indx+Len(ctemp(:indx2))
          End Do

       End If
    End Do

    Stop

  End Subroutine esdf_help

  ! 
  ! Stop execution due to an error cause by esdf
  !

  Subroutine esdf_die(string)
    Character(*), Intent(in) :: string
    Write (*,'(a,a)') ' ESDF ERROR: ',Trim(string)
    Write (*,'(a)') ' Stopping now'
    Stop    
  End Subroutine esdf_die

  ! 
  ! Warning due to an error cause by esdf
  !

  Subroutine esdf_warn(string)
    Character(*), Intent(in) :: string
    nwarns=nwarns+1
    warns(nwarns) = string
  End Subroutine esdf_warn

  !
  ! Dump the warnings to screen
  !

  Subroutine esdf_warnout
    Integer(I4B) :: i
    Do i=1,nwarns
       Write (*,*) 'ESDF WARNING: '//Trim(warns(i))
    End Do
  End Subroutine esdf_warnout

  !
  ! Deallocate the data arrays --- call this before re-initialising
  !

  Subroutine esdf_close
    Deallocate(llist,tlist,block_data)

!   added by JLM

    deallocate(warns,dump)
  End Subroutine esdf_close

  !
  ! Dump an input file which contains all set variables
  ! including defaults
  !

  Subroutine esdf_dump(filename)

    ! Local

    Integer(I4B) :: unit_dump,ierr,i,j,indx
    Character(*) :: filename
    Character(llength) :: cjunk

    ! Open the ESDF.esdf file

    Call esdf_newfile(unit_dump,filename,ierr)
    If(ierr.Eq.1) Call esdf_die('Unable to open main input file "ESDF.esdf"')

    indx = Maxval(Index(dump(1:ndmp),':'))

    Do i=1,ndmp
       j=Index(dump(i),':')
       If(j.Gt.0) Then
          cjunk = dump(i)(1:j-1)
          dump(i)=Trim(cjunk)//Repeat(' ',indx-j+1)//': '&
               & //Trim(Adjustl(dump(i)(j+1:)))//'#'
       Endif
    End Do

    indx = Maxval(Index(dump(1:ndmp),'#',back=.True.))

    Do i=1,ndmp
       j=Index(dump(i),'#',back=.True.)
       If(j.Gt.0) Then

          dump(i)=dump(i)(1:j-1)//Repeat(' ',indx-j+1)//'#'
       Endif
    End Do

    Do i=1,ndmp
       j=Index(dump(i),':')
       If(j.Gt.0) Then
          cjunk = dump(i)(1:j-1)
          Do j=1,numkw
             If(Index(cjunk,Trim(kw(j)%label)).Gt.0) Exit
          End Do
          Select Case(kw(j)%typ(1:1))
          Case('I') ; cjunk ='Integer'
          Case('S') ; cjunk ='Single Precision'
          Case('D') ; cjunk ='Double Precision'
          Case('P') ; cjunk ='Physical'
          Case('T') ; cjunk ='String'
          Case('E') ; cjunk ='Defined'
          Case('B') ; cjunk ='Block'
          Case('L') ; cjunk ='Boolean'   
          End Select
          indx=Index(kw(j)%dscrpt,'!*')
          dump(i)=Trim(dump(i))//Trim(kw(j)%dscrpt(1:indx-1))//' ('//Trim(Adjustl(cjunk))//')'
       Endif
    End Do

    Do i=1,ndmp
       Write(unit_dump,'(a)') Adjustl(dump(i))
    End Do

  End Subroutine esdf_dump

End Module esdf
 
