      subroutine read_ht ( lat, ht )

! read topography data for specified latitude

      implicit none

! parameters
      integer      lui ! logical unit for input
      integer      nlong ! number of longitudes
      parameter (lui = 99)
      parameter (nlong = 4320) ! 360 * 12

! arguments
      integer      ht(nlong) ! topographic height in metres (output)
      real         lat ! latitude (input)

! variables
      character    buf(2,nlong) ! buffer
      integer      i ! subscript
      integer      ios ! i/o status
      integer(kind=8) ival ! int. value
      integer      j ! subscript
      logical      new ! true 1st call
      save new
      data new / .true. /

! open file
      if (new) then
         write(6,*)"open file='topo2' lui=",lui
         open (lui, file='topo2',
     &   status='old', access='direct',
     &   form='unformatted', recl=nlong*2, iostat=ios)
         new = .false.
      end if

      i = nint((90.0 - lat) * 12.0 + 1.0)
      if ( i.lt.3 ) print *,'i=',i

      read (lui, rec=i, iostat=ios) buf

      do j = 1, nlong
        ival = 256 * ichar(buf(1,j)) + ichar(buf(2,j))
        if (ival .ge. X'8000') then ! sign extend
           write(6,*) "ERROR: Integer too big"
           stop
           !ival = ior( ival, X'FFFFFFFFFFFF0000' )
        end if
        ht(j) = ival
      end do

      end
