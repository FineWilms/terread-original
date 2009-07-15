      program readdem

      character*11 file
      character*1 ew,ns
      real lons,lone,lats,late
      real zs(4800*6000)
      real dl

      dl = 1.0/120.0
      nx=4800
      ny=6000

      call read1km("W020N40.DEM",nx,ny,zs)

      call ncdf_setup("dem.nc",idnc,3,nx,ny,1,20020101,0000)
      call nc2out(zs,nx,ny,1,1.,idnc,"zs","zs","m",0.,5000.)

      stop

      !do i = 1,9
      do i = 3,3
       lons=-180.+real(i-1)*40.
       lone=lons+real(nx-1)*dl
       ew="W"
       if ( nint(lons).gt.0 ) ew="E"

       !do j = 1,3
       do j = 1,1
        lats=90.-real(j-1)*50.
        late=lats-real(ny-1)*dl
        ns="N"
        if ( nint(lats).lt.0 ) ns="S"

        write(file,'(a1,i3.3,a1,i2.2,".DEM")')
     &               ew,abs(nint(lons)),ns,abs(nint(lats))

        write(6,*)"call read1km(file,nx,ny)"
        call read1km(file,nx,ny,zs)

       enddo ! j
      enddo ! i

      stop
      end
c=======================================================================
      subroutine read1km(file,nx,ny,zs)

      character*11 file
      integer zi(7200)
      integer*2 zi2(7200)
      !integer*2 zi2
      real zs(nx,ny)
      integer nrecl

      nrecl = nx*2
      nrecl = nx/2

      write(6,*)"read1km file,nx,ny=",file,nx,ny

       open(CONVERT='BIG_ENDIAN',unit=21,file=file,access='direct'
     &            ,form='unformatted',recl=nrecl)
!      open(21,file=file,status='old',form='unformatted'
!    &          ,access='direct',recl=2)
!    &          ,access='direct',recl=4)

!#######################################################################
      !do j=ny,1,-1
      do j=1,ny,1
        jr=ny+1-j
!#######################################################################

         read(21,rec=j) (zi2(i),i=1,nx)

         rmax=-1.e35
         rmin=+1.e35
         do i=1,nx
           zi2(i)=max(0,zi2(i))
           rmax=max(rmax,real(zi2(i)))
           rmin=min(rmin,real(zi2(i)))
           zs(i,jr)=max(0.,real(zi2(i)))
         enddo !i=1,nx
         !if(rmax.gt.-999.)write(6,*)"ny,jr,rmax,rmin=",ny,jr,rmax,rmin

      enddo ! j

 1001 continue

      close(21)

      call amap(zs, nx, ny, 'zs', 0., 0.)

      return
      end ! read1km
c=======================================================================
!     function int2_conv(in_dat)

!      integer*2 in_dat
!      integer*2 i1,i2
!      character ch1(2), ch2(2)
!      equivalence(i1,ch1)
!      equivalence(i2,ch2)

!      i1=in_dat

!      do k=1,2
!        ch2(k) = ch1(3-k)
!      enddo

!      int_conv = i2

!     end
!=======================================================================
!     function int_conv(in_dat)

!      integer in_dat
!      integer i1,i2
!      character ch1(4), ch2(4)
!      equivalence(i1,ch1)
!      equivalence(i2,ch2)

!      i1=in_dat

!      do k=1,4
!        ch2(k) = ch1(5-k)
!      enddo

!      int_conv = i2

!     end
c=======================================================================
      subroutine amap(w, im, jm, label, cinti, base)
c**********************************************************************c
c      a paper map shading routine for an 'a' grid                     c
c      w = data dimension im*jm                                        c
c      label = char*32 label                                           c
c      cinti = contour interval ( will pick one if 0 )                 c
c      base = base value for countours                                 c
c**********************************************************************c
      parameter ( ncol=120, nrow=80 )
      dimension w(1)
      character*1  blk, print(ncol), ctbl(20)
      character*(*) label
      logical oprt
      data blk/' '/, ctbl/'0',' ','1',' ','2',' ','3',' '
     .                   ,'4',' ','5',' ','6',' ','7'
     .                   ,' ','8',' ','9',' '/
      data c1/1./, p5/.5/, c0/0.0/, c1e5/1.e5/
      data spval/-1.e35/, spp/-9.e34/
      data oprt/.true./
c-----------------------------------------------------------------------
      if ( .not.oprt ) return
c-----------------------------------------------------------------------
c calculate addition grid paramters
      cint = cinti
      cncol = ncol-c1
      cnrow = nrow-c1
      imjm=im*jm
      rjs=1.
      rjf=jm
      ris=1.
      rif=im
c-----------------------------------------------------------------------
c determine number of valid points to use
c-----------------------------------------------------------------------
      ln = imjm
c-----------------------------------------------------------------------
c find maximum and minimum of field
c-----------------------------------------------------------------------
      wmin = -spval
      wmax =  spval
      do 45 n = 1 , ln
        if ( w(n) .lt. wmin ) then
           wmin = w(n)
           nmin = n
        elseif ( w(n) .gt. wmax ) then
           wmax = w(n)
           nmax = n
        endif
   45 continue
      jx=nmax/jm+1
      ix=nmax-(jx-1)*jm
      jn=nmin/jm+1
      in=nmin-(jn-1)*jm
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c change contour interval if number of contours > 40 or < 2
c or calculate contour interval if zero was specified
c-----------------------------------------------------------------------
      if ( cint .gt. 1.e-10 ) then
         ncont = ( wmax - wmin ) / cint
         if ( ncont.gt.40 .or. ncont.lt.2 ) then
            cint = ( wmax - wmin ) / 10.
         endif
      else
         cint = ( wmax - wmin ) / 10.
      endif
c-----------------------------------------------------------------------
c print heading with label, max, min, cint, and base
c-----------------------------------------------------------------------
      if ( wmax .ne. wmin ) then
         if ( max(wmax,abs(wmin)).ge.c1e5
     .        .or. min(wmax,abs(wmin)).lt.c1) then
            write(6,1) label, wmax, ix,jx, wmin, in,jn, cint, base
    1       format(1h0,a32,/,' max=',1p,e10.2,' at',2i6,' min='
     .           ,e10.2,' at',2i6,' cint=',e10.2,' base=',e10.2)
         else
            write(6,2) label, wmax, ix,jx, wmin, in,jn, cint, base
    2       format(1h0,a32,/,' max=',f10.2,' at',2i6,' min='
     .           ,f10.2,' at',2i6,' cint=',f10.2,' base=',f10.2)
         endif
      else
c-----------------------------------------------------------------------
c have a constant field
         write(6,3) label, wmax, nmax, wmin, nmin, cint, base
    3    format(1h0,a32,/,' max=',f10.2,' at',i6,' min=',f10.2
     .           ,' at',i6,' cint=',f10.2,' base=',f10.2)
         return
      endif
c-----------------------------------------------------------------------
c main loop for printing rows
c-----------------------------------------------------------------------
c assume (1,1) at lower left hand corner
      do 80 j=nrow, 1, -1
c fill in column printing array with blanks
        do 50 i = 1 , ncol
   50     print(i) = blk
c-----------------------------------------------------------------------
c main loop columns
c-----------------------------------------------------------------------
        do 60 i = 1 , ncol
c determine the x and y index for each point on the page
          x = ris+(rif-c1)*(i-c1)/cncol
          y = rjs+(rjf-c1)*(j-c1)/cnrow
          ibl = ifix(x)
          jbl = ifix(y)
c determine k index in 1 dimensional space for interpolation
          k  = ibl + im*(jbl-1)
c determine interpolation factors
          rf = x - ibl
          qf = y - jbl
c get the values of the lahm field at the diamond corner points.
          kbl=k
          kbr=k+1
          ktl=kbl+im
          ktr=kbr+im
          pbl=spval
          pbr=spval
          ptl=spval
          ptr=spval
          if ( kbl.ge.1 .and. kbl.le.ln ) pbl = w(kbl)
          if ( kbr.ge.1 .and. kbr.le.ln ) pbr = w(kbr)
          if ( ktl.ge.1 .and. ktl.le.ln ) ptl = w(ktl)
          if ( ktr.ge.1 .and. ktr.le.ln ) ptr = w(ktr)
c determine which points have valid data
          nspp=0
          if ( jbl.lt.1  .or. pbl.lt.spp) nspp=1
          if ( jbl.ge.jm .or. pbr.lt.spp) nspp=nspp+2
          if ( ibl.lt.1  .or. ptl.lt.spp) nspp=nspp+4
          if ( ibl.ge.im .or. ptr.lt.spp) nspp=nspp+8
c perform bilinear interpolation to the x-y gridpoint.
c all points valid
          if ( nspp.eq.0 ) then
             z = (c1-rf)*(c1-qf)*pbl + rf*(c1-qf)*pbr
     .          +(c1-rf)*    qf *ptl + rf*    qf *ptr
          else
c missing more than 1 point
             z = spval
             print(i) = '*'
             go to 60
          endif
c determine what number to print at this location
          miq = mod ( int ( abs(z-base)/cint ) , 20 )
          n   = miq+1
          if ( z.lt.base ) n = 20 - miq
          print(i) = ctbl(n)
c determine if near max/min point
          if ( kbl.eq.nmax .or. kbr.eq.nmax .or. ktl.eq.nmax
     .         .or. ktr.eq.nmax ) print(i) = 'x'
          if ( kbl.eq.nmin .or. kbr.eq.nmin .or. ktl.eq.nmin
     .         .or. ktr.eq.nmin ) print(i) = 'n'
c-----------------------------------------------------------------------
c end of column loop
c-----------------------------------------------------------------------
   60   continue
c print out this row
        write(6,70) (print(k),k=1,ncol)
   70   format (1x,131a1)
c-----------------------------------------------------------------------
c end of row loop
c-----------------------------------------------------------------------
   80 continue
      return
      end
c=======================================================================
