      program global_terrain

      include 'newmpar.h'
      include 'dates.h'   ! to pass ds from setxyz
      include 'indices.h'
      include 'parm.h'   ! to pass rlong0,rlat0,schmidt  to setxyz
      logical olam

c     various work arrays
      common/work/rmsk(il,jl),inum(il,jl),inumx(il,jl),zss(il,jl)
     &,add(il,jl),almsk(il,jl)
     &,tmax(il,jl),tmin(il,jl),tsd(il,jl)

      logical debug
      character*60 fileout
      character*9 formout

      integer mx,my
      integer nx,ny
      parameter (mx=4800,my=6000)
      integer i,j
      character*1 ns,ew
      character*11 file
      integer lons,lats

      data debug/.false./

      do i = 1,9
       lons=-140+(i-1)*40
       ew="E"
       if ( lons.gt.0 ) ew="W"
       write(6,*)"lons=",lons

       do j = 1,3
        lats=90-(j-1)*50 
        ns="N"
        if ( lats.lt.0 ) ns="S"
        write(6,*)"lons,lats=",lons,lats

        write(file,'(a1,i3.3,a1,i2.2,".DEM")') ew,abs(lons),ns,abs(lats)
        nx=4800
        ny=6000
        !call readglob(file,nx,ny,lons,lats,debug)
       enddo ! j

      enddo ! i

      call readglob("E140S10.DEM",nx,ny,140,-10,debug)

      do i = 1,6
       lons=-120+(i-1)*60
       ew="E"
       if ( lons.ge.0 ) ew="W"
       lats=-60
       write(6,*)"lons,lats=",lons,lats

       write(file,'(a1,i3.3,"S60.DEM")') ew,abs(lons)
       nx=7200
       ny=3600
       !call readglob(file,nx,ny,lons,lats,debug)
      enddo ! i

      stop
      end
c=======================================================================
      subroutine readglob(file,nx,ny,lons,lats,debug)

      character*11 file
      integer*2 zz(nx),iix,iin
      integer nrecl,lons,lats

      write(6,*)"file,nx,ny,nrecl=",file,nx,ny,nrecl
      nrecl = nx*2

      open(21,file=file,access='direct',form='unformatted',recl=nrecl)
      do j=1,ny
        !read(21,rec=1+ny-j)(zz(i),i=1,nx)
        read(21,rec=j)zz
        rix=-1.e35
        rin=+1.e35
        do i=1,nx
          rix=max(rix,real(zz(i)))
          rin=min(rin,real(zz(i)))
        enddo !i=1,nx
        if(rix.gt.0)write(6,*)"j,rix,rin=",j,rix,rin
      enddo

      close(21)

      return
      end
c=======================================================================
      subroutine amap(w, im, jm, label, cinti, base)
c**********************************************************************c
c      a paper map shading routine for an 'a' grid                     c
c      w = data dimension im*jm                                        c
c      label = char*32 label                                           c
c      cinti = contour interval ( will pick one if 0 )                 c
c      base = base value for countours                                 c
c**********************************************************************c
      parameter ( ncol=72, nrow=40 )
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
!     if ( cint .gt. 1.e-10 ) then
!        ncont = ( wmax - wmin ) / cint
!        if ( ncont.gt.40 .or. ncont.lt.2 ) then
!           cint = ( wmax - wmin ) / 10.
!        endif
!     else
!        cint = ( wmax - wmin ) / 10.
!     endif
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
