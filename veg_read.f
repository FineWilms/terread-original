      program veg

      include 'newmpar.h'

      parameter ( rtd=180.0/3.14159265 )
      parameter ( dl = 1.0/120.0 )
      parameter ( gslon = -180+dl/2. )
      parameter ( gslat = 90-dl/2. )
      parameter ( mx=43200, my=21600 )

      character*13 gfile
      parameter( gfile='gsib2_0ll.img')

      integer i,j
      byte zza(mx)
      integer iz(mx)
      real veg_mod(500*500)

! global 1 km data
      write(*,*) 'gfile=',gfile
      open(20,file=gfile,access='direct',form='unformatted',recl=mx)

      gridn = 1.0
      gridx = 1.0
      rlonn = -71.6
      rlonx = -71.0
      rlatx =  41.8
      rlatn =  41.4

      imin = nint((rlonn - gslon)/dl) + 1
      imax = nint((rlonx - gslon)/dl) + 1
      jmin = nint((gslat - rlatx)/dl) + 1
      jmax = nint((gslat - rlatn)/dl) + 1
                                                                                
      write(6,*)"grid n,x=",gridn,gridx
      write(6,*)"rlon n,x=",rlonn,rlonx
      write(6,*)"rlat n,x=",rlatn,rlatx
      write(6,*)"i    n,x=",imin,imax
      write(6,*)"j    n,x=",jmin,jmax

      im=imax-imin+1
      jm=jmax-jmin+1
      write(6,*)"im,jm=",im,jm

      do jg=jmin,jmax
        jmod=jmax-jg+1

        glat=gslat-(jg-1)*dl
!       write(6,*)"############## jg,glat=",jg,glat,jmax,my

        if (glat.gt.rlatn .and. glat.lt.rlatx) then

! read each row of relevant data from 1 km data set (row=jg)
!         write(6,*)"read(20,rec=jg) (zza(ig),ig=1,mx)",mx
          read(20,rec=jg) (zza(ig),ig=1,mx)

          do ig=imin,imax
            imod=1+ig-imin

            glon = gslon + (ig-1)*dl
            if (glon.gt.rlonn .and. glon.lt.rlonx) then
              if(zza(ig).ge.19)zza(ig)=0
              iq=imod+im*(jmod-1)
              veg_mod(iq)=zza(ig)
            endif

          enddo ! ig=1,mx

        endif ! (glat.gt.rlatn .and. glat.lt.rlatx) then

      enddo ! jg=jmin,jmax

      call amap(veg_mod, im, jm, 'veg', 1., 0.)

      stop
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
      parameter ( ncol=110, nrow=50 )
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
