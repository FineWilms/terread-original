      program global_terrain

      parameter (rtd=180.0/3.14159265 )
      parameter ( nter = 37 )

      include 'newmpar.h'
      include 'dates.h'   ! to pass ds from setxyz
      include 'map.h'   ! em
      include 'xyzinfo.h'   ! rlat,rlong
      include 'indices.h'
      include 'parm.h'   ! to pass rlong0,rlat0,schmidt  to setxyz
      logical olam

      include 'rwork.h'

      logical debug, do1km, okx, oky, do250
      character*60 fileout
      character*9 formout

      integer nx,ny
      integer i,j
      character*1 ns,ew
      character*11 file
      real lons,lats
      real lone,late
      character*8 terfil(nter)

      data debug/.true./
      data do1km/.true./
      data do250/.true./
      data idia/24/
      data jdia/72/

      namelist / topnml / ds, du, tanl, rnml, stl1, stl2, debug
     &  ,luout, fileout, olam, wbd, sbd, dlon, dlat
     &  ,idia, jdia ,rlong0, rlat0, schmidt
     &  ,do1km,id,jd,do250

      open ( unit=5,file='top.nml',status='unknown' )
      read ( 5,topnml, end=5 )
 5    write ( 6,topnml )

      call setxyz

      do n = 1,37
        read(88,*) terfil(n)
        write(6,*) terfil(n)
      enddo ! n = 1,37

      gridx=-1.e35
      rlatx=-1.e35
      rlonx=-1.e35
      gridn=+1.e35
      rlatn=+1.e35
      rlonn=+1.e35
      !imax=-999
      !jmax=-999
      !imin=+999
      !jmin=+999
      do j=1,jl
       do i=1,il
         n=i+(j-1)*il
         grid(i,j)=(ds/em(i,j))/1.e3 ! km
         rlond(i,j)=rlong(n)*rtd
         if(rlond(i,j).gt.180.)rlond(i,j)=rlond(i,j)-360.
         rlatd(i,j)=rlat (n)*rtd
         gridx=max(gridx,grid(i,j))
         gridn=min(gridn,grid(i,j))
         if ( grid(i,j).lt.5. ) then    ! use 1km when grid spc < 5 km
           rlatx=max(rlatx,rlatd(i,j))
           rlonx=max(rlonx,rlond(i,j))
           rlatn=min(rlatn,rlatd(i,j))
           rlonn=min(rlonn,rlond(i,j))
           !jmax=max(jmax,j)
           !imax=max(imax,i)
           !jmin=min(jmin,j)
           !imin=min(imin,i)
           !write(6,*)i,j,n,grid(i,j),rtd*rlong(n),rtd*rlat(n)
         endif ! grid < 20 km
       enddo
      enddo
      write(6,*)"grid(km) n,x=",gridn,gridx
      write(6,*)"rlon n,x=",rlonn,rlonx
      write(6,*)"rlat n,x=",rlatn,rlatx
      !write(6,*)"i    x,n=",imax,imin
      !write(6,*)"j    x,n=",jmax,jmin

      !write(97,'(30f6.1)') grid
      !write(98,'(30f6.2)') rlond
      !write(99,'(30f6.2)') rlatd

c     nested model topography (output)
      write(6,*) 'open',luout,fileout
      open(luout,file=fileout,form='formatted',status='unknown')

! netcdf file
      call ncdf_setup("zsm.nc",idnc,3,il,jl,1,20020101,0000)
      call nc2out(grid ,il,jl,1,1.,idnc,"grid","grid","km",0.,5000.)
      call nc2out(rlond,il,jl,1,1.,idnc,"lon","lon","degrees",0.,360.)
      call nc2out(rlatd,il,jl,1,1.,idnc,"lat","lat","degrees",-90.,90.)
      call ncsnc(idnc,ier)

c initialize min,max, and sd arrays
      do j=1,jl
       do i=1,il
        tmin(i,j)=99999.
        tmax(i,j)=-99999.
        tsd(i,j)=0.
        almsk(i,j)=0.
        zss(i,j)=0.
        inum(i,j)=0
       enddo
      enddo

!================================
      if ( do1km ) then
!================================

      do n = 1,37

        open(13,file=terfil(n),form="unformatted")
        read(13) nx,ny,lons,lats,dl
        write(6,*)"nx,ny,dl=", nx,ny,dl
        lone=lons+(nx-1)*dl
        late=lats+(ny-1)*dl

        okx = .false.
        okx = okx .or. (lons.ge.rlonn .and. lons.le.rlonx)
        okx = okx .or. (lone.ge.rlonn .and. lone.le.rlonx)
        okx = okx .or. (lons.le.rlonn .and. lone.ge.rlonx)
        oky = .false.
        oky = oky .or. (lats.ge.rlatn .and. lats.le.rlatx) ! start lat in area
        oky = oky .or. (late.ge.rlatn .and. late.le.rlatx) !  end  lat in area
        oky = oky .or. (lats.le.rlatn .and. late.ge.rlatx) ! grid encompasses hr reg
        oky = oky .or. (lats.ge.rlatn .and. late.le.rlatx) ! grid within hr reg

        write(6,'("lons,lone=",2f6.1," lats,late=",2f6.1,2l4)')
     &             lons,lone,lats,late,okx,oky

        if ( okx.and.oky ) then
          write(6,*)"call read250(nx,ny,lons,lats,dl,debugi,idia,jdia)"
          call read250(nx,ny,lons,lats,dl,debugi,idia,jdia)
        endif

      enddo ! n

      !write(6,*)"inum"
      !do j=49,96
      ! write(6,'(48i2)')(inum(i,j),i=1,48)
      !enddo
      !write(6,*)"rmsk"
      !do j=49,96
      ! write(6,'(48i2)')(nint(rmsk(i,j)),i=1,48)
      !enddo
      !write(6,*)"zss"
      !do j=49,96
      ! write(6,'(48i2)')(nint(zss(i,j)/10),i=1,48)
      !enddo

!================================
      endif ! ( do1km ) then
!================================

      write(6,*)"calling read10km"
      call read10km(debug)

      write(6,*)"inum"
      do j=49,96
       write(6,'(48i2)')(inum(i,j),i=1,48)
      enddo
      write(6,*)"rmsk"
      do j=49,96
       write(6,'(48i2)')(nint(rmsk(i,j)),i=1,48)
      enddo
      write(6,*)"zss"
      do j=49,96
       write(6,'(48i2)')(nint(zss(i,j)/10),i=1,48)
      enddo

      write(6,*)"now compute grid average values"
      numzer=0
      imin=9999
      jmin=9999
      imax=0
      jmax=0

      do j=1,jl
       do i=1,il
          !write(6,*)i,j,inum(i,j),zss(i,j)
          if ( inum(i,j).eq.0 ) then
! no input points found, assume sea point
             numzer=numzer+1
             imin=min(i,imin)
             imax=max(i,imax)
             jmin=min(j,jmin)
             jmax=max(j,jmax)
             zss(i,j) = 0.
             rmsk(i,j)= 0.
             tsd(i,j)= 0.
             tmin(i,j)= 0.
             tmax(i,j)= 0.
          else
! set land-sea mask (rmsk)
             rnum=1./inum(i,j)
             almsk(i,j) = almsk(i,j)*rnum

! if less than half of pnts are ocean, assume it is land pnt
             if ( almsk(i,j).lt. .5 ) then
! ocean point
                rmsk(i,j) = 0.
                zss(i,j) = -.01
             else ! almsk
! land point
                rmsk(i,j) = 1.
                zss(i,j) = zss(i,j)*rnum
             endif ! almsk

! compute sd
             tsd(i,j)=sqrt(abs(tsd(i,j)*rnum-zss(i,j)**2)) ! abs for rounding?
!            tsd(i,j)=sqrt(tsd(i,j)*rnum-zss(i,j)**2)
          end if

       enddo ! i
      enddo ! j

2     if(numzer.gt.0)then
!       this checks for non-data points, and inserts a neighbouring value
        !print *,'**** numzer= ',numzer
        numzer=0
        do j=1,jl
         do i=1,il
          iq=i+(j-1)*il
          inumx(i,j)=inum(i,j)
          if(inum(i,j).eq.0)then
            iq2=0
            if(inum(in(iq),1).ne.0)then
              iq2=in(iq)
            endif
            if(inum(ie(iq),1).ne.0)then
              iq2=ie(iq)
            endif
            if(inum(iw(iq),1).ne.0)then
              iq2=iw(iq)
            endif
            if(inum(is(iq),1).ne.0)then
              iq2=is(iq)
            endif
            if(iq2.ne.0)then
              zss(i,j)=zss(iq2,1)
              rmsk(i,j)=rmsk(iq2,1)
              tsd(i,j)=tsd(iq2,1)
              almsk(i,j)=almsk(iq2,1)
              inumx(i,j)=1
            else
              numzer=numzer+1
            endif    !  (iq2.ne.0)
          endif      !  (inum(i,j).eq.0)
         enddo
        enddo
        do j=1,jl
         do i=1,il
          inum(i,j)=inumx(i,j)
         enddo
        enddo
      endif         ! (numzer.gt.0)
      if(numzer.gt.0)go to 2

      write(luout,'(i3,i4,2f6.1,f6.3,f8.0,''  orog-mask-var'')')
     &                           il,jl,rlong0,rlat0,schmidt,ds

      call nc2out(zss  ,il,jl,1,1.,idnc,"zs","zs","m",-100.,30000.)
      call nc2out(rmsk ,il,jl,1,1.,idnc,"lsm","lsm","none",-1.,1.)
      call nc2out(tsd  ,il,jl,1,1.,idnc,"tsd","tsd","m",0.,30000.)
      call nc2out(tmax ,il,jl,1,1.,idnc,"zmax","zmax","m",-100.,30000.)
      call nc2out(tmin ,il,jl,1,1.,idnc,"zmin","zmin","m",-100.,30000.)

c     write out g*zs(il,jl) to formatted file
      do j=1,jl
        do i=1,il
          zss(i,j)=9.80616*zss(i,j)
        end do ! i=1,il
      end do ! j=1,jl
      ilout=min(il,30)
      write (formout,'(''(''i3''f7.0)'')')ilout   !  i.e. (<il>f7.0)
      write(luout,formout) zss

c     write out land/sea mask to formatted file
      write (formout,'(''(''i3''f4.1)'')')ilout   !  i.e. (<il>f4.1)
      write(luout,formout) rmsk

c     write out std.dev of top. to formatted file
      write (formout,'(''(''i3''f6.0)'')')ilout   !  i.e. (<il>f6.0)
      write(luout,formout) tsd
      print *,'zss, rmsk and tsd written to unit',luout,' for il=',il

      do j=1,jl
        do i=1,il
          zss(i,j)=real(inum(i,j))
        end do ! i=1,il
      end do ! j=1,jl
      call nc2out(zss,il,jl,1,1.,idnc,"inum","inum","none",0.,2000.)
      call ncsnc(idnc,ier)

      stop
      end
c=======================================================================
      subroutine read250(nx,ny,lons,lats,dl,debug,idia,jdia)

      integer*2 zz(3360*2800)
      real dl
      integer nrecl
      real lons,lats
      parameter ( zmin=-100 )

      logical debug,ok,ow

      data ow/.false./

      include 'newmpar.h'
      include 'rwork.h'

      write(6,*)"read250 nx,ny,lons,lats,dl=",nx,ny,lons,lats,dl,debug
      read(13) (zz(i),i=1,nx*ny)

!#######################################################################
      do jg=1,ny
!#######################################################################

        aglat=lats+real(jg-1)*dl
        if(ow)write(6,*)jg,aglat

        ok = (aglat.gt.rlatn.and.aglat.lt.rlatx)

        if(ok)then
          !write(6,*)"jg,lat,latn,latx,lons=",jg,aglat,rlatn,rlatx,lons
          if(ow)write(6,*)"zz=",(zz(n),n=1+(jg-1)*nx,nx+(jg-1)*nx)

!#######################################################################
         do ig=1,nx
!#######################################################################

           n=ig+(jg-1)*nx
           aglat=lats+real(jg-1)*dl
           aglon=lons+real(ig-1)*dl
           if(ow)write(6,*)ig,aglon

!-----------------------------------------------------------------------
           ok=(aglat.gt.rlatn.and.aglat.lt.rlatx)
           ok=ok .and. (aglon.gt.rlonn.and.aglon.lt.rlonx)
           if ( ok ) then
            if(ow) write(6,*)"ig,lon,lonn,lonx=",ig,aglon,rlonn,rlonx
!-----------------------------------------------------------------------

! compute model grid i,j
           call latltoij(aglon,aglat,alci,alcj,nface)  ! con-cubic/octagon
           lci = nint(alci)
           lcj = nint(alcj)
! convert to "double" (i,jg) notation
           lcj=lcj+nface*il
           if(ow)write(6,*)ig,jg,aglon,aglat,lci,lcj

! check to make sure within grid dimensions
           if(lci.gt.0.and.lci.le.il.and.lcj.gt.0.and.lcj.le.jl)then
             if( zz(n).lt.1 ) then
! ocean point
               amask = 0.
               zs=.0
             else  ! zs>=-1000
! land point
               amask = 1.
               zs=real(zz(n))
             end if  ! zs<-1000
! accumulate number of pnts in grid box
             inum(lci,lcj) = inum(lci,lcj) + 1
! accumulate topog. pnts
             zss(lci,lcj) = zss(lci,lcj) + zs
! accumulate lmask pnts
             almsk(lci,lcj) = almsk(lci,lcj) + amask
! find max/min topog. pnts in grid box
             tmax(lci,lcj)=max(tmax(lci,lcj),zs)
             tmin(lci,lcj)=min(tmin(lci,lcj),zs)
! sum of squares for sd calc.
             tsd(lci,lcj)=tsd(lci,lcj)+zs**2

             if ( debug ) then
                if ( lci.eq.idia .and. lcj.eq.jdia ) then
                 write(6,'("lci,j,ig,jg,zz,glo,gla,zss,lm, inum()="
     &   ,5i5,2f8.3,3f8.1,i6)') lci,lcj,ig,jg,zz(n),aglon,aglat
     &          ,zss(lci,lcj),almsk(lci,lcj), inum(lci,lcj)
                endif  ! selected points only
             endif  ! debug

           endif ! (lci.gt.0.and.lci.le.il.and.lcj.gt.0.and.lcj.le.jl)then

!-----------------------------------------------------------------------
           endif ! ( ok ) then
!-----------------------------------------------------------------------

!#######################################################################
         enddo ! i
        endif! (ok)then
!#######################################################################
        if ( mod(jg,100).eq.0 .and. ok ) then
          write(6,*)"ig,jg,lon,lat,lci,lcj=",ig,jg,aglon,aglat,lci,lcj
        endif
!#######################################################################
      enddo ! j
!#######################################################################

      close(21)

      return
      end ! read250
c=======================================================================
      subroutine read10km(debug)

      integer zz(7200)
      real dl
      real lons,lats
      parameter ( zmin=-100 )

      logical debug

      include 'newmpar.h'
      include 'rwork.h'

      nx=4320
      ny=2160
      lats=-90.
      lons=0.
      dl = 5./60.

!#######################################################################
      do j=1,ny
!#######################################################################
        aglat = lats + (j-1)*dl  ! 0.+...

        call read_ht ( aglat, zz )

!#######################################################################
        do i=1,nx
!#######################################################################

          aglon=(lons)+real(i-1)*dl

!-----------------------------------------------------------------------
          call latltoij(aglon,aglat,alci,alcj,nface)  ! con-cubic/octagon
          lci = nint(alci)
          lcj = nint(alcj)
! convert to "double" (i,j) notation
          lcj=lcj+nface*il
          !write(6,*)i,j,aglon,aglat,lci,lcj
          if ( inum(lci,lcj) .lt. 1 ) then

            if(lci.gt.0.and.lci.le.il.and.lcj.gt.0.and.lcj.le.jl)then
! all other points
! (fixup to 5 min data coastline)
             !if( zz(i).lt.-10 ) then
             if( zz(i).lt.0 ) then
! ocean point
               amask = 0.
               zs=0.
             else  ! zs>=-1000
! land point
               amask = 1.
               zs=real(zz(i))
             end if  ! zs<-1000
! accumulate topog. pnts
             zss(lci,lcj)=zss(lci,lcj)+zs
! accumulate lmask pnts
             almsk(lci,lcj) = almsk(lci,lcj) + amask
! accumulate number of  pnts in grid box
             inum(lci,lcj) = inum(lci,lcj) + 1
! find max/min topog. pnts in grid box
             tmax(lci,lcj)=max(tmax(lci,lcj),zs)
             tmin(lci,lcj)=min(tmin(lci,lcj),zs)
! sum of squares for sd calc.
             tsd(lci,lcj)=tsd(lci,lcj)+zs**2

             if ( debug ) then
                if ( lci.eq.idia .and. lcj.eq.jdia ) then
                  write(6,'("lci,lcj,i,zz(i),zs,zss(),almsk(), inum()="
     &           ,4i6,3f10.3,i6)') lci,lcj,i,zz(i),zs,zss(lci,lcj)
     &                       ,almsk(lci,lcj), inum(lci,lcj)
                endif  ! selected points only
             endif  ! debug

            endif ! (lci.gt.0.and.lci.le.il.and.lcj.gt.0.and.lcj.le.jl)then
          endif ! ( inum(i,j) .lt. 1 ) then

!-----------------------------------------------------------------------

!#######################################################################
        enddo ! i
!#######################################################################
        if (mod(j,100).eq.0) then
          write(6,*)"i,j,lon,lat,lci,lcj=",i,j,aglon,aglat,lci,lcj
         endif
!#######################################################################
      enddo ! j
!#######################################################################

      return
      end ! read10km
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
      subroutine ncdf_setup(ofile,idnc,ndim,il,jl,kl,kdate,ktime)

      include "netcdf.inc"   ! comment out on atmos

      character*(*) ofile
* netCDF id
      integer  idnc
      integer  ier
* dimension ids
      integer  dimil,dimjl,dimkl,dimtim
* variable ids
      integer idil,idjl,idkl,idnt
* input variables
      integer ndim,il,jl,kl
      character timorg*20
      character*3 month(12)
      data month/'jan','feb','mar','apr','may','jun'
     &          ,'jul','aug','sep','oct','nov','dec'/

      idnc = nccre ( ofile,ncclob,ier )
      print *,'###### create netcdf idnc,ndim,file=',idnc,ndim,ofile

c define dimensions
      print *,'define dimensions'
      dimil = ncddef(idnc,'longitude',il,ier)
      dimjl = ncddef(idnc,'latitude',jl,ier)
      if ( ndim.eq.4 ) then
        dimkl = ncddef(idnc,'dim3',kl,ier)
        dimtim = ncddef(idnc,'time',NCUNLIM,ier)
      elseif ( ndim.eq.3 ) then
        dimtim = ncddef(idnc,'time',NCUNLIM,ier)
      endif ! ( ndim.eq.4 ) then

c define variables
      print *,'define variables'
      idil = ncvdef(idnc,'longitude',NCFLOAT,1,dimil,ier)
      call ncaptc(idnc,idil,'point_spacing',NCCHAR,4,'even',ier)
      call ncaptc(idnc,idil,'units',NCCHAR,12,'degrees_east',ier)
      idjl = ncvdef(idnc,'latitude',NCFLOAT,1,dimjl,ier)
      call ncaptc(idnc,idjl,'point_spacing',NCCHAR,4,'even',ier)
      call ncaptc(idnc,idjl,'units',NCCHAR,13,'degrees_north',ier)

      if ( ndim.ge.3 ) then
        if ( ndim.eq.4 ) then
          idkl = ncvdef(idnc,'dim3',NCFLOAT,1,dimkl,ier)
          call ncaptc(idnc,idkl,'units',NCCHAR,5,'layer',ier)
          call ncaptc(idnc,idkl,'positive',NCCHAR,4,'down',ier)
        endif ! ( ndim.eq.4 ) then
        idnt = ncvdef(idnc,'time',NCFLOAT,1,dimtim,ier)
        call ncaptc(idnc,idnt,'units',NCCHAR,5,'hours',ier)
        call ncaptc(idnc,idnt,'point_spacing',NCCHAR,4,'even',ier)
        icy=kdate/10000
        icm=(kdate-icy*10000)/100
        icd=(kdate-icy*10000-icm*100)
        ich=ktime/100
        icmi=(ktime-ich*100)
        ics=0
        write(6,*) icy,icm,icd,ich,icmi,ics
        write(timorg,'(i2.2,"-",a3,"-",i4.4,3(":",i2.2))')
     &               icd,month(icm),icy,ich,icmi,ics
        print *,'timorg=',timorg
        call ncaptc(idnc,idnt,'time_origin',NCCHAR,20,timorg,ier)
      endif ! ( ndim.ge.3 ) then

      call ncendf(idnc,ier)

      return
      end
c=======================================================================
      subroutine nc2out(data,il,jl,nt,thr,idnc,sname,lname,units,dn,dx)

      logical debug
!     parameter ( debug=.false. )
      parameter ( debug=.true. )
      logical oscale
      parameter ( oscale=.false. )

c this program write out a 2-dim array of data with unlimited time dim.
c in netcdf format

      include "netcdf.inc"        ! comment out on atmos

      character*(*) sname
      character*(*) lname
      character*(*) units
* netCDF id
      integer  idnc
      integer  ier
* dimension ids
      integer  dimil,dimjl,dimkl,dimtim
* variable ids
      integer  idil,idjl,idkl,idnt
      integer  idvar
* variable shapes
      integer idims(3)
* corners and edge lengths
      integer corner(3),edges(3)
* data variables
      real data(il,jl)
* packing variables
      parameter (nx=4800,ny=6000)
      integer*2 ipack(nx*ny)
      integer*2 vrange(2)
      data vrange/-32500,32500/
* dim variables
      real xdim(nx),ydim(ny)

      save vrange

      write(6,'("*** nc2out *** idnc,nt,thr,sname=",2i4,f7.2,a6)')
     &          idnc,nt,thr,sname

      if ( debug ) then
        print *,'il,jl,nt,idnc=',il,jl,nt,idnc
        print *,'lname=',lname
        print *,'units=',units
        print *,'sx,sn=',dx,dn
      endif ! debug

      if ( il*jl .gt. nx*ny ) then
        print *,'il*jl=',il*jl,' too big in nc2out'
        stop
      endif

c***********************************************************************
      if ( nt.eq.1 ) then
c only do the following for the first call to each variable
c***********************************************************************

        call ncredf(idnc,ier)
        if ( debug ) print *,'ncredf ier=',ier

c get dimension ids
        dimil = ncdid(idnc,'longitude',ier)
        dimjl = ncdid(idnc,'latitude',ier)
        dimtim = ncdid(idnc,'time',ier)
        idims(1) = dimil
        idims(2) = dimjl
        idims(3) = dimtim
        if ( debug ) print *,'idims=',idims

c create variable using short name
        write(6,*)"thr=",thr
        if ( thr.ge.0. ) then
          ndims=3
        else
          ndims=2
        endif ! ( thr.ge.0. ) then
          if(oscale)then
            idvar = ncvdef (idnc,sname,ncshort,ndims,idims,ier)
          else
            idvar = ncvdef (idnc,sname,ncfloat,ndims,idims,ier)
          endif
        if(debug)write(6,*)'time idvar,sn,ndims=',idvar,sname,ndims

c give it a long name
        call ncaptc(idnc,idvar,'long_name',NCCHAR
     &             ,lngstr(lname),lname,ier)
        if ( debug ) print *,'lname ier=',ier

c give it units
        call ncaptc(idnc,idvar,'units',NCCHAR
     &             ,lngstr(units),units,ier)
        if ( debug ) print *,'units ier=',ier

        if(oscale)then
c give it a fill value
          call ncapt(idnc,idvar,'_FillValue',NCSHORT,1,vrange(1),ier)
          if ( debug ) print *,'fill ier=',ier
c define valid scaled max/min variables
          scalef= (dx-dn)/float(vrange(2)-vrange(1))
          addoff= dn-scalef*float(vrange(1))
          if ( debug ) then
            write(6,'("sx,sn,a,s=",1p,6e12.3)')dx,dn,addoff,scalef
          endif ! debug
          call ncapt(idnc,idvar,'add_offset',NCFLOAT,1,addoff,ier)
          call ncapt(idnc,idvar,'scale_factor',NCFLOAT,1,scalef,ier)
          call ncaptc(idnc,idvar,'FORTRAN_format',ncchar,5,'G11.4',ier)
          call ncapt(idnc,idvar,'valid_range',ncshort,2,vrange,ier)
        endif!(oscale)then

        if ( debug ) print *,'done attrib'

* leave define mode
        call ncendf(idnc,ier)

* store xdim
        do i=1,il
          xdim(i)=float(i)
        enddo
        idil = ncvid(idnc,'longitude',ier)
c       print *,'idil,xdim=',idil,(xdim(i),i=1,il)
        call ncvpt(idnc,idil,1,il,xdim,ier)

* store ydim
        do j=1,jl
          ydim(j)=float(j)
        enddo
        idjl = ncvid(idnc,'latitude',ier)
c       print *,'idjl,ydim=',idjl,(ydim(i),i=1,jl)
        call ncvpt(idnc,idjl,1,jl,ydim,ier)

c***********************************************************************
      endif ! nt=1
c***********************************************************************

      corner(1) = 1
      corner(2) = 1
      corner(3) = nt
      edges(1) = il
      edges(2) = jl
      edges(3) = 1
      if ( debug ) then
        print *,'corner=',corner
        print *,'edges=',edges
      endif ! debug

      write(6,*)"thr=",thr
      if ( thr.ge.0. ) then
c set time to number of hours since start = thr
        idnt = ncvid(idnc,'time',ier)
        if ( debug ) print *,'ncdiv idnt,ier=',idnt,ier
        call ncvpt1(idnc,idnt,nt,thr,ier)
        if ( debug ) print *,'ncvpt nt,ier=',nt,ier
      endif ! ( thr.ge.0. ) then

c find variable index
      idvar = ncvid(idnc,sname,ier)
      if ( debug ) then
        print *,'ncdiv sname,idvar,ier=',sname,idvar,ier
      endif ! debug

      if(oscale)then
c get scaling factors for this variable
        scalef= (dx-dn)/float(vrange(2)-vrange(1))
        addoff= dn-scalef*float(vrange(1))
        if ( debug ) write(6,'("sx,sn,a,s=",6f10.2)')dx,dn,addoff,scalef
c pack data into ipack
        adx=-1.e29
        adn= 1.e29
        idx= vrange(1)
        idn= vrange(2)
        do j=1,jl
          do i=1,il
           n=i+(j-1)*il
           adx=max(adx,data(i,j))
           adn=min(adn,data(i,j))
           ipack(n)=nint((data(i,j)-addoff)/scalef)
           idx=max(idx,(ipack(n)))
           idn=min(idn,(ipack(n)))
          end do ! i
        end do ! j
        if(adx.gt.dx)print *,'## actual dx .gt. set dx ##',adx,dx
        if(adn.lt.dn)print *,'## actual dn .lt. set dn ##',adn,dn
        if ( debug ) then
          print *,'adx,adn=',adx,adn
          print *,'idx,idn=',idx,idn
        endif ! debug
c put packed data into cdf file
        call ncvpt(idnc,idvar,corner,edges,ipack,ier)
      else!(oscale)then
        call ncvpt(idnc,idvar,corner,edges,data,ier)
      endif!(oscale)then
      if ( debug ) print *,'ier=',ier

      return
      end
c=======================================================================
      function lngstr( string )
      character*(*) string
      ilen = len(string)
c     print*,'string=',string
c     print*,'ilen=',ilen
      do 100 lngstr=ilen,1,-1
        if ( string(lngstr:lngstr) .ne. ' ' ) go to 99
  100 continue
      lngstr = 0
   99 continue
      return
      end
c=======================================================================
      include 'setxyz.f'
      include 'jimcc.f'
      include 'latltoij.f'
      include 'read_ht.f'
