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

! read data then close file
      read(13) (zz(i),i=1,nx*ny)
      close(13)

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
           if(lci.gt.0.and.lci.le.il.and.lcj.gt.0.and.lcj.le.jl) then
	   
             if( zz(n).lt.1 ) then
! ocean point if zz < 1
               amask = 0.
               zs=0.
             else
! land point zz >= 1
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

      return
      end ! read250
