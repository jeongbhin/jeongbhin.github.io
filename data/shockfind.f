      subroutine shockfind
c-----------------------------------------------------------------------
c  SHOCK DETECTOR using Entropy & Temperature & 3D div(u)
c        phi  = rho2 = xi*rrhomin
c        phi2 = T2   = tmax      
c        phi3 = vshock = Ms * csound     
c-----------------------------------------------------------------------
      include 'common.param'
      include 'common.WENO'
      include 'common.MHD'
c
c Local variables
      real*8 temp1(nx,nx,nx), tempb(nx,nx,nx), tempc(nx,nx,nx) 
      real*8 tempd(nx,nx,nx)
      real*8 bx,by,bz,EE,pg,vv2,bb2, pgjc, mcut, mjump
      real*8 pgr1,pgr2,pgjmp,xir,pi2
c
c 1D array for shock detection
      real*8 rho(nx),tp(nx),divu(nx),sjmp(nx),dtds(nx),pgr(nx)
      real*8 tem4, shkjmp, shockd, dum
      real*8 rrhomin, tmin, rhomax, tmax, xi
      real*8 amach, rtemp, bf, csound, vshock
      real*8 phi, phi2, phi3, ammax, phimax, phimax2, phimax3
c
      integer ishock(0:nx+1) 
      integer im1, ip1, iss, iq, iq1, ish, ish1, ish2
      integer i1, i2, i3, ifinal, ik, iknew, iflag, nflag, nsh
      integer i,j,k,ix,iy,iz, in, jn, kn, iw, coms
      integer iw1m,iw1p,iw2m,iw2p,iw3m,iw3p
c
c set tem4 = Minimum temperature for IGM
c temperature in Kelvin
      tem4 = 0.d0!1.0d04
      temp0 = 1.d0
c
      coms = 0 ! 0 for P jump, 1 for T jump

      pi2 = 4.D0*DATAN(1.D0)*2.d0

      mcut = 1.2d0  ! mach number cut
      mjump = 1.1d0 ! mach number for shock finding

c     T ratio for mjump
      shkjmp = (1.d0+(gam-1.d0)/2.d0*mjump**2)*(2.d0*gam/(gam-1.d0)*
     +          mjump**2-1.d0)/(mjump**2*(2.d0*gam/(gam-1.d0)+
     +          (gam-1.d0)/2.d0))
c     T jump condition
      shkjmp = 2.d0*(shkjmp-1.d0)/(shkjmp+1.d0)

c     pg ratio for mjump
      pgjc = 2.d0*gam*mjump**2/(gam+1.d0)-(gam-1.d0)/(gam+1.d0)

c     pg jump condition
      pgjc = 2.d0*(pgjc - 1.d0)/(pgjc + 1.d0)

c=======================================================================
c save varialbes at t^n and t^{n=1} for Shock source and CR transport 
c-----------------------------------------------------------------------
c denp, denn, vxp, vxn, bmagp, bmagn, temp, temn
c save

!$OMP PARALLEL DO DEFAULT(SHARED)
!$OMP+PRIVATE(ix,iy,iz)
      do iz=-1,nz+2
      do iy=-1,ny+2
      do ix=-1,nx+2
         denp(ix,iy,iz) =  denn(ix,iy,iz)
         vxp(ix,iy,iz)  = vxn(ix,iy,iz)
         vyp(ix,iy,iz)  = vyn(ix,iy,iz) 
         vzp(ix,iy,iz)  = vzn(ix,iy,iz)
         temp(ix,iy,iz) = temn(ix,iy,iz) 
         bmagp(ix,iy,iz) = bmagn(ix,iy,iz)
         divvp(ix,iy,iz) = divvn(ix,iy,iz)
         ishockp(ix,iy,iz) = ishockn(ix,iy,iz)
      enddo
      enddo
      enddo
!$OMP END PARALLEL DO

c
      nshp = nshn
      do i = 1, nshp
       ishxp(i)=ishxn(i)
       ishyp(i)=ishyn(i)
       ishzp(i)=ishzn(i)
       Msp(i) = Msn(i) 
       nep(i) = nen(i) 
       Tep(i) = Ten(i) 
       Vsp(i) = Vsn(i)
       Bfp(i) = Bfn(i)
      enddo
c
       iw=5
!$OMP PARALLEL DO DEFAULT(SHARED)
!$OMP+PRIVATE(ix,iy,iz,bx,by,bz,EE,BB2,vv2,pg)
      do iz=-1,nz+2
      do iy=-1,ny+2
      do ix=-1,nx+2
         denn(ix,iy,iz) = q(iw,1,ix,iy,iz)
         vxn(ix,iy,iz) = q(iw,2,ix,iy,iz)/denn(ix,iy,iz)
         vyn(ix,iy,iz) = q(iw,3,ix,iy,iz)/denn(ix,iy,iz)
         vzn(ix,iy,iz) = q(iw,4,ix,iy,iz)/denn(ix,iy,iz)
         bx  = q(iw,5,ix,iy,iz)
         by  = q(iw,6,ix,iy,iz)
         bz  = q(iw,7,ix,iy,iz)
         EE  = q(iw,8,ix,iy,iz)
         BB2 = bx**2+by**2+bz**2
         bmagn(ix,iy,iz) = sqrt(BB2)
         vv2 = vxn(ix,iy,iz)**2 + vyn(ix,iy,iz)**2 + vzn(ix,iy,iz)**2
         pg = (gam-1.d0)*(EE-0.5d0*(denn(ix,iy,iz)*vv2+BB2))
         temn(ix,iy,iz) = pg/denn(ix,iy,iz)*temp0
      enddo
      enddo
      enddo
!$OMP END PARALLEL DO

c
!$OMP PARALLEL DO DEFAULT(SHARED)
!$OMP+PRIVATE(ix,iy,iz)
      do iz=-1,nz+2
      do iy=-1,ny+2
      do ix=-1,nx+2

      if(ix.ge.1.and.iy.ge.1.and.iz.ge.1.and.
     >   ix.le.nx.and.iy.le.ny.and.iz.le.nz) then
       divvn(ix,iy,iz)=
     >   (vxn(ix+1,iy,iz)-vxn(ix-1,iy,iz))/dx/2.d0
     >  +(vyn(ix,iy+1,iz)-vyn(ix,iy-1,iz))/dy/2.d0
     >  +(vzn(ix,iy,iz+1)-vzn(ix,iy,iz-1))/dz/2.d0

        temp1(ix,iy,iz) = 0.0d0
        tempb(ix,iy,iz) = 0.0d0
        tempc(ix,iy,iz) = 0.0d0
        tempd(ix,iy,iz) = 0.0d0
      endif

        ishockn(ix,iy,iz)=0
      enddo
      enddo
      enddo
!$OMP END PARALLEL DO


c-----------------------------------------------------------------------
c Find Shock Quantities at t^{n+1}
c----- Z Pass ---------------------------------------------
C   Detecting shocks from jump condition

!$OMP PARALLEL DO DEFAULT(SHARED)
!$OMP+PRIVATE(ix,iy,iz,dum,shockd,im1,ip1,nsh,iss,ish1,ish2,ish,iq1,iq)
!$OMP+PRIVATE(i1,i2,rrhomin,rhomax,tmin,tmax,rtemp,bf,amach,xi)
!$OMP+PRIVATE(csound,vshock,phi,phi2,phi3,ishock,rho,tp,divu,dtds)
!$OMP+PRIVATE(pgr,pgr1,pgr2,pgjmp,xir,iw,iw1m,iw1p,iw2m,iw2p,iw3m,iw3p)

        DO ix = 1, nx
        DO iy = 1, ny
c
        do iz = 0, nz 
          ishock(iz)= 0
        enddo
        do iz= 1, nz 
         rho(iz)= denn(ix,iy,iz)
         tp(iz)  = temn(ix,iy,iz)
         divu(iz) = divvn(ix,iy,iz) 
         pgr(iz) = temn(ix,iy,iz)*denn(ix,iy,iz)/temp0
        enddo
c
        do iz = 1, nz 
         im1 = iz - 1
         ip1 = iz + 1
         if(im1.lt. 1) im1 = im1 + nz
         if(ip1.gt.nz) ip1 = ip1 - nz
         dum = abs( tp(ip1) - tp(im1) )

         sjmp(iz)= 2.0*dum / ( tp(ip1) + tp(im1) )
         shockd = max(divu(iz),shkjmp-sjmp(iz))

         pgjmp = 2.d0*abs(pgr(ip1)-pgr(im1))/(pgr(ip1)+pgr(im1))

         dtds(iz) = (log10(tp(ip1)) -log10(tp(im1)))
     +             *(log10(tp(ip1))-2./3.*log10(rho(ip1))
     +              -log10(tp(im1))+2./3.*log10(rho(im1)))

        IF ( coms.eq.1 )then
         if((shockd .lt.0.0d0).and.(dtds(iz).gt.0.0d0)) ishock(iz) = 1
        else
         if((max(divu(iz),pgjc-pgjmp).lt.0.d0)) ishock(iz) = 1
        endif

        enddo
        ishock(0) = ishock(nz)
        ishock(nz+1) = ishock(1)
c
        do iz = 1, nz 
         im1 = iz - 1
         ip1 = iz + 1
         if(im1.lt. 1) im1 = im1 + nz
         if(ip1.gt.nz) ip1 = ip1 - nz
         if(( ishock(iz).eq.0) .and. (ishock(im1).eq.1) .and.
     >      (ishock(ip1).eq.1)) then
         if(max(divu(iz),pgjc-pgjmp).lt.0.) ishock(iz)=1
         endif
        enddo
c
c find if div(i) is local maximum, => ishock=2
        do iz = 1, nz 
         if(ishock(iz).eq.1) then
         im1 = iz - 1
         ip1 = iz + 1
         if(im1.lt. 1) im1 = im1 + nz
         if(ip1.gt.nz) ip1 = ip1 - nz
            if( (abs(divu(iz)).ge.abs(divu(im1))).and.
     +          (abs(divu(iz)).ge.abs(divu(ip1)))) then
               ishock(iz) = 2
            endif
         endif
        enddo
c find the center of the shock

      IF ( coms.eq.0 )then
      do iw=1,nz
         if(ishock(iw).eq.1) then
            iw1m = iw-1
            iw1p = iw+1
            if(iw1m.lt. 1) iw1m = iw1m+nw
            if(iw1p.gt.nw) iw1p = iw1p-nw
c
            if((divu(iw).le.divu(iw1p)).and.(divu(iw).le.divu(iw1m)))
     +        ishock(iw) = 2
         endif
      enddo
c
      do iw=1,nw
         if(ishock(iw).eq.2) then
            iw2m = iw-2
            iw2p = iw+2
            if(iw2m.lt. 1) iw2m = iw2m+nw
            if(iw2p.gt.nw) iw2p = iw2p-nw
c
            if((ishock(iw2m).eq.2).and.(divu(iw2m).lt.divu(iw))) then
               ishock(iw) = 1
            endif
            if((ishock(iw2p).eq.2).and.(divu(iw2p).lt.divu(iw))) then
               ishock(iw) = 1
            endif
         endif
      enddo

      do iw=1,nw
         if(ishock(iw).eq.2) then
            iw3m = iw-3
            iw3p = iw+3
            if(iw3m.lt. 1) iw3m = iw3m+nw
            if(iw3p.gt.nw) iw3p = iw3p-nw
c
            if((ishock(iw3m).eq.2).and.(divu(iw3m).lt.divu(iw))) then
               ishock(iw) = 1
            endif
            if((ishock(iw3p).eq.2).and.(divu(iw3p).lt.divu(iw))) then
               ishock(iw) = 1
            endif
         endif
      enddo
      ENDIF


      if (coms.eq.1) then
      DO iz= 1, nz 
         if((ishock(iz-1).eq.0).and.(ishock(iz).eq.1)) then
            nsh = 1
            do iss= 1,19
              iq1 = iz + iss
               if (ishock(iq1).gt.0) then
                  nsh = nsh+1
               else
                  goto 111
               endif
            enddo
         else
            goto 112
         endif
 111     continue
         do iq1= iz, iz+nsh-1
            iq = iq1
            if(iq.gt.nz) iq = iq - nz
            if (ishock(iq1).eq.2) goto 112
         enddo
         if (mod(nsh,2).eq.1) then
            ish = iz + nsh/2
            if(ish.gt.nz) ish = ish-nz
            ishock(ish) = 2
         else
            ish1 = iz + nsh/2-1
            ish2 = iz + nsh/2
            if (abs(divu(ish1)).gt.abs(divu(ish2))) then
               ishock(ish1) = 2
            else
               ishock(ish2) = 2
            endif
         endif
 112    continue
      ENDDO 
      endif
c     shock center if ishock(iz)=2
c
        do iz= 1, nz 
        if(ishock(iz) .eq.2 )then
         ish = iz
c
         i1 = ish - 2
         if( ishock(ish-1). eq. 0) i1 = ish-1
         i2 = ish + 2
         if( ishock(ish+1). eq. 0) i2 = ish+1
         if(i1.lt.1)  i1 = i1 + nz
         if(i2.gt.nz) i2 = i2 - nz

         rrhomin=min(rho(i1), rho(i2))
         rhomax=max(rho(i1), rho(i2))
         tmin = min(tp(i1), tp(i2))
         tmin = max(tmin, tem4)
         tmax = max(tp(i1), tp(i2))
         tmax = max(tmax, tem4)
c 
c temperature, sound speed, shock seed in real units
         rtemp = tmax/tmin
         bf = 4.0d0* rtemp -3.5d0
         amach = (bf + sqrt(bf*bf+3.75d0))/2.5d0
         amach = sqrt(amach)

         if (coms.eq.0) then
         pgr1 = min(pgr(i1),pgr(i2))
         pgr2 = max(pgr(i1),pgr(i2))

         xir = pgr2/pgr1
         amach =  sqrt((gam+1.d0)/(2.d0*gam)
     +        *(xir + (gam-1.d0)/(gam+1.d0)))

         endif
        
         xi = 4.0d0*amach**2/(amach**2 + 3.0d0)
         csound = sqrt( gam* tmin * tem0)
         vshock = amach * csound
c
         phi = rrhomin*xi
         phi2 = tmax 
         phi3 = vshock
c
         if(amach.le.mcut) then 
           amach=0.0d0
           phi  =0.0d0
           phi2 =0.0d0
           phi3 =0.0d0
         endif
         temp1(ix,iy,iz) = amach
         tempb(ix,iy,iz) = phi
         tempc(ix,iy,iz) = phi2
         tempd(ix,iy,iz) = phi3

        endif
        enddo
c
        ENDDO 
        ENDDO 
c       print *, 'z-pass done'
!$OMP END PARALLEL DO
c -----------------------------------------------------
c----- X Pass ---------------------------------------------
c
!$OMP PARALLEL DO DEFAULT(SHARED)
!$OMP+PRIVATE(ix,iy,iz,dum,shockd,im1,ip1,nsh,iss,ish1,ish2,ish,iq1,iq)
!$OMP+PRIVATE(i1,i2,rrhomin,rhomax,tmin,tmax,rtemp,bf,amach,xi)
!$OMP+PRIVATE(csound,vshock,phi,phi2,phi3,ishock,rho,tp,divu,dtds)
!$OMP+PRIVATE(ifinal,ammax,phimax,phimax2,phimax3,iknew)
!$OMP+PRIVATE(pgr,pgr1,pgr2,pgjmp,xir,iw,iw1m,iw1p,iw2m,iw2p,iw3m,iw3p)
        DO iz = 1, nz 
        DO iy = 1, ny
c
        do ix = 0, nx 
          ishock(ix)= 0
        enddo
        do ix= 1, nx
         rho(ix)= denn(ix,iy,iz)
         tp(ix)  = temn(ix,iy,iz)
         divu(ix) = divvn(ix,iy,iz)
         pgr(ix) = temn(ix,iy,iz)*denn(ix,iy,iz)/temp0
        enddo
        do ix = 1, nx
         im1 = ix - 1
         ip1 = ix + 1
         if(im1.lt. 1) im1 = im1 + nx
         if(ip1.gt.nx) ip1 = ip1 - nx
         dum = abs( tp(ip1) - tp(im1) )
         sjmp(ix)= 2.0*dum / ( tp(ip1) + tp(im1) )
         shockd = max(divu(ix),shkjmp-sjmp(ix))
         dtds(ix) = (log10(tp(ip1)) -log10(tp(im1)))
     +             *(log10(tp(ip1))-2./3.*log10(rho(ip1))
     +              -log10(tp(im1))+2./3.*log10(rho(im1)))


         pgjmp = 2.d0*abs(pgr(ip1)-pgr(im1))/(pgr(ip1)+pgr(im1))
        IF ( coms.eq.1 )then
         if((shockd .lt.0.0d0).and.(dtds(ix).gt.0.0d0)) ishock(ix) = 1
        else
         if((max(divu(ix),pgjc-pgjmp).lt.0.d0)) ishock(ix) = 1
        endif
        enddo
        ishock(0) = ishock(nx)
        ishock(nx+1) = ishock(1)
        do ix = 1, nx 
         im1 = ix - 1
         ip1 = ix + 1
         if(im1.lt. 1) im1 = im1 + nx
         if(ip1.gt.nx) ip1 = ip1 - nx
         if((ishock(ix).eq.0).and.(ishock(im1).eq.1).and.
     >      (ishock(ip1).eq.1)) then
         if(max(divu(ix),pgjc-pgjmp).lt.0.) ishock(ix)=1
         endif
        enddo
c
c find if div(i) is local maximum
        do ix = 1, nx 
         if(ishock(ix).eq.1) then
         im1 = ix - 1
         ip1 = ix + 1
         if(im1.lt. 1) im1 = im1 + nx
         if(ip1.gt.nx) ip1 = ip1 - nx
            if( (abs(divu(ix)).ge.abs(divu(im1))).and.
     +          (abs(divu(ix)).ge.abs(divu(ip1)))) then
               ishock(ix) = 2
            endif
         endif
        enddo

      IF ( coms.eq.0 )then
      do iw=1,nx
         if(ishock(iw).eq.1) then
            iw1m = iw-1
            iw1p = iw+1
            if(iw1m.lt. 1) iw1m = iw1m+nw
            if(iw1p.gt.nw) iw1p = iw1p-nw
c
            if((divu(iw).le.divu(iw1p)).and.(divu(iw).le.divu(iw1m)))
     +        ishock(iw) = 2
         endif
      enddo
c
      do iw=1,nw
         if(ishock(iw).eq.2) then
            iw2m = iw-2
            iw2p = iw+2
            if(iw2m.lt. 1) iw2m = iw2m+nw
            if(iw2p.gt.nw) iw2p = iw2p-nw
c
            if((ishock(iw2m).eq.2).and.(divu(iw2m).lt.divu(iw))) then
               ishock(iw) = 1
            endif
            if((ishock(iw2p).eq.2).and.(divu(iw2p).lt.divu(iw))) then
               ishock(iw) = 1
            endif
         endif
      enddo

      do iw=1,nw
         if(ishock(iw).eq.2) then
            iw3m = iw-3
            iw3p = iw+3
            if(iw3m.lt. 1) iw3m = iw3m+nw
            if(iw3p.gt.nw) iw3p = iw3p-nw
c
            if((ishock(iw3m).eq.2).and.(divu(iw3m).lt.divu(iw))) then
               ishock(iw) = 1
            endif
            if((ishock(iw3p).eq.2).and.(divu(iw3p).lt.divu(iw))) then
               ishock(iw) = 1
            endif
         endif
      enddo
      ENDIF

      if (coms.eq.1) then
c find the center of the shock
      DO ix= 1, nx 
         if((ishock(ix-1).eq.0).and.(ishock(ix).eq.1)) then
            nsh = 1
            do iss= 1, 19
             iq1 = ix + iss 
             if(iq1.gt.nx) iq1 = iq1-nx
               if (ishock(iq1).gt.0) then
                  nsh = nsh+1
               else
                  goto 211
               endif
            enddo
         else
            goto 212
         endif
 211     continue
         do iq1= ix, ix+nsh-1
            iq = iq1
            if(iq.gt.nx) iq = iq - nx
            if (ishock(iq).eq.2) goto 212
         enddo
         if (mod(nsh,2).eq.1) then
            ish = ix +nsh/2
            if(ish.gt.nx) ish = ish-nx
            ishock(ish) = 2
         else
            ish1 = ix + nsh/2-1
            ish2 = ix + nsh/2
            if(ish1.gt.nx) ish1 = ish1 - nx
            if(ish2.gt.nx) ish2 = ish2 - nx
            if (abs(divu(ish1)).gt.abs(divu(ish2))) then
               ishock(ish1) = 2
            else
               ishock(ish2) = 2
            endif
         endif
 212   continue
      ENDDO
      endif
c
        do ix= 1, nx
        if(ishock(ix) .eq.2 )then
         ish = ix 
         i1 = ish - 2
         if( ishock(ish-1). eq. 0) i1 = ish-1
         i2 = ish + 2
         if( ishock(ish+1). eq. 0) i2 = ish+1
         if(i1.lt.1)  i1 = i1 + nx
         if(i2.gt.nx) i2 = i2 - nx
         rrhomin=min(rho(i1), rho(i2))
         rhomax=max(rho(i1), rho(i2))
         tmin = min(tp(i1), tp(i2))
         tmin = max(tmin, tem4)
         tmax = max(tp(i1), tp(i2))
         tmax = max(tmax, tem4)
c
         rtemp = tmax/tmin
         bf = 4.0d0* rtemp -3.5d0
         amach = (bf + sqrt(bf*bf+3.75d0))/2.5d0
         amach = sqrt(amach)

         if (coms.eq.0) then
         pgr1 = min(pgr(i1),pgr(i2))
         pgr2 = max(pgr(i1),pgr(i2))

         xir = pgr2/pgr1
         amach =  sqrt((gam+1.d0)/(2.d0*gam)
     +        *(xir + (gam-1.d0)/(gam+1.d0)))

         endif

         xi = 4.0d0*amach**2/(amach**2 + 3.0d0)
         csound = sqrt( gam* tmin * tem0)
         vshock = amach * csound
         phi = rrhomin*xi 
         phi2 = tmax 
         phi3 = vshock
c        print*, rtemp,bf

         if(amach.le.mcut) then 
           amach= 0.0d0
           phi  = 0.0d0
           phi2 = 0.0d0
           phi3 = 0.0d0
         endif
c
         temp1(ix,iy,iz) = max(temp1(ix,iy,iz),amach)
         tempb(ix,iy,iz) = max(tempb(ix,iy,iz),phi) 
         tempc(ix,iy,iz) = max(tempc(ix,iy,iz),phi2) 
         tempd(ix,iy,iz) = max(tempd(ix,iy,iz),phi3) 
        endif
        enddo
c
        ENDDO 
        ENDDO 
c       print *, 'x-pass done'
!$OMP END PARALLEL DO
c----------------------------------------------------------
c----- y Pass ---------------------------------------------
c
!$OMP PARALLEL DO DEFAULT(SHARED)
!$OMP+PRIVATE(ix,iy,iz,dum,shockd,im1,ip1,nsh,iss,ish1,ish2,ish,iq1,iq)
!$OMP+PRIVATE(i1,i2,rrhomin,rhomax,tmin,tmax,rtemp,bf,amach,xi)
!$OMP+PRIVATE(csound,vshock,phi,phi2,phi3,ishock,rho,tp,divu,dtds)
!$OMP+PRIVATE(ifinal,ammax,phimax,phimax2,phimax3,iknew)
!$OMP+PRIVATE(pgr,pgr1,pgr2,pgjmp,xir,iw,iw1m,iw1p,iw2m,iw2p,iw3m,iw3p)
        DO iz = 1, nz 
        DO ix = 1, nx
c
        do iy = 0, ny 
          ishock(iy)= 0
        enddo
        do iy= 1, ny
         rho(iy)= denn(ix,iy,iz)
         tp(iy)  = temn(ix,iy,iz)
         divu(iy) = divvn(ix,iy,iz)
         pgr(iy) = temn(ix,iy,iz)*denn(ix,iy,iz)/temp0
        enddo
c
        do iy = 1, ny
         im1 = iy - 1
         ip1 = iy + 1
         if(im1.lt. 1) im1 = im1 + ny
         if(ip1.gt.ny) ip1 = ip1 - ny
         dum = abs( tp(ip1) - tp(im1) )
         sjmp(iy)= 2.0*dum / ( tp(ip1) + tp(im1) )
         shockd = max(divu(iy),shkjmp-sjmp(iy))
         dtds(iy) = (log10(tp(ip1)) -log10(tp(im1)))
     +             *(log10(tp(ip1))-2./3.*log10(rho(ip1))
     +              -log10(tp(im1))+2./3.*log10(rho(im1)))

         pgjmp = 2.d0*abs(pgr(ip1)-pgr(im1))/(pgr(ip1)+pgr(im1))
        IF ( coms.eq.1 )then
         if((shockd .lt.0.0d0).and.(dtds(iy).gt.0.0d0)) ishock(iy) = 1
        else
         if((max(divu(iy),pgjc-pgjmp).lt.0.d0)) ishock(iy) = 1
        endif
        enddo
        ishock(0) = ishock(ny)
        ishock(ny+1) = ishock(1)
        do iy = 1, ny 
         im1 = iy - 1
         ip1 = iy + 1
         if(im1.lt. 1) im1 = im1 + ny
         if(ip1.gt.ny) ip1 = ip1 - ny
         if((ishock(iy).eq.0).and.(ishock(im1).eq.1).and.
     >      (ishock(ip1).eq.1)) then
         if(max(divu(iy),pgjc-pgjmp).lt.0.) ishock(iy)=1

         endif
        enddo
c find if div(i) is local maximum
        do iy = 1, ny 
         if(ishock(iy).eq.1) then
         im1 = iy - 1
         ip1 = iy + 1
         if(im1.lt. 1) im1 = im1 + ny
         if(ip1.gt.ny) ip1 = ip1 - ny
            if( (abs(divu(iy)).ge.abs(divu(im1))).and.
     +          (abs(divu(iy)).ge.abs(divu(ip1)))) then
               ishock(iy) = 2
            endif
         endif
        enddo

      IF ( coms.eq.0 )then
      do iw=1,ny
         if(ishock(iw).eq.1) then
            iw1m = iw-1
            iw1p = iw+1
            if(iw1m.lt. 1) iw1m = iw1m+nw
            if(iw1p.gt.nw) iw1p = iw1p-nw
c
            if((divu(iw).le.divu(iw1p)).and.(divu(iw).le.divu(iw1m)))
     +        ishock(iw) = 2
         endif
      enddo
c
      do iw=1,nw
         if(ishock(iw).eq.2) then
            iw2m = iw-2
            iw2p = iw+2
            if(iw2m.lt. 1) iw2m = iw2m+nw
            if(iw2p.gt.nw) iw2p = iw2p-nw
c
            if((ishock(iw2m).eq.2).and.(divu(iw2m).lt.divu(iw))) then
               ishock(iw) = 1
            endif
            if((ishock(iw2p).eq.2).and.(divu(iw2p).lt.divu(iw))) then
               ishock(iw) = 1
            endif
         endif
      enddo

      do iw=1,nw
         if(ishock(iw).eq.2) then
            iw3m = iw-3
            iw3p = iw+3
            if(iw3m.lt. 1) iw3m = iw3m+nw
            if(iw3p.gt.nw) iw3p = iw3p-nw
c
            if((ishock(iw3m).eq.2).and.(divu(iw3m).lt.divu(iw))) then
               ishock(iw) = 1
            endif
            if((ishock(iw3p).eq.2).and.(divu(iw3p).lt.divu(iw))) then
               ishock(iw) = 1
            endif
         endif
      enddo
      ENDIF

      if (coms.eq.1) then
c find the center of the shock
      DO iy= 1, ny 
         if((ishock(iy-1).eq.0).and.(ishock(iy).eq.1)) then
            nsh = 1
            do iss= 1, 19
             iq1 = iy + iss
             if(iq1.gt.ny) iq1 = iq1-ny
               if (ishock(iq1).gt.0) then
                  nsh = nsh+1
               else
                  goto 311
               endif
            enddo
         else
            goto 312
         endif
 311     continue
         do iq1= iy, iy+nsh-1
            iq = iq1
            if(iq.gt.ny) iq = iq - ny
            if (ishock(iq).eq.2) goto 312
         enddo
         if (mod(nsh,2).eq.1) then
            ish = iy +nsh/2
            if(ish.gt.ny) ish = ish-ny
            ishock(ish) = 2
         else
            ish1 = iy + nsh/2-1
            ish2 = iy + nsh/2
            if(ish1.gt.ny) ish1 = ish1 - ny
            if(ish2.gt.ny) ish2 = ish2 - ny
            if (abs(divu(ish1)).gt.abs(divu(ish2))) then
               ishock(ish1) = 2
            else
               ishock(ish2) = 2
            endif
         endif
 312   continue
      ENDDO
      endif

        do iy= 1, ny
        if(ishock(iy) .eq. 2 )then
         ish = iy 
         i1 = ish - 2
         if( ishock(ish-1). eq. 0) i1 = ish-1
         i2 = ish + 2
         if( ishock(ish+1). eq. 0) i2 = ish+1
         if(i1.lt.1)  i1 = i1 + ny
         if(i2.gt.ny) i2 = i2 - ny
         rrhomin=min(rho(i1), rho(i2))
         rhomax=max(rho(i1), rho(i2))
         tmin = min(tp(i1), tp(i2))
         tmin = max(tmin, tem4)
         tmax = max(tp(i1), tp(i2))
         tmax = max(tmax, tem4)

         rtemp = tmax/tmin
         bf = 4.0d0* rtemp -3.5d0
         amach = (bf + sqrt(bf*bf+3.75d0))/2.5d0
         amach = sqrt(amach)

         if (coms.eq.0) then
         pgr1 = min(pgr(i1),pgr(i2))
         pgr2 = max(pgr(i1),pgr(i2))

         xir = pgr2/pgr1
         amach =  sqrt((gam+1.d0)/(2.d0*gam)
     +        *(xir + (gam-1.d0)/(gam+1.d0)))

         endif

         xi = 4.0d0*amach**2/(amach**2 + 3.0d0)
         csound = sqrt( gam* tmin * tem0)
         vshock = amach * csound
         phi = rrhomin*xi 
         phi2 = tmax 
         phi3 = vshock
         if(amach.le.mcut) then 
           amach= 0.0d0
           phi  = 0.0d0
           phi2 = 0.0d0
           phi3 = 0.0d0
         endif
c
         temp1(ix,iy,iz) = max(temp1(ix,iy,iz),amach)
         tempb(ix,iy,iz) = max(tempb(ix,iy,iz),phi)
         tempc(ix,iy,iz) = max(tempc(ix,iy,iz),phi2)
         tempd(ix,iy,iz) = max(tempd(ix,iy,iz),phi3)
        endif
        enddo
        ENDDO 
        ENDDO 
c       print *, 'y-pass done'
!$OMP END PARALLEL DO

c detect isolated shocks and clean them up 	

!$OMP PARALLEL DO DEFAULT(SHARED)
!$OMP+PRIVATE(ix,iy,iz,nflag,kn,k,j,i,jn,in)
      do iz=1, nx
      do iy=1, ny
      do ix=1, nz
        if(temp1(ix,iy,iz).gt.0.0d0) then
          nflag=0
          do k = -2, 2
           kn = iz + k
            if(kn.lt.1) kn = kn + nz 
            if(kn.gt.nz) kn =  kn - nz 
           do j = -2, 2
            jn = iy + j
            if(jn.lt.1) jn = jn + ny
            if(jn.gt.ny) jn = jn - ny
            do i = -2, 2
             in = ix + i
             if(in.lt.1) in = in + nx
             if(in.gt.nx) in = in - nx
              if(temp1(in,jn,kn).gt.0.0d0) nflag=nflag+1
            enddo
           enddo
          enddo
          if(nflag.le.3) temp1(ix,iy,iz)=0.0d0
          if(nflag.le.3) tempb(ix,iy,iz)=0.0d0
          if(nflag.le.3) tempc(ix,iy,iz)=0.0d0
          if(nflag.le.3) tempd(ix,iy,iz)=0.0d0
        endif
      enddo
      enddo
      enddo
!$OMP END PARALLEL DO

c--------------------------------------------------------------------------
c
c store the indices and propoerties of shock zones with Ms> 1.5
      iflag=0

!$OMP PARALLEL DO DEFAULT(SHARED)
!$OMP+PRIVATE(ix,iy,iz)
      do iz=1,nz     
      do iy=1,ny
      do ix=1,nx
        temp_a(ix,iy,iz) = temp1(ix,iy,iz)
       if(temp1(ix,iy,iz) .gt. 1.5d0) then
        iflag = iflag + 1
        ishockn(ix,iy,iz)=1
        ishxn(iflag)=ix
        ishyn(iflag)=iy
        ishzn(iflag)=iz
       Msn(iflag)= temp1(ix,iy,iz)
       nen(iflag)= tempb(ix,iy,iz) 
       Ten(iflag)= tempc(ix,iy,iz) 
       Vsn(iflag)= tempd(ix,iy,iz)
       Bfn(iflag)= Bmagn(ix,iy,iz)
       endif
      enddo
      enddo
      enddo    
!$OMP END PARALLEL DO
      nshn = iflag 
c
c write the shock data 
c      do i = 1, nshn 
c       write(5,500) i,ishxn(i),ishyn(i),ishzn(i),Msn(i),nen(i),Ten(i),
c     > Vsn(i),Bfn(i)
c      enddo
c
500   format(4i8,1p5e15.6) 
c-----------------------------------------------------------------------
c  end
c-----------------------------------------------------------------------
      return
      end
