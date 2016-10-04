!-----------------------------------------------------------------------
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!     Driver for statistics computation
      subroutine stat_avg_3D

      !implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'TSTEP_DEF'
      include 'TSTEP'
      include 'INPUT_DEF'
      include 'INPUT'
      include 'STATS'

      if (ISTEP.eq.0) then
!     Initialize 3D statistics

         call stat_init3d       ! It could all be incorporated in statistics_2D?

         STAT_ATIME = 0. 
         STAT_TSTART = time
         
      else

         if (mod(ISTEP,stat_comp).eq.0) call stat_compute3d

         
         if (mod(ISTEP,stat_outp).eq.0) then
            
            call stat_mfo_outfld3D ! Output
            STAT_ATIME = 0.
            STAT_TSTART = time
            
         endif
         
         
      endif
      
      return
      end subroutine stat_avg_3D

!======================================================================
      
!     main interface for statistics initialisation
      subroutine stat_init3d()
      !implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'STATS'           ! Statistics specific variables
      include 'INPUT_DEF'
      include 'INPUT'           ! if3d

      call rzero(STAT,lx1*ly1*lz1*lelt*NSTAT3D)
      call rzero(STAT_TEMP,lx1*ly1*lz1*lelt)

      return
      end
c----------------------------------------------------------------------

!     compute statistics
      subroutine stat_compute3d

      !implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'SOLN_DEF'
      include 'SOLN'
      include 'TSTEP_DEF'
      include 'TSTEP'
      include 'STATS'           ! 2D statistics speciffic variables
      include 'INPUT_DEF'
      include 'INPUT'           ! if3d

      real dtime,alpha,beta,surf_mean
      integer ntot,ierr
      
!     Work arrays
      real wk1(lx1*ly1*lz1), wk2(lx1*ly1*lz1)
      real slvel(LX1,LY1,LZ1,LELT,3), tmpvel(LX1,LY1,LZ1,LELT,3)
      real dudx(LX1,LY1,LZ1,LELT,3) ! du/dx, du/dy and du/dz
      real dvdx(LX1,LY1,LZ1,LELT,3) ! dv/dx, dv/dy and dv/dz
      real dwdx(LX1,LY1,LZ1,LELT,3) ! dw/dx, dw/dy and dw/dz
      
!     Array to store pressure interpolated to the velocity mesh
      real p0(lx1,ly1,lz1,lelt)

!     Array to store the reference pressure, i.e., the pressure
!     averaged over the walls
      real pmean(lx1,ly1,lz1,lelt)
      
!     Number of points per processor in each of the 3D fields
      ntot=lx1*ly1*lz1*lelt

c     test   write(*,*) 'ntot=',ntot
      
!     Calculate time span of current statistical sample
      dtime=time-STAT_TSTART

!     Update total time over which the current stat file is averaged      
      STAT_ATIME=STAT_ATIME+dtime
c     test write(*,*) 'ts,t,dt,ta=',STAT_TSTART,time,dtime,STAT_ATIME

!     Time average is compuated as:
!     Accumulated=alpha*Accumulated+beta*New
!     Calculate alpha and beta
      beta=dtime/STAT_ATIME
      alpha=1.0-beta

c     test  write(*,*) 'alpha,beta=',alpha,beta
      
!     Map pressure to velocity mesh
      call mappr(p0,PR,wk1,wk2)

c     test if (nid.eq.0) write(*,*) 'po',p0(1,1,1,1)
      
!     Calculate reference pressure as wall-averaged pressure
!     Note the minus sign
      pmean = -surf_mean(p0,1,'W  ',ierr)
      
c     test if (nid.eq.0) write(*,*) 'pmean',pmean(1,1,1,1)
      
!     Add p0 and pmean (with minus sing),
!     This subtracts the reference pressure
      call cadd(p0,pmean,ntot)
      
c      if (nid.eq.2) write(*,*) 'pfinal',p0(1,1,1,1)

c!     Compute derivative tensor
c      call comp_derivat
      
!     Compute derivative tensor
      call user_stat_trnsv(tmpvel,dudx,dvdx,dwdx,slvel)

c      if (nid.eq.2) write(*,*) 'pfinal2',p0(1,1,1,1)
      
      
!     Test routines for deriative calculation
      
c      write(*,*) 'nid,dudx=',nid,DUIDXJ(1,1,1)
c      write(*,*) 'nid,dudxx=',nid,dudx(1,1,1,1,1)

c      write(*,*) 'nid,dudy=',nid,DUIDXJ(1,1,4)
c      write(*,*) 'nid,dudyy=',nid,dudx(1,1,1,1,2)

c      write(*,*) 'nid,dudz=',nid,DUIDXJ(1,1,7)
c      write(*,*) 'nid,dudzz=',nid,dudx(1,1,1,1,3)

c      write(*,*) 'nid,dvdx=',nid,DUIDXJ(1,1,8)
c      write(*,*) 'nid,dvdxx=',nid,dvdx(1,1,1,1,1)

c      write(*,*) 'nid,dvdy=',nid,DUIDXJ(1,1,2)
c      write(*,*) 'nid,dvdyy=',nid,dvdx(1,1,1,1,2)

c      write(*,*) 'nid,dvdz=',nid,DUIDXJ(1,1,5)
c      write(*,*) 'nid,dvdzz=',nid,dvdx(1,1,1,1,3)

c      write(*,*) 'nid,dwdx=',nid,DUIDXJ(1,1,6)
c      write(*,*) 'nid,dwdxx=',nid,dwdx(1,1,1,1,1)

c      write(*,*) 'nid,dwdy=',nid,DUIDXJ(1,1,9)
c      write(*,*) 'nid,dwdyy=',nid,dwdx(1,1,1,1,2)

c      write(*,*) 'nid,dwdz=',nid,DUIDXJ(1,1,3)
c      write(*,*) 'nid,dwdzz=',nid,dwdx(1,1,1,1,3)
      

      
      
!     Compute 3D statistics

c-----------------------------------------------------------------------
!     <u>t
c      call avgt(STAT(1,1),vx,alpha,beta,ntot)
!     <v>t
c      call avgt(STAT(1,2),vy,alpha,beta,ntot)
!     <w>t
c      call avgt(STAT(1,3),vz,alpha,beta,ntot)
!     <p>t
c      call avgt(STAT(1,4),p0,alpha,beta,ntot)
      
c      if (nid.eq.5) write(*,*) 'STAT,vx',STAT(1,1),vx(1,1,1,1)
c      if (nid.eq.3) write(*,*) 'STAT,vy',STAT(1,2),vy(1,1,1,1)
c      if (nid.eq.4) write(*,*) 'STAT,vz',STAT(1,3),vz(1,1,1,1)
c      if (nid.eq.7) write(*,*) 'STAT,p0',STAT(1,4),p0(1,1,1,1)

c      if(nid.eq.2) write(*,*) 'vz,p0=',vz(1,1,1,1),p0(1,1,1,1)
      
c-----------------------------------------------------------------------
      
!     <u>t
      call add2sxy(STAT(1,1),alpha,vx,beta,ntot)
!     <v>t
      call add2sxy(STAT(1,2),alpha,vy,beta,ntot)
!     <w>t
      call add2sxy(STAT(1,3),alpha,vz,beta,ntot)
!     <p>t
      call add2sxy(STAT(1,4),alpha,p0,beta,ntot)

c      if (nid.eq.2) write(*,*) 'pfinal3',p0(1,1,1,1)
      
      
c     test if (nid.eq.5) write(*,*) 'STAT,vx',STAT(1,1),vx(1,1,1,1)
c     test if (nid.eq.3) write(*,*) 'STAT,vy',STAT(1,2),vy(1,1,1,1)
c     test if (nid.eq.4) write(*,*) 'STAT,vz',STAT(1,3),vz(1,1,1,1)
c      if (nid.eq.7) write(*,*) 'STAT,p0',STAT(1,4),p0(1,1,1,1)
      
c-----------------------------------------------------------------------
      
!     uu
      call prod2(STAT_TEMP,vx,vx,ntot)
!     <uu>t
      call add2sxy(STAT(1,5),alpha,STAT_TEMP,beta,ntot)

!     vv
      call prod2(STAT_TEMP,vy,vy,ntot)
!     <vv>t
      call add2sxy(STAT(1,6),alpha,STAT_TEMP,beta,ntot)

!     ww
      call prod2(STAT_TEMP,vz,vz,ntot)
!     <ww>t
      call add2sxy(STAT(1,7),alpha,STAT_TEMP,beta,ntot)

!     pp
      call prod2(STAT_TEMP,p0,p0,ntot)
!     <pp>t
      call add2sxy(STAT(1,8),alpha,STAT_TEMP,beta,ntot)
      

c      test   write(*,*) 'nid,prod2,vx=',nid,STAT_TEMP(1),vx(1,1,1,1)

c       test    write(*,*) 'nid,prod2,p=',nid,STAT_TEMP(1),p0(1,1,1,1)
      
      
c      if (nid.eq.5) write(*,*) 'uu=',STAT(1,5)
c      if (nid.eq.3) write(*,*) 'vv=',STAT(1,6)
c      if (nid.eq.4) write(*,*) 'ww=',STAT(1,7)
c      if (nid.eq.7) write(*,*) 'pp=',STAT(1,8)
      
      
c-----------------------------------------------------------------------

!     uv
      call prod2(STAT_TEMP,vx,vy,ntot)
!     <uv>t
      call add2sxy(STAT(1,9),alpha,STAT_TEMP,beta,ntot)

!     vw
      call prod2(STAT_TEMP,vy,vz,ntot)
!     <vw>t
      call add2sxy(STAT(1,10),alpha,STAT_TEMP,beta,ntot)

!     uw
      call prod2(STAT_TEMP,vx,vz,ntot)
!     <uw>t
      call add2sxy(STAT(1,11),alpha,STAT_TEMP,beta,ntot)

c      if (nid.eq.2) write(*,*) 'uw,u,w',STAT_TEMP(1),vx(1,1,1,1),
c     &     vz(1,1,1,1)
      
      
c-----------------------------------------------------------------------

!     pu
      call prod2(STAT_TEMP,p0,vx,ntot)
!     <pu>t
      call add2sxy(STAT(1,12),alpha,STAT_TEMP,beta,ntot)

c      if (nid.eq.2) write(*,*) 'pu,p0,u',STAT_TEMP(1),p0(1,1,1,1),
c     &     vx(1,1,1,1)
      
      
!     pv
      call prod2(STAT_TEMP,p0,vy,ntot)
!     <pv>t
      call add2sxy(STAT(1,13),alpha,STAT_TEMP,beta,ntot)
      
!     pw
      call prod2(STAT_TEMP,p0,vz,ntot)
!     <pw>t
      call add2sxy(STAT(1,14),alpha,STAT_TEMP,beta,ntot)

c-----------------------------------------------------------------------

!     pdudx
      call prod2(STAT_TEMP,p0,dudx(1,1,1,1,1),ntot)
!     <pdudx>t
      call add2sxy(STAT(1,15),alpha,STAT_TEMP,beta,ntot)
      
c      if(nid.eq.2) write(*,*) 'prod,po,dudx=',STAT_TEMP(1),p0(1,1,1,1),
c     &     dudx(1,1,1,1,1)

!     pdudy
      call prod2(STAT_TEMP,p0,dudx(1,1,1,1,2),ntot)
!     <pdudy>t
      call add2sxy(STAT(1,16),alpha,STAT_TEMP,beta,ntot)

!     pdudz
      call prod2(STAT_TEMP,p0,dudx(1,1,1,1,3),ntot)
!     <pdudz>t
      call add2sxy(STAT(1,17),alpha,STAT_TEMP,beta,ntot)
      
c-----------------------------------------------------------------------

!     pdvdx
      call prod2(STAT_TEMP,p0,dvdx(1,1,1,1,1),ntot)
!     <pdvdx>t
      call add2sxy(STAT(1,18),alpha,STAT_TEMP,beta,ntot)

!     pdvdy
      call prod2(STAT_TEMP,p0,dvdx(1,1,1,1,2),ntot)
!     <pdvdy>t
      call add2sxy(STAT(1,19),alpha,STAT_TEMP,beta,ntot)

!     pdvdz
      call prod2(STAT_TEMP,p0,dvdx(1,1,1,1,3),ntot)
!     <pdvdz>t
      call add2sxy(STAT(1,20),alpha,STAT_TEMP,beta,ntot)

c-----------------------------------------------------------------------
      
!     pdwdx
      call prod2(STAT_TEMP,p0,dwdx(1,1,1,1,1),ntot)
!     <pdwdx>t
      call add2sxy(STAT(1,21),alpha,STAT_TEMP,beta,ntot)

!     pdwdy
      call prod2(STAT_TEMP,p0,dwdx(1,1,1,1,2),ntot)
!     <pdwdy>t
      call add2sxy(STAT(1,22),alpha,STAT_TEMP,beta,ntot)

!     pdwdz
      call prod2(STAT_TEMP,p0,dwdx(1,1,1,1,3),ntot)
!     <pdwdz>t
      call add2sxy(STAT(1,23),alpha,STAT_TEMP,beta,ntot)

c-----------------------------------------------------------------------
      
!     uuu
      call prod3(STAT_TEMP,vx,vx,vx,ntot)
!     <uuu>t
      call add2sxy(STAT(1,24),alpha,STAT_TEMP,beta,ntot)

c      if (nid.eq.2) write(*,*) 'stat,uuu,u=',STAT(1,24),STAT_TEMP(1),
c     &     vx(1,1,1,1)

!     vvv
      call prod3(STAT_TEMP,vy,vy,vy,ntot)
!     <vvv>t
      call add2sxy(STAT(1,25),alpha,STAT_TEMP,beta,ntot)

!     www
      call prod3(STAT_TEMP,vz,vz,vz,ntot)
!     <www>t
      call add2sxy(STAT(1,26),alpha,STAT_TEMP,beta,ntot)

!     ppp
      call prod3(STAT_TEMP,p0,p0,p0,ntot)
!     <ppp>t
      call add2sxy(STAT(1,27),alpha,STAT_TEMP,beta,ntot)

c-----------------------------------------------------------------------      




      
c     End of the routine is here
      
!     Assign starting time of next statistical sample
!     It does not depend  on whether we stored a stat file or not
      STAT_TSTART=time
      
      
      return
      end

c----------------------------------------------------------------------

!     Outpost of time-averaged fields
      subroutine stat_mfo_outfld3D()
      !implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'STATS'           ! Statistics specific variables
      include 'INPUT_DEF'
      include 'INPUT'           ! if3d
      include 'SOLN_DEF'
      include 'SOLN'
            
!     Fields to outpost: <u>t, <v>t, <w>t, <p>t, <uu>t
      call outpost(STAT(1,1),STAT(1,2),STAT(1,3),STAT(1,4),STAT(1,5),
     &     'st1')

!     Fields to outpost: <vv>t, <ww>t, <pp>t, <uv>t, <vw>t
      call outpost(STAT(1,6),STAT(1,7),STAT(1,8),STAT(1,9),STAT(1,10),
     &     'st2')

!     Fields to outpost: <uw>t, <pu>t, <pv>t, <pw>t, <pdudx>t
      call outpost(STAT(1,11),STAT(1,12),STAT(1,13),STAT(1,14),
     &     STAT(1,15),'st3')
      


      

      
      
      
      


!     Fields to outpost: <uu>t, <vv>t, <ww>t, <pp>t
c      call outpost(STAT(1,1),STAT(1,2),STAT(1,3),STAT(1,4),STAT(1,5)
c     &     ,'rms')

!     Fields to outpost: <uv>t, <vw>t, <uw>t
c      call outpost(STAT(1,9),STAT(1,10),STAT(1,11),vx,vy,'rm2')      

!     Fields to outpost: <pu>t, <pv>t, <pw>t
c      call outpost(STAT(1,12),STAT(1,13),STAT(1,14),vx,vy,'pui')
      
      
      return
      end
      
c----------------------------------------------------------------------
      
!     Time accumulation with alpha and beta
c      subroutine avgt(avg,f,alpha,beta,n)
c      !implicit none

c      include 'SIZE_DEF'
c      include 'SIZE'            !
c      include 'STATS'

c      integer k,n
c      real avg(n),f(n)
c      real alpha, beta
            
c      do k=1,n
c         avg(k) = alpha*avg(k)+beta*f(k)
c      enddo

c      return
c      end
!-----------------------------------------------------------------------

      subroutine prod2(avg,f,g,n)
      !implicit none

      integer k,n
      real avg(n),f(n),g(n)

      do k=1,n
         avg(k) = f(k)*g(k)
      enddo

      return
      end
      
c-----------------------------------------------------------------------

      subroutine prod3(avg,f,g,h,n)
      !implicit none

      integer k,n
      real avg(n),f(n),g(n),h(n)

      do k=1,n
         avg(k) = f(k)*g(k)*h(k)
      enddo

      return
      end

c-----------------------------------------------------------------------
      
      
c      subroutine comp_derivat
c      !implicit none
      
c      include 'SIZE_DEF'
c      include 'SIZE'
c      include 'STATS'           ! Statistics specific variables
c      include 'INPUT_DEF'
c      include 'INPUT'           ! if3d
c      include 'SOLN_DEF'
c      include 'SOLN'
c      include 'GEOM_DEF'
c      include 'GEOM'
c      include 'DXYZ_DEF'
c      include 'DXYZ'
      
c      integer e,n,nxyz,k

c      n    = nx1-1              ! Polynomial degree
c      nxyz = nx1*ny1*nz1

c      do e=1,nelv
c         call local_grad3(ur,us,ut,vx,N,e,dxm1,dxtm1)
c         call local_grad3(vr,vs,vt,vy,N,e,dxm1,dxtm1)
c         call local_grad3(wr,ws,wt,vz,N,e,dxm1,dxtm1)


c!     Derivative tensor computed by using the inverse of
c!     the Jacobian array jacmi
c         do k=1,nxyz
c            duidxj(k,e,1) = jacmi(k,e)*(ur(k)*rxm1(k,1,1,e)+
c     $           us(k)*sxm1(k,1,1,e)+
c     $           ut(k)*txm1(k,1,1,e))
c            duidxj(k,e,2) = jacmi(k,e)*(vr(k)*rym1(k,1,1,e)+
c     $           vs(k)*sym1(k,1,1,e)+
c     $           vt(k)*tym1(k,1,1,e))
c            duidxj(k,e,3) = jacmi(k,e)*(wr(k)*rzm1(k,1,1,e)+
c     $           ws(k)*szm1(k,1,1,e)+
c     $           wt(k)*tzm1(k,1,1,e))
c            duidxj(k,e,4) = jacmi(k,e)*(ur(k)*rym1(k,1,1,e)+
c     $           us(k)*sym1(k,1,1,e)+
c     $           ut(k)*tym1(k,1,1,e))
c            duidxj(k,e,5) = jacmi(k,e)*(vr(k)*rzm1(k,1,1,e)+
c     $           vs(k)*szm1(k,1,1,e)+
c     $           vt(k)*tzm1(k,1,1,e))
c            duidxj(k,e,6) = jacmi(k,e)*(wr(k)*rxm1(k,1,1,e)+
c     $           ws(k)*sxm1(k,1,1,e)+
c     $           wt(k)*txm1(k,1,1,e))
c            duidxj(k,e,7) = jacmi(k,e)*(ur(k)*rzm1(k,1,1,e)+
c     $           us(k)*szm1(k,1,1,e)+
c     $           ut(k)*tzm1(k,1,1,e))
c            duidxj(k,e,8) = jacmi(k,e)*(vr(k)*rxm1(k,1,1,e)+
c     $           vs(k)*sxm1(k,1,1,e)+
c     $           vt(k)*txm1(k,1,1,e))
c            duidxj(k,e,9) = jacmi(k,e)*(wr(k)*rym1(k,1,1,e)+
c     $           ws(k)*sym1(k,1,1,e)+
c     $           wt(k)*tym1(k,1,1,e))
c         enddo
c      enddo
c
c      return
c      end
c
c-----------------------------------------------------------------------
      
