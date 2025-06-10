c     Summed Faldt+corrections 12/20/03
      Program STRONG
      implicit double precision (a-h,o-z)
      character*8 carg
      integer nargs
      integer Z
      DOUBLE PRECISION  x(3)
      common/xc/x
      common/nucl/r_n,c_n,wc,ro_0,A_n
      common/kine/Q,E,q_L 
      common/pars/ifunc
      external F3

* ===============================================================      
*                       Defaults
* ===============================================================      
      E=5.0D0
      q_L=1.4822d0/E
      teta=0.d0
      
      pi = dacos(-1.d0)
      degr  = pi/180.d0

      Z=82  ! 208Pb 

* ================================================================
*           command line argument operation
* ================================================================
      nargs=iargc()
      do iarg=1,nargs
         call GetArg(iarg, carg)
         if ( carg(1:2) .eq. '-h' ) then
            call ErrMsg()
            stop
         elseif ( carg(1:2) .eq. '-e' ) then
            call GetArg(iarg+1, carg)
            read(carg,*)E
         elseif ( carg(1:2) .eq. '-t' ) then
            call GetArg(iarg+1, carg)
            read(carg,*)teta
         elseif ( carg(1:2) .eq. '-Z' ) then
            call GetArg(iarg+1, carg)
            read(carg,*)Z
         endif
      enddo

      call LoadTarget(Z)

      if (teta.eq.0) teta=1.d-10  ! ugly trick to calculate at theta_pi=0

        Q=10.D0*E*DSIN(teta*degr/2.d0)

        a=0.D0
        b=30.D0
	ifunc=1
        FS=DGMLT3(F3,a,b,10,6,x)
	ifunc=2
        FM=DGMLT3(F3,a,b,10,6,x)
	ifunc=3
	FF=DGMLT3(F3,a,b,10,6,x)
	ifunc=4
        FR=DGMLT3(F3,a,b,10,6,x)
        ifunc=5
        FI=DGMLT3(F3,a,b,10,6,x)
	
	WS=-2.d0*pi*ro_0*FS/Q
	WM=2.d0*pi*ro_0*FM  
        WF=-2.d0*pi*A_n*1.3*ro_0**2*FF /Q 
	WR=-2.d0*pi*A_n*1.3*ro_0**2*FR/Q 
	WI=-2.d0*pi*A_n*1.3*ro_0**2*FI/Q 
c        WW= (0.5403d0*(WM+WF-WR)+0.8415d0*WI)**2+
c     &  (-0.5403d0*WI+0.8415d0*(WM+WF-WR))**2  
c	W1=(10.d0*E*A_nucl*DSIN(teta*degr))**2*WW
c        W2=(10.d0*E*A_nucl*DSIN(teta*degr)*(WM+WF))**2

        xs = const(E) * A_n * A_n * sin(teta*degr) * sin(teta*degr) 
        xs = xs * WS * WS

 111    format(f5.1,6(f12.4))
	write(*,111) teta,WS,WM,WF,WR,WI,xs

      end

*======================================================================
*     subroutine const(E)
*======================================================================
      function const(E)
      implicit double precision (a-h,o-z)

      bn  = 1.48d0
      L2  = 100 * E*E   ! [ub]
      
      const = bn * L2

      return
      end
      
c-B---------------------------------------------------------------------
      subroutine F3(M,u3,F33,x)
      implicit double precision (a-h,o-z)
      external F2
      double precision u3(*),F33(*),x(*)
      common /pars/ifunc

        a=-10.D0
        b=10.D0

       do L=1,M
      	x(3)=u3(L)
       	F33(L)=DGMLT2(F2,a,b,10,6,x)
       enddo

       return
       end 
c x1-------------------------------------

      subroutine F2(M,u2,F22,x)
      implicit double precision (a-h,o-z)
      double precision u2(*),F22(*),x(*)
      external F1
      common /pars/ifunc
      
     
      do L=1,M
      x(2)=u2(L) 
       
      if (ifunc.eq.1) then
       a=0.d0
       b=1.d0
           
       elseif (ifunc.eq.2) then
       
       a=0.d0
       b=1.d0
      
      
      elseif (ifunc.eq.3) then
      a=x(2)
      b=10.d0
      
      
      elseif (ifunc.eq.4) then
      a=x(2)
      b=10.d0
      
      elseif (ifunc.eq.5) then
      a=x(2)
      b=10.d0
      
      endif
     
      F22(L)=DGMLT1(F1,a,b,10,6,x)

      enddo
      return
      end
c-x2---------------------------------------------
      subroutine F1(M,u1,F11,x)
      implicit double precision (a-h,o-z)
      DOUBLE PRECISION  F11(*),u1(*),x(*)
      common/kine/Q,E,q_L
      common/nucl/r_n,c_n,wc,ro_0,A_n
      common/pars/ifunc
      
      EXTERNAL Fro
      do  L=1,M   
        x(1)=u1(L)

        r1= DSQRT(x(3)**2+x(1)**2)
        r2=DSQRT(x(3)**2+x(2)**2)
        F11(L)=0.d0
	if(r1.lt.10.d0.OR.r.lt.10.d0) then 
	A1=DEXP((r1-r_n)/c_n)
        A2=DEXP((r2-r_n)/c_n)
c	AP1=A1/((1.d0+A1)**2*c_nucl*r1)
c	AP2=A2/((1.d0+A2)**2*c_nucl*r2)!w=0
	P1=2.d0*wc*r1/R_n**2-(1.d0+wc*r1**2/R_n**2)*A1/(c_n*(1.d0+A1))
        AP1=P1/(r1*(1.d0+A1))
        P2=2.d0*wc*r2/R_n**2-(1.d0+wc*r2**2/R_n**2)*A2/(c_n*(1.d0+A2))
        AP2=P2/(r2*(1.d0+A2))        
                AR=(1.D0+wc*r2**2/r_n**2)/(1.D0+A2)
                 ABS=DEXP(-1.3d0*A_n*ro_0*Fro(x(2)))         
	       
	       
		if (ifunc.eq.1) then
		  
                  F11(L)=DBESJ1(Q*x(3))*x(3)**2*AP2*ABS
	 	
		 elseif (ifunc.eq.2) then
		 
                F11(L)=DBESJ0(Q*x(3))*x(3)*AR*ABS
		
		elseif (ifunc.eq.3) then
		  
                F11(L)=DBESJ1(Q*x(3))*x(3)**2*AR*AP1*ABS
		
            elseif (ifunc.eq.4) then     !x(1)>x(2)
		  
	  F11(L)=DBESJ1(Q*x(3))*x(3)**2*AR*AP1*DCOS((x(2)-x(1))*q_L)*ABS
		  elseif (ifunc.eq.5) then
		  
         F11(L)=DBESJ1(Q*x(3))*x(3)**2*AR*AP1*DSIN((x(2)-x(1))*q_L)*ABS
	
              endif
              endif
             enddo
            return     
            end
	     
*-----------------------------------------------
      double precision FUNCTION Fro(x2)
      implicit double precision (a-h,o-z)

      double precision x(3)
      common/xc/x
      common/nucl/r_n,c_n,wc,ro_0,A_n
      EXTERNAL FG
       
      a=x(2)
      b=10.D0
      eps=0.001D0
      Fro=DGAUSS(FG,a,b,eps)        
      return
      end

*-----------------------------------------------
      DOUBLE PRECISION FUNCTION FG(xg)
      implicit double precision (a-h,o-z)
       
      double precision x(3)
      common/xc/x
      common/nucl/r_n,c_n,wc,ro_0,A_n
       
      r=DSQRT(x(3)**2+xg**2)
      FG=0.d0
      if(r.lt.15.d0) then
         
	FG=(1.D0+wc*r**2/r_n**2)/(1.D0+DEXP((r-r_n)/c_n))
        
      end if
      return
      end



*==========================================================================
*===                        subroutine LoadTarget(Z)                    ===
*==========================================================================

      subroutine LoadTarget(Z)
      implicit double precision (a-h,o-z)
      integer Z
      common/nucl/r_n,c_n,wc,ro_0,A_n

      if (Z.eq.82) then
     	A_n  = 208.d0			! for Pb
     	r_n  = 6.624d0
    	c_n  = 0.549d0
     	wc   = 0.
        ro_0 = 7.7d-4
      elseif (Z.eq.50) then
        A_n  = 50.d0			! for 120Sn
        r_n  = 5.315d0
        c_n  = 0.576d0
        wc   = 0.d0
 	ro_0 = 1.425e-3
      elseif (Z.eq.6) then
        A_n  = 12.d0			! for C12
        r_n  = 2.355d0
        c_n  = 0.5224d0
        wc   = -0.149d0
 	ro_0 = 1.519d-02
      else 
         print*,"LoadTarget: Error. No database for Z=",Z
         stop
      end if

      return
      end

*==========================================================================
*===                           subroutine ErrMsg()                      ===
*==========================================================================


      subroutine ErrMsg()

      print*,' '
      print*,' Moeller : Calculate Moeller Cross Section'
      print*,' '
      print*,'       -h          show this help'
      print*,'       -Z <Z>      Atomic number of nucleus[def]:82'
      print*,'       -e <energy> incident energy [GeV]'
      print*,'       -t <angle>  pion angle [deg.]'
      print*,' '
      print*,'   Example: strong -e 5.0 '
      print*,' '

      return
      end
  
