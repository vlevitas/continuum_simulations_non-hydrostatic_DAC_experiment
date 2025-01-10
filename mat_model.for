
      MODULE Zrsample	
	  
! Import 'm' values in FRIC
! No. of 'm' values cosidered along rad from analytical sol (upto r=0.25 mm)
        INTEGER, PARAMETER :: mp=26               	 

! Link FRIC and UMAT (use connect.txt)
! 'YieldMe' stores yield strength of mixture and con of omega-Zr for each finite element. It is used in FRIC.
        INTEGER, PARAMETER :: NodeNo=220             ! No. of material nodes along sample-anvil contact surface
        INTEGER, PARAMETER :: EleNO=3958             ! No. of material elements in sample
        REAL*8, DIMENSION(EleNO,2) :: YieldMe         

		Real*8, dimension(eleno,1) :: q_zero = 0.0
		
        INTEGER, PARAMETER :: ElNumTop=2485          ! First element on upper contact surface
        INTEGER, PARAMETER :: ElNumBot=2451          ! First element on bottom contact surface 			  
        INTEGER, PARAMETER :: ElTop1Down=2484        ! 1 element down from first element on upper contact surface			  

		REAL*8, DIMENSION(3000,1) :: thkrec			 ! '3000' is used as an idea for the number of time increments.
        REAL*8 ytop,ybot,ytop1,thkcor,thick
		real*8 thick1
		
		Integer xitermax,trackElas,ElMaxXiter
	  
		real*8 detFMax,detFeMax
		Integer ELdetFMax,ELdetFeMax
		
		real*8 detFMin,detFeMin
		Integer ELdetFMin,ELdetFeMin 
		
		REAL*8, DIMENSION(26,4) :: datarec=0.0d0	
		
! ! ! ! For generating files		
	    ! ! ! integer :: xmark
		
! For finding number of iterations in each increment		
		Integer xmarker
		integer iter_marker
		integer kinc_marker		
		
      END MODULE Zrsample    
	  
! Note: The variables that are defined in 'Module' are accessible to all subroutines.
! The module is called by "USE".  

C**********************************************************************

      SUBROUTINE UVARM(UVAR,DIRECT,T,TIME,DTIME,CMNAME,ORNAME,
     1 NUVARM,NOEL,NPT,LAYER,KSPT,KSTEP,KINC,NDI,NSHR,COORD,
     2 JMAC,JMATYP,MATLAYO,LACCFLA) 

      INCLUDE 'ABA_PARAM.INC'

      CHARACTER*80 CMNAME,ORNAME
      CHARACTER*3 FLGRAY(15)

      DIMENSION UVAR(NUVARM),DIRECT(3,3),T(3,3),TIME(2)
      DIMENSION ARRAY(15),JARRAY(15),JMAC(*),JMATYP(*),COORD(*)

C Stress tensor:
      CALL GETVRM('S',ARRAY,JARRAY,FLGRAY,JRCD,JMAC,JMATYP,
     1 MATLAYO,LACCFLA)

      UVAR(1) = ARRAY(1)       

C      WRITE(6,*) "UVARM", UVAR(1),'NOEL', NOEL,'NPT',NPT
           
      RETURN
      END  
	  
C********************************************************************** 	  

      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
	 
	 
c      INCLUDE 'ABA_PARAM.INC'	 
      
      USE Zrsample

	  IMPLICIT NONE

      CHARACTER*80 CMNAME
	  
      INTEGER NDI,NSHR,NTENS,NSTATV,NPROPS,NOEL,NPT,LAYER,
     1 KSPT,KSTEP,KINC
	 
      DOUBLE PRECISION STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,
     3 PROPS,COORDS,DROT,PNEWDT,CELENT,DFGRD0,DFGRD1  

      DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1 DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)

C**********************************************************************

! defining integers, parameters and variables

	  INTEGER M,N	  
	  
      PARAMETER (M=3,N=3)	 
	     
      INTEGER I,J,K,L,F,K1,K2,K3,K4,P3,P4 

! for dowhile loop

	  INTEGER GoOn,Xiter

! For diamond

      Real*8 c11,c12,c44,
     1     c111,c112,c123,c144,c155,c166,c456,
     2     c1111,c1112,c1122,c1123,c1144,c1155,
     3     c1255,c1266,c1456,c4444,c4455	 
	 
      Real*8 FT(3,3),UU(3,3),xlstr(3,3),pk2(3,3),
     1     Fpk2(3,3),Fpk2FT(3,3),sig(3,3),d2Hd2eta(3,3,3,3),
     2     XLc(3,3,3,3),xkir(3,3),XKc(3,3,3,3),Cc(3,3,3,3)    
	 
      Real*8 eta1,eta2,eta3,eta4,eta5,eta6,
     1     d2Hd11,d2Hd12,d2Hd13,d2Hd16,d2Hd22,d2Hd23,d2Hd26,
     2     d2Hd33,d2Hd36,d2Hd66
	 
! For diamond and sample
	 	 
	 real*8 xiden(3,3)	 
	 real*8 detf	 

! For sample
	 	 
	  Real*8 XELAM0,XEG0,Xlpar0,Xmpar0,Xnpar0,
     1     XELAM1,XEG1,Xlpar1,Xmpar1,Xnpar1,
     2     sya,syap,syb,sybp,tstr, 
     3     XpdH,Xpde,XprH,Xpre,xk,
     4     ELAM,EG,Xlp,Xmp,Xnp,xNU 
	  
	  Real*8 xUIne(3,3),xUIold(3,3),xi4(3,3,3,3),xUIneI(3,3),
     1     xFUIneI(3,3),XRe(3,3),Xue(3,3),XVe(3,3),FeFeT(3,3),XBe(3,3),
     2     BeBe(3,3),dI3dBe(3,3),dEdBe1(3,3),Stress1(3,3),Dev1(3,3),
     3     FeTFe(3,3),xLagEe(3,3),d2I3dBe2(3,3,3,3),d2EdBe2(3,3,3,3),
     4 	   xdJeInvdBe(3,3),dFdBe(3,3,3,3),YBe(3,3,3,3),CBe(3,3,3,3)
	 
	  Real*8 xc,xq,xlan,xdetUIne,BeI1,BeI2,XdetFesq,XdetFe, 
     1    pc1,pr1,sy0,sy1,sy,ydiff,check1,check2  
  
	  Real*8 FInv(3,3),xl(3,3)
	  
	  real*8 xDIne(3,3),XDp(3,3)
	 
	  Real*8 xEplastic(3,3) 
 
	  real*8 Delta1, Delta2, Delta3, xB, q0 
	  
	  real*8 DFGRD1T(3,3), xFTF(3,3), xEplastic2(3,3), xEplastic3(3,3)
	  
	  real*8 xslope_pdeq
C**********************************************************************	         
      
!=======================  For diamond  ============================
c achyut-see the material name	       
      IF (CMNAME.EQ.'DIAMOND') THEN	

! Elastic properties of diamond 

! Second order elastic constants  
	  c11=props(1)
	  c12=props(2)
	  c44=props(3)
      
! Third order elastic constants
      c111=props(4)     
      c112=props(5)
	  c123=props(6)
      c144=props(7)     
      c166=props(8)
      c155=props(9)
      c456=props(10)
	  
! Fourth order elastic constants
      c1111=props(11)      
	  c1112=props(12)
	  c1122=props(13)
      c1123=props(14)
      c1144=props(15)
      c1155=props(16)
      c1255=props(17)     
	  c1266=props(18)
      c1456=props(19)
      c4444=props(20)
      c4455=props(21) 

! Transpose of total deformation gradient                  	
	  call ktrans(DFGRD1,FT)

! U.U=F^T.F
      Call KMLT(FT,DFGRD1,UU)
	  
! Identity matrix	
	  call xidentity(xiden)	  	  

! Lagrangian strain
	  xlstr=0.5d0*(UU-XIDEN)

! Components of Lagrangian strain for stress evaluation
      eta1=xlstr(1,1)	! Radial strain
      eta2=xlstr(2,2)	! Normal strain
      eta3=xlstr(3,3)	! Azimuthal strain
      eta6=2.0d0*xlstr(1,2)	! Shear strain
! For axisymmetric case, xlstr(1,3)=xlstr(2,3)=0 
      eta4=0.0d0
      eta5=0.0d0
      	  
! Second Piola-Kirchoff stress		  
	  pk2=0.0d0	  
      pk2(1,1)=c11*eta1+c12*(eta2+eta3)	! Linear terms		
     +  +1.0d0/2.0d0*c111*eta1**2	! Quadratic terms
     +  +1.0d0/2.0d0*C112*(2.0d0*eta1*(eta2+eta3)+eta2**2+eta3**2)
     +  +c123*eta2*eta3+1.0d0/2.0d0*C144*eta4**2
     +  +1.0d0/2.0d0*c166*(eta5**2+eta6**2)	! Cubic terms
     +  +1.0d0/6.0d0*c1111*eta1**3+1.0d0/6.0d0*C1112*(eta2**3+eta3**3
     +  +3.0d0*eta1**2*(eta2+eta3))  	 
     +  +0.25d0*C1122*(2.0d0*eta1*eta2**2+2.0d0*eta1*eta3**2)
     +  +0.5d0*C1123*(eta2*eta3*(eta1+eta2+eta3)+eta1*eta2*eta3)
     +  +0.5d0*C1144*eta1*eta4**2
     +  +0.5d0*C1155*eta1*(eta6**2+eta5**2)
     +  +0.5d0*C1255*(eta2*(eta4**2+eta5**2)
     +  +eta3*(eta4**2+eta6**2))	 
     +  +0.5d0*C1266*(eta2*eta6**2+eta3*eta5**2)
     +  +C1456*(eta4*eta5*eta6)	 

      pk2(2,2)=c11*eta2+c12*(eta1+eta3)	! Linear terms 
     +  +0.5d0*c111*eta2**2	! Quadratic terms
     +  +0.5d0*C112*(2.0d0*eta2*(eta1+eta3)+eta1**2+eta3**2)
     +  +c123*eta1*eta3+0.5d0*C144*eta5**2
     +  +0.5d0*c166*(eta4**2+eta6**2)	 
     +  +1.0d0/6.0d0*C1111*eta2**3	! Cubic terms
     +  +1.0d0/6.0d0*C1112*(eta1**3+eta3**3
     +  +3.0d0*eta2**2*(eta1+eta3))	 
     +  +0.25d0*C1122*(2.0d0*eta2*eta1**2+2.0d0*eta2*eta3**2)
     +  +0.5d0*C1123*(eta1*eta3*(eta1+eta2+eta3)+eta1*eta2*eta3)
     +  +0.5d0*C1144*eta2*eta5**2
     +  +0.5d0*C1155*eta2*(eta4**2+eta6**2)
     +  +0.5d0*C1255*(eta1*(eta4**2+eta5**2)
     +  +eta3*(eta5**2+eta6**2))
     +  +0.5d0*C1266*(eta1*eta6**2+eta3*eta4**2)
     +  +C1456*(eta4*eta5*eta6)
	    
      pk2(3,3)=c11*eta3+c12*(eta1+eta2)	! Linear terms
     +  +0.5d0*c111*eta3**2	! Quadratic terms
     +  +0.5d0*C112*(2.0d0*eta3*(eta1+eta2)+eta1**2+eta2**2)
     +  +c123*eta1*eta2+0.5d0*C144*eta6**2
     +  +0.5d0*c166*(eta4**2+eta5**2)	
     +  +1.0d0/6.0d0*C1111*eta3**2	! Cubic terms
     +  +1.0d0/6.0d0*C1112*(eta1**3+eta2**3
     +  +3.0d0*eta3**2*(eta1+eta2))
     +  +0.25d0*C1122*(2.0d0*eta3*eta2**2+2.0d0*eta3*eta1**2)
     +  +0.5d0*C1123*(eta1*eta2*(eta1+eta2+eta3)+eta1*eta2*eta3)
     +  +0.5d0*C1144*eta3*eta6**2
     +  +0.5d0*C1155*eta3*(eta5**2+eta4**2)
     +  +0.5d0*C1255*(eta2*(eta5**2+eta6**2)
     +  +eta1*(eta4**2+eta6**2))
     +  +0.5d0*C1266*(eta2*eta4**2+eta1*eta5*2)
     +  +C1456*(eta4*eta5*eta6)

! Corrected for term corresponding to C1155 that is used by Mehdi. 
      pk2(1,2)=c44*eta6	! Linear terms
     +  +c456*eta4*eta5+c144*eta3*eta6+c166*(eta1+eta2)*eta6	! Quadratic terms
     +  +0.5d0*C1144*eta6*eta3**2	! Cubic terms
     +  +0.5d0*C1155*(eta6*eta1**2+eta6*eta2**2)
     +  +C1255*(eta2*eta3*eta6+eta1*eta3*eta6)+C1266*eta1*eta2*eta6
     +  +C1456*eta4*eta5*(eta1+eta2+eta3)+1.0d0/6.0d0*C4444*eta6**3
     +  +0.5d0*C4455*(eta4**2+eta5**2)*eta6

      pk2(2,1)=pk2(1,2)
      
! Cauchy stress (sig)  
	  Call KMLT(DFGRD1,pk2,Fpk2)
      Call KMLT(Fpk2,FT,Fpk2FT)	  
      CALL mdet(DFGRD1,detF)

	  if (detF.lt.0.0d0) then
	   print*,'Det of F is negative in diamond', detF
	   call stdb_abqerr(
     1 -3,'Det of F in diamond = %R is negative',0,detf,'')
	  endif
	  
	  sig=0.0d0
	  sig=1.0d0/detF*Fpk2FT

! Imposing symmetric constraint on Cauchy stress 
	  call xsymmetric2(sig)

! Imposing axisymmetric constraint  	  
	  call xaxicons(sig)
	  
	  stress=0.0d0
      DO I=1,3
         STRESS(I)=sig(I,I)
      END DO
      STRESS(4)=sig(1,2)	! Note how tensor gets stored in vector form in Abaqus
	  	  
! ntens=ndi + nshr.
! ntens=6 for CGAX4R element 
      IF(NTENS.GT.4) THEN
          STRESS(5)=0.0d0
          STRESS(6)=0.0d0
      END IF 	  

! Double derivative of Helmholtz energy wrt Lagrangian strain	        
      d2Hd11=c11
     +  +c111*eta1+c112*(eta2+eta3)
     +  +0.5d0*C1111*eta1**2
     +  +C1112*eta1*(eta2+eta3)+0.5d0*C1122*(eta2**2+eta3**2)
     +  +C1123*eta2*eta3+0.5d0*C1144*eta4**2
     +  +0.5d0*C1155*(eta5**2+eta6**2)	 
     
      d2Hd12=c12+c112*(eta1+eta2)+c123*eta3
     +  +0.5d0*C1112*(eta1**2+eta2**2)+C1122*eta1*eta2+
     +  +0.5d0*C1123*(eta3*(eta1+eta2+eta3)+eta2*eta3+eta1*eta3)
     +  +0.5d0*C1255*(eta4**2+eta5**2)+0.5d0*C1266*eta6**2
 
      d2Hd13=c12+c112*(eta1+eta3)+c123*eta2
     +  +0.5d0*C1112*(eta1**2+eta3**2)+C1122*eta1*eta3+
     +  +0.5d0*C1123*(eta2*(eta1+eta2+eta3)+eta2*eta3+eta1*eta2)
     +  +0.5d0*C1255*(eta4**2+eta6**2)+0.5d0*C1266*eta5**2	 
     
      d2Hd16=c166*eta6+C1155*eta1*eta6+C1255*eta3*eta6+
     +  +C1266*eta2*eta6+C1456*eta4*eta5 	 

! Please note that for axisymmetric case, d2Hd14=d2Hd15=0    
      d2Hd22=c11+c111*eta2+c112*(eta1+eta3)+0.5d0*C1111*eta2**2
     +  +C1112*eta2*(eta1+eta3)+0.5d0*C1122*(eta1**2+eta3**2)
     +  +C1123*eta1*eta3+0.5d0*C1144*eta5**2
     +  +0.5d0*C1155*(eta4**2+eta6**2)
    
      d2Hd23=c12+c112*(eta2+eta3)+c123*eta1
     +  +0.5d0*C1112*(eta2**2+eta3**2)+C1122*eta2*eta3+
     +  +0.5d0*C1123*(eta1*(eta1+eta2+eta3)+eta1*eta3+eta1*eta2)
     +  +0.5d0*C1255*(eta5**2+eta6**2)+0.5d0*C1266*eta4**2
     
      d2Hd26=c166*eta6+C1155*eta2*eta6+C1255*eta3*eta6
     +      +C1266*eta1*eta6+C1456*eta4*eta5
! Please note that for axisymmetric case, d2Hd24=d2Hd25=0 
	 
! Term corrected corresponding to C1112.      
      d2Hd33=c11+c111*eta3+c112*(eta1+eta2)+0.5d0*C1111*eta3**2
     +      +C1112*eta3*(eta1+eta2)+0.5d0*C1122*(eta1**2+eta2**2)
     +      +C1123*eta1*eta2+0.5d0*C1144*eta6**2
     +      +0.5d0*C1155*(eta4**2+eta5**2) 
     
      d2Hd36=c144*eta6+C1144*eta3*eta6+C1255*(eta2*eta6+eta1*eta6)
     1    +C1456*eta4*eta5
! Please note that for axisymmetric case, d2Hd34=d2Hd35=0 

! Term corrected corresponding to C1266.       
      d2Hd66=c44+c144*eta3+c166*(eta1+eta2)+0.5d0*C1144*eta3**2+
     1    0.5d0*C1155*(eta1**2+eta2**2)+C1255*eta3*(eta1+eta2)+
     2    C1266*eta1*eta2+0.5d0*C4444*eta6**2+     
     3    0.5d0*C4455*(eta4**2+eta5**2)
	 
	  d2Hd2eta=0.0d0
	  
      d2Hd2eta(1,1,1,1)=d2Hd11
      d2Hd2eta(1,1,2,2)=d2Hd12      
      d2Hd2eta(1,1,3,3)=d2Hd13      
      d2Hd2eta(1,1,1,2)=d2Hd16                  
      d2Hd2eta(2,2,2,2)=d2Hd22      
      d2Hd2eta(2,2,3,3)=d2Hd23    
      d2Hd2eta(2,2,1,2)=d2Hd26      
      d2Hd2eta(3,3,3,3)=d2Hd33       
      d2Hd2eta(3,3,1,2)=d2Hd36     
      d2Hd2eta(1,2,1,2)=d2Hd66     

! Major symmetries	  
      d2Hd2eta(2,2,1,1)=d2Hd12      
      d2Hd2eta(3,3,1,1)=d2Hd13      
      d2Hd2eta(1,2,1,1)=d2Hd16                      
      d2Hd2eta(3,3,2,2)=d2Hd23    
      d2Hd2eta(1,2,2,2)=d2Hd26             
      d2Hd2eta(1,2,3,3)=d2Hd36  

! Minor symmetries   
      d2Hd2eta(1,1,2,1)=d2Hd16                    
      d2Hd2eta(2,2,2,1)=d2Hd26          
      d2Hd2eta(3,3,2,1)=d2Hd36     
      d2Hd2eta(1,2,2,1)=d2Hd66 	  

! Minor symmetries              
      d2Hd2eta(2,1,1,2)=d2Hd66

! Minor symmetries   
      d2Hd2eta(2,1,1,1)=d2Hd16                       
      d2Hd2eta(2,1,2,2)=d2Hd26             
      d2Hd2eta(2,1,3,3)=d2Hd36

! Minor symmetries      
      d2Hd2eta(2,1,2,1)=d2Hd66
	  
! Tangent moduli
      XLc=0.0d0	  
      DO 1 I=1,M
       DO 1 J=1,M
        DO 1 K=1,M
         DO 1 L=1,M
          DO 1 K1=1,M
           DO 1 K2=1,M
            DO 1 K3=1,M
             DO 1 K4=1,M
          XLc(I,J,K,L)=XLc(I,J,K,L)+DFGRD1(I,K4)*DFGRD1(K4,K3)*
     1                 d2Hd2eta(K3,J,K,K2)*FT(K2,K1)*FT(K1,L)
1     CONTINUE

! Kirchoff stress = Cauchy stress*Jacobian
	  xkir=0.0d0
	  xkir=detF*sig
	  XKc=0.0d0
      DO 2 I=1,M
       DO 2 J=1,M
        DO 2 K=1,M
         DO 2 L=1,M
      XKc(I,J,K,L)=(xkir(I,K)*XIDEN(J,L)+xkir(J,K)*XIDEN(I,L)
     1     +xkir(I,L)*XIDEN(J,K)+xkir(J,L)*XIDEN(I,K))*0.25d0  
2     CONTINUE 
		
	  Cc=1.0d0/detF*(XLc+2.0d0*XKc)
      
! Imposing symmetries : This is used to avoid numerical errors	 
      DO 3 i=1,M
       DO 3 j=1,M
        DO 3 k=1,M
         DO 3 l=1,M
          Cc(k,l,i,j)=Cc(i,j,k,l)
          Cc(j,i,k,l)=Cc(i,j,k,l)
          Cc(i,j,l,k)=Cc(i,j,k,l)
3     CONTINUE

! Defining the DDSDDE matrix	  
      DDSDDE= 0.0d0
      DO 4 K1=1,3
       DO 4 K2=1,3
        DDSDDE(K1,K2)=(Cc(K1,K1,K2,K2)+Cc(K2,K2,K1,K1))/2.0d0
4     CONTINUE 
      
      DDSDDE(1,4)=(Cc(1,1,1,2) + Cc(1,2,1,1))/2.0d0
      DDSDDE(4,1)=DDSDDE(1,4)
      DDSDDE(2,4)=(Cc(2,2,1,2) + Cc(1,2,2,2))/2.0d0
      DDSDDE(4,2)=DDSDDE(2,4)
      DDSDDE(3,4)=(Cc(3,3,1,2) + Cc(1,2,3,3))/2.0d0
      DDSDDE(4,3)=DDSDDE(3,4)
      DDSDDE(4,4)=Cc(1,2,1,2)  
		
      DDSDDE(1,5)=(Cc(1,1,1,3) + Cc(1,1,1,3))/2.0d0
      DDSDDE(5,1)=DDSDDE(1,5)
      DDSDDE(2,5)=(Cc(2,2,1,3) + Cc(1,3,2,2))/2.0d0
      DDSDDE(5,2)=DDSDDE(2,5)
      DDSDDE(3,5)=(Cc(3,3,1,3) + Cc(1,3,3,3))/2.0d0
      DDSDDE(5,3)=DDSDDE(3,5)
      DDSDDE(4,5)=(Cc(1,2,1,3) + Cc(1,3,1,2))/2.0d0
      DDSDDE(5,4)=DDSDDE(4,5)
      DDSDDE(5,5)=Cc(1,3,1,3) 

      DDSDDE(1,6)=(Cc(1,1,2,3) + Cc(2,3,1,1))/2.0d0
      DDSDDE(6,1)=DDSDDE(1,6)
      DDSDDE(2,6)=(Cc(2,2,2,3) + Cc(2,3,2,2))/2.0d0
      DDSDDE(6,2)=DDSDDE(2,6)
      DDSDDE(3,6)=(Cc(3,3,2,3) + Cc(2,3,3,3))/2.0d0
      DDSDDE(6,3)=DDSDDE(3,6)
      DDSDDE(4,6)=(Cc(1,2,2,3) + Cc(2,3,1,2))/2.0d0
      DDSDDE(6,4)=DDSDDE(4,6)
      DDSDDE(5,6)=(Cc(1,3,2,3) + Cc(1,3,2,3))/2.0d0
      DDSDDE(6,5)=DDSDDE(5,6)    
      DDSDDE(6,6)=Cc(1,2,1,2)
    
	  endif ! END OF DIAMOND MATERIAL

!=======================  For sample  ============================ 
c achyut-see the material name	       
      IF (CMNAME.EQ.'ZIRCONIUM') THEN	

C**********This portion is for # itertions in each increment***********

	  if (kinc_marker .eq. kinc) then
		iter_marker = 0
	  endif

	  if (kinc_marker .ge. kinc .and. noel .eq. 1) then
		iter_marker = iter_marker + 1
		kinc_marker = kinc + 1
	  endif

C**********************************************************************
	  
	  if(noel.eq.1) then
		xitermax=0
		trackElas=0
		ElMaxXiter=0
		detFMax=0.0d0
		detFeMax=0.0d0
		ELdetFMax=0
		ELdetFeMax=0
		detFMin=100.0d0
		detFeMin=100.0d0
		ELdetFMin=0
		ELdetFeMin=0
	  endif	  

! Murnaghan properties of Alpha-Zr
      XELAM0=props(1)        ! this is Lamda
      XEG0=props(2)          ! G
      Xlpar0=props(3)		 ! l
      Xmpar0=props(4)		 ! m
      Xnpar0=props(5)		 ! n
      
! Murnaghan properties of Omega-Zr
      XELAM1=props(6)        
      XEG1=props(7)          
      Xlpar1=props(8)
      Xmpar1=props(9)
      Xnpar1=props(10)    
   
	  SyA=PROPS(11)    ! yield strength of alpha-Zr @ p=0 
      SyAP=PROPS(12)   ! linear pressure dependence of yield strength of alpha-Zr	  

      SyB=PROPS(13)    ! yield strength of omega-Zr @ p=0
      SyBP=PROPS(14)   ! linear pressure dependence of yield strength of omega-Zr
	  
	  tstr=props(41)      	! Transformation strain (volumetric)-not in percentage 
	  Xpde=props(42)
	  Xpdh=props(43)
	  XprH=props(44)
	  Xpre=props(45)
	  xk=props(46)
	  Delta1=props(47)
	  Delta2=props(48)
	  xB=props(49)

C**********************************************************************

! Cases 1, 2: From Improved constants
	  XELAM0 = 68.39D9
	  XEG0 = 35.89D9
	  Xlpar0 = -146.71D9
	  Xmpar0 = -119.10D9
	  Xnpar0 = -100.0D9
	  
	  XELAM1 = 75.41D9
	  XEG1 = 45.94D9
	  Xlpar1 = -91.28D9
	  Xmpar1 = -95.51D9
	  Xnpar1 = -4.0D9

	  Xpde=2.65D9
	  Xpdh = 5.4d9
	  
	  xk = 5.87d0
		
	  Delta1 = 0.0d0	
	  Delta2 = 1.0d0
	  
	  Delta3 = 0.0d0/10**(9)

	  xB = 0.0d0

C**********************************************************************
	  
! This is plastic right stretch tensor at the start of increment 
	  xUIne=0.0d0
      xUIne(1,1)=STATEV(1)
      xUIne(2,2)=STATEV(2)
      xUIne(3,3)=STATEV(3) 
      xUIne(1,2)=STATEV(4)
      xUIne(2,1)=xUIne(1,2)
	  
	  xUIold=xUIne

	  call xsymmetric2(xUIne)
	  call xaxicons(xUIne)

! 'xc' is concentration stored as state variable
! 'xc' is required when kinetic equation is modeled	  
 	  xc=statev(9)

! 'xq' is Odqvist parameter stored as state-variable
	  xq=statev(7)
	  
! Remember- All these state-variabes are stored for each finite element 
	  
! 'xlan' is plastic multiplier from the flow rule
      xlan=0.0d0    

! yieldme is initialized with zero at the start of every iteration.
! Important to note that UMAT functions in the increasing order of	elements.	  
	  if (noel.eq.1) then  
	   yieldme=0.0d0 
	  endif

C**********************************************************************
	  
      ELAM=(1.0d0-xc)*XELAM0+xc*XELAM1
      EG=(1.0d0-xc)*XEG0+xc*XEG1
      XLP=(1.0d0-xc)*Xlpar0+xc*Xlpar1
      XMP=(1.0d0-xc)*Xmpar0+xc*Xmpar1
      XNP=(1.0d0-xc)*Xnpar0+xc*Xnpar1
      
	  xNU=ELAM/(2.0d0*(ELAM+EG))

      IF (xNU .GT. 0.4d0) THEN
       PRINT*,'xNU is large',xNU
	   call stdb_abqerr(-3,'Poisson ratio = %R is large',0,xNU,'')
      ENDIF

! Initializing right stretch plastic tensor  
	  if (kinc.le.1) then
	   call xidentity(xUIne) 
	  endif

! Determinant of 'U_Ine'	  	  
	  CALL mdet(xUIne,xdetUIne)

! Nullity check for 'U_p' 	  
      IF (xdetUIne.Eq.0.0d0) THEN
	   PRINT*,'Det of U_Ine < 0',xdetUIne
	   call stdb_abqerr(
     1 -3,'Det of U_Ine in sample = %R is negative',0,xdetUIne,'')
      ENDIF

	  call xidentity(xiden) 
	  call xidentity4(xi4)
	  
! Inverse of 'U_Ine'            
      call xmatInv3D(xUIne,xUIneI) 

! Determinant of F 	  
      CALL mdet(DFGRD1,detF)	
	  
	  if (detF.le.0.0d0) then
	   print*,'Det of F is negative in sample',detF
	   call stdb_abqerr(
     1 -3,'Det of F in sample = %R is negative',0,detf,'')
	  endif

! F_e trial
      Call KMLT(DFGRD1,xUIneI,xFUIneI)

! Polar decomposition of F_e trial 
      Call skinem(xFUIneI,XRe,XUe)

! Left stretch tensor of F_e trial	  
      CALL gettheVe(xRe,XUe,XVe) 

! V_e.V_e = F_e.F_e^T                  
      Call KMLT(XVe,XVe,FeFeT)       

! 'XBe' above is trial elastic strain measure
	  xbe=0.0d0
	  xbe=0.5d0*(FeFeT-xiden)

! Imposing symmetric conditions
	  call xsymmetric2(XBe)
	  call xsymmetric2(XVe)

! Imposing axisymmetric conditions
	  call xaxicons(XBe)
	  call xaxicons(XVe)

! First invariant of Be	  
	  call xinvariant1(XBe,BeI1)	
	  
! Second invariant of Be	  
	  call xinvariant2(XBe,BeI2)	

! Another way to find invariants.
! But, it has more error than the method used.
	  ! ! ! call sinv(xbe,sinv1,sinv2,3,3)
	  
! BeBe=XBe*XBe
	  call kmlt(XBe,XBe,BeBe)

! Derivative of 3rd invariant of B_e wrt B_e.
	  dI3dBe=0.0d0
	  dI3dBe=BeBe-BeI1*xbe+BeI2*xiden
	  call xsymmetric2(dI3dBe)	  
	  
! Derivative of energy wrt B_e.	
	  dEdBe1=0.0d0 
      dEdBe1=ELAM*BeI1*XIDEN+2.0d0*EG*XBe
     1            +(XLP*BeI1**2-2.0d0*XMP*BeI2)*XIDEN
     2            +XNP*dI3dBe+2.0d0*XMP*BeI1*XBe 
	  call xsymmetric2(dEdBe1)	

! detFe: Determinant of F_e trial	  
	  CALL mdet(2.0d0*XBe+xiden,XdetFesq)
	  XdetFe=dsqrt(XdetFesq)

	  if (isNaN(XdetFe)) then
		PRINT*,'DetFe < 0 (elastic)',kinc,noel
		call xit
	  endif 
	  
! stress1=1/detFe*(2*Be+I)*dEdBe1
	  Stress1=0.0d0
      DO 101 I=1,3            
       DO 101 J=1,3			  
        DO 101 F=1,3
         stress1(I,J)=stress1(I,J)+1.0d0/XdetFe*(2.0d0*XBe(I,F)+
     +       XIDEN(I,F))*dEdBe1(F,J)
101   Continue

	  call xsymmetric2(Stress1)      
	  call xaxicons(Stress1)    

! Deviatoric stress  
      Call KDEVIA(stress1,dev1)	

! 'pc1' is Mises yield strength	
      Call KEFFP(dev1,pc1)		

! Pressure from Cauchy stress
	  call pressure(Stress1,pr1)

! Current yield strength in alpha phase
	  Sy0=SyA+SyAP*pr1	

! Current yield strength in omega phase	  
	  sy1=SyB+SyBP*pr1	

! Current yield strength in mixture  	  
      Sy=(1.0d0-xc)*sy0+xc*sy1	
	  
! Note yield strength is saturated with respect to plastic deformation  
       
      ydiff=pc1-Sy         	  
	  check1=0.0d0
	  check2=0.0d0
	  
      IF (ydiff.LE.0.0d0) THEN
	  check1=1.0d0
	  
	  trackElas=1+trackElas	

!============== Elasticity in sample  ======================== 

C**********************************************************************
	  
	  Call KMLT(xue,xue,FeTFe)
	  
	  xLagEe=0.0d0
  	  xLagEe=0.5d0*(FeTFe-XIDEN)

	  if (sy.lt.0.0d0 .or. xc.lt.-1.0D-5 .or. xq.lt.-1.0d-7) then
	   print*,'yield strength or concentration or Odqvist par <0:elas'
	   print*,sy,xc,xq,noel,kinc
	   call xit
	  endif

	  stress=0.0d0	
      DO I=1,3
       STRESS(I)=Stress1(I,I)
      END DO
      STRESS(4)=Stress1(1,2)

      IF(NTENS.GT.4) THEN
        STRESS(5)=0.0d0  !Stress1(1,3)
        STRESS(6)=0.0d0  !Stress1(2,3)
      END IF

	  yieldme(noel,1)=xc
      yieldme(NOEL,2)=sy

      STATEV(1)=xUIne(1,1)
      STATEV(2)=xUIne(2,2)
      STATEV(3)=xUIne(3,3) 
      STATEV(4)=xUIne(1,2)
      STATEV(5)=0.0d0	!xUIne(3,1)
      STATEV(6)=0.0d0	!xUIne(3,2)
	  
	  statev(7)=xq
	  statev(9)=xc
	  statev(10)=xlan
	  statev(11)=sy
	  statev(13)=pr1
	  	            
      K=14 
      DO 102 I=1,3
       DO 102 J=I,3
        K=K+1
        STATEV(K)=XBe(I,J)         ! Be-STATEV(15-20)
        STATEV(K+6)=XVe(I,J)       ! Ve-STATEV(21-26)
        STATEV(K+12)=Stress1(I,J)  ! Stress-STATEV(27-32)
102   CONTINUE
	  
	  statev(35)=xLagEe(1,1)
	  statev(36)=xLagEe(2,2)
	  statev(37)=xLagEe(3,3)
	  statev(38)=xLagEe(1,2)

	  statev(71)= COORDS(1) 

	  statev(80)= COORDS(2)	  

	  !!!statev(83) = xpde
	  
	  xDIne=0.0d0
	  XDp=0.0d0	  

! d2I3/d(Be)2	  
      d2I3dBe2=0.0d0            	  
      DO 103 I=1,3 
	   DO 103 J=1,3
		DO 103 l=1,3   
         DO 103 k=1,3
		   
	  d2I3dBe2(i,j,l,k)=0.5d0*(xbe(i,l)*xiden(j,k)+xbe(i,k)*xiden(j,l)+
     +    xbe(k,j)*xiden(i,l)+xbe(l,j)*xiden(i,k))-
     +    xbe(i,j)*xiden(l,k)-xbe(l,k)*xiden(i,j)+
     +    xiden(i,j)*xiden(l,k)*BeI1
     +    -0.5d0*BeI1*(xi4(i,j,l,k)+xi4(j,i,l,k))
103   CONTINUE	  

! d2(E)/d(Be)2
	  d2EdBe2=0.0d0
      DO 104 I=1,3 
	   DO 104 J=1,3
		DO 104 K=1,3   
         DO 104 L=1,3		   
            d2EdBe2(i,j,k,l)=elam*xi4(i,k,j,l)+EG*(xi4(i,j,k,l)+
     + 		xi4(j,i,k,l))+2.0d0*XLP*BeI1*xi4(i,k,j,l)+XMP*BeI1*
     + 		(xi4(i,j,k,l)+xi4(j,i,k,l))+2.0d0*XMP*xiden(k,l)*xbe(i,j)+
     +  	2.0d0*xmp*xiden(i,j)*xbe(k,l)-2.0d0*xmp*BeI1*xi4(i,k,j,l)+
     + 		xnp*d2I3dBe2(i,j,k,l)
104   CONTINUE

	  xdJeInvdBe=0.0d0
	  xdJeInvdBe=-(1.0d0/xDetFe**3)*
     + 		(xiden+2.0d0*(BeI1*xiden-xbe)+4.0d0*dI3dBe)
	 
! (df/dBe)  
	  dFdBe=0.0d0 
      DO 105 I=1,3 
	   DO 105 J=1,3
		DO 105 K=1,3   
         DO 105 L=1,3
		 
		  DO F=1,3		   
        dFdBe(i,j,k,l)=dFdBe(I,J,K,L)+
     + 	xdJeInvdBe(k,l)*(2.0d0*xbe(i,f)+xiden(i,f))*
     +dEdBe1(f,j)+1.0d0/xdetFe*(xi4(i,f,k,l)+xi4(f,i,k,l))*dEdBe1(f,j)	
     + 	+1.0d0/xdetFe*(2.0d0*xbe(i,f)+xiden(i,f))*d2EdBe2(f,j,k,l)	 
	      enddo
		  
105   CONTINUE	

! Y(Be)=(df/dBe).(2Be+I)
      YBe=0.0d0
	  CBe=0.0d0
      DO 106 I=1,3
       DO 106 J=1,3
        DO 106 K=1,3
         DO 106 L=1,3
		 
          DO F=1,3
           YBe(I,J,K,L)=YBe(I,J,K,L)+dFdBe(I,J,K,F)*
     +                 (2.0d0*XBe(F,L)+XIDEN(F,L))
          ENDDO 
		 
         CBe(I,J,K,L)=YBe(I,J,K,L)+stress1(I,J)*XIDEN(K,L)
106   CONTINUE
                   
! Imposing major and minor symmetries
      DO 107 I=1,3
       DO 107 J=1,3
        DO 107 K=1,3
         DO 107 L=1,3
          CBe(K,L,I,J)=CBe(I,J,K,L)
          CBe(J,I,K,L)=CBe(I,J,K,L)
          CBe(I,J,L,K)=CBe(I,J,K,L)
107   CONTINUE  
          
! Defining DDSDDE	  
      DDSDDE=0.0d0	  
      DO 108 K1=1,3
       DO 108 K2=1,3
        DDSDDE(K1,K2)=(CBe(K1,K1,K2,K2)+CBe(K2,K2,K1,K1))/2.0d0
108   CONTINUE
                     
      DDSDDE(1,4)=(CBe(1,1,1,2) + CBe(1,2,1,1))/2.0d0
      DDSDDE(4,1)=DDSDDE(1,4)
      DDSDDE(2,4)=(CBe(2,2,1,2) + CBe(1,2,2,2))/2.0d0
      DDSDDE(4,2)=DDSDDE(2,4)
      DDSDDE(3,4)=(CBe(3,3,1,2) + CBe(1,2,3,3))/2.0d0
      DDSDDE(4,3)=DDSDDE(3,4)
      DDSDDE(4,4)=CBe(1,2,1,2)  
		
      DDSDDE(1,5)=(CBe(1,1,1,3) + CBe(1,1,1,3))/2.0d0
      DDSDDE(5,1)=DDSDDE(1,5)
      DDSDDE(2,5)=(CBe(2,2,1,3) + CBe(1,3,2,2))/2.0d0
      DDSDDE(5,2)=DDSDDE(2,5)
      DDSDDE(3,5)=(CBe(3,3,1,3) + CBe(1,3,3,3))/2.0d0
      DDSDDE(5,3)=DDSDDE(3,5)
      DDSDDE(4,5)=(CBe(1,2,1,3) + CBe(1,3,1,2))/2.0d0
      DDSDDE(5,4)=DDSDDE(4,5)
      DDSDDE(5,5)=CBe(1,3,1,3) 

      DDSDDE(1,6)=(CBe(1,1,2,3) + CBe(2,3,1,1))/2.0d0
      DDSDDE(6,1)=DDSDDE(1,6)
      DDSDDE(2,6)=(CBe(2,2,2,3) + CBe(2,3,2,2))/2.0d0
      DDSDDE(6,2)=DDSDDE(2,6)
      DDSDDE(3,6)=(CBe(3,3,2,3) + CBe(2,3,3,3))/2.0d0
      DDSDDE(6,3)=DDSDDE(3,6)
      DDSDDE(4,6)=(CBe(1,2,2,3) + CBe(2,3,1,2))/2.0d0
      DDSDDE(6,4)=DDSDDE(4,6)
      DDSDDE(5,6)=(CBe(1,3,2,3) + CBe(1,3,2,3))/2.0d0
      DDSDDE(6,5)=DDSDDE(5,6)    
      DDSDDE(6,6)=CBe(1,2,1,2)
 	  
!============== Plasticity in sample  ======================== 

      ELSE

	  check2=1.0d0  
	  
! FINV:Inverse of F
! xL:velocity gradient       
      CALL xmatInv3D(DFGRD1,FINV)  
	  
	  xl=0.0d0  
      DO 201 I=1,3 
	   DO 201 J=1,3
	   
		DO K=1,3
		 xl(I,J)=xl(I,J)+(DFGRD1(I,K)-DFGRD0(I,K))*
     +              FInv(K,J)/DTIME
		enddo
		
201	  continue

	  xmarker=0

	  q0 = q_zero(noel, 1)

! Here, linear relationship between pde and q0
	  
	  xslope_pdeq = (2.55D9-2.65D9)/(0.52738-0.42806) 
	  
	  xpde = xslope_pdeq*(q0 - 0.42806) + 2.65D9
	  
	  if (xmarker .eq. 0) then 
	  
	  call NR4(noel, kinc,
     1  xl, dtime, 
     2 	statev(7), statev(9), statev(15), statev(16), statev(18),  	
     3	statev(20), xbe, stress1, xq, xlan, xc, xUIne, sy, pr1,
     4  XVe, xLagEe, xEplastic, CBe, DFGRD1, xUIold, xEplastic2, 
     4	xitermax, ElMaxXiter, ELdetFMax, ELdetFeMax, ELdetFMin, 
     5	ELdetFeMin, detFMax, detFeMax, detFMin, detFeMin, kstep,   
     6	ddsdde, iter_marker,     		
     7	XELAM0,XEG0,Xlpar0,Xmpar0,Xnpar0,
     8	XELAM1,XEG1,Xlpar1,Xmpar1,Xnpar1,
     9	tstr,sya,syap,syb,sybp,Xpdh,Xpde,xk,Delta1,Delta2,Delta3,xB,q0)
	 
	  endif

C**********************************************************************

	  if (xc.eq.0.0d0) then
		q_zero(noel, 1) = xq
	  endif

C**********************************************************************
	  
! Yield required for FRIC
	  YieldMe(NOEL,1)=xc
      YieldMe(NOEL,2)=sy      

! Update the stress tensor
      DO I=1,3               
         STRESS(I)=stress1(I,I)
      ENDDO
         STRESS(4)=stress1(1,2)
! Note element type is CGAX4R, dim(stress)is 6*1		  
      IF (NTENS.GT.4) THEN
         STRESS(5)=0.0d0
         STRESS(6)=0.0d0 
      ENDIF
	  	                     
      STATEV(1)=xUIne(1,1)
      STATEV(2)=xUIne(2,2)
      STATEV(3)=xUIne(3,3) 
      STATEV(4)=xUIne(1,2)
! xUIne(3,1)=xUIne(3,2)=0; Axisymmetric conditions
      STATEV(5)=0.0d0 
      STATEV(6)=0.0d0  
	  
	  STATEV(7)=Xq 
      STATEV(9)=Xc
      STATEV(10)=xlan  
	  
      Statev(11)=Sy  
	  statev(13)=pr1
	  	           
      K=14   !! starts with STATEV(15)
      DO I=1,3
       DO J=I,3
        K=K+1
        STATEV(K)=XBe(I,J)          ! Be  STATEV(15-20)
        STATEV(K+6)=XVe(I,J)        ! Ve  STATEV(21-26)
        STATEV(K+12)=stress1(I,J)   ! Stress  STATEV(27-32)
        ENDDO
      ENDDO

	  statev(35)=xLagEe(1,1)
	  statev(36)=xLagEe(2,2)
	  statev(37)=xLagEe(3,3)
	  statev(38)=xLagEe(1,2)

	  STATEV(67)=xEplastic(1,1)
      STATEV(68)=xEplastic(2,2)
      STATEV(69)=xEplastic(3,3) 
      STATEV(70)=xEplastic(1,2)	 

	  statev(71)= COORDS(1)	  

	  STATEV(72)=xEplastic2(1,1)
      STATEV(73)=xEplastic2(2,2)
      STATEV(74)=xEplastic2(3,3) 
      STATEV(75)=xEplastic2(1,2)	

	  call ktrans(DFGRD1,DFGRD1T)	
	  
	  call kmlt(DFGRD1T,DFGRD1,xFTF)
	  
	  xEplastic3 = 0.5d0*(xFTF-xiden)
	  
	  STATEV(76)=xEplastic3(1,1)
      STATEV(77)=xEplastic3(2,2)
      STATEV(78)=xEplastic3(3,3) 
      STATEV(79)=xEplastic3(1,2)
	  
	  statev(80)= COORDS(2)

	  statev(81) = q0
	  
	  statev(82) = xq - q0
	  
	  statev(83) = xpde
	  
	  if (xc.le.0.0d0) then
		statev(84) = 0.0d0
	  else
		statev(84) = xpde	  
	  endif 
	  
      ENDIF !! END OF PLASTIC PORTION OF SAMPLE   

C**********************************************************************

	  if (noel.eq.1) then
	   ybot=0.0d0
	   ytop=0.0d0
	   ytop1=0.0d0		
	  endif

! y-coordinate of first element on lower contact surface
      IF (NOEL.EQ.ElNumBot) THEN         
        ybot=COORDS(2)        
      ENDIF
	  
! y-coordinate of 1 element down from first element on upper contact surface	  
      IF (NOEL.EQ.ElTop1Down) THEN		 
       ytop1=COORDS(2)
      ENDIF

! y-coordinate of first element on upper contact surface                
      IF (NOEL.EQ.ElNumTop) THEN		 
       ytop=COORDS(2)
      ENDIF 	  
	  
      IF (NOEL.EQ.EleNo) THEN	
	   thkcor=dabs(ytop-ytop1)
	   thick=dabs(ytop-ybot)
	   thick1=thkcor+thick
	   thkrec(kinc,1)=thick1
      ENDIF

C**********************************************************************  
	  
      IF (NOEL.EQ.EleNO.and.npt.eq.1) THEN	  
	  PRINT*,KINC,dtime,trackElas,xitermax,ElMaxXiter,npt,
     1    detFMax,ELdetFMax,detFeMax,ELdetFeMax,
     2    detFMin,ELdetFMin,detFeMin,ELdetFeMin,kstep,iter_marker	  
      ENDIF	
	  
      IF (DTIME.LT.1.0D-7.and.kinc.ge.1) THEN
	   PRINT*, 'DTIME=', DTIME
	   print*,check1,check2,noel,kinc
	   call xit
	  ENDIF
  
c achyut-add the data files below
C**********************************************************************  

	  if (NOEL.EQ.EleNO.and.(thick+thkcor).le.0.165d0.and.
     +	statev(41).eq.0 .and. kinc.gt.1) then 	 
		statev(41)=1		
		datarec(1,1)=0.165d0
		datarec(1,2)=thick+thkcor
		datarec(1,3)=kinc
		datarec(1,4)=kinc/5	

	    call fprint(datarec) 
	 
	  endif
	  
	  if (NOEL.EQ.EleNO.and.(thick+thkcor).le.0.16477d0.and.
     +	statev(42).eq.0) then 	 
		statev(42)=1		
		datarec(2,1)=0.16477d0	
		datarec(2,2)=thick+thkcor
		datarec(2,3)=kinc	
		datarec(2,4)=kinc/5	

		call fprint(datarec) 
		
	  endif
	  
	  if (NOEL.EQ.EleNO.and.(thick+thkcor).le.0.163d0.and.
     +	statev(43).eq.0) then 	 
		statev(43)=1		
		datarec(3,1)=0.163d0	
		datarec(3,2)=thick+thkcor
		datarec(3,3)=kinc	
	    datarec(3,4)=kinc/5	
		
	    call fprint(datarec) 
		
	  endif
	  
	  if (NOEL.EQ.EleNO.and.(thick+thkcor).le.0.15432d0.and.
     +	statev(44).eq.0) then 	 
		statev(44)=1		
		datarec(4,1)=0.15432d0	
		datarec(4,2)=thick+thkcor
		datarec(4,3)=kinc	 
		datarec(4,4)=kinc/5	

		call fprint(datarec)	
		
	  endif

	  if (NOEL.EQ.EleNO.and.(thick+thkcor).le.0.14313d0.and.
     +	statev(45).eq.0) then 	 
		statev(45)=1		
		datarec(5,1)=0.14313d0	
		datarec(5,2)=thick+thkcor
		datarec(5,3)=kinc	
		datarec(5,4)=kinc/5	

		call fprint(datarec)
		
	  endif	  

	  if (NOEL.EQ.EleNO.and.(thick+thkcor).le.0.12844d0.and.
     +	statev(46).eq.0) then 	 
		statev(46)=1		
		datarec(6,1)=0.12844d0	
		datarec(6,2)=thick+thkcor
		datarec(6,3)=kinc	
		datarec(6,4)=kinc/5	

	    call fprint(datarec) 
				
	  endif	 	  

	  if (NOEL.EQ.EleNO.and.(thick+thkcor).le.0.11882d0.and.
     +	statev(47).eq.0) then 	 
		statev(47)=1		
		datarec(7,1)=0.11882d0	
		datarec(7,2)=thick+thkcor
		datarec(7,3)=kinc	
		datarec(7,4)=kinc/5	

	    call fprint(datarec) 
				
	  endif		  

	  if (NOEL.EQ.EleNO.and.(thick+thkcor).le.0.11217d0.and.
     +	statev(48).eq.0) then 	 
		statev(48)=1		
		datarec(8,1)=0.11217d0	
		datarec(8,2)=thick+thkcor
		datarec(8,3)=kinc
		datarec(8,4)=kinc/5	

	    call fprint(datarec) 
				
	  endif		  

	  if (NOEL.EQ.EleNO.and.(thick+thkcor).le.0.10518d0.and.
     +	statev(49).eq.0) then 	 
		statev(49)=1		
		datarec(9,1)=0.10518d0	
		datarec(9,2)=thick+thkcor
		datarec(9,3)=kinc	
		datarec(9,4)=kinc/5	

	    call fprint(datarec) 
				
	  endif	

	  if (NOEL.EQ.EleNO.and.(thick+thkcor).le.0.096852d0.and.
     +	statev(50).eq.0) then 	 
		statev(50)=1		
		datarec(10,1)=0.096852d0
		datarec(10,2)=thick+thkcor
		datarec(10,3)=kinc	
		datarec(10,4)=kinc/5	

	    call fprint(datarec) 
				
	  endif	  

	  if (NOEL.EQ.EleNO.and.(thick+thkcor).le.0.091421d0.and.
     +	statev(51).eq.0) then 	 
		statev(51)=1		
		datarec(11,1)=0.091421d0
		datarec(11,2)=thick+thkcor
		datarec(11,3)=kinc
		datarec(11,4)=kinc/5	

	    call fprint(datarec) 
				
	  endif		  

	  if (NOEL.EQ.EleNO.and.(thick+thkcor).le.0.087008d0.and.
     +	statev(52).eq.0) then 	 
		statev(52)=1		
		datarec(12,1)=0.087008d0
		datarec(12,2)=thick+thkcor
		datarec(12,3)=kinc
		datarec(12,4)=kinc/5	

	    call fprint(datarec) 
				
	  endif			  

	  if (NOEL.EQ.EleNO.and.(thick+thkcor).le.0.080639d0.and.
     +	statev(53).eq.0) then 	 
		statev(53)=1		
		datarec(13,1)=0.080639d0
		datarec(13,2)=thick+thkcor
		datarec(13,3)=kinc
		datarec(13,4)=kinc/5	

	    call fprint(datarec) 
				
	  endif		  

	  if (NOEL.EQ.EleNO.and.(thick+thkcor).le.0.076394d0.and.
     +	statev(54).eq.0) then 	 
		statev(54)=1		
		datarec(14,1)=0.076394d0
		datarec(14,2)=thick+thkcor
		datarec(14,3)=kinc
		datarec(14,4)=kinc/5	

	    call fprint(datarec) 
				
	  endif		  

	  if (NOEL.EQ.EleNO.and.(thick+thkcor).le.0.068753d0.and.
     +	statev(55).eq.0) then 	 
		statev(55)=1		
		datarec(15,1)=0.068753d0
		datarec(15,2)=thick+thkcor
		datarec(15,3)=kinc
		datarec(15,4)=kinc/5	

	    call fprint(datarec) 
					
	  endif		  

	  if (NOEL.EQ.EleNO.and.(thick+thkcor).le.0.058392d0.and.
     +	statev(56).eq.0) then 	 
		statev(56)=1		
		datarec(16,1)=0.058392d0
		datarec(16,2)=thick+thkcor
		datarec(16,3)=kinc
		datarec(16,4)=kinc/5	

	    call fprint(datarec) 
				
	  endif		  

	  if (NOEL.EQ.EleNO.and.(thick+thkcor).le.0.054014d0.and.
     +	statev(57).eq.0) then 	 
		statev(57)=1		
		datarec(17,1)=0.054014d0
		datarec(17,2)=thick+thkcor
		datarec(17,3)=kinc
		datarec(17,4)=kinc/5	

	    call fprint(datarec) 
				
	  endif	

	  if (NOEL.EQ.EleNO.and.(thick+thkcor).le.0.048949d0.and.
     +	statev(58).eq.0) then 	 
		statev(58)=1		
		datarec(18,1)=0.048949d0
		datarec(18,2)=thick+thkcor
		datarec(18,3)=kinc
		datarec(18,4)=kinc/5	

	    call fprint(datarec) 
				
	  endif	  

	  if (NOEL.EQ.EleNO.and.(thick+thkcor).le.0.046396d0.and.
     +	statev(59).eq.0) then 	 
		statev(59)=1		
		datarec(19,1)=0.046396d0
		datarec(19,2)=thick+thkcor
		datarec(19,3)=kinc
		datarec(19,4)=kinc/5	

	    call fprint(datarec) 
				
	  endif	  

	  if (NOEL.EQ.EleNO.and.(thick+thkcor).le.0.043092d0.and.
     +	statev(60).eq.0) then 	 
		statev(60)=1		
		datarec(20,1)=0.043092d0
		datarec(20,2)=thick+thkcor
		datarec(20,3)=kinc
		datarec(20,4)=kinc/5	

	    call fprint(datarec) 
				
	  endif	

	  if (NOEL.EQ.EleNO.and.(thick+thkcor).le.0.039124d0.and.
     +	statev(61).eq.0) then 	 
		statev(61)=1		
		datarec(21,1)=0.039124d0
		datarec(21,2)=thick+thkcor
		datarec(21,3)=kinc
		datarec(21,4)=kinc/5	

	    call fprint(datarec) 
				
	  endif	

	  if (NOEL.EQ.EleNO.and.(thick+thkcor).le.0.037096d0.and.
     +	statev(62).eq.0) then 	 
		statev(62)=1		
		datarec(22,1)=0.037096d0
		datarec(22,2)=thick+thkcor
		datarec(22,3)=kinc
		datarec(22,4)=kinc/5	

	    call fprint(datarec) 
				
	  endif	

	  if (NOEL.EQ.EleNO.and.(thick+thkcor).le.0.033383d0.and.
     +	statev(63).eq.0) then 	 
		statev(63)=1		
		datarec(23,1)=0.033383d0
		datarec(23,2)=thick+thkcor
		datarec(23,3)=kinc
		datarec(23,4)=kinc/5	

	    call fprint(datarec) 
				
	  endif	

	  if (NOEL.EQ.EleNO.and.(thick+thkcor).le.0.03107d0.and.
     +	statev(64).eq.0) then 	 
		statev(64)=1		
		datarec(24,1)=0.03107d0
		datarec(24,2)=thick+thkcor
		datarec(24,3)=kinc
		datarec(24,4)=kinc/5	

	    call fprint(datarec) 
				
	  endif	

	  if (NOEL.EQ.EleNO.and.(thick+thkcor).le.0.030553d0.and.
     +	statev(65).eq.0) then 	 
		statev(65)=1		
		datarec(25,1)=0.030553d0
		datarec(25,2)=thick+thkcor
		datarec(25,3)=kinc
		datarec(25,4)=kinc/5	

	    call fprint(datarec) 
				
	  endif

	  if (NOEL.EQ.EleNO.and.(thick+thkcor).le.0.028457d0.and.
     +	statev(66).eq.0) then 	 
		statev(66)=1		
		datarec(26,1)=0.028457d0
		datarec(26,2)=thick+thkcor
		datarec(26,3)=kinc
		datarec(26,4)=kinc/5	

	    call fprint(datarec) 
				
	  endif
	  
C**********************************************************************  

	  endif ! END OF SAMPLE MATERIAL 
	  
      RETURN
      END

C**********************************************************************


	  SUBROUTINE NR4(noel, kinc,
     1  xl, dtime, sv7, sv9, sv15, sv16,
     2 	sv18, sv20, xbe, stress1, xq, xlan, xc, xUIne, sy, pr1,
     3  XVe, xLagEe, xEplastic, CBe, DFGRD1, xUIold, xEplastic2,
     4	xitermax, ElMaxXiter, ELdetFMax, ELdetFeMax, ELdetFMin, 
     5	ELdetFeMin, detFMax, detFeMax, detFMin, detFeMin, kstep,   
     6	ddsdde, iter_marker,     	    	 
     7	XELAM0,XEG0,Xlpar0,Xmpar0,Xnpar0,
     8	XELAM1,XEG1,Xlpar1,Xmpar1,Xnpar1,
     9	tstr,sya,syap,syb,sybp,Xpdh,Xpde,xk,Delta1,Delta2,Delta3,xB,q0)
	 
	  implicit none

	  Integer xitermax,ElMaxXiter
	  
	  real*8 detFMax,detFeMax
	  Integer ELdetFMax,ELdetFeMax
		
	  real*8 detFMin,detFeMin
	  Integer ELdetFMin,ELdetFeMin 

	  real*8 sv7, sv9, sv15, sv16, sv18, sv20, dtime
	  real*8 ELAM, EG, XLP, XMP, XNP
	  real*8 XELAM0, XEG0, Xlpar0, Xmpar0, Xnpar0 
	  real*8 XELAM1, XEG1, Xlpar1, Xmpar1,Xnpar1
	  
	  real*8 xk, tstr
	  real*8 XpdH,Xpde,XprH,Xpre,sya,syap,syb,sybp	  
	  
	  INTEGER GoOn,Xiter, noel, kinc, i, j, k, l, f, kstep

	  real*8 xbe(3,3), Stress1(3,3), dev1(3,3), xiden(3,3)
	  real*8 xi4(3,3,3,3), BeBe(3,3), dI3dBe(3,3), dEdBe1(3,3) 
	  real*8 d2I3dBe2(3,3,3,3), d2EdBe2(3,3,3,3), xdJeInvdBe(3,3) 

	  Real*8 xc,xq,xlan,BeI1,BeI2,XdetFesq,XdetFe, 
     1    pc1,pr1,sy0,sy1,sy

! g1	  
	  Real*8 pc2,xsigybar
	  real*8 xcdot
	  
	  Real*8 xl(3,3),dev2(3,3),dg1dBe(3,3,3,3),
     +    dg1dS(3,3,3,3),dg1dlan(3,3),dg1dq(3,3),xdg1dc(3,3)
	 
	  real*8 dg1dBe1111,dg1dBe1122,dg1dBe1133,dg1dBe1112,
     1     dg1dBe2211,dg1dBe2222,dg1dBe2233,dg1dBe2212,
     2     dg1dBe3311,dg1dBe3322,dg1dBe3333,dg1dBe3312,
     3     dg1dBe1211,dg1dBe1222,dg1dBe1233,dg1dBe1212	 
	 	 	 
	  real*8 dg1dS1111,dg1dS1122,dg1dS1133,dg1dS1112,
     1     dg1dS2211,dg1dS2222,dg1dS2233,dg1dS2212,
     2     dg1dS3311,dg1dS3322,dg1dS3333,dg1dS3312,
     3     dg1dS1211,dg1dS1222,dg1dS1233,dg1dS1212
	  
	  real*8 dg1dlan11,dg1dlan22,dg1dlan33,dg1dlan12,
     1     dg1dq11,dg1dq22,dg1dq33,dg1dq12	
	 
      real*8 xdg1dc11,xdg1dc22,xdg1dc33,xdg1dc12

! g2	  
	  real*8 xlamdif,xGdif,xldif,xmdif,xndif	  
	  	  
	  real*8 dg2dBe(3,3,3,3),dg2dlan(3,3),dg2dS(3,3,3,3),
     +    dg2dq(3,3),xdg2dc(3,3)
	 
	  real*8 dg2dS1111,dg2dS1122,dg2dS1133,dg2dS1112,
     1     dg2dS2211,dg2dS2222,dg2dS2233,dg2dS2212,
     2     dg2dS3311,dg2dS3322,dg2dS3333,dg2dS3312,
     3     dg2dS1211,dg2dS1222,dg2dS1233,dg2dS1212

	  real*8 dg2dlan11,dg2dlan22,dg2dlan33,dg2dlan12,
     1     dg2dq11,dg2dq22,dg2dq33,dg2dq12 
	 
	  real*8 dg2dBe1111,dg2dBe1122,dg2dBe1133,dg2dBe1112,
     1     dg2dBe2211,dg2dBe2222,dg2dBe2233,dg2dBe2212,
     2     dg2dBe3311,dg2dBe3322,dg2dBe3333,dg2dBe3312,
     3     dg2dBe1211,dg2dBe1222,dg2dBe1233,dg2dBe1212
	
	  real*8 xdg2dc11,xdg2dc22,xdg2dc33,xdg2dc12	

! g3
      real*8 dg3dS(3,3),dg3dBe(3,3)	  
	  real*8 dg3dBe11,dg3dBe22,dg3dBe33,dg3dBe12,
     1       dg3dS11,dg3dS22,dg3dS33,dg3dS12
	  real*8 dg3dq,dg3dlan,xdg3dc
	  
! g4	  
	  real*8 dg4ds(3,3),dg4dBe(3,3)	  
	  real*8 dg4dBe11,dg4dBe22,dg4dBe33,dg4dBe12,
     1       dg4dS11,dg4dS22,dg4dS33,dg4dS12	 
      real*8 dg4dlan,dg4dq,xdg4dc,Xdg4_dlanda,Xdg4_dq

! g5	   
	  real*8 Xdg5dBe(3,3),Xdg5ds(3,3),dy1dybar(3,3),dy2dybar(3,3)	  
	  real*8 xdg5dBe11,xdg5dBe22,xdg5dBe33,xdg5dBe12,
     1    xdg5dS11,xdg5dS22,xdg5dS33,xdg5dS12,
     2    xdg5dlan,xdg5dc,xdg5dq,
     3    Xpd,Xpr
	 
	  real*8 XHevD,XHevR
	  
! for finding new variables	 
	  real*8 XJac(11,11),XJacI(11,11),Be_0(3,3),
     1   Xg1(3,3),Xg2(3,3),Xg(11,1),XVar(11,1),XVarN(11,1),
     2   yafter1(3,3),yafter2(3,3),yaf(11,1)	  
	  
	  real*8 yafm1	 

	  real*8 xeigvalvecBe(3),xeigvecmatBe(3,3),xeigvecmatBeT(3,3),
     1    xeigvalmatVe(3,3),xSVDBeVe(3,3),XVe_I(3,3),xReUIne(3,3),
     2 	  XUidot(3,3),XUiInv(3,3),xUidotUiInv(3,3),xUidotUiInvT(3,3),
     3    xReUidotUiInvSym(3,3),XReT(3,3),xDIne(3,3),XDp(3,3)

	  Real*8 xUtrans(3,3),xUtransInv(3,3),xUplastic(3,3)
	  Real*8 xU2plastic(3,3), xEplastic(3,3) 

	  real*8 xve(3,3), DFGRD1(3,3), XRe(3,3),xUIne(3,3)
	  real*8 xdetUIne, xUIold(3,3)

	  real*8 YBe(3,3,3,3),CBe(3,3,3,3)
	  Real*8 xndir(3,3),YBen(3,3),dphidsig(3,3),dphidSYBe(3,3)
	  real*8 hden, xLagEe(3,3), dFdBe(3,3,3,3), detF, xReTBe(3,3)	

	  real*8 Delta1, Delta2, Delta3, xB
	  
	  !!!logical :: violated
	  !!!real*8 penalty_param
	  !!!integer :: update_interval

	  real*8 xjacM(10,10),XgM(10,1),xvarM(10,1),XVarNM(10,1),yafM(10,1)
	  real*8 XJacIM(10,10)

	  real*8 xEplastic2(3,3), xF_P(3,3), xF_PT(3,3), xF_pTF_p(3,3) 
	  real*8 xdetUplastic, xdetReUIne
	  
C***************************************************************************

	  real*8 xdG2,xdL2,XPJ6,xndi,xmm
	  real*8 xfdc2(3,3),V1eV1e(3,3),xFI(3,3,3,3),dfidt(3,3),xnkm(3,3),
     1  beita(3,3),Fplast(3,3,3,3),bmft(3,3,3,3), ddsdde(6,6)
	 
	  real*8 xIone1, xac
	  integer M, N, K1, K2	 
	  
	  integer iter_marker

	  real*8 dfdc(3,3),Zfac1(3,3),Zfac2(3,3),Zfac3(3,3),Zfac(3,3)
	  real*8 Afac,Hden_new

	  real*8 q0
	  
C***************************************************************************
	  
	  !!!update_interval = 5
	  !!!penalty_param = 0.001

	  call xidentity(xiden)	 

	  call xidentity4(xi4)

      CALL mdet(DFGRD1,detF)
	  
      GoOn=1.0d0
      Xiter=0 

! Using the goto statement
	  if (xc.ge.1.0d0) then
		xc = 1.0d0
		goto 167
	  endif
	  
      DO WHILE (GoOn.GT.0.0d0)

! Gradients for the Newton Raphson	  		  
      Call KDEVIA(stress1,dev1)	  
	  Call KEFFP(dev1,pc1)		
	  pc2=2.0d0/3.0d0*(pc1**2)
	  dev2=dev1/DSQRT(pc2)   
	  call pressure(Stress1,pr1)  
	  
C**********************************************************************
! g1		  	

! dg1dBe  
	  xcdot=(xc-sv9)/dtime	
	  
	  Sy0=SyA+SyAP*pr1
	  sy1=SyB+SyBP*pr1
      Sy=(1.0d0-xc)*sy0+xc*sy1  
	  
	  ! ! ! xsigybar=xc*sy0+(1.0d0-xc)*sy1	 
	  xsigybar=xc*SyA+(1.0d0-xc)*SyB	! To make pressure-independent yield strength in kinetic eqn
              
	  dg1dBe=0.0d0		
      DO 202 I=1,3
       DO 202 J=1,3
        DO 202 K=1,3    
         DO 202 L=1,3

	  dg1dBe(i,j,k,l)=(xi4(i,j,k,l)+xi4(j,i,k,l))/2.0d0-
     1  (xl(I,K)*xiden(J,L)+xl(I,L)*xiden(J,k)+ 
     2  xl(j,k)*xiden(i,l)+xl(j,l)*xiden(i,k))*DTIME*0.5d0
     3  +DTIME*(xiden(i,k)*(xlan*dev2(l,j)+tstr*xiden(l,j)*xcdot/3.0d0)
     4  +xiden(i,l)*(xlan*dev2(k,j)+tstr*xiden(k,j)*xcdot/3.0d0))
202	  continue  

	  call xsymmetric4_1(dg1dBe)
	  call xsymmetric4_2(dg1dBe)  

	  dg1dBe1111=dg1dBe(1,1,1,1)
	  dg1dBe1122=dg1dBe(1,1,2,2)
	  dg1dBe1133=dg1dBe(1,1,3,3)
	  dg1dBe1112=dg1dBe(1,1,1,2)+dg1dBe(1,1,2,1)
	  
	  dg1dBe2211=dg1dBe(2,2,1,1)
	  dg1dBe2222=dg1dBe(2,2,2,2)
	  dg1dBe2233=dg1dBe(2,2,3,3)
	  dg1dBe2212=dg1dBe(2,2,1,2)+dg1dBe(2,2,2,1)
	  
	  dg1dBe3311=dg1dBe(3,3,1,1)
	  dg1dBe3322=dg1dBe(3,3,2,2)
	  dg1dBe3333=dg1dBe(3,3,3,3)
	  dg1dBe3312=dg1dBe(3,3,1,2)+dg1dBe(3,3,2,1)

	  dg1dBe1211=dg1dBe(1,2,1,1)
	  dg1dBe1222=dg1dBe(1,2,2,2)
	  dg1dBe1233=dg1dBe(1,2,3,3)
	  dg1dBe1212=dg1dBe(1,2,1,2)+dg1dBe(1,2,2,1)

! dg1dlan
      dg1dlan=0.0d0   
      DO 203 I=1,3
       DO 203 J=1,3
        DO 203 L=1,3
         dg1dlan(I,J)=dg1dlan(I,J)+DTIME*
     +               (2.0d0*XBe(I,L)+xiden(I,L))*dev2(L,J)     
203   CONTINUE	  		  
	  call xsymmetric2(dg1dlan)	  
 
	  dg1dlan11=dg1dlan(1,1)
	  dg1dlan22=dg1dlan(2,2)
	  dg1dlan33=dg1dlan(3,3)
	  dg1dlan12=dg1dlan(1,2)          

! dg1dS
      dg1dS=0.0D0     
      DO 204 I=1,3
       DO 204 J=1,3
        DO 204 K=1,3    
         DO 204 L=1,3   
		 
		  DO F=1,3
	       dg1dS(I,J,K,L)=dg1dS(I,J,K,L)+
     +      DTIME*xlan*(2.0d0*XBe(I,F)+xiden(I,F))/DSQRT(pc2)*
     +      (-1.0d0/3.0d0*xi4(F,K,J,L)+0.5d0*xi4(F,J,K,L)
     +      +0.5d0*xi4(J,F,K,L)-Dev1(F,J)*Dev1(L,K)/pc2) 	 
		  end do
		  
204   CONTINUE	

	  call xsymmetric4_1(dg1dS)
	  call xsymmetric4_2(dg1dS)

	  dg1dS1111=dg1dS(1,1,1,1)
	  dg1dS1122=dg1dS(1,1,2,2)
	  dg1dS1133=dg1dS(1,1,3,3)
	  dg1dS1112=dg1dS(1,1,1,2)+dg1dS(1,1,2,1)
	  
	  dg1dS2211=dg1dS(2,2,1,1)
	  dg1dS2222=dg1dS(2,2,2,2)
	  dg1dS2233=dg1dS(2,2,3,3)
	  dg1dS2212=dg1dS(2,2,1,2)+dg1dS(2,2,2,1)
	  
	  dg1dS3311=dg1dS(3,3,1,1)
	  dg1dS3322=dg1dS(3,3,2,2)
	  dg1dS3333=dg1dS(3,3,3,3)
	  dg1dS3312=dg1dS(3,3,1,2)+dg1dS(3,3,2,1)
	  
	  dg1dS1211=dg1dS(1,2,1,1)
	  dg1dS1222=dg1dS(1,2,2,2)
	  dg1dS1233=dg1dS(1,2,3,3)
	  dg1dS1212=dg1dS(1,2,1,2)+dg1dS(1,2,2,1)

! dg1dq
	  dg1dq=0.0d0
	  dg1dq11=0.0d0;dg1dq22=0.0d0;dg1dq33=0.0d0;dg1dq12=0.0d0
	  call xsymmetric2(dg1dq)
	  
! dg1dc
	  xdg1dc=0.0d0 
      xdg1dc=(2.0d0*XBe+xiden)*tstr/3.0d0	  	  
	  call xsymmetric2(xdg1dc)
	
	  xdg1dc11=xdg1dc(1,1)
	  xdg1dc22=xdg1dc(2,2)
	  xdg1dc33=xdg1dc(3,3)
	  xdg1dc12=xdg1dc(1,2)

C********************************************************************** 
! g2  

! dg2dlan
	  dg2dlan=0.0d0
	  dg2dlan11=0.0d0;dg2dlan22=0.0d0;dg2dlan33=0.0d0;dg2dlan12=0.0d0
	  call xsymmetric2(dg2dlan)

! dg2dq
	  dg2dq=0.0d0	  
	  dg2dq11=0.0d0;dg2dq22=0.0d0;dg2dq33=0.0d0;dg2dq12=0.0d0
	  call xsymmetric2(dg2dq)

! dg2dS
	  dg2dS=0.0d0
      DO 205 I=1,3 
	   DO 205 J=1,3
		DO 205 K=1,3   
         DO 205 L=1,3
          dg2dS(I,J,K,L)=0.5d0*(xi4(I,J,K,L)+xi4(j,i,K,L))
205   CONTINUE	
	  
	  call xsymmetric4_1(dg2dS)
	  call xsymmetric4_2(dg2dS)
	  
! Components for the Jacobian
	  dg2dS1111=dg2dS(1,1,1,1)
	  dg2dS1122=dg2dS(1,1,2,2)
	  dg2dS1133=dg2dS(1,1,3,3)
	  dg2dS1112=dg2dS(1,1,1,2)+dg2dS(1,1,2,1)
	  
	  dg2dS2211=dg2dS(2,2,1,1)
	  dg2dS2222=dg2dS(2,2,2,2)
	  dg2dS2233=dg2dS(2,2,3,3)
	  dg2dS2212=dg2dS(2,2,1,2)+dg2dS(2,2,2,1)
	  
	  dg2dS3311=dg2dS(3,3,1,1)
	  dg2dS3322=dg2dS(3,3,2,2)
	  dg2dS3333=dg2dS(3,3,3,3)
	  dg2dS3312=dg2dS(3,3,1,2)+dg2dS(3,3,2,1)
	  
	  dg2dS1211=dg2dS(1,2,1,1)
	  dg2dS1222=dg2dS(1,2,2,2)
	  dg2dS1233=dg2dS(1,2,3,3)
	  dg2dS1212=dg2dS(1,2,1,2)+dg2dS(1,2,2,1)
	
! dg2dBe
	  call xinvariant1(XBe,BeI1)
	  call xinvariant2(XBe,BeI2)
	  
	  call kmlt(XBe,XBe,BeBe)

      ELAM=(1.0d0-xc)*XELAM0+xc*XELAM1
      EG=(1.0d0-xc)*XEG0+xc*XEG1
      XLP=(1.0d0-xc)*Xlpar0+xc*Xlpar1
      XMP=(1.0d0-xc)*Xmpar0+xc*Xmpar1
      XNP=(1.0d0-xc)*Xnpar0+xc*Xnpar1

	  dI3dBe=0.0d0 
	  dI3dBe=BeBe-BeI1*xbe+BeI2*xiden
	  call xsymmetric2(dI3dBe)	  
	  	  
	  dEdBe1=0.0d0 
      dEdBe1=ELAM*BeI1*XIDEN+2.0d0*EG*XBe
     1            +(XLP*BeI1**2-2.0d0*XMP*BeI2)*XIDEN
     2            +XNP*dI3dBe+2.0d0*XMP*BeI1*XBe 
	  call xsymmetric2(dEdBe1)

	  d2I3dBe2=0.0d0
      DO 206 I=1,3 
	   DO 206 J=1,3
		DO 206 K=1,3   
         DO 206 L=1,3
		   
	  d2I3dBe2(i,j,k,l)=0.5d0*(xbe(i,k)*xiden(j,l)+xbe(i,l)*xiden(j,k)+
     +    xbe(l,j)*xiden(i,k)+xbe(k,j)*xiden(i,l))-
     +    xbe(i,j)*xiden(k,l)-xbe(k,l)*xiden(i,j)+
     +    xiden(i,j)*xiden(k,l)*BeI1
     +    -0.5d0*BeI1*(xi4(i,j,k,l)+xi4(j,i,k,l))
206   CONTINUE	

	  call xsymmetric4_1(d2I3dBe2)
	  call xsymmetric4_2(d2I3dBe2)

	  d2EdBe2=0.0d0
      DO 207 I=1,3 
	   DO 207 J=1,3
		DO 207 K=1,3   
         DO 207 L=1,3
		   
            d2EdBe2(i,j,k,l)=elam*xi4(i,k,j,l)+EG*(xi4(i,j,k,l)+
     + 		xi4(j,i,k,l))+2.0d0*XLP*BeI1*xi4(i,k,j,l)+XMP*BeI1*
     + 		(xi4(i,j,k,l)+xi4(j,i,k,l))+2.0d0*XMP*xiden(k,l)*xbe(i,j)+
     +  	2.0d0*xmp*xiden(i,j)*xbe(k,l)-2.0d0*xmp*BeI1*xi4(i,k,j,l)+
     + 		xnp*d2I3dBe2(i,j,k,l)	 
207   CONTINUE	

	  call xsymmetric4_1(d2EdBe2)
	  call xsymmetric4_2(d2EdBe2)
	  
	  CALL mdet(2.0d0*XBe+xiden,XdetFesq)
	  XdetFe=dsqrt(XdetFesq)	  

	  if (isNaN(XdetFe)) then
		PRINT*,'DetFe < 0 (plastic_1_1)',kinc,noel,xiter,xc,sv9
		call xit
	  endif 
	  
	  xdJeInvdBe=0.0d0
	  xdJeInvdBe=-(1.0d0/xDetFe**3)*
     +    (xiden+2.0d0*(BeI1*xiden-xbe)+4.0d0*dI3dBe)	
	  call xsymmetric2(xdJeInvdBe)

	  dg2dBe=0.0d0
      DO 208 I=1,3 
	   DO 208 J=1,3
		DO 208 K=1,3   
         DO 208 L=1,3
		 
		  DO F=1,3		   
        dg2dBe(i,j,k,l)=dg2dBe(i,j,k,l)-
     + 	xdJeInvdBe(k,l)*(2.0d0*xbe(i,f)+xiden(i,f))*dEdBe1(f,j)	
     +	-1.0d0/xdetFe*(xi4(i,f,k,l)+xi4(f,i,k,l))*dEdBe1(f,j)	
     + 	-1.0d0/xdetFe*(2.0d0*xbe(i,f)+xiden(i,f))*d2EdBe2(f,j,k,l)	 
		  end do
		  
208   CONTINUE

	  call xsymmetric4_1(dg2dBe)
	  call xsymmetric4_2(dg2dBe)
	  
! Components for the Jacobian
	  dg2dBe1111=dg2dBe(1,1,1,1)
	  dg2dBe1122=dg2dBe(1,1,2,2)
	  dg2dBe1133=dg2dBe(1,1,3,3)
	  dg2dBe1112=dg2dBe(1,1,1,2)+dg2dBe(1,1,2,1)
	  
	  dg2dBe2211=dg2dBe(2,2,1,1)
	  dg2dBe2222=dg2dBe(2,2,2,2)
	  dg2dBe2233=dg2dBe(2,2,3,3)
	  dg2dBe2212=dg2dBe(2,2,1,2)+dg2dBe(2,2,2,1)

	  dg2dBe3311=dg2dBe(3,3,1,1)
	  dg2dBe3322=dg2dBe(3,3,2,2)
	  dg2dBe3333=dg2dBe(3,3,3,3)
	  dg2dBe3312=dg2dBe(3,3,1,2)+dg2dBe(3,3,2,1)

	  dg2dBe1211=dg2dBe(1,2,1,1)
	  dg2dBe1222=dg2dBe(1,2,2,2)
	  dg2dBe1233=dg2dBe(1,2,3,3)
	  dg2dBe1212=dg2dBe(1,2,1,2)+dg2dBe(1,2,2,1)	          
		  
! dg2dc
	  xlamdif=XELAM1-XELAM0
	  xGdif=XEG1-XEG0
	  xldif=Xlpar1-Xlpar0
	  xmdif=Xmpar1-Xmpar0
	  xndif=Xnpar1-Xnpar0
	  
	  xdg2dc=0.0D0  
      DO 209 I=1,3 
	   DO 209 J=1,3
	   
		DO K=1,3 	  		  
      xdg2dc(I,J)=xdg2dc(I,J)-(1.0/XdetFe)*(2.0*XBe(I,K)+xiden(I,K))*
     +(xlamdif*BeI1*xiden(K,J)+2.0d0*xGdif*XBe(K,J)+
     +xldif*(BeI1**2)*xiden(K,J)+
     +2.0d0*xmdif*(BeI1*XBe(K,J)-BeI2*xiden(K,J))+ xndif*dI3dBe(K,J))	 
		end do
		
209   CONTINUE   
	  call xsymmetric2(xdg2dc)
	  
	  xdg2dc11=xdg2dc(1,1)
	  xdg2dc22=xdg2dc(2,2)
	  xdg2dc33=xdg2dc(3,3)
	  xdg2dc12=xdg2dc(1,2)

C**********************************************************************
! g3

! dg3dBe
	  dg3dBe=0.0d0
	  dg3dBe11=0.0d0;dg3dBe22=0.0d0;dg3dBe33=0.0d0;dg3dBe12=0.0d0
	  call xsymmetric2(dg3dBe)

! dg3dS	  
      dg3dS=0.0d0
	  dg3dS11=0.0d0;dg3dS22=0.0d0;dg3dS33=0.0d0;dg3dS12=0.0d0
	  call xsymmetric2(dg3dS)

! dg3dlan  
      dg3dlan=-((2.0d0/3.0d0)**0.5)*DTIME

! dg3dq	  
      dg3dq=1.0d0

! dg3dc	  
	  xdg3dc=0.0d0

C**********************************************************************
! g4

! dg4dBe
	  dg4dBe=0.0d0
	  dg4dBe11=0.0d0;dg4dBe22=0.0d0;dg4dBe33=0.0d0;dg4dBe12=0.0d0
	  call xsymmetric2(dg4dBe)

! dg4dlan 	  
	  dg4dlan=0.0d0

! dg4dq	  
	  dg4dq=0.0d0

! dg4dS	  
	  dg4ds=0.0d0
	  dg4dS=(1.5d0**0.5)*dev2 + 
     +	1.0d0/3.0d0*xiden*((1.0d0-xc)*SyAP+xc*SyBP)	  
	  call xsymmetric2(dg4dS)	
	  
	  dg4dS11=dg4dS(1,1)
	  dg4dS22=dg4dS(2,2)
	  dg4dS33=dg4dS(3,3)
	  dg4dS12=dg4dS(1,2)+dg4dS(2,1)

! dg4dc
	  xdg4dc=sy0-sy1	  

C**********************************************************************
! g5

! dg5dBe
	  xdg5dBe=0.0d0 
	  xdg5dBe11=0.0d0;xdg5dBe22=0.0d0;xdg5dBe33=0.0d0;xdg5dBe12=0.0d0;
	  call xsymmetric2(xdg5dBe)
	  
! dg5dlan	  
	  xdg5dlan=0.0D0  

! dg5dS	
      IF (pr1.Gt.Xpde) THEN
       XHevD=1.0d0
      ELSE
       XHevD=0.0d0
      ENDIF
      Xpd=(pr1-Xpde)/(XpdH-Xpde)
	  	  
      IF (pr1.gt.Xpre) THEN
       XHevR=0.0d0 ! This should be changed if reverse PT is there.
      ELSE
       XHevR=0.0d0
      ENDIF
      Xpr=(pr1-Xpre)/(XprH-Xpre)	

	  dy2dybar=-SyBP/(3.0d0*xsigybar)*xiden+
     +  (sy1/xsigybar**2)*(1.0d0/3.0d0)*(xc*SyAP+(1.0d0-xc)*SyBP)*xiden
	 
	  dy1dybar=-SyAP/(3.0d0*xsigybar)*xiden+
     +  (sy0/xsigybar**2)*(1.0d0/3.0d0)*(xc*SyAP+(1.0d0-xc)*SyBP)*xiden
 
	  xdg5ds=0.0d0
	 
      xdg5ds=xiden*xk*(1.0d0+xB*xq)/3.0d0*(Xq-Sv7)*syB/xsigybar*XHevD*
     +  (1.0+Delta1-Delta2*xc)*	  
     +  (Xpd*Delta3 + (1.0+Delta3*pr1)/(XpdH-Xpde))
	 
	  call xsymmetric2(xdg5ds)	 
	  
	  xdg5dS11=xdg5ds(1,1);xdg5dS22=xdg5ds(2,2);xdg5dS33=xdg5ds(3,3)
	  xdg5dS12=xdg5ds(1,2)+xdg5ds(2,1)
 
! dg5dq  
	 
      xdg5dq=-xk*(1.0d0+xB*xq)/xsigybar*
     +	(1.0d0+Delta1-Delta2*Xc)*(1.0d0+Delta3*pr1)*Xpd*XHevD*syB - 	 
     +	(1.0+Delta1-Delta2*Xc)*(1.0+Delta3*pr1)*Xpd*XHevD*syB*(Xq-Sv7)*
     +	xk*(xB)/xsigybar	 
	 
! dg5dc	 
	 
      xdg5dc=1.0d0 + xk*(1.0d0+xB*xq)*(Xq-Sv7)/xsigybar*
     + (1.0d0+Delta3*pr1)*Xpd*XHevD*syB*(Delta2 + 
     + (1.0d0+Delta1-Delta2*xc)*
     + (syA-syB)/xsigybar)
	 
C**********************************************************************  
! Initializing the Jacobian matrix 		   
      XJac=0.0d0             

! g1
	  xjac(1,1)=dg1dBe1111;xjac(1,2)=dg1dBe1122
	  xjac(1,3)=dg1dBe1133;xjac(1,7)=dg1dBe1112
	  
	  xjac(2,1)=dg1dBe2211;xjac(2,2)=dg1dBe2222
	  xjac(2,3)=dg1dBe2233;xjac(2,7)=dg1dBe2212
	  
	  xjac(3,1)=dg1dBe3311;xjac(3,2)=dg1dBe3322
	  xjac(3,3)=dg1dBe3333;xjac(3,7)=dg1dBe3312
	  	  
	  xjac(7,1)=dg1dBe1211;xjac(7,2)=dg1dBe1222
	  xjac(7,3)=dg1dBe1233;xjac(7,7)=dg1dBe1212;
	  
	  xjac(1,4)=dg1dS1111;xjac(1,5)=dg1dS1122;xjac(1,6)=dg1dS1133
	  xjac(1,8)=dg1dS1112
	  
	  xjac(2,4)=dg1dS2211;xjac(2,5)=dg1dS2222;xjac(2,6)=dg1dS2233
	  xjac(2,8)=dg1dS2212
	  
	  xjac(3,4)=dg1dS3311;xjac(3,5)=dg1dS3322;xjac(3,6)=dg1dS3333
	  xjac(3,8)=dg1dS3312
	  
	  xjac(7,4)=dg1dS1211;xjac(7,5)=dg1dS1222;xjac(7,6)=dg1dS1233	  
	  xjac(7,8)=dg1dS1212
	  
	  xjac(1,9)=dg1dq11;xjac(2,9)=dg1dq22;xjac(3,9)=dg1dq33
	  xjac(7,9)=dg1dq12
	  
	  xjac(1,10)=dg1dlan11;xjac(2,10)=dg1dlan22
	  xjac(3,10)=dg1dlan33;xjac(7,10)=dg1dlan12
	  
	  xjac(1,11)=xdg1dc11;xjac(2,11)=xdg1dc22;xjac(3,11)=xdg1dc33
	  xjac(7,11)=xdg1dc12	  

! g2	  	 
	  xjac(4,1)=dg2dBe1111;xjac(4,2)=dg2dBe1122
	  xjac(4,3)=dg2dBe1133;xjac(4,7)=dg2dBe1112

	  xjac(5,1)=dg2dBe2211;xjac(5,2)=dg2dBe2222
	  xjac(5,3)=dg2dBe2233;xjac(5,7)=dg2dBe2212

	  xjac(6,1)=dg2dBe3311;xjac(6,2)=dg2dBe3322
	  xjac(6,3)=dg2dBe3333;xjac(6,7)=dg2dBe3312		  
	  
	  xjac(8,1)=dg2dBe1211;xjac(8,2)=dg2dBe1222
	  xjac(8,3)=dg2dBe1233;xjac(8,7)=dg2dBe1212
 
	  xjac(4,4)=dg2dS1111;xjac(4,5)=dg2dS1122;xjac(4,6)=dg2dS1133
	  xjac(4,8)=dg2dS1112
	  
	  xjac(5,4)=dg2dS2211;xjac(5,5)=dg2dS2222;xjac(5,6)=dg2dS2233
	  xjac(5,8)=dg2dS2212

	  xjac(6,4)=dg2dS3311;xjac(6,5)=dg2dS3322;xjac(6,6)=dg2dS3333
	  xjac(6,8)=dg2dS3312

	  xjac(8,4)=dg2dS1211;xjac(8,5)=dg2dS1222;xjac(8,6)=dg2dS1233
	  xjac(8,8)=dg2dS1212
	  
	  xjac(4,9)=dg2dq11;xjac(5,9)=dg2dq22;xjac(6,9)=dg2dq33
	  xjac(8,9)=dg2dq12
	  
	  xjac(4,10)=dg2dlan11;xjac(5,10)=dg2dlan22
	  xjac(6,10)=dg2dlan33;xjac(8,10)=dg2dlan12
	  
	  xjac(4,11)=xdg2dc11;xjac(5,11)=xdg2dc22;xjac(6,11)=xdg2dc33
	  xjac(8,11)=xdg2dc12

! g3	  
	  xjac(9,1)=dg3dBe11;xjac(9,2)=dg3dBe22;xjac(9,3)=dg3dBe33
	  xjac(9,4)=dg3dS11;xjac(9,5)=dg3dS22;xjac(9,6)=dg3dS33
	  xjac(9,7)=dg3dBe12;xjac(9,8)=dg3dS12;xjac(9,9)=dg3dq
	  xjac(9,10)=dg3dlan;xjac(9,11)=xdg3dc;

! g4
	  xjac(10,1)=dg4dBe11;xjac(10,2)=dg4dBe22;xjac(10,3)=dg4dBe33
	  xjac(10,4)=dg4dS11;xjac(10,5)=dg4dS22;xjac(10,6)=dg4dS33
	  xjac(10,7)=dg4dBe12;xjac(10,8)=dg4dS12;xjac(10,9)=dg4dq
	  xjac(10,10)=dg4dlan;xjac(10,11)=xdg4dc;

! g5	  
	  xjac(11,1)=xdg5dBe11;xjac(11,2)=xdg5dBe22;xjac(11,3)=xdg5dBe33
	  xjac(11,4)=xdg5dS11;xjac(11,5)=xdg5dS22;xjac(11,6)=xdg5dS33
	  xjac(11,7)=xdg5dBe12;xjac(11,8)=xdg5dS12;xjac(11,9)=xdg5dq
	  xjac(11,10)=xdg5dlan;xjac(11,11)=xdg5dc;  

C**********************************************************************	  

! The converged solution Be from previous time-step  
	  Be_0=0.0d0
	  be_0(1,1) = sv15; be_0(2,2) = sv18; be_0(3,3) = sv20
	  be_0(1,2) = sv16; be_0(2,1) = be_0(1,2) 

! Old functional values

! g1 
	  xg1=0.0d0
      DO 211 I=1,3
       DO 211 J=1,3
	   
        DO L=1,3
			Xg1(I,J)=Xg1(I,J)-(xl(I,L)*XBe(L,J)+XBe(I,L)*xl(J,L))
     +  +(2.0d0*XBe(I,L)+XIDEN(I,L))*(xlan*dev2(L,J)
     +  +tstr*XIDEN(L,J)*xcdot/3.0d0)
		enddo
			
        Xg1(I,J)=Xg1(I,J)*DTIME+XBe(I,J)-Be_0(I,J)-DTIME*(xl(I,J)
     +                 +xl(J,I))/2.0d0
211   CONTINUE

! g2   	  
	  Xg2=0.0d0
      DO 212 I=1,3
       DO 212 J=1,3
        
		DO L=1,3
         Xg2(I,J)=Xg2(I,J)-(2.0d0*xbe(i,l)+XIDEN(I,L))*dEdBe1(L,J)
		enddo
		Xg2(I,J)=Xg2(I,J)/XdetFe + stress1(I,J)
		
212     CONTINUE
      
! Xg(1:11,1)
      DO I=1,3
         Xg(I,1)=Xg1(I,I)
		 Xg(I+3,1)=Xg2(I,I)
      ENDDO
	  
      Xg(7,1)=Xg1(1,2)
      Xg(8,1)=Xg2(1,2)
      
! g3
      Xg(9,1)=Xq-Sv7-((2.0d0/3.0d0)**0.5)*xlan*DTIME     	  
          
! g4   
      Xg(10,1)=pc1-Sy  
	  
! g5

      Xg(11,1)=(Xc-Sv9)-xk*(1.0d0+xB*xq)*(Xq-Sv7)/xsigybar*
     +  (1.0d0+Delta1-Delta2*Xc)*(1.0d0+Delta3*pr1)*Xpd*XHevD*syB 
	 	    
! Code portion for new variables using Newton Raphson 		  
      CALL matrixinv(XJac,XJacI,11)

! Initializing the solution vector for strain components with 
! trial elastic strain tensor & stress components with 
! trial stress tensor and also concentration,odqvist and lagrange parameter. 	  

! Old variables
	  DO I=1,3
       XVar(I,1)=XBe(I,I)
	   XVar(I+3,1)=Stress1(I,I)
      ENDDO		  

      XVar(7,1)=XBe(1,2)	  
      XVar(8,1)=Stress1(1,2)	  
      XVar(9,1)=Xq 
      XVar(10,1)=xlan
	  xvar(11,1)=xc	 

! New variables
	  xvarn=0.0d0
      DO I=1,11
        DO J=1,11
          XVarN(I,1)=XVarN(I,1)+XJacI(I,J)*Xg(J,1)
	    enddo
        XVarN(I,1)=XVar(I,1)-XVarN(I,1)
	  enddo

! Using the goto statement
	  if (XVarN(11,1).gt.1.0d0) then
		xc = 1.0d0
		goto 167
	  endif
	  
	  xbe=0.0d0
	  stress1=0.0d0	  

      DO I=1,3
       XBe(I,I)=XVarN(I,1)
	   Stress1(I,I)=XVarN(I+3,1)
      ENDDO
      XBe(1,2)=XVarN(7,1)
      XBe(2,1)=XBe(1,2)

      Stress1(1,2)=XVarN(8,1)
      Stress1(2,1)=Stress1(1,2)
	  
      Xq=XVarN(9,1)
      xlan=XVarN(10,1)
      Xc=XVarN(11,1)

! Using new variables         
      Call KDEVIA(stress1,dev1)	  
	  Call KEFFP(dev1,pc1)	  
	  pc2=2.0d0/3.0d0*(pc1**2)
	  dev2=dev1/DSQRT(pc2) 
	  call pressure(Stress1,pr1)

	  xcdot=(xc-sv9)/dtime

	  CALL mdet(2.0d0*XBe+xiden,XdetFesq)
	  XdetFe=dsqrt(XdetFesq)

	  if (isNaN(XdetFe)) then
		PRINT*,'DetFe < 0 (plastic_1_2)',kinc,noel,xiter,xc,sv9
		call xit
	  endif 
	  
	  call xinvariant1(XBe,BeI1)
	  call xinvariant2(XBe,BeI2)
	  
	  call kmlt(XBe,XBe,BeBe)

      ELAM=(1.0d0-xc)*XELAM0+xc*XELAM1
      EG=(1.0d0-xc)*XEG0+xc*XEG1
      XLP=(1.0d0-xc)*Xlpar0+xc*Xlpar1
      XMP=(1.0d0-xc)*Xmpar0+xc*Xmpar1
      XNP=(1.0d0-xc)*Xnpar0+xc*Xnpar1
  
	  sy0=sya+syap*pr1
	  sy1=syb+sybp*pr1
      sy=(1.0d0-xc)*sy0+xc*sy1 
	  
	  ! ! ! xsigybar=xc*sy0+(1.0d0-xc)*sy1	 
	  xsigybar=xc*SyA+(1.0d0-xc)*SyB	! To make pressure-independent yield strength in kinetic eqn	  

! new addition starts
      IF (pr1.Gt.Xpde) THEN
       XHevD=1.0d0
      ELSE
       XHevD=0.0d0
      ENDIF
!	  
      IF (pr1.gt.Xpre) THEN
       XHevR=0.0d0 ! This should be changed if reverse PT is there.
      ELSE
       XHevR=0.0d0
      ENDIF
! new addition ends
	  
	  Xpd=(pr1-Xpde)/(XpdH-Xpde)
	  Xpr=(pr1-Xpre)/(XprH-Xpre)

	  dI3dBe=0.0d0
	  dI3dBe=BeBe-BeI1*xbe+BeI2*xiden
	  call xsymmetric2(dI3dBe)	  
	  	  
	  dEdBe1=0.0d0 
      dEdBe1=ELAM*BeI1*XIDEN+2.0d0*EG*XBe
     1            +(XLP*BeI1**2-2.0d0*XMP*BeI2)*XIDEN
     2            +XNP*dI3dBe+2.0d0*XMP*BeI1*XBe 
	  call xsymmetric2(dEdBe1)

C**********************************************************************
! New functional values
	  
! g1 
	  yafter1=0.0d0
      DO 214 I=1,3
       DO 214 J=1,3
	   
        DO L=1,3
	  yafter1(I,J)=yafter1(I,J)-(xl(I,L)*XBe(L,J)+XBe(l,i)*xl(J,L))
     +  +(2.0d0*XBe(I,L)+xiden(I,L))*(xlan*dev2(L,J)
     +  +tstr*xiden(L,J)*xcdot/3.0d0)
		enddo		
		
        yafter1(I,J)=yafter1(I,J)*DTIME+XBe(I,J)-Be_0(I,J)
     +  -DTIME*(xl(I,J)+xl(J,I))/2.0d0
	  
214   continue
      
! g2   	  
	  yafter2=0.0d0
      DO 215 I=1,3
       DO 215 J=1,3
	   
        DO L=1,3		
	     yafter2(I,J)=yafter2(I,J)-
     +  (2.0d0*XBe(I,L)+xiden(I,L))*dEdBe1(L,J)
		enddo	
		
        yafter2(I,J)=yafter2(I,J)/XdetFe+Stress1(I,J)
215   continue

! yaf(1:11,1)
	  yaf=0.0d0
      DO I=1,3
       yaf(I,1)=yafter1(I,I)
	   yaf(I+3,1)=yafter2(I,I)
      ENDDO
	  
      yaf(7,1)=yafter1(1,2)
      yaf(8,1)=yafter2(1,2)
      
! g3
      yaf(9,1)=Xq-Sv7-((2.0d0/3.0d0)**0.5)*xlan*DTIME         
          
! g4	    
      yaf(10,1)=pc1-sy 
	  
! g5	 
      yaf(11,1)=(Xc-Sv9)-xk*(1.0d0+xB*xq)*(Xq-Sv7)/xsigybar*
     +  (1.0d0+Delta1-Delta2*Xc)*(1.0d0+Delta3*pr1)*Xpd*XHevD*syB		 
	 
C**********************************************************************   

	  yafm1=dsqrt(yaf(1,1)**2+yaf(2,1)**2+yaf(3,1)**2+yaf(4,1)**2+
     1    yaf(5,1)**2+yaf(6,1)**2+yaf(7,1)**2+yaf(8,1)**2+
     2    yaf(9,1)**2+yaf(10,1)**2+yaf(11,1)**2) 	 

      !!!IF (yafm1.LE.0.0001d0 .and. .not. violated) THEN
      IF (yafm1.LE.1.0d-4) THEN
         GoOn=-1	
		 ! ! ! print*, XdetFe, detF
		 
		if ((1.0+Delta3*pr1)*
     1		(1.0-Delta2*xc)*(1.0+xB*xq)*XHevD.lt.0.0d0)then
		   print*,'No reverse PT allowed',xc,sv9,noel,kinc,xq,q0
		   print*,'Parameters are not correct, plastic: dc/dq<1.0'
		   call xit
		endif		 
		 
      ENDIF
	  
      IF (Xiter.GT.25000000) THEN	
	   print*,'hello'	
       PRINT*, 'Jacobian not converge',kinc,NOEL,Xiter,yaf,yafm1,dtime,
     1    detF,detFMax,ELdetFMax,detFeMax,ELdetFeMax,	
     2    detFMin,ELdetFMin,detFeMin,ELdetFeMin,kstep	   
       CALL XIT			   
      ENDIF
             
      Xiter=Xiter+1 
      
      ENDDO		! While loop ends

C**********************************************************************

	  if (sy.lt.0.0d0.or.xc.lt.-1.0D-5.or.xq.lt.-1.0D-6) then
	  print*,'yield strength or concentration or Odqvist par <0:plas'
		print*,sy,xc,xq,noel,kinc,q0
		call xit
	  endif   

	  if (xiter.gt.xitermax) then
	    xitermax=xiter
		ElMaxXiter=NOEL
	  endif	 	

	  if (detF.gt.detFMax) then
	    detFMax=detF
		ELdetFMax=NOEL
	  endif
	  
	  if (XdetFe.gt.detFeMax) then
	    detFeMax=XdetFe
		ELdetFeMax=NOEL
	  endif
	  
	  if (detF.lt.detFMin) then
	    detFMin=detF
		ELdetFMin=NOEL
	  endif
	  
	  if (XdetFe.lt.detFeMin) then
	    detFeMin=XdetFe
		ELdetFeMin=NOEL
	  endif

	  
!	From Be, calculating Ve using SVD.	  

	  call spectral(XBe,xeigvalvecBe,xeigvecmatBe)

	  xeigvalmatVe=0.0d0
	  xeigvalmatVe(1,1)=dsqrt(2.0d0*xeigvalvecBe(1)+1.0d0)
	  xeigvalmatVe(2,2)=dsqrt(2.0d0*xeigvalvecBe(2)+1.0d0)
	  xeigvalmatVe(3,3)=dsqrt(2.0d0*xeigvalvecBe(3)+1.0d0)
	  
	  if (isNaN(xeigvalmatVe(1,1))) then
	   print*,'Error in sign of Ve(1,1)',kinc,noel
	   call xit
	  endif 

	  if (isNaN(xeigvalmatVe(2,2))) then
	   print*,'Error in sign of Ve(2,2)',kinc,noel
	   call xit
	  endif 	  

	  if (isNaN(xeigvalmatVe(3,3))) then
	   print*,'Error in sign of Ve(3,3)',kinc,noel
	   call xit
	  endif 
	  
	  call ktrans(xeigvecmatBe,xeigvecmatBeT)
	  
	  call kmlt(xeigvecmatBe,xeigvalmatVe,xSVDBeVe)

	  call kmlt(xSVDBeVe,xeigvecmatBeT,xVe)
	  
	  CALL xmatInv3D(XVe,XVe_I)

	  call kmlt(XVe_I,DFGRD1,xReUIne)
	  
      Call skinem(xReUIne,XRe,xUIne)
	  
	  CALL mdet(xUIne,xdetUIne)

	  xUtrans = xiden*(xdetUIne**(1.0d0/3.0d0))
	  
	  call xmatInv3D(xUtrans,xUtransInv)
	  
	  call kmlt(xUIne,xUtransInv,xUplastic)


! Extra portion	  
	  CALL mdet(xUplastic,xdetUplastic)
	  
	  if (xdetUplastic .lt. 0.99 .or. xdetUplastic .gt. 1.01) then
		print*,'Error in determinant of U_p'
		call xit
	  endif
	  
	  CALL mdet(xReUIne,xdetReUIne)
	  
	  xF_P = xReUIne/xdetReUIne
	  
	  call ktrans(xF_P,xF_PT)
	  
	  call kmlt(xF_PT, xF_P, xF_pTF_p)

	  xEplastic2 = 0.5d0*(xF_pTF_p-xiden)
	  
!
	  
	  call kmlt(xUplastic, xUplastic, xU2plastic)
	  
	  xEplastic=0.5d0*(xU2plastic-xiden)
	  
	  XUidot=(xUIne-xUIold)/dtime
	  
	  call xmatInv3D(xUIne,XUiInv)

	  call kmlt(XUidot,XUiInv,xUidotUiInv)
	  
	  call ktrans(xUidotUiInv,xUidotUiInvT)
	  call ktrans(XRe,XReT)
	  
	  call kmlt(XRe,0.5d0*(xUidotUiInv+xUidotUiInvT),xReUidotUiInvSym) 
	  call kmlt(xReUidotUiInvSym,XReT,xDIne)
	 
	  XDp=XDIne-xiden*tstr*xcdot/3.0d0
	                      
! Adding for Lagrangian elastic strains 
	  call kmlt(xReT,xbe,xReTBe)
	  call kmlt(xReTBe,XRe,xLagEe)

C********************************************************************** 
	  
	  d2I3dBe2=0.0d0
      DO 701 I=1,3 
	   DO 701 J=1,3
		DO 701 K=1,3   
         DO 701 L=1,3
		   
	  d2I3dBe2(i,j,k,l)=0.5d0*(xbe(i,k)*xiden(j,l)+xbe(i,l)*xiden(j,k)+
     +    xbe(l,j)*xiden(i,k)+xbe(k,j)*xiden(i,l))-
     +    xbe(i,j)*xiden(k,l)-xbe(k,l)*xiden(i,j)+
     +    xiden(i,j)*xiden(k,l)*BeI1
     +    -0.5d0*BeI1*(xi4(i,j,k,l)+xi4(j,i,k,l))
701   CONTINUE	

	  call xsymmetric4_1(d2I3dBe2)
	  call xsymmetric4_2(d2I3dBe2)

	  d2EdBe2=0.0d0
      DO 702 I=1,3 
	   DO 702 J=1,3
		DO 702 K=1,3   
         DO 702 L=1,3
		   
            d2EdBe2(i,j,k,l)=elam*xi4(i,k,j,l)+EG*(xi4(i,j,k,l)+
     + 		xi4(j,i,k,l))+2.0d0*XLP*BeI1*xi4(i,k,j,l)+XMP*BeI1*
     + 		(xi4(i,j,k,l)+xi4(j,i,k,l))+2.0d0*XMP*xiden(k,l)*xbe(i,j)+
     +  	2.0d0*xmp*xiden(i,j)*xbe(k,l)-2.0d0*xmp*BeI1*xi4(i,k,j,l)+
     + 		xnp*d2I3dBe2(i,j,k,l)	 
702   CONTINUE	

	  call xsymmetric4_1(d2EdBe2)
	  call xsymmetric4_2(d2EdBe2)

	  xdJeInvdBe=0.0d0
	  xdJeInvdBe=-(1.0d0/xDetFe**3)*
     + 		(xiden+2.0d0*(BeI1*xiden-xbe)+4.0d0*dI3dBe)

! (df/dBe)  
	  dFdBe=0.0d0 
      DO 703 I=1,3 
	   DO 703 J=1,3
		DO 703 K=1,3   
         DO 703 L=1,3
		 
		  DO F=1,3		   
        dFdBe(i,j,k,l)=dFdBe(I,J,K,L)+
     + 	xdJeInvdBe(k,l)*(2.0d0*xbe(i,f)+xiden(i,f))*
     +  dEdBe1(f,j)+1.0/xdetFe*(xi4(i,f,k,l)+xi4(f,i,k,l))*dEdBe1(f,j)	
     + 	+ 1.0d0/xdetFe*(2.0d0*xbe(i,f)+xiden(i,f))*d2EdBe2(f,j,k,l)	 

	      enddo
		  
703   CONTINUE

! Y(Be)=(df/dBe).(2Be+I)
	  YBe=0.0d0
      DO 704 I=1,3
       DO 704 J=1,3
        DO 704 K=1,3
         DO 704 L=1,3        
          DO 704 F=1,3
           YBe(I,J,K,L)=YBe(I,J,K,L)+dFdBe(I,J,K,F)*
     +                 (2.0d0*XBe(F,L)+XIDEN(F,L))
704   continue   

	  ! ! ! xndir=(1.5**0.5)*dev1/sy
	  xndir=(1.5d0**0.5d0)*dev2 + 
     +	1.0d0/3.0d0*xiden*((1.0d0-xc)*SyAP+xc*SyBP)
	 
	  YBen=0.0d0
	  do 705 i=1,3
	   do 705 j=1,3
	    do 705 k=1,3
	     do 705 l=1,3
		  YBen(i,j)=YBen(i,j)+YBe(i,j,k,l)*xndir(l,k)
705   continue 

	  dphidsig=0.0d0
	  dphidsig=(1.5d0**0.5)*dev2 +
     +	1.0d0/3.0d0*xiden*((1.0d0-xc)*SyAP+xc*SyBP) 	  
	  call xsymmetric2(dphidsig)		
		
	  dphidSYBe=0.0d0	
	  do 706 i=1,3
	   do 706 j=1,3
	    do 706 k=1,3
	     do 706 l=1,3	  
		  dphidSYBe(i,j)=dphidSYBe(i,j)+dphidsig(k,l)*YBe(l,k,i,j)
706   continue
	  
	  Hden=0.0d0
	  do 707 i=1,3
	   do 707 j=1,3
		Hden=hden+dphidsig(i,j)*YBen(j,i)
707   continue


C**********************************************************************

	  dFdc = 0.0D0  
      DO 122 I=1,3 
	   DO 122 J=1,3
	   
		DO K=1,3 	  		  
       dFdc(I,J)=dFdc(I,J)+(1.0d0/XdetFe)*(2.0d0*XBe(I,K)+xiden(I,K))*
     +(xlamdif*BeI1*xiden(K,J)+2.0d0*xGdif*XBe(K,J)+
     +xldif*(BeI1**2.0)*xiden(K,J)+
     +2.0d0*xmdif*(BeI1*XBe(K,J)-BeI2*xiden(K,J))+ xndif*dI3dBe(K,J))	 
		end do
		
122   CONTINUE   
	  call xsymmetric2(dFdc)

! Afac=dc/dq
	  Afac=xk*(1.0d0-xc)*Xpd*XHevD*sy1/xsigybar

	  Zfac1 = -1.0d0*YBen
	  
	  Zfac2 = 0.0d0
	  do 123 i=1,3
	   do 123 j=1,3
	    do 123 k=1,3
	     do 123 l=1,3
		  Zfac2(i,j)=Zfac2(i,j)+YBe(i,j,k,l)*xiden(l,k)	  
123   continue 
	  Zfac2 = -1.0d0*Zfac2*tstr*Afac*dsqrt(2.0d0/3.0D0)
	
	  Zfac3 = dFdc*Afac*dsqrt(2.0d0/3.0D0)	

	  Zfac = Zfac1 + Zfac2 + Zfac3 

	  Hden_new=0.0d0
	  do 124 i=1,3
	   do 124 j=1,3
		Hden_new=hden_new+dphidsig(i,j)*Zfac(j,i)+
     1	(sy0-sy1)*Afac*dsqrt(2.0d0/3.0D0)
124   continue
	  
	  CBe=0.0d0
      DO 708 I=1,3
       DO 708 J=1,3
        DO 708 K=1,3
         DO 708 L=1,3 
          CBe(I,J,K,L)=YBe(I,J,K,L)+stress1(I,J)*XIDEN(K,L) 
     1    - yBen(i,j)*dphidSYBe(k,l)/Hden
708   continue

C**********************************************************************

      DO 709 I=1,3
       DO 709 J=1,3
        DO 709 K=1,3
         DO 709 L=1,3
          CBe(K,L,I,J)=CBe(I,J,K,L)
          CBe(J,I,K,L)=CBe(I,J,K,L)
          CBe(I,J,L,K)=CBe(I,J,K,L)
709   continue

C**********************************************************************

      DDSDDE= 0.0d0	  
      DO 710 K1=1,3
       DO 710 K2=1,3
        DDSDDE(K1,K2)=(CBe(K1,K1,K2,K2)+CBe(K2,K2,K1,K1))/2.0d0
710   continue
                     
      DDSDDE(1,4)=(CBe(1,1,1,2) + CBe(1,2,1,1))/2.0d0
      DDSDDE(4,1)=DDSDDE(1,4)
      DDSDDE(2,4)=(CBe(2,2,1,2) + CBe(1,2,2,2))/2.0d0
      DDSDDE(4,2)=DDSDDE(2,4)
      DDSDDE(3,4)=(CBe(3,3,1,2) + CBe(1,2,3,3))/2.0d0
      DDSDDE(4,3)=DDSDDE(3,4)
      DDSDDE(4,4)=CBe(1,2,1,2)  
		
      DDSDDE(1,5)=(CBe(1,1,1,3) + CBe(1,1,1,3))/2.0d0
      DDSDDE(5,1)=DDSDDE(1,5)
      DDSDDE(2,5)=(CBe(2,2,1,3) + CBe(1,3,2,2))/2.0d0
      DDSDDE(5,2)=DDSDDE(2,5)
      DDSDDE(3,5)=(CBe(3,3,1,3) + CBe(1,3,3,3))/2.0d0
      DDSDDE(5,3)=DDSDDE(3,5)
      DDSDDE(4,5)=(CBe(1,2,1,3) + CBe(1,3,1,2))/2.0d0
      DDSDDE(5,4)=DDSDDE(4,5)
      DDSDDE(5,5)=CBe(1,3,1,3) 

      DDSDDE(1,6)=(CBe(1,1,2,3) + CBe(2,3,1,1))/2.0d0
      DDSDDE(6,1)=DDSDDE(1,6)
      DDSDDE(2,6)=(CBe(2,2,2,3) + CBe(2,3,2,2))/2.0d0
      DDSDDE(6,2)=DDSDDE(2,6)
      DDSDDE(3,6)=(CBe(3,3,2,3) + CBe(2,3,3,3))/2.0d0
      DDSDDE(6,3)=DDSDDE(3,6)
      DDSDDE(4,6)=(CBe(1,2,2,3) + CBe(2,3,1,2))/2.0d0
      DDSDDE(6,4)=DDSDDE(4,6)
      DDSDDE(5,6)=(CBe(1,3,2,3) + CBe(1,3,2,3))/2.0d0
      DDSDDE(6,5)=DDSDDE(5,6)    
      DDSDDE(6,6)=CBe(1,2,1,2)

C**********************************************************************

	  return




167	  DO WHILE (GoOn.GT.0.0d0)

! Gradients for the Newton Raphson	  		  
      Call KDEVIA(stress1,dev1)	  
	  Call KEFFP(dev1,pc1)		
	  pc2=2.0d0/3.0d0*(pc1**2)
	  dev2=dev1/DSQRT(pc2)   
	  call pressure(Stress1,pr1)  
	  
C**********************************************************************
! g1		  	

! dg1dBe  
	  xcdot=(xc-sv9)/dtime	
	  
	  Sy0=SyA+SyAP*pr1
	  sy1=SyB+SyBP*pr1
      Sy=(1.0d0-xc)*sy0+xc*sy1  
	  xsigybar=xc*sy0+(1.0d0-xc)*sy1	      
              
	  dg1dBe=0.0d0		
      DO 902 I=1,3
       DO 902 J=1,3
        DO 902 K=1,3    
         DO 902 L=1,3

	  dg1dBe(i,j,k,l)=(xi4(i,j,k,l)+xi4(j,i,k,l))/2.0d0-
     1  (xl(I,K)*xiden(J,L)+xl(I,L)*xiden(J,k)+ 
     2  xl(j,k)*xiden(i,l)+xl(j,l)*xiden(i,k))*DTIME*0.5d0
     3  +DTIME*(xiden(i,k)*(xlan*dev2(l,j)+tstr*xiden(l,j)*xcdot/3.0d0)
     4  +xiden(i,l)*(xlan*dev2(k,j)+tstr*xiden(k,j)*xcdot/3.0d0))
902	  continue  

	  call xsymmetric4_1(dg1dBe)
	  call xsymmetric4_2(dg1dBe)  

	  dg1dBe1111=dg1dBe(1,1,1,1)
	  dg1dBe1122=dg1dBe(1,1,2,2)
	  dg1dBe1133=dg1dBe(1,1,3,3)
	  dg1dBe1112=dg1dBe(1,1,1,2)+dg1dBe(1,1,2,1)
	  
	  dg1dBe2211=dg1dBe(2,2,1,1)
	  dg1dBe2222=dg1dBe(2,2,2,2)
	  dg1dBe2233=dg1dBe(2,2,3,3)
	  dg1dBe2212=dg1dBe(2,2,1,2)+dg1dBe(2,2,2,1)
	  
	  dg1dBe3311=dg1dBe(3,3,1,1)
	  dg1dBe3322=dg1dBe(3,3,2,2)
	  dg1dBe3333=dg1dBe(3,3,3,3)
	  dg1dBe3312=dg1dBe(3,3,1,2)+dg1dBe(3,3,2,1)

	  dg1dBe1211=dg1dBe(1,2,1,1)
	  dg1dBe1222=dg1dBe(1,2,2,2)
	  dg1dBe1233=dg1dBe(1,2,3,3)
	  dg1dBe1212=dg1dBe(1,2,1,2)+dg1dBe(1,2,2,1)

! dg1dlan
      dg1dlan=0.0d0   
      DO 903 I=1,3
       DO 903 J=1,3
        DO 903 L=1,3
         dg1dlan(I,J)=dg1dlan(I,J)+DTIME*
     +               (2.0d0*XBe(I,L)+xiden(I,L))*dev2(L,J)     
903   CONTINUE	  		  
	  call xsymmetric2(dg1dlan)	  
 
	  dg1dlan11=dg1dlan(1,1)
	  dg1dlan22=dg1dlan(2,2)
	  dg1dlan33=dg1dlan(3,3)
	  dg1dlan12=dg1dlan(1,2)          

! dg1dS
      dg1dS=0.0D0     
      DO 904 I=1,3
       DO 904 J=1,3
        DO 904 K=1,3    
         DO 904 L=1,3   
		 
		  DO F=1,3
	       dg1dS(I,J,K,L)=dg1dS(I,J,K,L)+
     +      DTIME*xlan*(2.0d0*XBe(I,F)+xiden(I,F))/DSQRT(pc2)*
     +      (-1.0d0/3.0d0*xi4(F,K,J,L)+0.5d0*xi4(F,J,K,L)
     +      +0.5d0*xi4(J,F,K,L)-Dev1(F,J)*Dev1(L,K)/pc2) 	 
		  end do
		  
904   CONTINUE	

	  call xsymmetric4_1(dg1dS)
	  call xsymmetric4_2(dg1dS)

	  dg1dS1111=dg1dS(1,1,1,1)
	  dg1dS1122=dg1dS(1,1,2,2)
	  dg1dS1133=dg1dS(1,1,3,3)
	  dg1dS1112=dg1dS(1,1,1,2)+dg1dS(1,1,2,1)
	  
	  dg1dS2211=dg1dS(2,2,1,1)
	  dg1dS2222=dg1dS(2,2,2,2)
	  dg1dS2233=dg1dS(2,2,3,3)
	  dg1dS2212=dg1dS(2,2,1,2)+dg1dS(2,2,2,1)
	  
	  dg1dS3311=dg1dS(3,3,1,1)
	  dg1dS3322=dg1dS(3,3,2,2)
	  dg1dS3333=dg1dS(3,3,3,3)
	  dg1dS3312=dg1dS(3,3,1,2)+dg1dS(3,3,2,1)
	  
	  dg1dS1211=dg1dS(1,2,1,1)
	  dg1dS1222=dg1dS(1,2,2,2)
	  dg1dS1233=dg1dS(1,2,3,3)
	  dg1dS1212=dg1dS(1,2,1,2)+dg1dS(1,2,2,1)

! dg1dq
	  dg1dq=0.0d0
	  dg1dq11=0.0d0;dg1dq22=0.0d0;dg1dq33=0.0d0;dg1dq12=0.0d0
	  call xsymmetric2(dg1dq)
	  
! dg1dc
	  xdg1dc=0.0d0	 
      xdg1dc=(2.0d0*XBe+xiden)*tstr/3.0d0	  	  
	  call xsymmetric2(xdg1dc)
	
	  xdg1dc11=xdg1dc(1,1)
	  xdg1dc22=xdg1dc(2,2)
	  xdg1dc33=xdg1dc(3,3)
	  xdg1dc12=xdg1dc(1,2)

C********************************************************************** 
! g2  

! dg2dlan
	  dg2dlan=0.0d0
	  dg2dlan11=0.0d0;dg2dlan22=0.0d0;dg2dlan33=0.0d0;dg2dlan12=0.0d0
	  call xsymmetric2(dg2dlan)

! dg2dq
	  dg2dq=0.0d0	  
	  dg2dq11=0.0d0;dg2dq22=0.0d0;dg2dq33=0.0d0;dg2dq12=0.0d0
	  call xsymmetric2(dg2dq)

! dg2dS
	  dg2dS=0.0d0
      DO 905 I=1,3 
	   DO 905 J=1,3
		DO 905 K=1,3   
         DO 905 L=1,3
          dg2dS(I,J,K,L)=0.5d0*(xi4(I,J,K,L)+xi4(j,i,K,L))
905   CONTINUE	
	  
	  call xsymmetric4_1(dg2dS)
	  call xsymmetric4_2(dg2dS)
	  
! Components for the Jacobian
	  dg2dS1111=dg2dS(1,1,1,1)
	  dg2dS1122=dg2dS(1,1,2,2)
	  dg2dS1133=dg2dS(1,1,3,3)
	  dg2dS1112=dg2dS(1,1,1,2)+dg2dS(1,1,2,1)
	  
	  dg2dS2211=dg2dS(2,2,1,1)
	  dg2dS2222=dg2dS(2,2,2,2)
	  dg2dS2233=dg2dS(2,2,3,3)
	  dg2dS2212=dg2dS(2,2,1,2)+dg2dS(2,2,2,1)
	  
	  dg2dS3311=dg2dS(3,3,1,1)
	  dg2dS3322=dg2dS(3,3,2,2)
	  dg2dS3333=dg2dS(3,3,3,3)
	  dg2dS3312=dg2dS(3,3,1,2)+dg2dS(3,3,2,1)
	  
	  dg2dS1211=dg2dS(1,2,1,1)
	  dg2dS1222=dg2dS(1,2,2,2)
	  dg2dS1233=dg2dS(1,2,3,3)
	  dg2dS1212=dg2dS(1,2,1,2)+dg2dS(1,2,2,1)
	
! dg2dBe
	  call xinvariant1(XBe,BeI1)
	  call xinvariant2(XBe,BeI2)
	  
	  call kmlt(XBe,XBe,BeBe)

      ELAM=(1.0d0-xc)*XELAM0+xc*XELAM1
      EG=(1.0d0-xc)*XEG0+xc*XEG1
      XLP=(1.0d0-xc)*Xlpar0+xc*Xlpar1
      XMP=(1.0d0-xc)*Xmpar0+xc*Xmpar1
      XNP=(1.0d0-xc)*Xnpar0+xc*Xnpar1

	  dI3dBe=0.0d0 
	  dI3dBe=BeBe-BeI1*xbe+BeI2*xiden
	  call xsymmetric2(dI3dBe)	  
	  	  
	  dEdBe1=0.0d0 
      dEdBe1=ELAM*BeI1*XIDEN+2.0d0*EG*XBe
     1            +(XLP*BeI1**2-2.0d0*XMP*BeI2)*XIDEN
     2            +XNP*dI3dBe+2.0d0*XMP*BeI1*XBe 
	  call xsymmetric2(dEdBe1)

	  d2I3dBe2=0.0d0
      DO 906 I=1,3 
	   DO 906 J=1,3
		DO 906 K=1,3   
         DO 906 L=1,3
		   
	  d2I3dBe2(i,j,k,l)=0.5d0*(xbe(i,k)*xiden(j,l)+xbe(i,l)*xiden(j,k)+
     +    xbe(l,j)*xiden(i,k)+xbe(k,j)*xiden(i,l))-
     +    xbe(i,j)*xiden(k,l)-xbe(k,l)*xiden(i,j)+
     +    xiden(i,j)*xiden(k,l)*BeI1
     +    -0.5d0*BeI1*(xi4(i,j,k,l)+xi4(j,i,k,l))
906   CONTINUE	

	  call xsymmetric4_1(d2I3dBe2)
	  call xsymmetric4_2(d2I3dBe2)

	  d2EdBe2=0.0d0
      DO 907 I=1,3 
	   DO 907 J=1,3
		DO 907 K=1,3   
         DO 907 L=1,3
		   
            d2EdBe2(i,j,k,l)=elam*xi4(i,k,j,l)+EG*(xi4(i,j,k,l)+
     + 		xi4(j,i,k,l))+2.0d0*XLP*BeI1*xi4(i,k,j,l)+XMP*BeI1*
     + 		(xi4(i,j,k,l)+xi4(j,i,k,l))+2.0d0*XMP*xiden(k,l)*xbe(i,j)+
     +  	2.0d0*xmp*xiden(i,j)*xbe(k,l)-2.0d0*xmp*BeI1*xi4(i,k,j,l)+
     + 		xnp*d2I3dBe2(i,j,k,l)	 
907   CONTINUE	

	  call xsymmetric4_1(d2EdBe2)
	  call xsymmetric4_2(d2EdBe2)
	  
	  CALL mdet(2.0d0*XBe+xiden,XdetFesq)
	  XdetFe=dsqrt(XdetFesq)	  

	  if (isNaN(XdetFe)) then
		PRINT*,'DetFe < 0 (plastic_2_1)',kinc,noel,xiter,xc,sv9
		call xit
	  endif 
 
	  xdJeInvdBe=0.0d0
	  xdJeInvdBe=-(1.0d0/xDetFe**3)*
     +    (xiden+2.0d0*(BeI1*xiden-xbe)+4.0d0*dI3dBe)	
	  call xsymmetric2(xdJeInvdBe)

	  dg2dBe=0.0d0
      DO 908 I=1,3 
	   DO 908 J=1,3
		DO 908 K=1,3   
         DO 908 L=1,3
		 
		  DO F=1,3		   
        dg2dBe(i,j,k,l)=dg2dBe(i,j,k,l)-
     + 	xdJeInvdBe(k,l)*(2.0d0*xbe(i,f)+xiden(i,f))*dEdBe1(f,j)	
     +	-1.0d0/xdetFe*(xi4(i,f,k,l)+xi4(f,i,k,l))*dEdBe1(f,j)	
     + 	-1.0d0/xdetFe*(2.0d0*xbe(i,f)+xiden(i,f))*d2EdBe2(f,j,k,l)	 
		  end do
		  
908   CONTINUE

	  call xsymmetric4_1(dg2dBe)
	  call xsymmetric4_2(dg2dBe)
	  
! Components for the Jacobian
	  dg2dBe1111=dg2dBe(1,1,1,1)
	  dg2dBe1122=dg2dBe(1,1,2,2)
	  dg2dBe1133=dg2dBe(1,1,3,3)
	  dg2dBe1112=dg2dBe(1,1,1,2)+dg2dBe(1,1,2,1)
	  
	  dg2dBe2211=dg2dBe(2,2,1,1)
	  dg2dBe2222=dg2dBe(2,2,2,2)
	  dg2dBe2233=dg2dBe(2,2,3,3)
	  dg2dBe2212=dg2dBe(2,2,1,2)+dg2dBe(2,2,2,1)

	  dg2dBe3311=dg2dBe(3,3,1,1)
	  dg2dBe3322=dg2dBe(3,3,2,2)
	  dg2dBe3333=dg2dBe(3,3,3,3)
	  dg2dBe3312=dg2dBe(3,3,1,2)+dg2dBe(3,3,2,1)

	  dg2dBe1211=dg2dBe(1,2,1,1)
	  dg2dBe1222=dg2dBe(1,2,2,2)
	  dg2dBe1233=dg2dBe(1,2,3,3)
	  dg2dBe1212=dg2dBe(1,2,1,2)+dg2dBe(1,2,2,1)	          
		  
! dg2dc
	  xlamdif=XELAM1-XELAM0
	  xGdif=XEG1-XEG0
	  xldif=Xlpar1-Xlpar0
	  xmdif=Xmpar1-Xmpar0
	  xndif=Xnpar1-Xnpar0
	  
	  xdg2dc=0.0D0  
      DO 909 I=1,3 
	   DO 909 J=1,3
	   
		DO K=1,3 	  		  
       xdg2dc(I,J)=xdg2dc(I,J)-(1.0/XdetFe)*(2.0*XBe(I,K)+xiden(I,K))*
     +(xlamdif*BeI1*xiden(K,J)+2.0d0*xGdif*XBe(K,J)+
     +xldif*(BeI1**2)*xiden(K,J)+
     +2.0d0*xmdif*(BeI1*XBe(K,J)-BeI2*xiden(K,J))+ xndif*dI3dBe(K,J))	 
		end do
		
909   CONTINUE   
	  call xsymmetric2(xdg2dc)
	  
	  xdg2dc11=xdg2dc(1,1)
	  xdg2dc22=xdg2dc(2,2)
	  xdg2dc33=xdg2dc(3,3)
	  xdg2dc12=xdg2dc(1,2)

C**********************************************************************
! g3

! dg3dBe
	  dg3dBe=0.0d0
	  dg3dBe11=0.0d0;dg3dBe22=0.0d0;dg3dBe33=0.0d0;dg3dBe12=0.0d0
	  call xsymmetric2(dg3dBe)

! dg3dS	  
      dg3dS=0.0d0
	  dg3dS11=0.0d0;dg3dS22=0.0d0;dg3dS33=0.0d0;dg3dS12=0.0d0
	  call xsymmetric2(dg3dS)

! dg3dlan  
      dg3dlan=-((2.0d0/3.0d0)**0.5)*DTIME

! dg3dq	  
      dg3dq=1.0d0

! dg3dc	  
	  xdg3dc=0.0d0

C**********************************************************************
! g4

! dg4dBe
	  dg4dBe=0.0d0
	  dg4dBe11=0.0d0;dg4dBe22=0.0d0;dg4dBe33=0.0d0;dg4dBe12=0.0d0
	  call xsymmetric2(dg4dBe)

! dg4dlan 	  
	  dg4dlan=0.0d0

! dg4dq	  
	  dg4dq=0.0d0

! dg4dS	  
	  dg4ds=0.0d0
	  dg4dS=(1.5**0.5)*dev2+1.0d0/3.0d0*xiden*((1.0d0-xc)*SyAP+xc*SyBP) 	
	  call xsymmetric2(dg4dS)	
	  
	  dg4dS11=dg4dS(1,1)
	  dg4dS22=dg4dS(2,2)
	  dg4dS33=dg4dS(3,3)
	  dg4dS12=dg4dS(1,2)+dg4dS(2,1)

! dg4dc
	  xdg4dc=sy0-sy1	  

C**********************************************************************
! g5

! dg5dBe
	  xdg5dBe=0.0d0 
	  xdg5dBe11=0.0d0;xdg5dBe22=0.0d0;xdg5dBe33=0.0d0;xdg5dBe12=0.0d0
	  call xsymmetric2(xdg5dBe)
	  
! dg5dlan	  
	  xdg5dlan=0.0D0  

! dg5dS	
      IF (pr1.Gt.Xpde) THEN
       XHevD=1.0d0
      ELSE
       XHevD=0.0d0
      ENDIF
      Xpd=(pr1-Xpde)/(XpdH-Xpde)
	  	  
      IF (pr1.gt.Xpre) THEN
       XHevR=0.0d0 ! This should be changed if reverse PT is there.
      ELSE
       XHevR=0.0d0
      ENDIF
      Xpr=(pr1-Xpre)/(XprH-Xpre)	

	  dy2dybar=-SyBP/(3.0d0*xsigybar)*xiden+
     +  (sy1/xsigybar**2)*(1.0d0/3.0d0)*(xc*SyAP+(1.0d0-xc)*SyBP)*xiden
	 
	  dy1dybar=-SyAP/(3.0d0*xsigybar)*xiden+
     +  (sy0/xsigybar**2)*(1.0d0/3.0d0)*(xc*SyAP+(1.0d0-xc)*SyBP)*xiden
 
	  xdg5ds=0.0d0
	 
      xdg5ds=xiden*xk*(1.0d0+xB*xq)/3.0d0*(Xq-Sv7)/xsigybar*
     +  (1.0d0+Delta1-Delta2*xc)/(XpdH-Xpde)*sy1*XHevD- 
     +  xk*(1.0d0+xB*xq)*(xq-sv7)*(1.0d0+Delta1-Delta2*xc)* 
     +  Xpd*XHevD*dy2dybar	
	  
	  call xsymmetric2(xdg5ds)	 
	  
	  xdg5dS11=xdg5ds(1,1);xdg5dS22=xdg5ds(2,2);xdg5dS33=xdg5ds(3,3)
	  xdg5dS12=xdg5ds(1,2)+xdg5ds(2,1)
 
! dg5dq 

      xdg5dq=-xk*(1.0d0+xB*xq)/xsigybar*
     +	(1.0d0+Delta1-Delta2*Xc)*Xpd*XHevD*sy1 - 	 
     +	(1.0d0+Delta1-Delta2*Xc)*Xpd*XHevD*sy1*(Xq-Sv7)*
     +	xk*(xB)/xsigybar

! dg5dc	  

      xdg5dc=1.0d0 + xk*(1.0d0+xB*xq)*(Xq-Sv7)/xsigybar*
     + ((Delta2*Xpd*XHevD*sy1) + 
     + ((1.0d0+Delta1-Delta2*xc)*Xpd*XHevD*sy1)*
     + (sy0-sy1)/xsigybar)

C********************************************************************** 
 
! Initializing the Jacobian matrix 		   
      XJac=0.0d0             

! g1
	  xjac(1,1)=dg1dBe1111;xjac(1,2)=dg1dBe1122
	  xjac(1,3)=dg1dBe1133;xjac(1,7)=dg1dBe1112
	  
	  xjac(2,1)=dg1dBe2211;xjac(2,2)=dg1dBe2222
	  xjac(2,3)=dg1dBe2233;xjac(2,7)=dg1dBe2212
	  
	  xjac(3,1)=dg1dBe3311;xjac(3,2)=dg1dBe3322
	  xjac(3,3)=dg1dBe3333;xjac(3,7)=dg1dBe3312
	  	  
	  xjac(7,1)=dg1dBe1211;xjac(7,2)=dg1dBe1222
	  xjac(7,3)=dg1dBe1233;xjac(7,7)=dg1dBe1212;
	  
	  xjac(1,4)=dg1dS1111;xjac(1,5)=dg1dS1122;xjac(1,6)=dg1dS1133
	  xjac(1,8)=dg1dS1112
	  
	  xjac(2,4)=dg1dS2211;xjac(2,5)=dg1dS2222;xjac(2,6)=dg1dS2233
	  xjac(2,8)=dg1dS2212
	  
	  xjac(3,4)=dg1dS3311;xjac(3,5)=dg1dS3322;xjac(3,6)=dg1dS3333
	  xjac(3,8)=dg1dS3312
	  
	  xjac(7,4)=dg1dS1211;xjac(7,5)=dg1dS1222;xjac(7,6)=dg1dS1233	  
	  xjac(7,8)=dg1dS1212
	  
	  xjac(1,9)=dg1dq11;xjac(2,9)=dg1dq22;xjac(3,9)=dg1dq33
	  xjac(7,9)=dg1dq12
	  
	  xjac(1,10)=dg1dlan11;xjac(2,10)=dg1dlan22
	  xjac(3,10)=dg1dlan33;xjac(7,10)=dg1dlan12
	  
	  xjac(1,11)=xdg1dc11;xjac(2,11)=xdg1dc22;xjac(3,11)=xdg1dc33
	  xjac(7,11)=xdg1dc12	  

! g2	  	 
	  xjac(4,1)=dg2dBe1111;xjac(4,2)=dg2dBe1122
	  xjac(4,3)=dg2dBe1133;xjac(4,7)=dg2dBe1112

	  xjac(5,1)=dg2dBe2211;xjac(5,2)=dg2dBe2222
	  xjac(5,3)=dg2dBe2233;xjac(5,7)=dg2dBe2212

	  xjac(6,1)=dg2dBe3311;xjac(6,2)=dg2dBe3322
	  xjac(6,3)=dg2dBe3333;xjac(6,7)=dg2dBe3312		  
	  
	  xjac(8,1)=dg2dBe1211;xjac(8,2)=dg2dBe1222
	  xjac(8,3)=dg2dBe1233;xjac(8,7)=dg2dBe1212
 
	  xjac(4,4)=dg2dS1111;xjac(4,5)=dg2dS1122;xjac(4,6)=dg2dS1133
	  xjac(4,8)=dg2dS1112
	  
	  xjac(5,4)=dg2dS2211;xjac(5,5)=dg2dS2222;xjac(5,6)=dg2dS2233
	  xjac(5,8)=dg2dS2212

	  xjac(6,4)=dg2dS3311;xjac(6,5)=dg2dS3322;xjac(6,6)=dg2dS3333
	  xjac(6,8)=dg2dS3312

	  xjac(8,4)=dg2dS1211;xjac(8,5)=dg2dS1222;xjac(8,6)=dg2dS1233
	  xjac(8,8)=dg2dS1212
	  
	  xjac(4,9)=dg2dq11;xjac(5,9)=dg2dq22;xjac(6,9)=dg2dq33
	  xjac(8,9)=dg2dq12
	  
	  xjac(4,10)=dg2dlan11;xjac(5,10)=dg2dlan22
	  xjac(6,10)=dg2dlan33;xjac(8,10)=dg2dlan12
	  
	  xjac(4,11)=xdg2dc11;xjac(5,11)=xdg2dc22;xjac(6,11)=xdg2dc33
	  xjac(8,11)=xdg2dc12

! g3	  
	  xjac(9,1)=dg3dBe11;xjac(9,2)=dg3dBe22;xjac(9,3)=dg3dBe33
	  xjac(9,4)=dg3dS11;xjac(9,5)=dg3dS22;xjac(9,6)=dg3dS33
	  xjac(9,7)=dg3dBe12;xjac(9,8)=dg3dS12;xjac(9,9)=dg3dq
	  xjac(9,10)=dg3dlan;xjac(9,11)=xdg3dc;

! g4
	  xjac(10,1)=dg4dBe11;xjac(10,2)=dg4dBe22;xjac(10,3)=dg4dBe33
	  xjac(10,4)=dg4dS11;xjac(10,5)=dg4dS22;xjac(10,6)=dg4dS33
	  xjac(10,7)=dg4dBe12;xjac(10,8)=dg4dS12;xjac(10,9)=dg4dq
	  xjac(10,10)=dg4dlan;xjac(10,11)=xdg4dc;

! g5	  
	  xjac(11,1)=xdg5dBe11;xjac(11,2)=xdg5dBe22;xjac(11,3)=xdg5dBe33
	  xjac(11,4)=xdg5dS11;xjac(11,5)=xdg5dS22;xjac(11,6)=xdg5dS33
	  xjac(11,7)=xdg5dBe12;xjac(11,8)=xdg5dS12;xjac(11,9)=xdg5dq
	  xjac(11,10)=xdg5dlan;xjac(11,11)=xdg5dc;  

	  do i=1,10
		do j=1,10
		  xjacM(i,j) = xjac(i,j)
	  	enddo
	  enddo
	  
C**********************************************************************
  
	  Be_0=0.0d0
	  be_0(1,1) = sv15; be_0(2,2) = sv18; be_0(3,3) = sv20
	  be_0(1,2) = sv16; be_0(2,1) = be_0(1,2) 

! g1 
	  xg1=0.0d0
      DO 911 I=1,3
       DO 911 J=1,3
	   
        DO L=1,3
			Xg1(I,J)=Xg1(I,J)-(xl(I,L)*XBe(L,J)+XBe(I,L)*xl(J,L))
     +  +(2.0d0*XBe(I,L)+XIDEN(I,L))*(xlan*dev2(L,J)
     +  +tstr*XIDEN(L,J)*xcdot/3.0d0)
		enddo
			
        Xg1(I,J)=Xg1(I,J)*DTIME+XBe(I,J)-Be_0(I,J)-DTIME*(xl(I,J)
     +                 +xl(J,I))/2.0d0
911   CONTINUE

! g2   	  
	  Xg2=0.0d0
      DO 912 I=1,3
       DO 912 J=1,3
        
		DO L=1,3
         Xg2(I,J)=Xg2(I,J)-(2.0d0*xbe(i,l)+XIDEN(I,L))*dEdBe1(L,J)
		enddo
		Xg2(I,J)=Xg2(I,J)/XdetFe + stress1(I,J)
		
912     CONTINUE
      
! Xg(1:11,1)
      DO I=1,3
         Xg(I,1)=Xg1(I,I)
		 Xg(I+3,1)=Xg2(I,I)
      ENDDO
	  
      Xg(7,1)=Xg1(1,2)
      Xg(8,1)=Xg2(1,2)
      
! g3
      Xg(9,1)=Xq-Sv7-((2.0d0/3.0d0)**0.5)*xlan*DTIME     	  
          
! g4   
      Xg(10,1)=pc1-Sy  
	  
! g5

      Xg(11,1)=(Xc-Sv9)-xk*(1.0d0+xB*xq)*(Xq-Sv7)/xsigybar*
     +  (1.0d0+Delta1-Delta2*Xc)*Xpd*XHevD*sy1 	

	  do i=1,10
		xgm(i,1) = Xg(i,1)
	  enddo
    
! Code portion for new variables using Newton Raphson 		  
      !!!CALL matrixinv(XJac,XJacI,11)
	  
	  CALL matrixinv(XJacM,XJacIM,10)

! Initializing the solution vector for strain components with 
! trial elastic strain tensor & stress components with 
! trial stress tensor and also concentration,odqvist and lagrange parameter. 	  

! Old variables
	  DO I=1,3
       XVar(I,1)=XBe(I,I)
	   XVar(I+3,1)=Stress1(I,I)
      ENDDO		  

      XVar(7,1)=XBe(1,2)	  
      XVar(8,1)=Stress1(1,2)	  
      XVar(9,1)=Xq 
      XVar(10,1)=xlan
	  xvar(11,1)=xc	 

	  do i=1,10
		xvarM(i,1) = Xvar(i,1)
	  enddo

	
! New variables

	  xvarnM=0.0d0
      DO I=1,10
        DO J=1,10
          XVarNM(I,1)=XVarNM(I,1)+XJacIM(I,J)*XgM(J,1)
	    enddo
        XVarNM(I,1)=XVarM(I,1)-XVarNM(I,1)
	  enddo
	  
	  xbe=0.0d0
	  stress1=0.0d0	  

      DO I=1,3
       XBe(I,I)=XVarNM(I,1)
	   Stress1(I,I)=XVarNM(I+3,1)
      ENDDO
      XBe(1,2)=XVarNM(7,1)
      XBe(2,1)=XBe(1,2)

      Stress1(1,2)=XVarNM(8,1)
      Stress1(2,1)=Stress1(1,2)
	  
      Xq=XVarNM(9,1)
      xlan=XVarNM(10,1)


! Using new variables         
      Call KDEVIA(stress1,dev1)	  
	  Call KEFFP(dev1,pc1)	  
	  pc2=2.0d0/3.0d0*(pc1**2)
	  dev2=dev1/DSQRT(pc2) 
	  call pressure(Stress1,pr1)

	  xcdot=(xc-sv9)/dtime

	  CALL mdet(2.0d0*XBe+xiden,XdetFesq)
	  XdetFe=dsqrt(XdetFesq)

	  if (isNaN(XdetFe)) then
		PRINT*,'DetFe < 0 (plastic_2_2)',kinc,noel,xiter,xc,sv9
		call xit
	  endif 
	  
	  call xinvariant1(XBe,BeI1)
	  call xinvariant2(XBe,BeI2)
	  
	  call kmlt(XBe,XBe,BeBe)

      ELAM=(1.0d0-xc)*XELAM0+xc*XELAM1
      EG=(1.0d0-xc)*XEG0+xc*XEG1
      XLP=(1.0d0-xc)*Xlpar0+xc*Xlpar1
      XMP=(1.0d0-xc)*Xmpar0+xc*Xmpar1
      XNP=(1.0d0-xc)*Xnpar0+xc*Xnpar1
  
	  sy0=sya+syap*pr1
	  sy1=syb+sybp*pr1
      sy=(1.0d0-xc)*sy0+xc*sy1 
	  xsigybar=xc*sy0+(1.0d0-xc)*sy1
	  
	  Xpd=(pr1-Xpde)/(XpdH-Xpde)
	  Xpr=(pr1-Xpre)/(XprH-Xpre)

	  dI3dBe=0.0d0
	  dI3dBe=BeBe-BeI1*xbe+BeI2*xiden
	  call xsymmetric2(dI3dBe)	  
	  	  
	  dEdBe1=0.0d0 
      dEdBe1=ELAM*BeI1*XIDEN+2.0d0*EG*XBe
     1            +(XLP*BeI1**2-2.0d0*XMP*BeI2)*XIDEN
     2            +XNP*dI3dBe+2.0d0*XMP*BeI1*XBe 
	  call xsymmetric2(dEdBe1)

C**********************************************************************
! New functional values
! g1 
	  yafter1=0.0d0
      DO 914 I=1,3
       DO 914 J=1,3
	   
        DO L=1,3
	  yafter1(I,J)=yafter1(I,J)-(xl(I,L)*XBe(L,J)+XBe(l,i)*xl(J,L))
     +  +(2.0d0*XBe(I,L)+xiden(I,L))*(xlan*dev2(L,J)
     +  +tstr*xiden(L,J)*xcdot/3.0d0)
		enddo		
		
        yafter1(I,J)=yafter1(I,J)*DTIME+XBe(I,J)-Be_0(I,J)
     +  -DTIME*(xl(I,J)+xl(J,I))/2.0d0
	  
914   continue
      
! g2   	  
	  yafter2=0.0d0
      DO 915 I=1,3
       DO 915 J=1,3
	   
        DO L=1,3		
	     yafter2(I,J)=yafter2(I,J)-
     +  (2.0d0*XBe(I,L)+xiden(I,L))*dEdBe1(L,J)
		enddo	
		
        yafter2(I,J)=yafter2(I,J)/XdetFe+Stress1(I,J)
915   continue

! yaf(1:11,1)
	  yaf=0.0d0
      DO I=1,3
       yaf(I,1)=yafter1(I,I)
	   yaf(I+3,1)=yafter2(I,I)
      ENDDO
	  
      yaf(7,1)=yafter1(1,2)
      yaf(8,1)=yafter2(1,2)
      
! g3
      yaf(9,1)=Xq-Sv7-((2.0d0/3.0d0)**0.5)*xlan*DTIME         
          
! g4	    
      yaf(10,1)=pc1-sy 
	  
! g5
      yaf(11,1)=(Xc-Sv9)-xk*(1.0d0+xB*xq)*(Xq-Sv7)/xsigybar*
     +  (1.0d0+Delta1-Delta2*Xc)*Xpd*XHevD*sy1

	  do i=1,10
		yafM(i,1) = yaf(i,1)
	  enddo
   
C**********************************************************************    

	  yafm1=dsqrt(yaf(1,1)**2+yaf(2,1)**2+yaf(3,1)**2+yaf(4,1)**2+
     1    yaf(5,1)**2+yaf(6,1)**2+yaf(7,1)**2+yaf(8,1)**2+
     2    yaf(9,1)**2+yaf(10,1)**2 + 0) 

      !!!IF (yafm1.LE.0.0001d0 .and. .not. violated) THEN
      IF (yafm1.LE.1.0d-4) THEN
        GoOn=-1	 
      ENDIF

      IF (Xiter.GT.25000000) THEN	
	   print*,'hello'	  
       PRINT*, 'Jacobian not converge',kinc,NOEL,Xiter,yaf,yafm1,dtime,
     1    detF,detFMax,ELdetFMax,detFeMax,ELdetFeMax,	
     2    detFMin,ELdetFMin,detFeMin,ELdetFeMin,kstep	
       CALL XIT			   
      ENDIF
             
      Xiter=Xiter+1 
      
      ENDDO		! While loop ends

C**********************************************************************

	  if (sy.lt.0.0d0.or.xc.lt.-1.0D-5.or.xq.lt.-1.0D-6) then
	  print*,'yield strength or concentration or Odqvist par <0:plas'
		print*,sy,xc,xq,noel,kinc
		call xit
	  endif   

	  if (xiter.gt.xitermax) then
	    xitermax=xiter
		ElMaxXiter=NOEL
	  endif	 	

	  if (detF.gt.detFMax) then
	    detFMax=detF
		ELdetFMax=NOEL
	  endif
	  
	  if (XdetFe.gt.detFeMax) then
	    detFeMax=XdetFe
		ELdetFeMax=NOEL
	  endif
	  
	  if (detF.lt.detFMin) then
	    detFMin=detF
		ELdetFMin=NOEL
	  endif
	  
	  if (XdetFe.lt.detFeMin) then
	    detFeMin=XdetFe
		ELdetFeMin=NOEL
	  endif

	  
!	From Be, calculating Ve using SVD.	  

	  call spectral(XBe,xeigvalvecBe,xeigvecmatBe)

	  xeigvalmatVe=0.0d0
	  xeigvalmatVe(1,1)=dsqrt(2.0d0*xeigvalvecBe(1)+1.0d0)
	  xeigvalmatVe(2,2)=dsqrt(2.0d0*xeigvalvecBe(2)+1.0d0)
	  xeigvalmatVe(3,3)=dsqrt(2.0d0*xeigvalvecBe(3)+1.0d0)
	  
	  if (isNaN(xeigvalmatVe(1,1))) then
	   print*,'Error in sign of Ve(1,1)',kinc,noel
	   call xit
	  endif 

	  if (isNaN(xeigvalmatVe(2,2))) then
	   print*,'Error in sign of Ve(2,2)',kinc,noel
	   call xit
	  endif 	  

	  if (isNaN(xeigvalmatVe(3,3))) then
	   print*,'Error in sign of Ve(3,3)',kinc,noel
	   call xit
	  endif 
	  
	  call ktrans(xeigvecmatBe,xeigvecmatBeT)
	  
	  call kmlt(xeigvecmatBe,xeigvalmatVe,xSVDBeVe)

	  call kmlt(xSVDBeVe,xeigvecmatBeT,xVe)
	  
	  CALL xmatInv3D(XVe,XVe_I)

	  call kmlt(XVe_I,DFGRD1,xReUIne)
	  
      Call skinem(xReUIne,XRe,xUIne)
	  
	  CALL mdet(xUIne,xdetUIne)

	  xUtrans = xiden*(xdetUIne**(1.0d0/3.0d0))
	  
	  call xmatInv3D(xUtrans,xUtransInv)
	  
	  call kmlt(xUIne,xUtransInv,xUplastic)
	  
	  call kmlt(xUplastic, xUplastic, xU2plastic)
	  
	  xEplastic=0.5d0*(xU2plastic-xiden)
	  
	  XUidot=(xUIne-xUIold)/dtime
	  
	  call xmatInv3D(xUIne,XUiInv)

	  call kmlt(XUidot,XUiInv,xUidotUiInv)
	  
	  call ktrans(xUidotUiInv,xUidotUiInvT)
	  call ktrans(XRe,XReT)
	  
	  call kmlt(XRe,0.5d0*(xUidotUiInv+xUidotUiInvT),xReUidotUiInvSym) 
	  call kmlt(xReUidotUiInvSym,XReT,xDIne)
	 
	  XDp=XDIne-xiden*tstr*xcdot/3.0d0
	                      
! Adding for Lagrangian elastic strains 
	  call kmlt(xReT,xbe,xReTBe)
	  call kmlt(xReTBe,XRe,xLagEe)

C********************************************************************** 
	  
	  d2I3dBe2=0.0d0
      DO 801 I=1,3 
	   DO 801 J=1,3
		DO 801 K=1,3   
         DO 801 L=1,3
		   
	  d2I3dBe2(i,j,k,l)=0.5d0*(xbe(i,k)*xiden(j,l)+xbe(i,l)*xiden(j,k)+
     +    xbe(l,j)*xiden(i,k)+xbe(k,j)*xiden(i,l))-
     +    xbe(i,j)*xiden(k,l)-xbe(k,l)*xiden(i,j)+
     +    xiden(i,j)*xiden(k,l)*BeI1
     +    -0.5d0*BeI1*(xi4(i,j,k,l)+xi4(j,i,k,l))
801   CONTINUE	

	  call xsymmetric4_1(d2I3dBe2)
	  call xsymmetric4_2(d2I3dBe2)

	  d2EdBe2=0.0d0
      DO 802 I=1,3 
	   DO 802 J=1,3
		DO 802 K=1,3   
         DO 802 L=1,3
		   
            d2EdBe2(i,j,k,l)=elam*xi4(i,k,j,l)+EG*(xi4(i,j,k,l)+
     + 		xi4(j,i,k,l))+2.0d0*XLP*BeI1*xi4(i,k,j,l)+XMP*BeI1*
     + 		(xi4(i,j,k,l)+xi4(j,i,k,l))+2.0d0*XMP*xiden(k,l)*xbe(i,j)+
     +  	2.0d0*xmp*xiden(i,j)*xbe(k,l)-2.0d0*xmp*BeI1*xi4(i,k,j,l)+
     + 		xnp*d2I3dBe2(i,j,k,l)	 
802   CONTINUE	

	  call xsymmetric4_1(d2EdBe2)
	  call xsymmetric4_2(d2EdBe2)

	  xdJeInvdBe=0.0d0
	  xdJeInvdBe=-(1.0d0/xDetFe**3)*
     + 		(xiden+2.0d0*(BeI1*xiden-xbe)+4.0d0*dI3dBe)

! (df/dBe)  
	  dFdBe=0.0d0 
      DO 803 I=1,3 
	   DO 803 J=1,3
		DO 803 K=1,3   
         DO 803 L=1,3
		 
		  DO F=1,3		   
        dFdBe(i,j,k,l)=dFdBe(I,J,K,L)+
     + 	xdJeInvdBe(k,l)*(2.0d0*xbe(i,f)+xiden(i,f))*
     +  dEdBe1(f,j)+1.0/xdetFe*(xi4(i,f,k,l)+xi4(f,i,k,l))*dEdBe1(f,j)	
     + 	+ 1.0d0/xdetFe*(2.0d0*xbe(i,f)+xiden(i,f))*d2EdBe2(f,j,k,l)	 

	      enddo
		  
803   CONTINUE

! Y(Be)=(df/dBe).(2Be+I)
	  YBe=0.0d0
      DO 804 I=1,3
       DO 804 J=1,3
        DO 804 K=1,3
         DO 804 L=1,3        
          DO 804 F=1,3
           YBe(I,J,K,L)=YBe(I,J,K,L)+dFdBe(I,J,K,F)*
     +                 (2.0d0*XBe(F,L)+XIDEN(F,L))
804   continue   

	  ! ! ! xndir=(1.5**0.5)*dev1/sy
	  xndir=(1.5d0**0.5d0)*dev2 + 
     +	1.0d0/3.0d0*xiden*((1.0d0-xc)*SyAP+xc*SyBP)
	 
	  YBen=0.0d0
	  do 805 i=1,3
	   do 805 j=1,3
	    do 805 k=1,3
	     do 805 l=1,3
		  YBen(i,j)=YBen(i,j)+YBe(i,j,k,l)*xndir(l,k)
805   continue 
	  
	  ! ! ! dphidsig=sqrt(3.0/2.0)*dev2
	  ! ! ! dphidsig=3.0/2.0*dev1/sy
	  ! ! ! dphidsig=(1.5**0.5)*dev1/sy

	  dphidsig=0.0d0
	  dphidsig=(1.5d0**0.5)*dev2 +
     +	1.0d0/3.0d0*xiden*((1.0d0-xc)*SyAP+xc*SyBP) 	  
	  call xsymmetric2(dphidsig)		
		
	  dphidSYBe=0.0d0	
	  do 806 i=1,3
	   do 806 j=1,3
	    do 806 k=1,3
	     do 806 l=1,3	  
		  dphidSYBe(i,j)=dphidSYBe(i,j)+dphidsig(k,l)*YBe(l,k,i,j)
806   continue
	  
	  Hden=0.0d0
	  do 807 i=1,3
	   do 807 j=1,3
		Hden=hden+dphidsig(i,j)*YBen(j,i)
807   continue


C**********************************************************************

	  dFdc = 0.0D0  
      DO 822 I=1,3 
	   DO 822 J=1,3
	   
		DO K=1,3 	  		  
       dFdc(I,J)=dFdc(I,J)+(1.0d0/XdetFe)*(2.0d0*XBe(I,K)+xiden(I,K))*
     +(xlamdif*BeI1*xiden(K,J)+2.0d0*xGdif*XBe(K,J)+
     +xldif*(BeI1**2.0)*xiden(K,J)+
     +2.0d0*xmdif*(BeI1*XBe(K,J)-BeI2*xiden(K,J))+ xndif*dI3dBe(K,J))	 
		end do
		
822   CONTINUE   
	  call xsymmetric2(dFdc)

! Afac=dc/dq
	  Afac=xk*(1.0d0-xc)*Xpd*XHevD*sy1/xsigybar

	  Zfac1 = -1.0d0*YBen
	  
	  Zfac2 = 0.0d0
	  do 823 i=1,3
	   do 823 j=1,3
	    do 823 k=1,3
	     do 823 l=1,3
		  Zfac2(i,j)=Zfac2(i,j)+YBe(i,j,k,l)*xiden(l,k)	  
823   continue 
	  Zfac2 = -1.0d0*Zfac2*tstr*Afac*dsqrt(2.0d0/3.0D0)
	
	  Zfac3 = dFdc*Afac*dsqrt(2.0d0/3.0D0)	

	  Zfac = Zfac1 + Zfac2 + Zfac3 

	  Hden_new=0.0d0
	  do 824 i=1,3
	   do 824 j=1,3
		Hden_new=hden_new+dphidsig(i,j)*Zfac(j,i)+
     1	(sy0-sy1)*Afac*dsqrt(2.0d0/3.0D0)
824   continue

C**********************************************************************

	  ! ! ! CBe=0.0d0
      ! ! ! DO 808 I=1,3
       ! ! ! DO 808 J=1,3
        ! ! ! DO 808 K=1,3
         ! ! ! DO 808 L=1,3 
          ! ! ! CBe(I,J,K,L)=YBe(I,J,K,L)+stress1(I,J)*XIDEN(K,L) 
     ! ! ! 1    + Zfac(i,j)*-1.0d0*dphidSYBe(k,l)/Hden_new 
! ! ! 808   continue
	  
	  CBe=0.0d0
      DO 808 I=1,3
       DO 808 J=1,3
        DO 808 K=1,3
         DO 808 L=1,3 
          CBe(I,J,K,L)=YBe(I,J,K,L)+stress1(I,J)*XIDEN(K,L) 
     1    - yBen(i,j)*dphidSYBe(k,l)/Hden
808   continue

C**********************************************************************

      DO 809 I=1,3
       DO 809 J=1,3
        DO 809 K=1,3
         DO 809 L=1,3
          CBe(K,L,I,J)=CBe(I,J,K,L)
          CBe(J,I,K,L)=CBe(I,J,K,L)
          CBe(I,J,L,K)=CBe(I,J,K,L)
809   continue

      DDSDDE= 0.0d0	  
      DO 810 K1=1,3
       DO 810 K2=1,3
        DDSDDE(K1,K2)=(CBe(K1,K1,K2,K2)+CBe(K2,K2,K1,K1))/2.0d0
810   continue
                     
      DDSDDE(1,4)=(CBe(1,1,1,2) + CBe(1,2,1,1))/2.0d0
      DDSDDE(4,1)=DDSDDE(1,4)
      DDSDDE(2,4)=(CBe(2,2,1,2) + CBe(1,2,2,2))/2.0d0
      DDSDDE(4,2)=DDSDDE(2,4)
      DDSDDE(3,4)=(CBe(3,3,1,2) + CBe(1,2,3,3))/2.0d0
      DDSDDE(4,3)=DDSDDE(3,4)
      DDSDDE(4,4)=CBe(1,2,1,2)  
		
      DDSDDE(1,5)=(CBe(1,1,1,3) + CBe(1,1,1,3))/2.0d0
      DDSDDE(5,1)=DDSDDE(1,5)
      DDSDDE(2,5)=(CBe(2,2,1,3) + CBe(1,3,2,2))/2.0d0
      DDSDDE(5,2)=DDSDDE(2,5)
      DDSDDE(3,5)=(CBe(3,3,1,3) + CBe(1,3,3,3))/2.0d0
      DDSDDE(5,3)=DDSDDE(3,5)
      DDSDDE(4,5)=(CBe(1,2,1,3) + CBe(1,3,1,2))/2.0d0
      DDSDDE(5,4)=DDSDDE(4,5)
      DDSDDE(5,5)=CBe(1,3,1,3) 

      DDSDDE(1,6)=(CBe(1,1,2,3) + CBe(2,3,1,1))/2.0d0
      DDSDDE(6,1)=DDSDDE(1,6)
      DDSDDE(2,6)=(CBe(2,2,2,3) + CBe(2,3,2,2))/2.0d0
      DDSDDE(6,2)=DDSDDE(2,6)
      DDSDDE(3,6)=(CBe(3,3,2,3) + CBe(2,3,3,3))/2.0d0
      DDSDDE(6,3)=DDSDDE(3,6)
      DDSDDE(4,6)=(CBe(1,2,2,3) + CBe(2,3,1,2))/2.0d0
      DDSDDE(6,4)=DDSDDE(4,6)
      DDSDDE(5,6)=(CBe(1,3,2,3) + CBe(1,3,2,3))/2.0d0
      DDSDDE(6,5)=DDSDDE(5,6)    
      DDSDDE(6,6)=CBe(1,2,1,2)

C***************************************************************************	  
	  
	  return

	  end subroutine nr4	  
	  









***********************************************
**           UTILITY    SUBROUTINES          **
***********************************************

**************************************
** 		Matrix multiplication		**
**************************************

      SUBROUTINE KMLT(DM1,DM2,DM)
     
	  IMPLICIT NONE
	      
	  Integer i,j,k
      REAL*8 DM1(3,3),DM2(3,3),DM(3,3)
	
	  DM=0.0
      DO I=1,3
       DO J=1,3
        DO K=1,3 
         DM(I,J)=DM(I,J)+DM1(I,K)*DM2(K,J)
        enddo
	   enddo
	  enddo 
      
      RETURN
      END


*****************************************************
**   DEVIATORIC STRESS CALCULATION    **
*****************************************************

      SUBROUTINE KDEVIA(SIG,DEV)
     
	  IMPLICIT NONE

      Real*8 SIG(3,3),XIDEN(3,3),DEV(3,3)
      REAL*8 trace

	  call xidentity(xiden)	
	  
	  call xinvariant1(sig,trace)
    
	  dev=sig-(1.0/3.0)*trace*xiden
	  
      RETURN
      END


***************************************
**          EFFECTIVE STRESS          *
**   (CONTRACTED MATRIX CALCULATION)  *
***************************************

      SUBROUTINE KEFFP(xmatrix,vonMises)
    
	  IMPLICIT NONE

      INTEGER I,J
      Real*8 xmatrix(3,3)
      REAL*8 ContracSum, vonMises

      ContracSum=0.0D0
      DO I=1,3
       DO J=1,3
       ContracSum=ContracSum+xmatrix(I,J)*xmatrix(J,I)
	   enddo
	  enddo 
      vonMises=DSQRT((3.0/2.0)*ContracSum)
      
      RETURN
      END


***************************************
**         right polar decomposition         *
***************************************
       subroutine skinem(F,R,U)
      !
      ! This subroutine performs the right polar decomposition
      !  F = RU of the deformation gradient F into a rotation
      !  R and the right stretch tensor U.  The logarithmic 
      !  strain E = ln(U) is also returned.
      !
      !	F(3,3):       the deformation gradient; input
      !	detF:         the determinant of F; detF > 0
      !	R(3,3):       the rotation matrix; output
      !	U(3,3):       the right stretch tensor; output
      !	Uinv(3,3):    the inverse of U
      !	C(3,3):       the right Cauchy-Green tensor
      !	omega(3):     the squares of the principal stretches
      ! Ueigval(3):   the principal stretches
      !	eigvec(3,3):  matrix of eigenvectors of U
      !	E(3,3):       the logarithmic strain tensor; output
      ! istat:        success flag, istat=0 for a failed attempt; output
      !
	  IMPLICIT NONE      
      !
      integer istat,line
      !
      real*8 F(3,3),C(3,3),omega(3),Ueigval(3),eigvec(3,3),
     +  U(3,3),E(3,3),Uinv(3,3),R(3,3),detF
	 
	  INTEGER I,J
     

      !!	Store the identity matrix in R, U, and Uinv     
      DO I=1,3
          DO J=1,3
              R(I,J)=0.0D0
              R(I,I)=0.0D0
              U(I,J)=0.0D0
              U(I,I)=0.0D0
              Uinv(I,J)=0.0D0
              Uinv(I,I)=0.0D0
              E(I,J)=0.0D0
          ENDDO
      ENDDO
      

      ! Check if the determinant of F is greater than zero.
      !  If not, then print a diagnostic and cut back the 
      !  time increment.
      !
      call mdet(F,detF)
      if (detF.le.0.d0) then
        PRINT*, 'problem in kinematics, the',
     +       ' determinant of F is not greater than 0'
        istat = 0
        return
      end if
      

      ! Calculate the right Cauchy-Green tensor C
      !
      C = matmul(transpose(F),F)
      
 
      ! Calculate the eigenvalues and eigenvectors of C
      !
      call spectral(C,omega,eigvec)
      

      ! Calculate the principal values of U and E
      !
      Ueigval(1) = dsqrt(omega(1))
      Ueigval(2) = dsqrt(omega(2))
      Ueigval(3) = dsqrt(omega(3))
      !
      U(1,1) = Ueigval(1)
      U(2,2) = Ueigval(2)
      U(3,3) = Ueigval(3)
      !
      E(1,1) = dlog(Ueigval(1))
      E(2,2) = dlog(Ueigval(2))
      E(3,3) = dlog(Ueigval(3))
      

      ! Calculate the complete tensors U and E
      !
      U = matmul(matmul(eigvec,U),transpose(eigvec))
      E = matmul(matmul(eigvec,E),transpose(eigvec))
      

      ! Calculate Uinv
      !
      call xmatInv3D(U,Uinv)
      

      ! calculate R
      !
      R = matmul(F,Uinv)
      

      return
      end 


*****************************************************************************
**     The following subroutines calculate the spectral
**      decomposition of a symmetric 3 by 3 matrix
*****************************************************************************

      subroutine spectral(A,D,V)
      !
      ! This subroutine calculates the eigenvalues and eigenvectors of
      !  a symmetric 3 by 3 matrix A.
      !
      ! The output consists of a vector D containing the three
      !  eigenvalues in ascending order, and a matrix V whose
      !  columns contain the corresponding eigenvectors.
      !
      implicit none
      !
      integer np,nrot,i,j,istat
      parameter(np=3)
      !
      real*8 D(3),V(3,3),A(3,3),E(3,3)


      DO I=1,3
          DO J=1,3
              E(I,J)=A(I,J)
          ENDDO
      ENDDO
      
      call XJacobi(E,3,np,D,V,nrot,istat)
      call eigsrt(D,V,3,np)
	

      return
      end    
	
C****************************************************************************

      subroutine XJacobi(A,n,np,D,V,nrot)
      !
      ! Computes all eigenvalues and eigenvectors of a real symmetric
      !  matrix A, which is of size n by n, stored in a physical
      !  np by np array.  On output, elements of A above the diagonal
      !  are destroyed, but the diagonal and sub-diagonal are unchanged
      !  and give full information about the original symmetric matrix.
      !  Vector D returns the eigenvalues of A in its first n elements.
      !  V is a matrix with the same logical and physical dimensions as
      !  A whose columns contain, upon output, the normalized
      !  eigenvectors of A.  nrot returns the number of XJacobi rotation
      !  which were required.
      !
      ! This subroutine is taken from 'Numerical Recipes.'
      !
      implicit none
      !
      integer ip,iq,n,nmax,np,nrot,i,j,istat
      parameter (nmax=100)
      !
      real*8 A(np,np),D(np),V(np,np),B(nmax),Z(nmax),
     +  sm,tresh,G,T,H,theta,S,C,tau
     
      
      ! Initialize V to the identity matrix
      DO I=1,3
          DO J=1,3
              V(I,J)=0.0D0
              V(I,I)=1.0D0
          ENDDO
      ENDDO
      
      
      ! Initialize B and D to the diagonal of A, and Z to zero.
      !  The vector Z will accumulate terms of the form T*A_PQ as
      !  in equation (11.1.14)
      !
      do ip = 1,n
	    B(ip) = A(ip,ip)
	    D(ip) = B(ip)
	    Z(ip) = 0.d0
      end do
      
      
      ! Begin iteration
      !
      nrot = 0
      do i=1,200
          !
          ! Sum off-diagonal elements
          !
          sm = 0.d0
          do ip=1,n-1
            do iq=ip+1,n
	      sm = sm + dabs(A(ip,iq))
            end do
          end do
          !
          ! If sm = 0., then return.  This is the normal return,
          !  which relies on quadratic convergence to machine
          !  underflow.
          !
          if (sm.eq.0.d0) return
          !
          ! In the first three sweeps carry out the PQ rotation only if
          !  |A_PQ| > tresh, where tresh is some threshold value,
          !  see equation (11.1.25).  Thereafter tresh = 0.
          !
          if (i.lt.4) then
            tresh = 0.2d0*sm/n**2
          else
            tresh = 0.d0
          end if
          !
          do ip=1,n-1
            do iq=ip+1,n
              G = 100.d0*dabs(A(ip,iq))
              !
              ! After four sweeps, skip the rotation if the 
              !  off-diagonal element is small.
              !
	      if ((i.gt.4).and.(dabs(D(ip))+G.eq.dabs(D(ip)))
     +            .and.(dabs(D(iq))+G.eq.dabs(D(iq)))) then
                A(ip,iq) = 0.d0
              else if (dabs(A(ip,iq)).gt.tresh) then
                H = D(iq) - D(ip)
                if (dabs(H)+G.eq.dabs(H)) then
                  !
                  ! T = 1./(2.*theta), equation (11.1.10)
                  !
	          T =A(ip,iq)/H
	        else
	          theta = 0.5d0*H/A(ip,iq)
	          T =1.d0/(dabs(theta)+dsqrt(1.d0+theta**2.d0))
	          if (theta.lt.0.d0) T = -T
	        end if
	        C = 1.d0/dsqrt(1.d0 + T**2.d0)
	        S = T*C
	        tau = S/(1.d0 + C)
	        H = T*A(ip,iq)
	        Z(ip) = Z(ip) - H
	        Z(iq) = Z(iq) + H
	        D(ip) = D(ip) - H
	        D(iq) = D(iq) + H
	        A(ip,iq) = 0.d0
                !
                ! Case of rotations 1 <= J < P
		!		
	        do j=1,ip-1
	          G = A(j,ip)
	          H = A(j,iq)
	          A(j,ip) = G - S*(H + G*tau)
	          A(j,iq) = H + S*(G - H*tau)
	        end do
                !
                ! Case of rotations P < J < Q
                !
	        do j=ip+1,iq-1
	          G = A(ip,j)
	          H = A(j,iq)
	          A(ip,j) = G - S*(H + G*tau)
	          A(j,iq) = H + S*(G - H*tau)
	        end do
                !
                ! Case of rotations Q < J <= N
                !
	        do j=iq+1,n
                  G = A(ip,j)
	          H = A(iq,j)
	          A(ip,j) = G - S*(H + G*tau)
	          A(iq,j) = H + S*(G - H*tau)
	        end do
	        do j = 1,n
	          G = V(j,ip)
	          H = V(j,iq)
	          V(j,ip) = G - S*(H + G*tau)
	          V(j,iq) = H + S*(G - H*tau)
	        end do
	        nrot = nrot + 1
              end if
	    end do
	  end do
          !
          ! Update D with the sum of T*A_PQ, and reinitialize Z
          !
	  do ip=1,n
	    B(ip) = B(ip) + Z(ip)
	    D(ip) = B(ip)
	    Z(ip) = 0.d0
	  end do
	end do


      ! If the algorithm has reached this stage, then there
      !  are too many sweeps.  Print a diagnostic and cut the 
      !  time increment.
      !
      write (*,'(/1X,A/)') '50 iterations in XJacobi should never happen'
      istat = 0
      

      return
      end 
	
C****************************************************************************

      subroutine eigsrt(D,V,n,np)
      !
      ! Given the eigenvalues D and eigenvectors V as output from
      !  XJacobi, this subroutine sorts the eigenvales into ascending
      !  order and rearranges the colmns of V accordingly.
      !
      ! The subroutine is taken from 'Numerical Recipes.'
      !
      implicit none
      !
      integer n,np,i,j,k
      !
      real*8 D(np),V(np,np),P
      

      do i=1,n-1
	k = i
	P = D(i)
	do j=i+1,n
	  if (D(j).ge.P) then
	    k = j
	    P = D(j)
	  end if
	end do
	if (k.ne.i) then
	  D(k) = D(i)
	  D(i) = P
	  do j=1,n
	    P = V(j,i)
	    V(j,i) = V(j,k)
	    V(j,k) = P
	  end do
  	end if
      end do
      

      return
      end subroutine eigsrt


C********************************
***     Matrix inverse		***
C********************************

      subroutine xmatInv3D(A,A_inv)
      !
      ! Returns A_inv, the inverse and det_A, the determinant
      ! Note that the det is of the original matrix, not the
      ! inverse
      !
      implicit none
      !
      integer istat
      !
      real*8 A(3,3),A_inv(3,3),det_A,det_A_inv


      istat = 1
      
      det_A = A(1,1)*(A(2,2)*A(3,3) - A(3,2)*A(2,3)) -
     +        A(2,1)*(A(1,2)*A(3,3) - A(3,2)*A(1,3)) +
     +        A(3,1)*(A(1,2)*A(2,3) - A(2,2)*A(1,3))
        
      det_A_inv = 1.d0/det_A
        
      A_inv(1,1) = det_A_inv*(A(2,2)*A(3,3)-A(3,2)*A(2,3))
      A_inv(1,2) = det_A_inv*(A(3,2)*A(1,3)-A(1,2)*A(3,3))
      A_inv(1,3) = det_A_inv*(A(1,2)*A(2,3)-A(2,2)*A(1,3))
      A_inv(2,1) = det_A_inv*(A(3,1)*A(2,3)-A(2,1)*A(3,3))
      A_inv(2,2) = det_A_inv*(A(1,1)*A(3,3)-A(3,1)*A(1,3))
      A_inv(2,3) = det_A_inv*(A(2,1)*A(1,3)-A(1,1)*A(2,3))
      A_inv(3,1) = det_A_inv*(A(2,1)*A(3,2)-A(3,1)*A(2,2))
      A_inv(3,2) = det_A_inv*(A(3,1)*A(1,2)-A(1,1)*A(3,2))
      A_inv(3,3) = det_A_inv*(A(1,1)*A(2,2)-A(2,1)*A(1,2))
      

      return
      end subroutine xmatInv3D


C********************************
***     Matrix determinant		***
C********************************

      subroutine mdet(A,det)
      !
      ! This subroutine calculates the determinant
      ! of a 3 by 3 matrix [A]
      !
      implicit none
      !
      real*8  A(3,3),det


      det = A(1,1)*A(2,2)*A(3,3) 
     +	  + A(1,2)*A(2,3)*A(3,1)
     +	  + A(1,3)*A(2,1)*A(3,2)
     +	  - A(3,1)*A(2,2)*A(1,3)
     +	  - A(3,2)*A(2,3)*A(1,1)
     +	  - A(3,3)*A(2,1)*A(1,2)

      return
      end subroutine mdet


C********************************
***     Matrix inverse		***
C******************************** 

      SUBROUTINE matrixinv(A,AINV,N)
C         N: dimension of square matrix
C         A: square matrix of dimension NxN to be inverted
C      AINV: the inverted matrix
      IMPLICIT REAL*8(A-H,O-Z)
	  
	  INTEGER I,J,K,N		  
	  
      DIMENSION A(N,N),AINV(N,N),B(N,2*N)
      
C MAKE AUGMENTED MATRIX
      DO I=1,N
          DO J=1,N
              B(I,J)=0.0D0
              B(I,J+N)=0.0D0
              
              B(I,J)=A(I,J)
              IF (I.EQ.J) THEN
                  B(I,J+N)=1.0D0
              ENDIF
          ENDDO
      ENDDO
      
      DO I=1,N
C CHOOSE THE LEFTMOST NON-ZERO ELEMENT AS PIVOT   
          DO J=1,N
              IF (DABS(B(I,J)).GT.0) THEN
                  PIVOT=B(I,J)
                  EXIT
              ENDIF
          ENDDO
C STEP 1: Change the chosen pivot into "1" by dividing
C the pivot's row by the pivot number 
          DO J=1,2*N
              B(I,J)=B(I,J)/PIVOT
          ENDDO
          PIVOT=B(I,I)    !UPDATE PIVOT VALUE
C STEP 2: Change the remainder of the pivot's COLUMN into 0's
C by adding to each row a suitable multiple of the PIVOT ROW  
          DO K=1,N ! ROW
              IF (K.NE.I) THEN
                  XNUM=B(K,I)/PIVOT  !SAME COLUMN WITH THE CURRENT PIVOT
                  DO J=1,2*N ! COL
                      B(K,J)=B(K,J)-XNUM*B(I,J)
                  ENDDO
              ENDIF
          ENDDO
      ENDDO

          
C PREPARE THE FINAL INVERTED MATRIX
      DO I=1,N
          DO J=1,N
              AINV(I,J)=B(I,J+N)
          ENDDO
      ENDDO
          
          
      RETURN
      END


C********************************
***     Matrix determinant		***
C********************************
             
      SUBROUTINE getdet(A,MM,N,B)
      
	  IMPLICIT NONE
	  
      INTEGER I,J,K,L,N,NN,MM
	  
      REAL*8 ELEM(N,N),A(MM,MM)
      REAL*8 TEMP,B,M

      LOGICAL DETEXISTS
      
      DO I=1,N
          DO J=1,N
              ELEM(I,J)=A(I,J)
          ENDDO
      ENDDO
      
      DETEXISTS= .TRUE.
      L=1
      DO K=1,N-1
          IF (DABS(ELEM(K,K)).LE.1.0D-20) THEN
              DETEXISTS= .FALSE.
              DO I=K+1,N
                  IF (ELEM(I,K).NE.0.0) THEN
                      DO J=1,N
                          TEMP=ELEM(I,J)
                          ELEM(I,J)=ELEM(K,J)
                          ELEM(K,J)=TEMP
                      ENDDO
                      DETEXISTS= .TRUE.
                      L=-L
                      EXIT
                  ENDIF
              ENDDO
              IF (DETEXISTS .EQV. .FALSE.) THEN
                  B=0
                  RETURN
              ENDIF
          ENDIF
          DO J=K+1,N
              M=ELEM(J,K)/ELEM(K,K)
              DO I=K+1,N
                  ELEM(J,I)=ELEM(J,I)-M*ELEM(K,I)
              ENDDO
          ENDDO
      ENDDO
      B=L
      DO I=1,N
          B=B*ELEM(I,I)
      ENDDO
      
      RETURN
      END


C********************************
***     Left stretch tensor		***
C********************************

      SUBROUTINE gettheVe(Re1,Ue1,Ve1)
      
	  IMPLICIT NONE

	  INTEGER I,J,M,N
	  
	  PARAMETER (M=3,N=3)
	  
	  Real*8 Re1(M,N),Ve1(M,N),X1(M,N),X2(M,N),X3(M,N),Ue1(M,N)
     
      CALL KMLT(Re1,Ue1,X1)
C  x1=Re*Ue1
      
      CALL KTRANS(Re1,X2)
C  x2 is the transposed Re    
      
      CALL KMLT(X1,X2,X3)
      
      DO I=1,M
       DO J=1,N
        Ve1(I,J)=X3(I,J)
	   enddo
	  enddo 

      RETURN

      END


C********************************
***     Matrix transposition		***
C********************************

      SUBROUTINE KTRANS(ORIGN,TRAN)

	  IMPLICIT NONE
	  
	  INTEGER I,J,M,N
	  Real*8 ORIGN,TRAN
	  
      PARAMETER (M=3,N=3)
      DIMENSION ORIGN(M,N),TRAN(M,N)

      DO I=1,M
       DO J=1,N
       TRAN(J,I)=ORIGN(I,J)
       enddo
	  enddo 
	  
      RETURN
      END
	  
	  
C********************************
***     Matrix symmterizing		***
C********************************

      SUBROUTINE xsymmetric2(xmatrix)  
	  
	  IMPLICIT NONE
	  
	  INTEGER M,N
	  Real*8 xmatrix
	  
      PARAMETER (M=3,N=3)
	  
      DIMENSION xmatrix(M,N)

	  xmatrix(1,2)=(xmatrix(1,2)+xmatrix(2,1))/2.0
	  xmatrix(2,1)=xmatrix(1,2)
	  xmatrix(1,3)=(xmatrix(1,3)+xmatrix(3,1))/2.0
	  xmatrix(3,1)=xmatrix(1,3)
	  xmatrix(2,3)=(xmatrix(2,3)+xmatrix(3,2))/2.0
	  xmatrix(3,2)=xmatrix(2,3)

      RETURN
      END	  


C********************************
***     Axisymmetric constraints	***
C********************************

      SUBROUTINE xaxicons(xmatrix) 
	  
	  IMPLICIT NONE
	  
	  INTEGER M,N
	  Real*8 xmatrix
	  
      PARAMETER (M=3,N=3)
	  
      DIMENSION xmatrix(M,N)

	  xmatrix(1,3)=0.0
	  xmatrix(3,1)=0.0
      xmatrix(2,3)=0.0
	  xmatrix(3,2)=0.0

      RETURN
      END


C********************************
***     Identity matrix 		***
C********************************

      SUBROUTINE xidentity(xmatrix) 

	  IMPLICIT NONE
	  
	  INTEGER M,N
	  Real*8 xmatrix
	  
      PARAMETER (M=3,N=3)
	  
      DIMENSION xmatrix(M,N)

	  xmatrix=0.0
	  xmatrix(1,1)=1.0;xmatrix(2,2)=1.0;xmatrix(3,3)=1.0
	  
      RETURN
      END	


C********************************
***     Fourth order identity matrix	***
C********************************

      SUBROUTINE xidentity4(xmatrix) 
	  
	  IMPLICIT NONE
	  
	  INTEGER I,J,K,L
	  Real*8 xmatrix	
	 	  
      DIMENSION xmatrix(3,3,3,3)
	   
	  xmatrix=0.0
      DO I=1,3
       DO J=1,3
        Do K=1,3    
         DO L=1,3
		 
          IF (I.EQ.K .AND. J.EQ.L) then
		  xmatrix(I,J,K,L)=1.0
		  endif
		  
		 end do
		end do		
	   end do
	  end do 
	  
      RETURN
      END	


C********************************
***     Matrix first invariant		***
C********************************

      SUBROUTINE xinvariant1(xmatrix,xinvar1) 
	  
	  IMPLICIT NONE
	  
	  INTEGER M,N
	  REAL*8  xinvar1,xmatrix
	  
      PARAMETER (M=3,N=3)
	  
      DIMENSION xmatrix(M,N)
	  
	  xinvar1=xmatrix(1,1)+xmatrix(2,2)+xmatrix(3,3)
	  
      RETURN
      END	


C********************************
***     Matrix second invariant		***
C********************************

      SUBROUTINE xinvariant2(xmatrix,xinvar2) 
	  
	  IMPLICIT NONE
	  
	  INTEGER M,N
	  REAL*8  xinvar2,xmatrix
	  
      PARAMETER (M=3,N=3)
	  
      DIMENSION xmatrix(M,N)

	  xinvar2=xmatrix(1,1)*xmatrix(2,2)+xmatrix(2,2)*xmatrix(3,3)+
     +       xmatrix(3,3)*xmatrix(1,1)-xmatrix(1,2)*xmatrix(2,1)-
     +       xmatrix(2,3)*xmatrix(3,2)-xmatrix(1,3)*xmatrix(3,1)	  
      RETURN
      END	
	  
	  
C********************************
***     Pressure evaluation		***
C********************************	

      SUBROUTINE pressure(xmatrix,xpressure) 
	  
	  IMPLICIT NONE
	  
	  INTEGER M,N
	  REAL*8  xpressure,xmatrix
	  
      PARAMETER (M=3,N=3)
	  
      DIMENSION xmatrix(M,N)

	  xpressure=-1.0/3.0*(xmatrix(1,1)+xmatrix(2,2)+xmatrix(3,3))
	  
      RETURN
      END  
	  
	  
C********************************
***     Imposing minor symmetries		***
C********************************

      SUBROUTINE xsymmetric4_1(xmatrix)  
	  
	  IMPLICIT NONE
	      
	  Integer i,j,k,l
	  REAL*8  xmatrix
	  
	  DIMENSION xmatrix(3,3,3,3)

      DO I=1,3
       DO J=I,3
        DO K=1,3    
         DO L=1,3
          xmatrix(I,J,K,L)=xmatrix(J,I,K,L)
		 end do 
		end do 
	   end do 
	  end do 

      RETURN
      END	  
	 
	 
C********************************
***     Imposing minor symmetries		***
C********************************

      SUBROUTINE xsymmetric4_2(xmatrix) 
	  
	  IMPLICIT NONE
	      
	  Integer i,j,k,l
	  REAL*8  xmatrix
	  
	  DIMENSION xmatrix(3,3,3,3)

      DO I=1,3
       DO J=1,3
        DO K=1,3    
         DO L=K,3
          xmatrix(I,J,K,L)=xmatrix(I,J,L,K)
		 end do
		end do
	   end do 
	  end do 

      RETURN
      END	 
	  
	  
C********************************
***     FRIC code		***
C******************************** 

      SUBROUTINE FRIC(LM,TAU,DDTDDG,DDTDDP,DSLIP,SED,SFD,
     1    DDTDDT,PNEWDT,STATEV,DGAM,TAULM,PRESS,DPRESS,DDPDDH,SLIP,
     2    KSTEP,KINC,TIME,DTIME,NOEL,CINAME,SLNAME,MSNAME,NPT,NODE,
     3    NPATCH,COORDS,RCOORD,DROT,TEMP,PREDEF,NFDIR,MCRD,NPRED,
     4    NSTATV,CHRLNGTH,PROPS,NPROPS)
            
      USE Zrsample
      
      INCLUDE 'ABA_PARAM.INC'
                   
      CHARACTER*80 CINAME,SLNAME,MSNAME
      DIMENSION TAU(NFDIR),DDTDDG(NFDIR,NFDIR),DDTDDP(NFDIR),
     1    DSLIP(NFDIR),DDTDDT(NFDIR,2),STATEV(*),DGAM(NFDIR),
     2    TAULM(NFDIR),SLIP(NFDIR),TIME(2),COORDS(MCRD),
     3    RCOORD(MCRD),DROT(2,2),TEMP(2),PREDEF(2,*),PROPS(NPROPS)
                    
      Real*8 GAMMA(2),thk1(26,1)
	 
      INTEGER I
	  
	  PARAMETER(ZERO=0.0D0)
	  
      Real*8 con,ysm,xmum,dxmumdp,dyshdp,
     1     thk2,xm,yield,taueqv,gameqv,gcrit,dgsleq,
     2     taucrit1,taucrit2,stiff1,stiff2
	 
	  Real*8 xmu0,xmu1,xmup0,xmup1,ysp0,ysp1
	  Real*8 flag_10
            
      ! ! ! INTEGER, DIMENSION(:,:), allocatable :: connect   
	  	  
	  INTEGER, DIMENSION(nodeno,3) :: connect	  
	  
      REAL*8, DIMENSION(mp,27) :: xm1
      REAL*8, DIMENSION(mp,1) :: xm2
	  
! xm1 is 'm' values imported from Matlab for all the load cases
! xm2 is 'm' values based on current thickness at the center from UMAT
	  
      LOGICAL :: file_exists          

      xmu0=props(1)	! 'mu' for alpha phase @ p=0
	  xmu1=props(2)	! 'mu' for omega phase @ p=0
	  
	  xmup0=props(3)	! contact pressure dependence of 'mu' for alpha phase
      xmup1=props(4)	! contact pressure dependence of 'mu' for omega phase
	  
	  ysp0=props(5)	! Tensorial pressure dependence of yield strength in tension of alpha phase used in UMAT
      ysp1=props(6)	! Tensorial pressure dependence of yield strength in tension of omega phase used in UMAT
	    
! "connect.txt" contains surface nodes and corresponding element connectivities.
! "yieldme" has 2 columns. Ist column contains concentration of omega phase and 2nd column contains
! pressurized yield strength (in tension) of the corresponding finite element from the UMAT.

! Reading the connect.txt and storing in 'connect' variable	
! The nodes in connect.txt are in increasing order
  
      ! ! ! ALLOCATE(connect(NodeNo,3))
      OPEN (30,file='C:\Krishan5\Non-Kinetics\connect.txt')
        DO I=1,NodeNo 
          READ(30,*) connect(I,:)
        ENDDO
      CLOSE(30)	  

! Assigning yield strength to the nodal point  
! From the UMAT, the yield strength and concentration are at integration points of finite elements 
! So, using the following method, we are assigning yield strength and concentration to surface nodes
	 
      DO I=1, NodeNo
	  
        IF (NODE .EQ. connect(I,1)) THEN
! Concentration of omega phase 		  
         con = (yieldme(connect(I,2),1)+YieldMe(connect(I,3),1))/2.0d0
! Pressurized (tensorial) yield strength of mixture in tension			 
         ysm = (yieldme(connect(I,2),2)+yieldme(connect(I,3),2))/2.0d0	       		 
         EXIT
			 
        ENDIF
      ENDDO

      IF (ysm.LE.ZERO.or.con.lt.-1.0D-2.or.con.gt.1.0d0) THEN
		print*,ysm,con,node,kinc
		print*,'error in yield strength or concentration from FRIC'
		call xit  
      ENDIF

! Thicknesses for various load cases from Krishan's experiment	  
	  thk1(1,1)=props(7)
	  thk1(2,1)=props(8)
	  thk1(3,1)=props(9)
	  thk1(4,1)=props(10)
	  thk1(5,1)=props(11)
	  thk1(6,1)=props(12)
	  thk1(7,1)=props(13)
	  thk1(8,1)=props(14)
	  thk1(9,1)=props(15)
	  thk1(10,1)=props(16)
	  thk1(11,1)=props(17)
	  thk1(12,1)=props(18)
	  thk1(13,1)=props(19)
	  thk1(14,1)=props(20)
	  thk1(15,1)=props(21)
! Please note thickness data 15 is removed as concentration data 15 is not available
	  thk1(16,1)=props(22)
	  thk1(17,1)=props(23)
	  thk1(18,1)=props(24)
	  thk1(19,1)=props(25)
	  thk1(20,1)=props(26)
	  thk1(21,1)=props(27)
	  thk1(22,1)=props(28)
	  thk1(23,1)=props(29)
	  thk1(24,1)=props(30)
	  thk1(25,1)=props(31)
	  thk1(26,1)=props(32)
	  
! Please note: below thicknesses are for various load cases from DAC experiments	 

	  thk1(1,1)=0.165d0
	  thk1(2,1)=0.16477d0	  
	  thk1(3,1)=0.16300d0
	  thk1(4,1)=0.15432d0
	  thk1(5,1)=0.14313d0
	  thk1(6,1)=0.12844d0
	  thk1(7,1)=0.11882d0
	  thk1(8,1)=0.11217d0
	  thk1(9,1)=0.10518d0
	  thk1(10,1)=0.096852d0
	  thk1(11,1)=0.091421d0
	  thk1(12,1)=0.087008d0
	  thk1(13,1)=0.080639d0
	  thk1(14,1)=0.076394d0
	  thk1(15,1)=0.068753d0
	  thk1(16,1)=0.058392d0
	  thk1(17,1)=0.054014d0
	  thk1(18,1)=0.048949d0
	  thk1(19,1)=0.046396d0
	  thk1(20,1)=0.043092d0
	  thk1(21,1)=0.039124d0	  
	  thk1(22,1)=0.037096d0
	  thk1(23,1)=0.033383d0
	  thk1(24,1)=0.03107d0
	  thk1(25,1)=0.030553d0
	  thk1(26,1)=0.028457d0	   
	  
	  thk2=thick+thkcor
C****************************************************************************
! Here, it is the source of uncertainty. Please take care
C****************************************************************************

	  if (thk2.ge.thk1(3,1)) then
		xmu0=1.5d0;xmu1=1.5d0
		
	  ! ! ! elseif (thk2.ge.thk1(4,1)) then
		! ! ! xmu0=1.0d0;xmu1=1.0d0
		
	  elseif (thk2.ge.thk1(7,1)) then
		xmu0=0.20d0;xmu1=0.20d0  
	  elseif (thk2.ge.thk1(10,1)) then
		xmu0=0.30d0;xmu1=0.30d0
	  elseif (thk2.ge.thk1(11,1)) then
		xmu0=0.5d0;xmu1=0.5d0

	  else 
		xmu0=0.5d0;xmu1=0.5d0
	  endif
 
  
C****************************************************************************
	  
! 'xmum' is Poisson's ratio of mixture at current contact pressure and concentration      
	  xmum = (1.0-con)*xmu0+con*xmu1+(xmup0*(1.0-con)+xmup1*con)*press

! Contact pressure derivative of 'mu' of mixture
	  dxmumdp = xmup0*(1-con)+xmup1*con

! Contact pressure derivative of yield strength in shear of mixture	  
	  dyshdp=1.0d0/3.0d0*(1.0d0-con)*ysp0/(3.D0**0.5)+
     +	1.0d0/3.0d0*con*ysp1/(3.D0**0.5)
	 
! sqrt(3) in the above line is to get the shear strength pressure coefficient from yield strength
!  pressure coefficient in tension or compression.	     
	  
! mp is number of radial points for each load case
! '26' is the number of load cases +1 ('+1' is for various radii)
	  
      OPEN (31,file='C:\Krishan5\MatlabFiles2\msol_AN.txt')   ! This is 'm' that we read from 'm.txt' obtained from Matlab
        DO I=1,mp 
          READ(31,*) xm1(I,:)
        ENDDO
      CLOSE(31)     
	  
! Interpolating for correct radial distribution of 'm' based on thickness of sample at centre from experiments and UMAT.
	  if (thk2.ge.thk1(1,1)) then
		xm2(:,1) = xm1(:,2)
		
	  elseif (thk2.lt.thk1(26,1)) then
		xm2(:,1) = xm1(:,27)
			
	  else
		do i=2,26		
		  if (thk2.ge.thk1(i,1)) then			
			xm2(:,1)=xm1(:,i)+(xm1(:,i+1)-xm1(:,i))*
     +	    (thk2-thk1(i-1,1))/(thk1(i,1)-thk1(i-1,1))			  
		    exit
		  endif
		enddo
	  endif	  

! Interpolating for correct 'm' based on 'x' coordinate of node.	  
      if (COORDS(1).ge.xm1(mp,1)) then
        xm=xm2(mp,1)
!!!		xm=1  
	  elseif (coords(1).eq.xm1(1,1)) then			 
		xm=xm2(1,1)
		  
      else
        DO I=2,mp
		   if (coords(1).le.xm1(i,1)) then
			xm=xm2(i-1,1)+(xm2(i,1)-xm2(i-1,1))*
     +	(coords(1)-xm1(i-1,1))/(xm1(i,1)-xm1(i-1,1))			
		   exit
		   endif
		enddo
	  endif
		 
! Implementation of 3-D isotropic Coulomb friction law using elastic stick formulation

! Variables used:
! GCRIT = Critical elastic slide
! GAMMA(1),GAMMA(2) are elastic predictor slides in directions '1' and '2' respectively.
! STATEV(1),STATEV(2) are elastic slides at the start of increment in directions '1' and '2'. 

! CHECK IF PRESSURE IS NON-POSITIVE

	  IF (PRESS.LE.ZERO .OR. LM.EQ.2) THEN
		STATEV(1) = ZERO
		STATEV(2) = ZERO
		GAMMA(1) = DGAM(1)
		GAMMA(2) = DGAM(2)
		TAU(1) = ZERO
		TAU(2) = ZERO
		DSLIP(1) = ZERO
		DSLIP(2) = ZERO	
		
	  RETURN
	  ELSE
	  
	  LM = 0
	  GCRIT = 0.005d0*CHRLNGTH
	  YIELD = xM*ysm/(3.D0**0.5) 
! 'yield' is yield strength in shear
	  
! Computing critical stresses and artificial stiffnesses
	  TAUCRIT1 = xmum*PRESS
	  TAUCRIT2 = YIELD
	  
	  if (coords(1).le.0.25d0) then
		taucrit1=2.0d0*taucrit2	
	  endif  
	  
c Achyut- see the condition below	  
c      taucrit1=2.0d0*taucrit2
	  
	  if (taucrit1.lt.0.0d0 .or. taucrit2 .lt. 0.0d0) then
		print*,'Critical stresses are negative in FRIC'
		call xit
	  endif
	  	  
	  STIFF1 = TAUCRIT1/GCRIT
	  STIFF2 = TAUCRIT2/GCRIT	  
	  end if

! Computing elastic slide predictors
		GAMMA(1) = STATEV(1) + DGAM(1)
		GAMMA(2) = STATEV(2) + DGAM(2)
		
	  if (taucrit1 .le. taucrit2) then

		TAU(1) = STIFF1*GAMMA(1)
		TAU(2) = STIFF1*GAMMA(2)
		TAUEQV = dSQRT(TAU(1)**2 + TAU(2)**2)

		IF (TAUEQV .LT. TAUCRIT1) THEN	! Elastic stick regime
			DDTDDG(1,1) = STIFF1
			DDTDDG(1,2) = ZERO
			DDTDDG(2,1) = ZERO
			DDTDDG(2,2) = STIFF1
			DDTDDP(1) = (xmum+press*dxmumdp)*GAMMA(1)/GCRIT
			DDTDDP(2) = (xmum+press*dxmumdp)*GAMMA(2)/GCRIT
			DSLIP(1) = ZERO
			DSLIP(2) = ZERO
			STATEV(1) = GAMMA(1)
			STATEV(2) = GAMMA(2)				
			
		ELSE							! Sliding regime
			GAMEQV = SQRT(GAMMA(1)**2 + GAMMA(2)**2)
			DDTDDG(1,1) = TAUCRIT1/GAMEQV*(1-(GAMMA(1)/GAMEQV)**2)
	  DDTDDG(1,2)=-TAUCRIT1/GAMEQV*(GAMMA(1)/GAMEQV)*(GAMMA(2)/GAMEQV)
			DDTDDG(2,1) = DDTDDG(1,2)
			DDTDDG(2,2) = TAUCRIT1/GAMEQV*(1-(GAMMA(2)/GAMEQV)**2)
			DDTDDP(1) = (xmum+press*dxmumdp)*GAMMA(1)/GAMEQV
			DDTDDP(2) = (xmum+press*dxmumdp)*GAMMA(2)/GAMEQV
			TAU(1) = GAMMA(1)/GAMEQV*TAUCRIT1
			TAU(2) = GAMMA(2)/GAMEQV*TAUCRIT1

! Computing elastic slide and plastic slip
C			STATEV(1) = GAMMA(1)*GCRIT/GAMEQV
C 			STATEV(2) = GAMMA(2)*GCRIT/GAMEQV
! I don't think that above line is correct because gamma(1) is elastic predictor and not the
! elastic slide at the end of increment. I am wrong! Above is correct too!				
			STATEV(1) = TAU(1)*GCRIT/TAUCRIT1
			STATEV(2) = TAU(2)*GCRIT/TAUCRIT1		
			DGSLEQ = GAMEQV - GCRIT
C			DSLIP(1) = GAMMA(1)*DGSLEQ/GAMEQV
C			DSLIP(2) = GAMMA(2)*DGSLEQ/GAMEQV
! Using the same arguments, the above statment is also not correct. The below is correct.			
			DSLIP(1) = TAU(1)/TAUCRIT1*DGSLEQ
			DSLIP(2) = TAU(2)/TAUCRIT1*DGSLEQ
					
			
			flag_10 = DGSLEQ + GCRIT
			if (flag_10.le.zero) then
			 print*,node,kinc
			 print*,'From Fric!'
			 call xit
			endif
! If the flag_10 activates the 'xit' function, then we have to introduce the negative sign in the 
! plastic regime as mentioned in Fric_Formulation document.

! Now I think the idea about sign is redundant and not completely correct.
			
		ENDIF
		
	  else
	
		TAU(1) = STIFF2*GAMMA(1)
		TAU(2) = STIFF2*GAMMA(2)
		TAUEQV = dSQRT(TAU(1)**2 + TAU(2)**2)

		! ! ! IF (TAUEQV .LT. TAUCRIT2 .and. kinc.le.135) THEN	! Elastic stick regime
		IF (TAUEQV .LT. TAUCRIT2) THEN	! Elastic stick regime

     ! ! ! 1  (kinc.le.3 .or. kinc.eq.23 .or. kinc.eq.24 		
     ! ! ! 2  .or. kinc.eq.36 .or. kinc.eq.37 .or. kinc.eq.38 
     ! ! ! 3  .or. kinc.eq.39 .or. kinc.eq.40 .or. kinc.eq.54 
     ! ! ! 4  .or. kinc.eq.55 .or. kinc.eq.56 .or. kinc.eq.57	
     ! ! ! 5  .or. kinc.eq.58 .or. kinc.eq.59 .or. kinc.eq.60
     ! ! ! 6  .or. kinc.eq.61 .or. kinc.eq.129 .or. kinc.eq.130 )) THEN

			DDTDDG(1,1) = STIFF2
			DDTDDG(1,2) = ZERO
			DDTDDG(2,1) = ZERO
			DDTDDG(2,2) = STIFF2
			DDTDDP(1) = GAMMA(1)/GCRIT*dyshdp
			DDTDDP(2) = GAMMA(2)/GCRIT*dyshdp
			DSLIP(1) = ZERO
			DSLIP(2) = ZERO
			STATEV(1) = GAMMA(1)
			STATEV(2) = GAMMA(2)	
			
		ELSE									! Sliding regime
			GAMEQV = SQRT(GAMMA(1)**2 + GAMMA(2)**2)
			DDTDDG(1,1) = TAUCRIT2/GAMEQV*(1-(GAMMA(1)/GAMEQV)**2)
	  DDTDDG(1,2)=-TAUCRIT2/GAMEQV*(GAMMA(1)/GAMEQV)*(GAMMA(2)/GAMEQV)
			DDTDDG(2,1) = DDTDDG(1,2)
			DDTDDG(2,2) = TAUCRIT2/GAMEQV*(1-(GAMMA(2)/GAMEQV)**2)
			DDTDDP(1) = (GAMMA(1)/GAMEQV)*dyshdp
			DDTDDP(2) = (GAMMA(2)/GAMEQV)*dyshdp
			TAU(1) = GAMMA(1)*TAUCRIT2/GAMEQV
			TAU(2) = GAMMA(2)*TAUCRIT2/GAMEQV
			
! Computing elastic slide and plastic slip
C			STATEV(1) = GAMMA(1)*GCRIT/GAMEQV
C 			STATEV(2) = GAMMA(2)*GCRIT/GAMEQV
! I don't think that above line is correct because gamma(1) is elastic predictor and not the
! elastic slide at the end of increment. I am wrong! Above is correct too!		
			STATEV(1) = TAU(1)*GCRIT/TAUCRIT2
			STATEV(2) = TAU(2)*GCRIT/TAUCRIT2		
			DGSLEQ = GAMEQV - GCRIT
C			DSLIP(1) = GAMMA(1)*DGSLEQ/GAMEQV
C			DSLIP(2) = GAMMA(2)*DGSLEQ/GAMEQV
! Using the same arguments, the above statment is also not correct. I am wrong! Above is correct too!			
			DSLIP(1) = TAU(1)/TAUCRIT2*DGSLEQ
			DSLIP(2) = TAU(2)/TAUCRIT2*DGSLEQ
			

		ENDIF
	  END IF	


	  RETURN
	  END  

C****************************************************************************

      SUBROUTINE fprint(datarec) 
	 
      INCLUDE 'ABA_PARAM.INC'
	  
	  Real*8 datarec(26,3)	  
	 	  
	  OPEN(10,file = 'C:\Krishan5\Kinetics\Kinetics.txt',	
     1 	  	status='replace',ACTION='READWRITE')
	 
        DO I=1,26 
          WRITE(10,'(F,F,F,F)') 
     1 	  	datarec(I,1),datarec(I,2),datarec(I,3),datarec(I,4)
        ENDDO	 
		  
	  CLOSE(10)

	  RETURN
	  END 	

C****************************************************************************
	  
