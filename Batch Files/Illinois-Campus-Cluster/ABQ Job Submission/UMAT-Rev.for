!-------------------------------------------------------------------------
!---User Subroutine for nonlinear granular material, MEPDG model----------
!---------------By: Hao Wang - Modified By: Izak M. Said------------------
!---User Subroutine for stress softening soils, Bilinear model------------
!------------------------By: Izak M. Said---------------------------------
!-------------------------------------------------------------------------
!
	SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     &RPL,DDSDDT,DRPLDE,DRPLDT,
     &STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     &NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     &CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
!
	INCLUDE 'ABA_PARAM.INC'
!
	CHARACTER*80 CMNAME
	DIMENSION STRESS(NTENS),STATEV(NSTATV),
     &DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     &STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     &PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)
!
! Programmer defined variables below
!
	REAL*8:: oct_shear,EMODH,EMODS,EMOD
	REAL*8:: v13,v21,v12,tt
	REAL*8:: MRMAX, K1, K2, K3, K4, MRMIN, MSTRESS
	REAL*8:: ANU, E, NUM1, NUM2, NUM3
	REAL*8:: D,rr
	REAL*8:: THETA=0.D0
	REAL*8, DIMENSION(3) :: PS=0
	REAL*8, DIMENSION(3,3) :: AN=0
	REAL*8, DIMENSION(3) :: PSTEST=0
	REAL*8, DIMENSION(3,3) :: ANTEST=0
	REAL*8, DIMENSION(NTENS) :: MODSTRESS
	REAL*8, DIMENSION(6) :: MODSTRESSTEST
	REAL*8, DIMENSION(NTENS) :: SIGMASTRESS
	INTEGER::LSTR
	PARAMETER(ZERO=0.D0, ONE=1.D0, TWO=2.D0, THREE=3.D0, SIX=6.D0,
     &ENUMAX=.4999D0, NEWTON=10, TOLER=1.0D-6, GRAVITY=9.81)
!
	IF (NDI/=3) THEN
	WRITE (7, *) 'THIS UMAT MAY ONLY BE USED FOR 3D and AXISYMMETRIC
     &ELEMENTS WITH THREE DIRECT STRESS COMPONENTS'
	CALL XIT
	END IF
!
	IF (CMNAME(1:4) .EQ. 'P2') THEN
	   CALL Base(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     &RPL,DDSDDT,DRPLDE,DRPLDT,
     &STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     &NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     &CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
	ELSE IF(CMNAME(1:4) .EQ. 'P3') THEN
	   CALL Subgrade(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     &RPL,DDSDDT,DRPLDE,DRPLDT,
     &STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     &NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     &CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
	END IF
	END SUBROUTINE UMAT
!-----------------------------------------------------------------------
! If Stress Dependent - Hardening - Loop for 3D
!-----------------------------------------------------------------------
	SUBROUTINE Base(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     &RPL,DDSDDT,DRPLDE,DRPLDT,
     &STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     &NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     &CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
!
	INCLUDE 'ABA_PARAM.INC'
!
	CHARACTER*80 CMNAME
	DIMENSION STRESS(NTENS),STATEV(NSTATV),
     &DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     &STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     &PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)
!
! Programmer defined variables below
!
	REAL*8:: oct_shear,EMODH,EMODS,EMOD
	REAL*8:: v13,v21,v12,tt
	REAL*8:: MRMAX, K1, K2, K3, K4, MRMIN, MSTRESS
	REAL*8:: ANU, E, NUM1, NUM2, NUM3
	REAL*8:: D,rr
	REAL*8:: THETA=0.D0
	REAL*8, DIMENSION(3) :: PS=0
	REAL*8, DIMENSION(3,3) :: AN=0
	REAL*8, DIMENSION(3) :: PSTEST=0
	REAL*8, DIMENSION(3,3) :: ANTEST=0
	REAL*8, DIMENSION(NTENS) :: MODSTRESS
	REAL*8, DIMENSION(6) :: MODSTRESSTEST
	REAL*8, DIMENSION(NTENS) :: SIGMASTRESS
	INTEGER::LSTR
	PARAMETER(ZERO=0.D0, ONE=1.D0, TWO=2.D0, THREE=3.D0, SIX=6.D0,
     &ENUMAX=.4999D0, NEWTON=10, TOLER=1.0D-6, GRAVITY=9.81)
!
!----------------------------------------------------------------------
! UMAT for anisotropic elsticity with principal stresses
! Cannot be used for plane stress
!----------------------------------------------------------------------
!---------------Stress-Depedent Input Paramters -----------------------
!----------------------------------------------------------------------
! PROPS(1) - Material model - 1
! PROPS(2) - NU12
! PROPS(3) - K1
! PROPS(4) - K2
! PROPS(5) - K3
! PROPS(6) - NU31
! PROPS(7) - 1 for principal stress / 0 for S11/S22/S33 to calculate EMOD
! PROPS(8) - Minimum EMODH (MPa)
! PROPS(9) - Minimum EMODS (MPa)
! PROPS(10) - Minimum EMOD (MPa)
! PROPS(11) - EMODH for first increment (zero applied stress) (MPa)
! PROPS(12) - EMODS for first increment (zero applied stress) (MPa)
! PROPS(13) - EMOD for first increment (zero applied stress) (MPa)
! PROPS(14) - k4
! PROPS(15) - k5
! PROPS(16) - k6
! PROPS(17) - k7
! PROPS(18) - k8
! PROPS(19) - k9
! PROPS(20) - Depth of pavement 
! PROPS(21) - Density of Base Materials
!-----------------------------------------------------------------------
! Change sign on stress values to reflect geomech principal of
! compression = positive
! set stress to 1 kPa if zero or tensile. These values are used in the EMOD
! calcs only
! Units of stress are in MPa
! Modified values stored in SIGMASTRESS array which has the same dimension
! as STRESS
!
!-----------------------------------------------------------------------
! If Stress Dependent - Hardening - Loop for 3D
!-----------------------------------------------------------------------
	D=PROPS(20)
	rr=PROPS(21)
	IF (NSHR==3) THEN
	SIGMASTRESS=STRESS*-1
	MODSTRESS=STRESS*-1
	MODSTRESSTEST=STRESS*-1
	MODSTRESS(1)=MODSTRESS(1)+((D-COORDS(2))*1.0*rr*GRAVITY)
	MODSTRESS(2)=MODSTRESS(2)+((D-COORDS(2))*rr*GRAVITY)
	MODSTRESS(3)=MODSTRESS(3)+((D-COORDS(2))*1.0*rr*GRAVITY)
	SIGMASTRESS(1)=SIGMASTRESS(1)+((D-COORDS(2))*1.0*rr*GRAVITY)
	SIGMASTRESS(2)=SIGMASTRESS(2)+((D-COORDS(2))*rr*GRAVITY)
	SIGMASTRESS(3)=SIGMASTRESS(3)+((D-COORDS(2))*1.0*rr*GRAVITY)
!	
	DO K1=1, NTENS
	IF (SIGMASTRESS(K1)<0) THEN
	SIGMASTRESS(K1)=0.001
	END IF
	END DO
!
!Calculate principal stresses if PROPS(7)=1
	IF (PROPS(7)==1) THEN
	LSTR=1
	CALL SPRIND(MODSTRESSTEST,PSTEST,ANTEST,LSTR,3,3)
	IF (PSTEST(1)==ZERO .AND. PSTEST(2)==ZERO .AND. 
     &PSTEST(3)==ZERO) THEN
	K2=0
	ELSE IF (PSTEST(1)>=PSTEST(2) .AND. PSTEST(1)>=PSTEST(3)) THEN
	K2=1
	ELSE IF (PSTEST(2)>=PSTEST(1) .AND. PSTEST(2)>=PSTEST(3)) THEN
	K2=2
	ELSE IF (PSTEST(3)>=PSTEST(1) .AND. PSTEST(3)>=PSTEST(2)) THEN
	K2=3
	ELSE
	K2=0
	END IF
	IF (K2/=0) THEN
	THETA=ACOS(ANTEST(K2,3))
	END IF
!
	PSTEST=ZERO
	ANTEST=ZERO
	CALL SPRIND(MODSTRESS,PSTEST,ANTEST,LSTR,3,3)
	DO K1=1, NDI
	IF (PSTEST(K1)<0) THEN
	SIGMASTRESS(K1)=0.001
	ELSE
	SIGMASTRESS(K1)=PSTEST(K1)
	END IF
	END DO
	oct_shear=SQRT(((SIGMASTRESS(1)-SIGMASTRESS(2))*(SIGMASTRESS(1)
     &-SIGMASTRESS(2))+(SIGMASTRESS(2)-SIGMASTRESS(3))*(SIGMASTRESS(2)
     &-SIGMASTRESS(3))+(SIGMASTRESS(1)-SIGMASTRESS(3))*(SIGMASTRESS(1)
     &-SIGMASTRESS(3)))/9)
	END IF
! end of principal stress loop
!
	IF (PROPS(7)==0) THEN
	oct_shear=SQRT(((SIGMASTRESS(1)-SIGMASTRESS(2))*(SIGMASTRESS(1)
     &-SIGMASTRESS(2))+(SIGMASTRESS(2)-SIGMASTRESS(3))*(SIGMASTRESS(2)
     &-SIGMASTRESS(3))+(SIGMASTRESS(1)-SIGMASTRESS(3))*(SIGMASTRESS(1)
     &-SIGMASTRESS(3)))/9+2/3*(SIGMASTRESS(4)*SIGMASTRESS(4)+
     &SIGMASTRESS(5)*SIGMASTRESS(5)+SIGMASTRESS(6)*SIGMASTRESS(6)))
	END IF
!
	IF (PROPS(7)==ONE) THEN
	THETA=ZERO
	END IF
!
	STATEV(1)=SIGMASTRESS(1)
	STATEV(2)=SIGMASTRESS(2)
	STATEV(3)=SIGMASTRESS(3)
!
	IF (oct_shear<=0.001) THEN
	oct_shear=0.001
	END IF
!
	STATEV(4)=oct_shear
!
!
!
	IF (PROPS(1)==1) THEN
! Uzan Model for Granular material
! stresses are defined by abaqus in MPa, formula based on kPay
! orig formula as declared by Uzan, second formula adj for stress magnitude
!EMOD=k1*100*(((sig1+2*sig2)/100)**k2)*(((sig2-sig1)/100)*(2**0.5)/3)**k3
	EMOD=PROPS(3)*0.100*(((SIGMASTRESS(1)+SIGMASTRESS(2)+SIGMASTRESS(3))
     &*10)**PROPS(4))*((oct_shear/0.10+ONE)**PROPS(5))
	EMODH=PROPS(14)*0.100*(((SIGMASTRESS(1)+SIGMASTRESS(2)+SIGMASTRESS(3))
     &*10)**PROPS(15))*((oct_shear/0.10+ONE)**PROPS(16))
	EMODS=PROPS(17)*0.100*(((SIGMASTRESS(1)+SIGMASTRESS(2)+SIGMASTRESS(3))
     &*10)**PROPS(18))*((oct_shear/0.10+ONE)**PROPS(19))
	EMODMODIFY=EMOD/(1+(SIN(THETA)*SIN(THETA)-0.25*SIN(2*THETA))
     &*(PROPS(6)-1))
!
	IF (TIME(2)==0) THEN
	EMODH=PROPS(11)
	EMODS=PROPS(12)
	EMOD=PROPS(13)
	END IF
!
	IF (EMODH<PROPS(8)) THEN
	EMODH=PROPS(8)
	END IF
	IF (EMODS<PROPS(9)) THEN
	EMODS=PROPS(9)
	END IF
	IF (EMOD<PROPS(10)) THEN
	EMOD=PROPS(10)
	END IF
!
	END IF
	STATEV(5)=EMODH
	STATEV(6)=EMODS
	STATEV(7)=EMOD
	STATEV(8)=THETA
	STATEV(9)=EMODMODIFY
! ELASTIC PROPERTIES
	v13=PROPS(2)
	v21=PROPS(6)
	v12=EMODH/EMOD*v21
	tt=1/(1-v13*v13-2*v12*v21-2*v13*v12*v21)
! ELASTIC STIFFNESS
	DDSDDE=0
	DDSDDE(1,1)=EMODH*(1-v12*v21)*tt
	DDSDDE(3,3)=DDSDDE(1,1)
	DDSDDE(2,2)=EMOD*(1-v13*v13)*tt
	DDSDDE(5,5)=EMODH/(2*(1+v13))
	DDSDDE(4,4)=EMODS
	DDSDDE(6,6)=DDSDDE(4,4)
	DDSDDE(1,3)=EMODH*(v13+v21*v12)*tt
	DDSDDE(3,1)=DDSDDE(1,3)
	DDSDDE(1,2)=EMODH*(v21+v13*v21)*tt
	DDSDDE(2,3)=DDSDDE(1,2)
	DDSDDE(2,1)=DDSDDE(1,2)
	DDSDDE(3,2)=DDSDDE(1,2)
! CALCULATE STRESS
	DO K1=1, NTENS
	DO K2=1, NTENS
	STRESS(K2)=STRESS(K2)+DDSDDE(K2, K1)*DSTRAN(K1)
	END DO
	END DO
	END IF
!*********************************************************************
! end of 3D loop
!*********************************************************************
!*********************************************************************
	IF (TIME(1)>0) THEN
	PNEWDT=(ONE-TIME(1)-DTIME)/DTIME
	IF (PNEWDT<ONE) THEN
	PNEWDT=ONE
	END IF
	END IF
!
	RETURN
	END
!-----------------------------------------------------------------------
! If Stress Dependent - Softening - Loop for 3D
!-----------------------------------------------------------------------
	SUBROUTINE Subgrade(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     &RPL,DDSDDT,DRPLDE,DRPLDT,
     &STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     &NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     &CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
!
	INCLUDE 'ABA_PARAM.INC'
!
	CHARACTER*80 CMNAME
	DIMENSION STRESS(NTENS),STATEV(NSTATV),
     &DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     &STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     &PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)
!
! Programmer defined variables below
!
	REAL*8:: oct_shear,EMODH,EMODS,EMOD
	REAL*8:: v13,v21,v12,tt
	REAL*8:: MRMAX, K1, K2, K3, K4, MRMIN, MSTRESS
	REAL*8:: ANU, E, NUM1, NUM2, NUM3
	REAL*8:: D,rr
	REAL*8:: THETA=0.D0
	REAL*8, DIMENSION(3) :: PS=0
	REAL*8, DIMENSION(3,3) :: AN=0
	REAL*8, DIMENSION(3) :: PSTEST=0
	REAL*8, DIMENSION(3,3) :: ANTEST=0
	REAL*8, DIMENSION(NTENS) :: MODSTRESS
	REAL*8, DIMENSION(6) :: MODSTRESSTEST
	REAL*8, DIMENSION(NTENS) :: SIGMASTRESS
	INTEGER::LSTR
	PARAMETER(ZERO=0.D0, ONE=1.D0, TWO=2.D0, THREE=3.D0, SIX=6.D0,
     &ENUMAX=.4999D0, NEWTON=10, TOLER=1.0D-6, GRAVITY=9.81)
!
!-----------------------------------------------------------------------
!---------------Stress Softening Input Paramters -----------------------
!-----------------------------------------------------------------------
! PROPS(1) - Material model - 2
! PROPS(2) - MR Max
! PROPS(3) - K1
! PROPS(4) - K2
! PROPS(5) - K3
! PROPS(6) - K4
! PROPS(7) - MR Min
! PROPS(8) - Poisson
! PROPS(9) - Density of Subbase
! PROPS(10) - Depth of pavement
!----------------------------------------------------------------------
! Change sign on stress values to reflect geomech principal of
! compression = positive
! set stress to 1 kPa if zero or tensile. These values are used in the EMOD
! calcs only
! Units of stress are in MPa
! Modified values stored in SIGMASTRESS array which has the same dimension
! as STRESS
!-----------------------------------------------------------------------
	rr=PROPS(9)
	D=PROPS(10)
	IF (NSHR==3) THEN
	SIGMASTRESS=STRESS*-1
	MODSTRESS=STRESS*-1
	MODSTRESSTEST=STRESS*-1
	MODSTRESS(1)=MODSTRESS(1)+((D-COORDS(2))*0.6*rr*GRAVITY)
	MODSTRESS(2)=MODSTRESS(2)+((D-COORDS(2))*rr*GRAVITY)
	MODSTRESS(3)=MODSTRESS(3)+((D-COORDS(2))*0.6*rr*GRAVITY)
	SIGMASTRESS(1)=SIGMASTRESS(1)+((D-COORDS(2))*0.6*rr*GRAVITY)
	SIGMASTRESS(2)=SIGMASTRESS(2)+((D-COORDS(2))*rr*GRAVITY)
	SIGMASTRESS(3)=SIGMASTRESS(3)+((D-COORDS(2))*0.6*rr*GRAVITY)
!
!	Calculate principal stresses
!
	PSTEST=ZERO
	ANTEST=ZERO
    	LSTR=1
	CALL SPRIND(MODSTRESS,PSTEST,ANTEST,LSTR,3,3)
	DO K1=1, NDI
	IF (PSTEST(K1)<0) THEN
	SIGMASTRESS(K1)=10**-19
	ELSE
	SIGMASTRESS(K1)=PSTEST(K1)
	END IF
	END DO
!
	NUM1 = SIGMASTRESS(1)
	NUM2 = SIGMASTRESS(2)
	NUM3 = SIGMASTRESS(3)
	IF (SIGMASTRESS(1)>=SIGMASTRESS(2) .AND. SIGMASTRESS(1)>=SIGMASTRESS(3)) THEN
	SIGMASTRESS(2)= NUM1
	SIGMASTRESS(1)= NUM2
	SIGMASTRESS(3)= NUM3
	ELSE IF (SIGMASTRESS(2)>=SIGMASTRESS(1) .AND. SIGMASTRESS(2)>=SIGMASTRESS(3)) THEN
	SIGMASTRESS(1)= NUM1
	SIGMASTRESS(2)= NUM2
	SIGMASTRESS(3)= NUM3
	ELSE IF (SIGMASTRESS(3)>=SIGMASTRESS(1) .AND. SIGMASTRESS(3)>=SIGMASTRESS(2)) THEN
	SIGMASTRESS(1)= NUM1
	SIGMASTRESS(2)= NUM3
	SIGMASTRESS(3)= NUM2
	END IF
!
	MRMAX=PROPS(2)
	K1=PROPS(3)
	K2=PROPS(4)
	K3=PROPS(5)
	K4=PROPS(6)
	MRMIN=PROPS(7)
	MSTRESS = SIGMASTRESS(2) - (0.5*(SIGMASTRESS(1)+SIGMASTRESS(3)))
	IF (MSTRESS.LT.0.0206842800079188) THEN
	E=MRMAX
	END IF
	IF ((MSTRESS.GE.0.0206842800079188) .AND. (MSTRESS.LT.K4)) THEN
	E=MRMAX+(K1*(MSTRESS-0.0206842800079188))
	END IF
	IF ((MSTRESS.GE.K4) .AND. (MSTRESS.LT.0.186158520071269)) THEN
	E=K3+(K2*(MSTRESS-K4))
	END IF
	IF (MSTRESS.GE.0.186158520071269) THEN
	E=MRMIN
	END IF
!
! ELASTIC PROPERTIES
	ANU=PROPS(8)
	ALAMBDA=E/(ONE+ANU)/(ONE-TWO*ANU)
	BLAMBDA=(ONE-ANU)
    	CLAMBDA=(ONE-TWO*ANU)
! ELASTIC STIFFNESS
	DO I=1,NTENS
	  DO J=1,NTENS
	    DDSDDE(I,J)=0.0D0
	  ENDDO
	ENDDO
	DDSDDE(1,1)=(ALAMBDA*BLAMBDA)
	DDSDDE(2,2)=(ALAMBDA*BLAMBDA)
	DDSDDE(3,3)=(ALAMBDA*BLAMBDA)
	DDSDDE(4,4)=(ALAMBDA*CLAMBDA)/2.0
	DDSDDE(5,5)=(ALAMBDA*CLAMBDA)/2.0
	DDSDDE(6,6)=(ALAMBDA*CLAMBDA)/2.0
	DDSDDE(1,2)=(ALAMBDA*ANU)
	DDSDDE(1,3)=(ALAMBDA*ANU)
	DDSDDE(2,3)=(ALAMBDA*ANU)
	DDSDDE(2,1)=(ALAMBDA*ANU)
	DDSDDE(3,1)=(ALAMBDA*ANU)
	DDSDDE(3,2)=(ALAMBDA*ANU)
! CALCULATE STRESS
    	DO I=1,NTENS
	  DO J=1,NTENS
	    STRESS(I)=STRESS(I)+DDSDDE(I,J)*DSTRAN(J)
	  ENDDO
	ENDDO
	END IF
!*********************************************************************
! end of 3D loop
!*********************************************************************
!*********************************************************************
	IF (TIME(1)>0) THEN
	PNEWDT=(ONE-TIME(1)-DTIME)/DTIME
	IF (PNEWDT<ONE) THEN
	PNEWDT=ONE
	END IF
	END IF
!
	RETURN
	END