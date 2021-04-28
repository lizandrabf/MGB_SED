!@ *********************************************************************************************
!@ Criado por Diogo Buarque
!@ Data: Abr de 2011
!@
!@ Atualizado: Abr 2011
!@
!@ >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> SUBROTINA SED_CT >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!@
!@ - ESTA ROTINA CALCULA A CAPACIDADE DE TRANSPORTE DAS PARTÍCULAS DE SEDIMENTO NO TRECHO DE RIO
!@   POR DIFERENTES FORMULAÇÕES
!@
!@ *********************************************************************************************

SUBROUTINE SED_CT(Uat, VmTREC, SfTREC, HmJ)

USE VARS_MAIN
USE SED_VARS

IMPLICIT NONE

INTEGER i
REAL Uat, VmTREC, SfTREC, HmJ
!!!!!!!!!!!!YANG!!!!!!!!!!!!!       Hugo Fagundes 07/08/17        
REAL REX(3), UcWs(3)
REAL FMY(3), FNY(3), LOGCTS, LOGCTSaux(3)
REAL aux

!!!!!!!!!!!!MEYER-PETER E MULLER!!!!!!!!!!!!!       Hugo Fagundes 07/08/17        
REAL GSS !Peso espécífico do sedimento (t/m3) submerso (Carvalho, 2008)
REAL Kr,Ks
REAL D90 !D90[m] médio dissertação Livia Meneghel
REAL G !Peso especifico da água (t/m3) (Carvalho, 2008)
REAL QLT !parte da descarga líquida que influencia no leito (L/s)... verificar pubs.usgs.gov/wri/1989/4026/report.pdf

!!!!!!!!!!!!ENGELUND E HANSEN!!!!!!!!!!!!!       Hugo Fagundes 07/08/17        
REAL, ALLOCATABLE:: VEH(:) !Velocidade média de jusante do trecho em ft/s
REAL GEH !Peso específico do sedimento em lbf/ft³
REAL SG !Peso específico do sedimento/ Peso específico da água
REAL HEH !Profundidade média de jusante do trecho em lb
REAL grEH !gravidade em ft/s
REAL, ALLOCATABLE:: BEH(:) !Largura do rio em ft

!!!!!!!!!!!!ACKERS AND WHITE!!!!!!!!!!!!!       Hugo Fagundes 07/08/17        
REAL UXAW, RhAW, DXAW !velocidade de atrito, raio hidraulico, diametro adimensional
REAL D35
REAL alfaAW, nAW, AAW, mAW, CA, Fgr, Ggr !coeficientes da equação Dissertação Meneghel (2012)

!!!!!!!!!!!!!!!!!!!SCHOKLITSCH!!!!!!!!!!!!!!!!!!!       Hugo Fagundes 07/08/17        
REAL QoS

!!!!!Alocar variáveis!!!!!       Hugo Fagundes 07/08/17        
IF (xyz==1) THEN
    ALLOCATE (VEH(1173))
    ALLOCATE (BEH(1173))
    xyz=2
ENDIF


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!ROTTNER!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       Hugo Fagundes 07/08/17        
    IF (FLAG_TC==3) THEN
        VEH(IC)=0.0
        BEH(IC)=0.0
        VEH(IC)=3.2808*QJ2(IC)/(BRIO(IC)*HmJ)
        HEH=3.2808*HmJ
        grEH= 32.2175 !!!!SED_VARS??????
        SG= 2.65 !!!!SED_VARS??????
        GEH = 168.70 !!!!SED_VARS??????
        BEH(IC)=BRIO(IC)*3.2808
        
        DO i= 1,3
            CTS(IC,i)=BEH(IC)*GEH*(((SG-1)*grEH*(HEH**3))**0.5)*(((VEH(ic)/(((SG-1)*grEH*HEH)**0.5))*(0.667*((DMP(i)*3.2808/HEH)**(2/3))+0.14)-0.778*((DMP(i)*3.2808/HEH)**(2/3)))**3) !lb/s
            CTS(IC,i)=CTS(IC,i)*0.0454/(10**3) !TON/s
            CTS(IC,i)=CTS(IC,i)/QJ2(IC) !(TON/M3)   
        ENDDO
    ENDIF


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!SCHOKLITSCH!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       Hugo Fagundes 07/08/17        
    IF (FLAG_TC==4) THEN
        DO i= 1,3
            QoS=(0.00532*DMP(i)*39.37)/(SfTREC**(4/3))
            CTS(IC,i)=(86.7*(SfTREC**1.5)*(35.31*QJ2(IC)-BRIO(IC)*3.2808*QoS))/((0.008*39.37)**0.5) !lb/s
            CTS(IC,i)=CTS(IC,i)*0.0454/(10**3) !TON/s
            CTS(IC,i)=CTS(IC,i)/QJ2(IC) !(TON/M3)   
        ENDDO
    ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ACKERS E WHITE para AREIAS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       Hugo Fagundes 07/08/17        
    IF (FLAG_TC==5) THEN
        SG= 2.65
        alfaAW=10
        D35=0.0003
    
    
        RhAW=(BRIO(IC)*HmJ)/(BRIO(IC)+2*HmJ)
        UXAW= (9.81*RhAW*SfTREC)**0.5
        DXAW=D35*(((9.82*(SG-1))/(VISC**2))**(1/3))
        nAW=1-0.56*log10(DXAW)
        AAW=0.23/(DXAW**0.5)+0.14
        mAW=9.66/DXAW+1.34
        CA=10**(2.86*log10(DXAW)-(log10(DXAW)**2)-3.53)
    
        Fgr=(UXAW/(9.82*D35*(SG-1)))*((VmTREC/(5.657*log10(alfaAW*HmJ/D35)))**(1-nAW))
        Ggr=CA*((Fgr/AAW-1)**mAW)
    
    
        CTS(IC,1)=(10**6)*Ggr*SG*D35/(HmJ*((UXAW/VmTREC)**nAW)) !ppm
    
        write(*,*) CA, Ggr
    
        CTS(IC,1)=CTS(IC,1)/(10**6) !(TON/M3)
        CTS(IC,2)=CTS(IC,1)
        CTS(IC,3)=CTS(IC,1)
    ENDIF


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ENGELUND E HANSEN!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       Hugo Fagundes 07/08/17        
    IF (FLAG_TC==6) THEN
        VEH(IC)=0.0
        BEH(IC)=0.0
        VEH(IC)=3.2808*QJ2(IC)/(BRIO(IC)*HmJ)
        HEH=3.2808*HmJ
        grEH= 32.2175 !!!!SED_VARS??????
        SG= 2.65 !!!!SED_VARS??????
        GEH = 168.70 !!!!SED_VARS??????
        BEH(IC)=BRIO(IC)*3.2808
    
        DO i= 1,3
    
            CTS(IC,i)=(BEH(IC)*0.05*GEH*(VEH(IC)**2)*(HEH**0.5)*(SfTREC**1.5))/((DMP(i)*3.2808)*grEH*(1.65**2)) !lb/s - Carvalho (2008)
            CTS(IC,i)=CTS(IC,i)*0.0454/1000 !TON/s
            !!CTS(IC,i)=(BRIO(IC)*0.05*26500*(VmTREC**2)*(10000**1.5)*(HmJ**1.5)*(SfTREC**1.5))/((DMP(i))*(9.82**1.5)*(1.65**0.5)*(16500**1.5)) !kg/s - Reservoir Sedimentatio Handbook
            !!CTS(IC,i)=CTS(IC,i)/1000 !TON/s
            CTS(IC,i)=CTS(IC,i)/QJ2(IC) !(TON/M3)
    
    
        ENDDO
    ENDIF
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! MEYER-PETER E MULLER !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       Hugo Fagundes 07/08/17        
    IF (FLAG_TC==1) THEN
        Kr=26/(D90**(1/6))
        Ks=1/RUGMAN(IC)
        GSS=1.65 !!!!SED_VARS??????
        D90=0.00157  !!!!SED_VARS??????
        G=1 !!!!SED_VARS??????
        DO i= 1,3
    
            QLT=QJ2(IC)*0.6
    
            CTS(IC,i)=BRIO(IC)*((((G*(QLT/QJ2(IC))*((Ks/Kr)**1.5))*HmJ*SfTREC-0.047*GSS*DMP(i))/((G/9.82)**(1/3)))**1.5) !(TON/s)
            CTS(IC,i)=CTS(IC,i)/QJ2(IC) !(TON/M3)
  
        ENDDO
    ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! YANG !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
IF (FLAG_TC==2) THEN
    DO i= 1,3
        !@ Reynolds de atrito - mínimo igual a 1.2 (VERIFICAR OPÇÕES MELHORES)
        REX(i) = max(Uat*DMP(i)/VISC,1.2)

        !@ Estimativa da Velocidade Crítica sobre Velocidade de Queda (Uc/Ws)
	    IF (REX(i) < 70.) THEN
		    UcWs(i) = 0.66 + 2.5/(log10(REX(i)) - 0.06)
	    ELSE
		    UcWs(i) = 2.05
	    ENDIF

        !@ Parâmetros M e N da equação de capacidade de transporte de YANG
	    IF (DMP(i) <= 0.002) THEN
		    FMY(i) = 5.435 - 0.286*log10(WSP(i)*DMP(i)/VISC) - 0.457*log10(Uat/WSP(i))
		    FNY(i) = 1.799 - 0.409*log10(WSP(i)*DMP(i)/VISC) - 0.314*log10(Uat/WSP(i))
	    ELSE
		    FMY(i) = 6.681 - 0.633*log10(WSP(i)*DMP(i)/VISC) - 4.816*log10(Uat/WSP(i))
		    FNY(i) = 2.874 - 0.305*log10(WSP(i)*DMP(i)/VISC) - 0.282*log10(Uat/WSP(i))
	    ENDIF

        aux = VmTREC*SfTREC/WSP(i) - UcWs(i)*SfTREC
        LOGCTS = 0.0
        CTS(IC,i) = 0.0
        IF (aux > 0.) THEN
            LOGCTS = FMY(i) + FNY(i)*log10( aux )
            LOGCTSaux(i) = LOGCTS
            !@ Potencial de transporte do trecho de rio 
	        CTS(IC,i) = 10.**(LOGCTS) !@ (PPM)
            !@ Potencial de transporte do trecho de rio
	        CTS(IC,i) = CTS(IC,i)/(10.**(6.)) !@ (TON/M3)
            !@ Capacidade de transporte do trecho de rio
	        CTS(IC,i) = CTS(IC,i)*FracS(IC,i) !@ (TON/M3)
        ENDIF
  
    ENDDO
ENDIF
!IF (IC == 3834) THEN
!write(*,*) 'IT, IC, Sf  = ', IT, IC, SfTREC
!write(*,*) 'DMP, Visc   = ', DMP(1), VISC
!write(*,*) 'REX, WcWs   = ', REX(1), UcWs(1)
!write(*,*) 'FMY, FNY    = ', FMY(1), FNY(1)
!write(*,*) 'VERIF, Log  = ', aux, LOGCTSaux(1)
!write(*,*) 'CTS*, Frac  = ', CTS(IC,1), FracS(IC,1)
!pause
!ENDIF


RETURN
END SUBROUTINE