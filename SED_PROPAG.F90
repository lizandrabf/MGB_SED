!@ *********************************************************************************************
!@ Criado por Diogo Buarque
!@ Data: Março de 2013
!@
!@ Atualizado: Mar 2013
!@
!@ >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> SUBROTINA SED_REDE >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!@
!@ - SUBROTINA PARA PROPAGAR A CARGA DE SEDIMENTO AO LONGO DA DEDE DE DRENAGEM
!@
!@
!@ ---------------------------------------------------------------------------------------------
!@  VARIÁVEIS
!@ Descrição das variáveis locais:
!@
!@ ---------------------------------------------------------------------------------------------
!@
!@ *********************************************************************************************

SUBROUTINE SED_PROPAG(HmJ, VmJ, SfTREC, UatJ, Qlim, QtFL, VFL1X, VFL2X, HFLSD, DEPaux2)

!@ **********************************
USE VARS_MAIN
USE SED_VARS
USE VARS_INERC
USE VARS_CALIB
!@ **********************************

IMPLICIT NONE

!@ **********************************
!@ Declaração de variáveis:
REAL HmJ, VmJ, RHmJ, SfTREC, UatJ, Qlim
REAL LIMe, DeT, LIMd
REAL PONDt
REAL EQ1(3), EQ2(3), EQ3(3), EQ4(3), EQ5(3), EQ6(3), DEPaux2(NC,3)
REAL(8) DEPaux, EROSaux, DELTA(NC,3),cargatrec(NC,3),lamb,densiSED 
INTEGER i, JUS
REAL LIMeaux(3), LIMdaux(3)
REAL(8) iENTRAaux(3), iSAIaux(3)
REAL QtFL, HFLSD, VFL1X, VFL2X   !@ DCB_HD_Sed
REAL GFL(3), DFLaux(3)          !@ DCB_HD_Sed
!@ **********************************

!@ INICIALIZANDO VARIÁVEIS
PONDt  = 1.0   !@ Ponderador temporal da discretização por diferenças finitas
GFL    = 0.0
DFLaux = 0.0
lamb = 0.4 !hof 19.11.2020
densiSED= 1.5 !hof 19.11.2020

!@ #########################################################################################
!@ -------------------    CELULA COM RIO    ------------------------------------------------
!@ #########################################################################################
IF(NSUBT(IC).GT.0)THEN !((NSUBT(IC).GT.0).OR.(HdFLAG0>0))THEN
    !@ *****************************************************************************************
    !@ CONCENTRAÇÃO DE SEDIMENTOS NA SEÇÃO DE JUSANTE AO FINAL DO PASSO DE TEMPO(TON/M3)
    !@ -----------------------------------------------------------------------------------------
    !@  Equação de transporte: d(AC)/dt + d(AuC)/dx = ql
    !@ PARTICULAS DE AREIA NÃO INTERAGEM COM PLANÍCIES

    CSM2(IC,1) = 0.0
    IF (QM2(IC) > Qlim) CSM2(IC,1) = (CARGM(IC,1)+DEPaux2(IC,1))/(QM2(IC)*DTP) !Hugo Fagundes 19/08/2019
    CSJ2(IC,1) = 0.0
    
    !************************* Hugo Fagundes e Diogo Buarque 18.11.2020 *************************!
    !Equação de Exner assumindo que a variação da área do depósito (Delta Ab, equação 29. tese Diogo Burque)
    !de um dia para outro é desprezível.
    IF (QJ2(IC) > Qlim) THEN  !@ Só tem concentração a jusante se tem vazão
        EQ1(1) = PONDt*QM2(IC)*CSM2(IC,1) 
        EQ2(1) = (1.0-PONDt)*(QJ1(IC)*CSJ1(IC,1)-QM1(IC)*CSM1(IC,1)) 
        EQ3(1) = VolTREC1(IC)*(CSJ1(IC,1))/DTP
        EQ4(1) = VolTREC2(IC)/DTP + PONDt*QJ2(IC)
        CSJ2(IC,1) = (EQ1(1) - EQ2(1) + EQ3(1) + PSC(IC,2)/DTP)/EQ4(1)      !@ TON/M3
        
        !EQ1(1) = PONDt*QM2(IC)*CSM2(IC,1) 
        !EQ2(1) = (1.0-PONDt)*(QJ1(IC)*CSJ1(IC,1)-QM1(IC)*CSM1(IC,1)) 
        !EQ3(1) = PONDt*QJ2(IC)
        !CSJ2(IC,1) = (EQ1(1) - EQ2(1) + PSC(IC,2)/DTP)/EQ3(1)      !@ TON/M3
         
        IF (CSJ2(IC,1)<0.0) THEN !HOF 25/09/2020
            CSJ2(IC,1) =0.0;
        ENDIF
    ENDIF
    !@ -----------------------------------------------------------------------------------------
    !@ PARTÍCULAS DE SILTE E ARGILA SÃO TROCADAS ENTRE RIO E PLANÍCIES
    DO i= 2,3
    
        !@ Concentração na seção de montante (ton/m3)
        CSM2(IC,i) = 0.0
        IF (QM2(IC) > Qlim) CSM2(IC,i) = CARGM(IC,i)/(QM2(IC)*DTP)

        !@ Concentração na seção de jusante (ton/m3)
        CSJ2(IC,i) = 0.0

        IF (QJ2(IC) > Qlim .AND. QtFL >= 0.0) THEN  !@ Só tem concentração a jusante se tem vazão
            EQ1(i) = PONDt*QM2(IC)*(CSM2(IC,i)) 
            EQ2(i) = (1.0-PONDt)*(QJ1(IC)*(CSJ1(IC,i))-QM1(IC)*(CSM1(IC,i))) 
            EQ3(i) = (VolTREC1(IC)*(CSJ1(IC,i))/DTP) 
            EQ4(i) = (VolTREC2(IC)/DTP + PONDt*QJ2(IC))
            EQ5(i) = 0.5*QtFL*(0.5*(CSM2(IC,i)+CSJ1(IC,i)))   !HOF 
            EQ6(i) = (0.5*QtFL)              !@ DCB_HD_Sed    
            CSJ2(IC,i) = ABS((EQ1(i) - EQ2(i) + EQ3(i) + (PSC(IC,i+1)/DTP) - EQ5(i))/(EQ4(i) + EQ6(i)))  !@ TON/M3 
            
            IF (CSJ2(IC,i)<0.0) THEN !HOF 25/09/2020
                CSJ2(IC,i) =0.0;
            ENDIF
            
            GFL(i) = max(CFL1(IC,i)*VFL1(IC) + 0.5*QtFL*((CSJ2(IC,i)) + (CSM2(IC,i)))*DTP,0.0)   !@ TON/DIA 
            LIMd = 0.0
      
            IF(HFLSD>0.) THEN
                LIMd = min(WSP(i)*DTP/HFLSD, 1.0)
                DFLaux(i)  = GFL(i)*LIMd                !Carga que deposita na planície
            ELSEIF (CFL1(IC,i)*VFL1(IC)>0.0)  THEN      !Hugo Fagundes 08/06/2020
                DFLaux(i)  = GFL(i)
            ELSE                                        !Hugo Fagundes 08/06/2020
                DFLaux(i)  = 0.0
            ENDIF

            CFL2(IC,i) = 0.0
            IF(VFL2X>0.0) CFL2(IC,i) = max((GFL(i) - DFLaux(i))/VFL2X,0.0)
            DFL(IC,i)  = DFL(IC,i) + DFLaux(i)
 
        ELSEIF (QJ2(IC) > Qlim .AND. QtFL < 0.0) THEN
            EQ1(i) = PONDt*QM2(IC)*(CSM2(IC,i)) 
            EQ2(i) = (1.0-PONDt)*(QJ1(IC)*(CSJ1(IC,i))-QM1(IC)*(CSM1(IC,i))) 
            EQ3(i) = VolTREC1(IC)*(CSJ1(IC,i))/DTP 
            EQ4(i) = VolTREC2(IC)/DTP + PONDt*QJ2(IC)
            EQ5(i) = 0.5*min( QtFL*CFL1(IC,i) , CFL1(IC,i)*VFL1(IC)/DTP )  !@ DCB_HD_Sed Limita este termo, pois pode não haver volume suficiente na planície
            EQ6(i) = 0.0                   !@ DCB_HD_Sed
            CSJ2(IC,i) = ((EQ1(i) - EQ2(i) + EQ3(i) + (PSC(IC,i+1)/DTP) - EQ5(i))/(EQ4(i) + EQ6(i)))  !@ TON/M3     
            
            IF (CSJ2(IC,i)<0.0) THEN !HOF 25/09/2020
                CSJ2(IC,i) =0.0;
            ENDIF
            
            GFL(i) = max(CFL1(IC,i)*VFL1(IC) + QtFL*CFL1(IC,i)*DTP,0.0)  !@ TON/DIA
            LIMd = 0.0
            
            IF(HFLSD>0.) THEN
                LIMd = min(WSP(i)*DTP/HFLSD, 1.0)
                DFLaux(i)  = GFL(i)*LIMd                !Carga que deposita na planície     
            ELSEIF (CFL1(IC,i)*VFL1(IC)>0.0)  THEN      !Hugo Fagundes 08/06/2020
                DFLaux(i)  = CFL1(IC,i)*VFL1(IC)
            ELSE                                        !Hugo Fagundes 08/06/2020
                DFLaux(i)  = 0.0
            ENDIF
            
            IF(VFL2X>0.0) CFL2(IC,i) = max((GFL(i) - DFLaux(i))/VFL2X,0.0)
            DFL(IC,i)  = DFL(IC,i) + DFLaux(i)
       
        ENDIF

     
        !@ Acumulado da carga que entra na mini-bacia
        IF(QtFL>0.) inFL(iSEDaux,i)   = inFL(iSEDaux,i)  + 0.5*QtFL*(CSJ2(IC,i) + CSM2(IC,i))*DTP   !@ (ton)
        !@ Acumulado da vazão que entra na planície
        IF(QtFL>0.) BalQFL(iSEDaux,1)  = BalQFL(iSEDaux,1) + QtFL
        !@ Acumulado da carga que sai na mini-bacia
        IF(QtFL<0.) outFL(iSEDaux,i)  = outFL(iSEDaux,i) + QtFL*CFL1(IC,i)*DTP                      !@ (ton)
        !@ Acumulado da vazão que sai na planície
        IF(QtFL<0.) BalQFL(iSEDaux,2)  = BalQFL(iSEDaux,2) + QtFL
    
    ENDDO

    !@ CAPACIDADE DE TRANSPORTE DO TRECHO DE RIO EM TON/M3 (FORMULA DE YANG)
    !@ -----------------------------------------------------------------------------------------
    FracS(IC,:) = 0.0
    CTS(IC,:) = 0.0
    IF (QJ2(IC) > Qlim) THEN !@ Só tem capacidade a jusante se tem vazão                
        IF (CSJ2(IC,1) > 0.0) THEN
            FracS(IC,1) = CSJ2(IC,1)/sum(CSJ2(IC,:)) !@ porcentagem da concentração de areia
            FracS(IC,2) = CSJ2(IC,2)/sum(CSJ2(IC,:)) !@ porcentagem da concentração de silte
            FracS(IC,3) = CSJ2(IC,3)/sum(CSJ2(IC,:)) !@ porcentagem da concentração de argila
        ENDIF
        !@ Capacidade de transporte
        CALL SED_CT(UatJ, VmJ, SfTREC, HmJ)
    ENDIF
    !@ *****************************************************************************************

    !@ *****************************************************************************************
    !@ VERIFICAÇÃO DE TENDÊNCIA A EROSÃO OU DEPOSIÇÃO
    !@ -----------------------------------------------------------------------------------------
    LIMdaux     = 0.0   !@ IMPRESSÃO - Limitador da deposição no rio no passo de tempo
    LIMeaux     = 0.0   !@ IMPRESSÃO - Limitador da erosão no rio no passo de tempo
    iENTRAaux   = 0.0   !@ IMPRESSÃO - Carga que entra no rio no passo de tempo (ton/dia)
    iSAIaux     = 0.0   !@ IMPRESSÃO - Carga que sai do rio no passo de tempo (ton/dia)
    DO i= 1,3   !@ INICIO DO LOOP DAS PARTICULAS
        DEPaux  = 0.0
        EROSaux = 0.0
        LIMd    = 0.0
        LIMe    = 0.0
        DEPaux2 (IC,i) = 0.0   !Hugo Fagundes 19/08/2019
        !DELTA   = 0.0       !Hugo Fagundes 19/08/2019
        !cargatrec = 0.0     !Hugo Fagundes 19/08/2019
        
        IF (QJ2(IC) > Qlim) THEN   !@ Só há passagem de sedimentos se há vazão de jusante!
            
            !@ DEPOSIÇÃO ---------------------------------------------------------------------
            IF (CSJ2(IC,i) >= CTS(IC,i)) THEN
            !@ limitador temporal da deposição (igual ao HEC-RAS 4.0) - PARA D > 0.062 mm
                IF (DMP(i) <= 0.000062) THEN        !@ limite máximo silte (Wu, 2008)
                    CSJ2(IC,i) = CSJ2(IC,i)         !@ Toda carga de finos passa para jusante
                ELSEIF (DMP(i) <= 0.000125) THEN    !@ limite máximo areia muito fina (Wu, 2008)
                    DeT = Hmj        ! Distância efetiva
                    LIMd = min(WSP(i)*DTP/DeT, 1.0)
                    DEPaux = (CSJ2(IC,i) - CTS(IC,i))*(VolTREC2(IC)+QJ2(IC)*DTP)*LIMd   !@ (TON) !HUGO FAGUNDES 16/08/2019
                    DEPaux2(IC,i) = (CSJ2(IC,i) - CTS(IC,i))*(VolTREC2(IC)+QJ2(IC)*DTP)*(1.0-LIMd)  !HUGO FAGUNDES 16/08/2019
                    DEPTREC(IC,i) = DEPTREC(IC,i) + DEPaux !@ (TON)
                    CSJ2(IC,i) = CTS(IC,i) !@ (TON/M3)                              !Hugo Fagundes 16/08/2019
                    CSdep(IC,i)= CSdep(IC,i) + DEPaux  !@(TON)                            Hugo Fagundes 13/08/19
                ELSEIF (DMP(i) <= 0.00025) THEN     !@ limite máximo areia fina (Wu, 2008)
                    DeT = HmJ/2.5    ! Distância efetiva
                    LIMd = min(WSP(i)*DTP/DeT, 1.0)
                    DEPaux = (CSJ2(IC,i) - CTS(IC,i))*(VolTREC2(IC)+QJ2(IC)*DTP)*LIMd  !HUGO FAGUNDES 16/08/2019
                    DEPaux2(IC,i) = (CSJ2(IC,i) - CTS(IC,i))*(VolTREC2(IC)+QJ2(IC)*DTP)*(1.0-LIMd)  !HUGO FAGUNDES 16/08/2019
                    DEPTREC(IC,i) = DEPTREC(IC,i) + DEPaux  !@ (TON)
                    CSJ2(IC,i) = CTS(IC,i) !@ (TON/M3)                              !Hugo Fagundes 16/08/2019
                    CSdep(IC,i)= CSdep(IC,i) + DEPaux!@(TON)                            Hugo Fagundes 13/08/19
                ELSE
                    DeT = HmJ/11.24  ! Distância efetiva
                    LIMd = min(WSP(i)*DTP/DeT, 1.0)
                    DEPaux = (CSJ2(IC,i) - CTS(IC,i))*(VolTREC2(IC)+QJ2(IC)*DTP)*LIMd  !HUGO FAGUNDES 16/08/2019
                    DEPaux2(IC,i) = (CSJ2(IC,i) - CTS(IC,i))*(VolTREC2(IC)+QJ2(IC)*DTP)*(1.0-LIMd)  !HUGO FAGUNDES 16/08/2019
                    DEPTREC(IC,i) = DEPTREC(IC,i) + DEPaux  !@ (TON)
                    CSJ2(IC,i) = CTS(IC,i) !@ (TON/M3)                              !Hugo Fagundes 16/08/2019
                    CSdep(IC,i)= CSdep(IC,i) + DEPaux!@(TON)                        !Hugo Fagundes 13/08/19   
                ENDIF
                LIMdaux(i) = LIMd
                IF(i==1) THEN
                    HRIO1(IC)=HRIO2(IC) !hof 19.11.2020
                    HRIO2(IC)=HRIO2(IC) - (DEPaux/densiSED)/(SRIO(IC)*BRIO(IC)) !hof 19.11.2020
                ENDIF
            !@ EROSÃO INICIO ----------------------------------------------------------------------- Hugo Fagundes 12/08/19
            ELSE
            !@ limitador temporal da erosão (igual ao HEC-RAS 4.0) - igual para todas as partículas
                IF (DMP(i) <= 0.000062) THEN    !@ silte médio (Wu, 2008)
                CSJ2(IC,i) = CSJ2(IC,i)
                ELSE
                !*********Sem Limitador*********!
                !    LIMe = min( (1.0 + exp(-1.0)) - exp(-SRIO(IC)*1000./(30.*HmJ)), 1.0)
                !    !EROSaux = (CTS(IC,i) - CSJ2(IC,i))*QJ2(IC)*DTP*LIMe
                !    EROSaux = (CTS(IC,i) - CSJ2(IC,i))*(VolTREC2(IC) + QJ2(IC)*DTP)*LIMe     !@ (TON)  !HUGO FAGUNDES 16/08/2019  
	            !    CSJ2(IC,i) = CSJ2(IC,i) + EROSaux/(VolTREC2(IC)+QJ2(IC)*DTP)
	            !    EROSTREC(IC,i) = EROSTREC(IC,i) + EROSaux
	            !    !CSJ2(IC,i) = (CSJ2(IC,i)*QJ2(IC)*DTP + EROSaux)/(QJ2(IC)*DTP)
                !ENDIF
                !*********Sem Limitador********!
                
                !*********Com Limitador*********!
                    LIMe = min( (1.0 + exp(-1.0)) - exp(-SRIO(IC)*1000./(30.*HmJ)), 1.0)
                    EROSaux = (CTS(IC,i) - CSJ2(IC,i))*(VolTREC2(IC) + QJ2(IC)*DTP)*LIMe     !@ (TON)  !HUGO FAGUNDES 16/08/2019  
                    IF (EROSaux <= CSdep(IC,i)) THEN
                        CSJ2(IC,i) = CSJ2(IC,i) + EROSaux/(VolTREC2(IC)+QJ2(IC)*DTP)     !@ (TON/M3)    !HUGO FAGUNDES 16/08/2019  
                        CSdep(IC,i) = CSdep(IC,i) - EROSaux                 !@ (TON)
                        EROSTREC(IC,i) = EROSTREC(IC,i) + EROSaux           !@ (TON)             
                    ELSE
                        CSJ2(IC,i) = CSJ2(IC,i) + CSdep(IC,i)/(VolTREC2(IC)+QJ2(IC)*DTP)     !@ (TON/M3)    !HUGO FAGUNDES 16/08/2019  
                        EROSaux = CSdep(IC,i)                               !@ (TON)
                        EROSTREC(IC,i) = EROSTREC(IC,i) + EROSaux           !@ (TON)
                        CSdep(IC,i)= 0.0                                    !@ (TON/M3)
                    ENDIF  
                !*********Com Limitador*********!
                    
                ENDIF         
                LIMeaux(i) = LIMe
                IF(i==1) THEN
                    HRIO1(IC)=HRIO2(IC) !hof 19.11.2020
                    HRIO2(IC)=HRIO2(IC) + (EROSaux/densiSED)/(SRIO(IC)*BRIO(IC)) !hof 19.11.2020
                ENDIF
            !@ EROSÃO FIM ----------------------------------------------------------------------- Hugo Fagundes 12/08/19   

            ENDIF
        ELSE
            !@ Tudo vira deposição se não há vazão de saída!
            DEPaux = PSC(IC,i+1) + CSM2(IC,i)*abs(QM2(IC))*DTP! !@DCB_sed set2012
            DEPaux = DEPaux + CSJ1(IC,i)*VolTREC1(IC)   !Hugo Fagundes 10/08/2019
            
            !*********Com Limitador*********!
            CSdep(IC,i)= CSdep(IC,i) + DEPaux!@(TON)    !Hugo Fagundes 10/08/2019      
            !*********Com Limitador*********!
            
            DEPTREC(IC,i) = DEPTREC(IC,i) + DEPaux
            CSJ2(IC,i) = 0.0                                     !Hugo Fagundes 13/08/2019
        ENDIF
        !@ *******************************************************************************

        !@ *******************************************************************************
        !@ BALANÇO DE SEDIMENTO NO TRECHO DE RIO
        !@ -------------------------------------------------------------------------------
        !@ Descarga sólida na seção de jusante (ton/s)
        QSJ2(IC,i) = 0.0                                        !@DCB_sed set2012
        IF (QJ2(IC) > Qlim) QSJ2(IC,i) = CSJ2(IC,i)*QJ2(IC)    !@DCB_sed set2012
        
        JUS=CELJUS(IC)
        
        IF (JUS>0) CARGM(JUS,i) = CARGM(JUS,i) + QSJ2(IC,i)*DTP !@ (ton/dia)
        !@ *******************************************************************************

        !@ armazena carga de entrada na minibacia (ton/dia)
        iENTRAaux(i)      = CARGM(IC,i) + PSC(IC,i+1) + EROSaux 
        !@ acumula carga de entrada na minibacia (ton)
        iENTRA(iSEDaux,i) = iENTRA(iSEDaux,i) + iENTRAaux(i)
        !@ armazena carga de saída na minibacia (ton/dia)
        iSAIaux(i)        = CSJ2(IC,i)*QJ2(IC)*DTP + DEPaux 
        !@ acumula carga de saída na minibacia (ton)
        iSAI(iSEDaux,i)   = iSAI(iSEDaux,i) + iSAIaux(i)    
        
    ENDDO   !@ FIM DO LOOP DAS PARTICULAS

    !@ Salva acumulado para cada minibacia com calculo de sedimento
    DEPT(iSEDaux,:) = DEPTREC(IC,:)
    EROT(iSEDaux,:) = EROSTREC(IC,:)


!@ #########################################################################################
!@ -------------------    CELULA SEM RIO    ------------------------------------------------
!@ #########################################################################################       
ELSE
    !@ *****************************************************************************************
    !@ CONCENTRAÇÃO DE SEDIMENTOS NA SEÇÃO DE JUSANTE AO FINAL DO PASSO DE TEMPO(TON/M3)
    !@ -----------------------------------------------------------------------------------------
    DO i = 1,3
        !@ A carga da minibacia vai compor a concentração de jusante
        CSJ2(IC,i) = 0.0
        IF (QJ2(IC) > Qlim) CSJ2(IC,i) = PSC(IC,i+1)/(QJ2(IC)*DTP) + DEPaux2(IC,i)   ! ton/m3
    ENDDO

    !@ *****************************************************************************************

    !-------------------------------------------------------------------------------------------
    FracS(IC,:) = 0.0
    CTS(IC,:) = 0.0
    IF (QJ2(IC) > Qlim) THEN !@ Só tem capacidade a jusante se tem vazão
        IF (CSJ2(IC,1) > 0.0) THEN
            FracS(IC,1) = CSJ2(IC,1)/sum(CSJ2(IC,:)) !@ porcentagem da concentração de areia
            FracS(IC,2) = CSJ2(IC,2)/sum(CSJ2(IC,:)) !@ porcentagem da concentração de silte
            FracS(IC,3) = CSJ2(IC,3)/sum(CSJ2(IC,:)) !@ porcentagem da concentração de argila
        ENDIF
        !@ Capacidade de transporte
        CALL SED_CT(UatJ, VmJ, SfTREC, HmJ)
    ENDIF
    !@ *****************************************************************************************
        
    !@ *****************************************************************************************
    !@ VERIFICAÇÃO DE TENDÊNCIA A EROSÃO OU DEPOSIÇÃO
    !@ -----------------------------------------------------------------------------------------
    
    LIMdaux     = 0.0   !@ IMPRESSÃO
    LIMeaux     = 0.0   !@ IMPRESSÃO
    iENTRAaux   = 0.0   !@ IMPRESSÃO
    iSAIaux     = 0.0   !@ IMPRESSÃO

    DO i= 1,3   !@ INICIO DO LOOP DAS PARTICULAS
        DEPaux  = 0.0
        EROSaux = 0.0
        LIMd    = 0.0
        LIMe    = 0.0
        DEPaux2 (IC,i) = 0.0   !Hugo Fagundes 19/08/2019
        !DELTA   = 0.0   !Hugo Fagundes 19/08/2019
        !cargatrec = 0.0 !Hugo Fagundes 19/08/2019
        IF (QJ2(IC) > Qlim) THEN   !@ Só há passagem de sedimentos se há vazão de jusante!
            
            !@ DEPOSIÇÃO ---------------------------------------------------------------------
            IF (CSJ2(IC,i) >= CTS(IC,i)) THEN
            !@ limitador temporal da deposição (igual ao HEC-RAS 4.0) - PARA D > 0.062 mm
                IF (DMP(i) <= 0.000062) THEN        !@ limite máximo silte (Wu, 2008)
                    CSJ2(IC,i) = CSJ2(IC,i)         !@ Toda carga de finos passa para jusante
                ELSEIF (DMP(i) <= 0.000125) THEN    !@ limite máximo areia muito fina (Wu, 2008)
                    DeT = Hmj        ! Distância efetiva
                    LIMd = min(WSP(i)*DTP/DeT, 1.0)
                    DEPaux = (CSJ2(IC,i) - CTS(IC,i))*QJ2(IC)*DTP*LIMd  !@ (TON) DIOGO BUARQUE
                    DEPaux2(IC,i) = (CSJ2(IC,i) - CTS(IC,i))*QJ2(IC)*DTP*(1.0-LIMd)  !HUGO FAGUNDES 16/08/2019
                    DEPTREC(IC,i) = DEPTREC(IC,i) + DEPaux  !@ (TON)
                    CSJ2(IC,i) = CTS(IC,i) !@ (TON/M3)                              !Hugo Fagundes 16/08/2019
                    CSdep(IC,i)= CSdep(IC,i) + DEPaux!@(TON)                        !Hugo Fagundes 13/08/19   
                ELSEIF (DMP(i) <= 0.00025) THEN     !@ limite máximo areia fina (Wu, 2008)
                    DeT = HmJ/2.5    ! Distância efetiva
                    LIMd = min(WSP(i)*DTP/DeT, 1.0)
                    DEPaux = (CSJ2(IC,i) - CTS(IC,i))*QJ2(IC)*DTP*LIMd  !@ (TON) DIOGO BUARQUE
                    DEPaux2(IC,i) = (CSJ2(IC,i) - CTS(IC,i))*QJ2(IC)*DTP*(1.0-LIMd)  !HUGO FAGUNDES 16/08/2019
                    DEPTREC(IC,i) = DEPTREC(IC,i) + DEPaux  !@ (TON)
                    CSJ2(IC,i) = CTS(IC,i) !@ (TON/M3)                              !Hugo Fagundes 16/08/2019
                    CSdep(IC,i)= CSdep(IC,i) + DEPaux!@(TON)                        !Hugo Fagundes 13/08/19   
                ELSE
                    DeT = HmJ/11.24  ! Distância efetiva
                    LIMd = min(WSP(i)*DTP/DeT, 1.0)
                    DEPaux = (CSJ2(IC,i) - CTS(IC,i))*QJ2(IC)*DTP*LIMd  !@ (TON) DIOGO BUARQUE
                    DEPaux2(IC,i) = (CSJ2(IC,i) - CTS(IC,i))*QJ2(IC)*DTP*(1.0-LIMd)  !HUGO FAGUNDES 16/08/2019
                    DEPTREC(IC,i) = DEPTREC(IC,i) + DEPaux  !@ (TON)
                    CSJ2(IC,i) = CTS(IC,i) !@ (TON/M3)                              !Hugo Fagundes 16/08/2019
                    CSdep(IC,i)= CSdep(IC,i) + DEPaux!@(TON)                        !Hugo Fagundes 13/08/19   
                ENDIF
                LIMdaux(i) = LIMd
        
            !@ EROSÃO INICIO ----------------------------------------------------------------------- Hugo Fagundes 12/08/19
            ELSE
            !@ limitador temporal da erosão (igual ao HEC-RAS 4.0) - igual para todas as partículas
                IF (DMP(i) <= 0.000062) THEN    !@ silte médio (Wu, 2008)
                    CSJ2(IC,i) = CSJ2(IC,i)
                ELSE
                
                !************Sem Limitador**********
                !    LIMe = min( (1.0 + exp(-1.0)) - exp(-SRIO(IC)*1000./(30.*HmJ)), 1.0)
                !    EROSaux = (CTS(IC,i) - CSJ2(IC,i))*QJ2(IC)*DTP*LIMe
	               ! EROSTREC(IC,i) = EROSTREC(IC,i) + EROSaux
	               ! CSJ2(IC,i) = (CSJ2(IC,i)*QJ2(IC)*DTP + EROSaux)/(QJ2(IC)*DTP)
                !ENDIF
                !************Sem Limitador**********  
                
                !************Com Limitador (CSdep)**********
                    LIMe = min( (1.0 + exp(-1.0)) - exp(-SRIO(IC)*1000./(30.*HmJ)), 1.0)
                    EROSaux = (CTS(IC,i) - CSJ2(IC,i))*QJ2(IC)*DTP*LIMe     !@ (TON)
                    IF (EROSaux <= CSdep(IC,i)) THEN
                        CSJ2(IC,i) = CSJ2(IC,i) + EROSaux/(QJ2(IC)*DTP)     !@ (TON/M3)    !HUGO FAGUNDES 16/08/2019  
                        CSdep(IC,i) = CSdep(IC,i) - EROSaux                 !@ (TON)
                        EROSTREC(IC,i) = EROSTREC(IC,i) + EROSaux           !@ (TON)             
                    ELSE
                        CSJ2(IC,i) = CSJ2(IC,i) + CSdep(IC,i)/(QJ2(IC)*DTP)     !@ (TON/M3)    !HUGO FAGUNDES 16/08/2019  
                        EROSaux = CSdep(IC,i)                               !@ (TON)
                        EROSTREC(IC,i) = EROSTREC(IC,i) + EROSaux           !@ (TON)
                        CSdep(IC,i)= 0.0                                    !@ (TON/M3)
                    ENDIF  
                !************Com Limitador (CSdep)**********
                ENDIF         
                LIMeaux(i) = LIMe
            
            !@ EROSÃO FIM ----------------------------------------------------------------------- Hugo Fagundes 12/08/19   
                
            !@ EROSÃO -----------------------------------------------------------------------

            ENDIF
        ELSE
            !@ Tudo vira deposição se não há vazão de saída!
            DEPaux = PSC(IC,i+1)
            
            !************Com Limitador (CSdep)**********
            CSdep(IC,i)= CSdep(IC,i) + DEPaux!@(TON)    !Hugo Fagundes 10/08/2019    
            !************Com Limitador (CSdep)**********
            
            DEPTREC(IC,i) = DEPTREC(IC,i) + DEPaux
            CSJ2(IC,i) = 0.0                                     !Hugo Fagundes 13/08/2019
        ENDIF
        !@ *******************************************************************************    
        
        !@ *******************************************************************************
        !@ BALANÇO DE SEDIMENTO NO TRECHO DE RIO
        !@ -------------------------------------------------------------------------------
        !@ Descarga sólida na seção de jusante (ton/s)
        QSJ2(IC,i) = 0.0                                        !@DCB_sed set2012
        IF (QJ2(IC) > Qlim) QSJ2(IC,i) = CSJ2(IC,i)*QJ2(IC)    !@DCB_sed set2012
        
        JUS=CELJUS(IC)
        
        IF (JUS>0) CARGM(JUS,i) = CARGM(JUS,i) + QSJ2(IC,i)*DTP !@ (ton/dia)
        !@ *******************************************************************************
        
        !@ armazena carga de entrada na minibacia (ton/dia)
        iENTRAaux(i)      = PSC(IC,i+1) + EROSaux
        !@ acumula carga de entrada na minibacia (ton)
        iENTRA(iSEDaux,i) = iENTRA(iSEDaux,i) + iENTRAaux(i)
        !@ armazena carga de saída na minibacia (ton/dia)
        iSAIaux(i)        = CSJ2(IC,i)*QJ2(IC)*DTP + DEPaux
        !@ acumula carga de saída na minibacia (ton)
        iSAI(iSEDaux,i)   = iSAI(iSEDaux,i) + iSAIaux(i)
        
    ENDDO   !@ FIM DO LOOP DAS PARTICULAS
		    
    !@ Salva acumulado para cada minibacia com calculo de sedimento
    DEPT(iSEDaux,:) = DEPTREC(IC,:)
    EROT(iSEDaux,:) = EROSTREC(IC,:)

ENDIF
       
RETURN
END SUBROUTINE