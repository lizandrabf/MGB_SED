!@ *********************************************************************************************
!@ Criado por Diogo Buarque
!@ Data: Abril de 2011
!@
!@ Atualizado: Abr 2011
!@
!@ >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> SUBROTINA SED_REDE >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!@
!@ - SUBROTINA PARA PREPARAÇÃO DAS VARIÁVEIS DO ESCOAMENTO E SEDIMENTOS PARA PROPAGAÇÃO DA CERGA
!@   DE SEDIMENTOS AO LONGO DA DEDE DE DRENAGEM
!@
!@ *********************************************************************************************

SUBROUTINE SED_REDE

!@ **********************************
USE VARS_MAIN
USE SED_VARS
USE VARS_INERC
USE VARS_CALIB !Hugo Fagundes 03/10/2019
!USE PAR1_MOD    !@ DCB_HD_Sed
!USE TIME_MOD    !@ DCB_HD_Sed
!USE PP_MOD      !@ DCB_HD_Sed
!@ **********************************

IMPLICIT NONE

!@ **********************************
!@ Declaração de variáveis:
REAL HmJx, VmJx, RHmJx, SfTRECx, UatJx, Qlimx
REAL(8) iENTRAaux(3), iSAIaux(3)
REAL, ALLOCATABLE:: DEPaux2(:,:)
INTEGER SEDTR, i, iX1, iX2      !@ DCB_HD_Sed
REAL QtFL, AtFL, Vplan, Cotaplan
REAL*8 H1             !@ DCB_HD_Sed
REAL*8 FINT
REAL HFLSD                  !@ DCB_HD_Sed

!ALLOCATE(CSJ2aux(nSEDmini,3))
REAL(4) CSJ2aux(nSEDmini,3)
REAL VTOTALaux(nSEDmini,2)  !@ DCB_HD_Sed
!@ **********************************
ALLOCATE (DEPaux2(NC,3))

!@ *****************************************************************************************
!@ EM CADA INTERVALO DE TEMPO O QUE ERA i+1 VIRA i
!@ -----------------------------------------------------------------------------------------
DO IC=1,NC
	CSM1(IC,:)      = CSM2(IC,:)
	CSJ1(IC,:)      = CSJ2(IC,:)
	VolTREC1(IC)    = VolTREC2(IC)
	CFL1(IC,:)      = CFL2(IC,:)
ENDDO
VFL1 = VFL2 !@ DCB_HD_Sed
VFL2 = 0.0  !@ DCB_HD_Sed
!@ *****************************************************************************************

!@ INICIALIZANDO VARIÁVEIS
iSEDaux = 0     !@ Contador de minibacia com cálculo de sedimentos
CARGM   = 0.0   !@ Carga a montante do trecho
CSJ2    = 0.0
Qlimx   = 0.01
CFL2    = 0.0
VTOTAL  = 0.0       !@ DCB_HD_Sed
DEPaux2 = 0.0 !Hugo Fagundes 19/08/2019
!QTROCaux    = 0.0   !@ DCB_HD_Sed
!QTROC       = 0.0   !@ DCB_HD_Sed

!    ******************** Hugo Fagundes 31.07.19 ******************** 
!!Vinicius -> Inseri esta parte para considerar as vazões de montante, já que não existe QM1, QM2 e QJ1 no inercial             
If(hdFLAG0>0)then
    DO IC=1,NC
        !Update both upstream and downstream flow at previous time step 
        QM1(IC)=QM2(IC)        
        
        !Number of upstream catchments of IC
        Nentradas = MINIMONT(iC,1)            
        !Sum of upstream IC flows
        if(Nentradas==0)then
            QM2(iC)=0.0
        else
            QM2(iC) = SUM(Q2fl(MINIMONT(iC,2:1+Nentradas)))
        endif
    ENDDO
ENDIF
!    ******************** Hugo Fagundes 31.07.19 ******************** 
DO IC = 1,NC    !@ INICIO DO LOOP DAS MINIBACIAS
    IB=IBAC(IC)
    
    !@ CALCULA APENAS SUB-BACIAS SELECIONADAS
	IF(IB>SUBfim .OR. IB<SUBini)CYCLE 
	!IF(IB>82 .OR. IB<57)CYCLE   !MADEIRA LBF 02/2019

	iSEDaux = iSEDaux + 1 !@ Contador de minibacia com cálculo de sedimentos
    
    QtFL = 0.0  !@ Vazão de troca total na minibacia
    AtFL = 0.0  !@ Área total da planície de inundação na minibacia
    HFLSD  = 0.0  !@ Profundidade média da planície
    !H1   = 0.0  !@ DCB_HD_Sed Profundidade média no trecho de rio

    !@ #########################################################################################
    ! ------------------------------------------------------------------------------------------
    !@ ##########################      TRECHOS SIMULADOS COM MC       ##########################
    ! ------------------------------------------------------------------------------------------
    !@ #########################################################################################
	If(hdFLAG0==0)then

        !@ *****************************************************************************************
	    !@ CARACTERÍSTICAS HIDRÁULICAS DA SEÇÃO DE JUSANTE
	    !@ -----------------------------------------------------------------------------------------
        !@ Declividade de atrito no trecho (mínima de 0.01 m/km para evitar problemas com o
        !@ calculo da CT por YANG, pois valores menores geram CT=0 para vazões maiores
        !@ que 5000 m3/s)
	    SfTRECx = max(DECL(IC),0.00001) 
    	
        !@ Profundidade média na seção de jusante por Manning assumindo Raio = HmTREC (m)
        HmJx = ((RUGMAN(IC)*QJ2(IC))/(BRIO(IC)*SfTRECx**0.5))**(3.0/5.0)

        !@ Volume médio de água no trecho
        VolTREC2(IC) = HmJx*BRIO(IC)*SRIO(IC)*1000.
        VTOTAL(IC,1) = VolTREC2(IC) !@ DCB_HD_Sed

        !@ Velocidade média na seção de jusante (m/s)
        VmJx = 0.0
        IF (QJ2(IC) > Qlimx) VmJx = QJ2(IC)/(BRIO(IC)*HmJx)

        !@ Raio Hidráulico médio da seção de jusante
        RHmJx = HmJx

        !@ Velocidade de atrito na seção de jusante (m^2/s)
        UatJx = sqrt(9.81*RHmJx*SfTRECx)
	    !@ *****************************************************************************************
        
        !@ QtFL, VFL1(IC), VFL2(IC), HFL são nulos na propagação por MC, pois não há planície
        QtFL        = 0.0
        VFL1(IC)    = 0.0
        VFL2(IC)    = 0.0
        HFLSD       = 0.0
        CALL SED_PROPAG(HmJx, VmJx, SfTRECx, UatJx, Qlimx, QtFL, VFL1(IC), VFL2(IC), HFLSD, DEPaux2)

    !@ ##################################################################################################################
    ! -------------------------------------------------------------------------------------------------------------------
    !@ ##########################     TRECHOS SIMULADOS COM INERCIAL  ###################################################
    ! -------------------------------------------------------------------------------------------------------------------
    !@ ##################################################################################################################
	ELSE

        !@ *****************************************************************************************
	    !@ CARACTERÍSTICAS HIDRÁULICAS DA SEÇÃO DE JUSANTE
	    !@ -----------------------------------------------------------------------------------------
        !@ Declividade de atrito no trecho (mínima de 0.01 m/km para evitar problemas com o
        !@ calculo da CT por YANG, pois valores menores geram CT=0 para vazões maiores
        !@ que 5000 m3/s)

    !    ******************** Hugo Fagundes 31.07.19 ********************	    
        !@ Profundidade d'água na seção de jusante

        
         ! Bottom level and water level of IC catchment:
        z1=ZTAB(1,iC)
        y1=Hfl(iC)+z1                
        
        !Calcula largura média de fluxo
        iCJus = CELJUS(iC)
          if(iCJus /= -1)then            
                bflow=DBLE(BRIO(iC)) + DBLE(BRIO(iCJus))
                bflow=bflow/2.                                                 
                !Computes downstream water elevation
                z2=ZTAB(1,iCJus)
                y2=Hfl(iCJus)+z2
       
        endif
        
            ! Calculates the hflow variable:
            hflow=max(y2,y1)-max(z2,z1)
            hflow=max(hflow,0.0)

        !Vazão por unidade de largura
        q0=Q2fl(iC)/bflow ! in m2/s    
        
        !River manning coefficient
        xMan=nMan(iC)
        
        HmJx = hflow !Hfl(iC) ---------->mesma coisa

        !@ Declividade de atrito no trecho (mínima de 0.01 m/km para evitar problemas com o
        !@ calculo da CT por YANG, pois valores menores geram CT=0 para vazões maiores
        !@ que 5000 m3/s)
        !@ Declividade da linha de energia na seção de jusante !LBF 06/2019	    
        !SfTRECx = max(((xMan*xMan*abs(q0)*q0)/(hflow**(10.0/3.0))),0.00001)       
        SfTRECx = max(((xMan*xMan*abs(q0)*q0)/(HmJx**(10.0/3.0))),0.00001)       
!    ******************** Hugo Fagundes 31.07.19 ********************	 


        !@ Volume médio de água no trecho
        VolTREC2(IC) = HmJx*BRIO(iC)*SRIO(IC)*1000.
        VTOTAL(IC,1) = VolTREC2(IC) !@ DCB_HD_Sed

        !@ Velocidade média na seção de jusante (m/s)
        VmJx = 0.0
        IF (QJ2(IC) > Qlimx) VmJx = QJ2(IC)/(BRIO(IC)*HmJx)

        !@ Raio Hidráulico médio da seção de jusante
        !RHmJx = 2.0*HmJx + BRIO(IC)
        !RHmJx = HmJx     !  Hugo Fagundes 31.07.19
        RHmJx = (HmJx*BRIO(IC))/(2.0*HmJx + BRIO(IC))!  Hugo Fagundes 31.07.19!

        !@ Velocidade de atrito na seção de jusante (m^2/s)
        UatJx = sqrt(9.81*RHmJx*SfTRECx)
    
		!@ Área da planície de inundação:
		!AtFL = (Area2(iC)*1000000.)-(BRIO(iC)*SRIO(IC)*1000.)  
		AtFL = max(0.,(Area2(iC)*1000000.)-(BRIO(iC)*SRIO(IC)*1000.))  !Vinicius
        
        !@ Calculo do volume na planicie:
        VFL2(IC) = max(0.,Vol2(ic)-VTAB(3,iC)) !LBF
        !VFL2(IC) = max(0.,Vol2(ic)-VTAB(2,iC)) !Vinicius

        !@ Vazão de troca:
        QtFL = (VFL2(iC)-VFL1(iC))/DTP 
        
        !@ Volume de água na planícies
!        VFL2(IC) = max(VFL1(IC) + QtFL*DTP,0.0)    !@ deixar volume igual ao balanço que entra e sai
        VTOTAL(IC,2) = VFL2(IC) !@ DCB_HD_Sed
       
        !@ Cálculo da profundidade média da planície
        HFLSD = VFL2(IC)/AtFL 
        
        IF (ICALIB==3) THEN
            HmJx_SC(IC,IT)=HmJx     !Hugo Fagundes 03/10/2019
            SfTRECx_SC(IC,IT)=SfTRECx     !Hugo Fagundes 03/10/2019
            VOLTREC2_SC(IC,IT)=VOLTREC2(IC)     !Hugo Fagundes 03/10/2019
            VFL2_SC(IC,IT)=VFL2(IC)            !Hugo Fagundes 03/10/2019
            HFLSD_SC(IC,IT)=HFLSD            !Hugo Fagundes 03/10/2019
        ENDIF
     !######################################################################################################
       
!        IF (IC == 6717) WRITE(*,*) 'IC, IT, VFL1, VFL2 2 = ', IC, IT, VFL1(IC), VFL2(IC)
!        IF (IC == 6717) WRITE(*,*) 'HFL, WSsilt, WSarg   = ', HFL, WSP(2), WSP(3)
!        IF (IC == 6717) PAUSE

!       !@ Para desconsiderar planícies, comentar as próximas 4 linhas - em relação aos sedimentos.
!        QtFL        = 0.0
!        VFL1(IC)    = 0.0
!        VFL2(IC)    = 0.0
!        HFLSD       = 0.0
         CALL SED_PROPAG(HmJx, VmJx, SfTRECx, UatJx, Qlimx, QtFL, VFL1(IC), VFL2(IC), HFLSD, DEPaux2)
        
	ENDIF
    !@ #########################################################################################
    

    !@ *************************************
	!@ ORDENANDO MINIBACIAS DAS SUBBACIAS DE INTERESSE
    CSJ2aux(iSEDaux,1) = CSJ2(IC,1) 
    CSJ2aux(iSEDaux,2) = CSJ2(IC,2)
    CSJ2aux(iSEDaux,3) = CSJ2(IC,3)

!    QTROCaux(iSEDaux)  = QTROC(IC)
    VTOTALaux(iSEDaux,1) = VTOTAL(IC,1)
    VTOTALaux(iSEDaux,2) = VTOTAL(IC,2)
!    VTOTALaux(iSEDaux,1) = CFL2(IC,2)*(10.**(6.)) + CFL2(IC,3)*(10.**(6.))
    !@ *************************************

ENDDO   !@ FIM DO LOOP DAS MINIBACIAS

!    ******************** Hugo Fagundes 31.07.19 ********************
!Vinicius
If(hdFLAG0>0)then
    DO IC=1,NC
        !Update downstream flow at previous time step 
        QJ1(IC)=QJ2(IC)        
    ENDDO
ENDIF
!    ******************** Hugo Fagundes 31.07.19 ********************
IF(ICALIB==0)THEN
    !@ GRAVA CARGAS DE SEDIMENTO DAS MINIBACIAS QUE CHEGAM A REDE DE DRENAGEM
    WRITE(RIOCAR,'(3I10,<iSEDaux>F15.3)')  IDIA, IMES, IANO, (CSJ2aux(IC,1)*(10.**(6.)),IC=1,iSEDaux)   !@ CONC_RIO_areia.txt  (mg/L)
    WRITE(RIOCSI,'(3I10,<iSEDaux>F15.3)')  IDIA, IMES, IANO, (CSJ2aux(IC,2)*(10.**(6.)),IC=1,iSEDaux)   !@ CONC_RIO_silte.txt  (mg/L)
    WRITE(RIOCCL,'(3I10,<iSEDaux>F15.3)')  IDIA, IMES, IANO, (CSJ2aux(IC,3)*(10.**(6.)),IC=1,iSEDaux)   !@ CONC_RIO_argila.txt (mg/L)
    WRITE(RIODEP,'(3I10,<iSEDaux>F15.3)')  IDIA, IMES, IANO, (CSdep(IC,1),IC=1,iSEDaux)   !@ DEPdin_RIO_.txt (ton/dia)        ! Hugo Fagundes 10/08/17
    !WRITE(RIOERO,'(3I10,<iSEDaux>F15.3)')  IDIA, IMES, IANO, ( --- ,IC=1,iSEDaux)  !@ DEP_RIO.txt (ton/dia)         ! Hugo Fagundes 10/08/17
    !WRITE(VAZFLP,'(3I10,<iSEDaux>F15.3)')  IDIA, IMES, IANO, (QTROCaux(IC),IC=1,iSEDaux)  !@ DCB_HD_Sed
!    WRITE(PROFUNDIDADE,'(3I10,<iSEDaux>F17.3)')   IDIA, IMES, IANO, (Hfl(iC),IC=1,iSEDaux) !LBF 01/09/2019
!    WRITE(NIVELEL,'(3I10,<iSEDaux>F17.3)')   IDIA, IMES, IANO, (Yfl(iC),IC=1,iSEDaux) !LBF 01/09/2019
!    DEALLOCATE(CSJ2aux)


!     IF(IC==4074) THEN
!       WRITE(*,*) 'CHECK', IT, CSJ2aux(IC,2)  
   ENDIF 
   
DEALLOCATE (DEPaux2)

RETURN
END SUBROUTINE