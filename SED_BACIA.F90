    !*********************************************************************************
    !
    !  SUBROUTINE SED_BACIA estimates sediment yield by MUSLE equation for each catchment
    !   
    !---------------------------------------------------------------------------------
    !  Discussion:
    ! 
    !    This subroutine estimates sediment yield by MUSLE equation for each catchment
	!
    !    - MUSLE Equation (Williams, 1975): Sed = α*[(Qsup*qpico*A)^(β)]*K*C*P*LS*Rgros
    !    where:
    !    Sed [t/day] is the sediment load resulting from soil erosion, 
    !    Q_sup [mm/ha] is the surface runoff volume, 
    !    q_pico [m³/s] is the peak flow rate, 
    !    A is the surface area [km²], 
    !    α and β are the adjustment coefficients (which can be calibrated ), whose values originally estimated by Williams (1975) were 11.8 and 0.56, Rspectively, 
    !    K [0,013.t.m².h./m³.t.cm] is the soil erodibility factor, 
    !    C [-] is the cover and management factor, 
    !    P [-] is the conservation practice factor, 
    !    LS [-] is  the topographic factor,
    !    Rgros [-] exp(-0.053*Mrocha) - Rock factor
    !    Mrocha = Rock percentage in the first soil layers
    !
	!	 SED_BACIA is called inside MODELO.
	!
	!	 Saves global variable: 
	!
    !
    !  	Usage:
    !
    !    * CALDAT subroutine is called in this subroutine
    !
    !    where
    !
    !    * IDINI + IT - 1,IMES,IDIA,IANO arguments are passed in this subroutine
    !
    !    uses modules and functions
    !
    !    * module     VARS_MAIN   in      VARS_MAIN.f90
    !    * module     SED_VARS    in      SED_VARS.f90
    !
    !	 opens and Reads
    !
    !    no files are opened or read
    !
    !    creates
    !
    !    * Sediment yield estimated by MUSLE is write in SED_MINI.txt, SED_MINI_areia.txt, SED_MINI_argila.txt, and SED_MINI_silte.txt
    !    
    !   This Subroutine was created by Diogo Costa Buarque at 2010, September.
    !   
    !---------------------------------------------------------------------------------
    !  Licensing:
    !
    !    This code is distributed under the ...
    !
    !  Version/Modified: 
    !
    !    DIOGO COSTA BUARQUE, 09/2010
    !
    !  Authors:
    !
    !    Original fortran version by Walter Collischonn
    !    Present fortran version by:
    !    * Walter Collischonn
    !    * Rodrigo Cauduro Dias de Paiva
    !    * Diogo da Costa Buarque
    !    * Paulo Pontes Rogenes
    !    * Mino  Viana Sorribas
    !    * Fernando Mainardi Fan
    !    * Juan Martin Bravo 
    !    * Hugo de Oliveira Fagundes
    !
    !  Main Reference:
    !
    !    Walter Collischonn,
    !    Modelo de Grandes Bacias - Thesis
    !    Porto Alegre, 2001
    !    ISBN: XXXXXXXXXXX,
    !
    !---------------------------------------------------------------------------------
    !  Variables and Parameters:
    !  
    !
    ! TKSed     = retardos do escoamento superficial (s)
    ! DSUP    = lâmina do escoamento superficial (mm)
    ! VSUP    = volume do escoamento supeficial na minibacia (m3)
    ! sbtFLAG = flag para indicar minibacia com substituição de dados
    ! IB,IBAC = código da subbacia da minibacia
    ! IU,NU   = usos das minibacias
    ! PUSO    = porcentagem de uso em cada minibacia
    ! SLX     = carga e sedimentos remanescente na minibacia no passo de tempo
    ! SLXU    = carga e sedimentos remanescente por bloco da minibacia no passo de tempo
    ! PSCaux  = aporte das frações de sedimentos que chegam à drenagem das minibacias das subbacias selecionadas
    ! QSUPaux = vazão superficial que chega à drenagem apenas das minibacias das subbacias selecionadas 2013_02_06 (antes QCEaux)
    ! QCELaux = vazão que chega à drenagem apenas das minibacias das subbacias selecionadas 2013_02_06
    ! FRAC    = porcentagem de sedimento em cada bloco da minibacia
    
    ! SLC     = Total soil loss from the catchment (ton)
    ! SLU     = Total soil loss from each HRU of catchment (ton)
    ! QSAUX   = Storage, for surface runoff of each catchment, the delay time (s), total volume (m³) and blades for each HRU
    ! NPXU    = Pixels number for each HRU

    !---------------------------------------------------------------------------------	
    

SUBROUTINE SED_BACIA

USE SED_VARS
USE VARS_MAIN
USE VARS_CALIB !Hugo Fagundes 03/10/2019

IMPLICIT NONE

!@ RETARDOS DO ESCOAMENTO SUPERFICIAL (s)
REAL TKSed
!@ LAMINAS DO ESCOAMENTO SUPERFICIAL (mm)
REAL DSUP

REAL(4) QSUPaux(NC),FRAC(NU),QPICaux(NC)
REAL(4) QCELaux(NC) !@ 2013_02_06
INTEGER iaux, iUSO,i
    
    !ABAIXO ESTÃO DEFINIDOS OS VALORES DE ALFA E BETA PARA CADA SUBBACIA (5) RELACIONADA AS SUBBACIAS DA PARTE DO MODELO HIDROLOGICO (17).
    !ESSA ASSOCIAÇÃO DEVE SER FEITA PARA CADA CASO DE SIMULAÇÃO
IF (icalib==0) THEN     !HUGO FAGUNDES 21/09/17
    
!----------------------------MUSLE PARAMETERS--------------------    
ALFsed(1:NB)=11.8
BETsed(1:NB)=0.56
AuxTKS(1:NB)=1.0


ENDIF                   !HUGO FAGUNDES 21/09/17

    ALLOCATE(PSCaux(NC,4))
    ALLOCATE(DSUPa(NC)) !Hugo Fagundes 19/09/2019


PSC = 0.0
PSCaux = 0.0
DSUPa = 0.0

        !@ QSAUX armazena
            !@      TKSed da minibacia IC
            !@	    VSUP da minibacia IC
            !@	    DSUP bloco 1 da minibacia IC
            !@	    DSUP bloco 2 da minibacia IC
            !@	    DSUP bloco 3 da minibacia IC
            !@	    	*
            !@	    	*
            !@	    	*
            !@	    DSUP bloco N da minibacia IC
!----------------------
iaux = 0    !@ contador de minibacia com cálculo de sedimentos
!----------------------


!VERIFICA A DATA CORRESPONDE O DIA JULIANO
CALL CALDAT(IDINI + IT - 1,IMES,IDIA,IANO)

QPICaux = 0.0
DO IC=1,NC	!@ LOOP DAS MINIBACIAS
	
	IF (sbtFLAG(iC)==1) CYCLE !@ Nao computa minibacias a montante de minibacias com substituicao de dados
	!@ Neste caso, para sedimentos, ainda terá que ser fornecida uma concentração de sedimentos

    IB=IBAC(IC)
    
!@ DCB 30/04/1011 ###########################################################################
	IF(IB>SUBfim.OR.IB<SUBini)CYCLE

!@ DCB 30/04/1011 ###########################################################################

	iaux = iaux + 1 !@ contador de minibacia com cálculo de sedimentos
	
	DO IU=1,NU	!@ LOOP DOS USOS

		APIX = 0.0
		FCTE = 0.0
		QPIC = 0.0

		IF (NPXU(IC,IU) .LE. 0.0) THEN  !@ NAO TEM PIXEL(S) COM ESTE USO NESTA MINIBACIA
			CYCLE !@ PASSA PARA O PROXIMO USO
		ENDIF

		IF(PUSO(IC,IU).LT.0.0001)THEN !@ NAO TEM ESTE USO NESTA MINIBACIA
			CYCLE !@ PASSA PARA O PROXIMO USO
		ENDIF
		

		DSUP = QSAUX(IC,IU+2)
        IF (ICALIB==3) THEN
            DSUP_SC(IC,IT,IU)=DSUP !Hugo Fagundes 03/10/2019
        ENDIF
        
        
        !!!DSUPa(IC)=DSUPa(IC)+DSUP/NU    !Calcula e armazena o DSUP médio para cada minibacia HUGO FAGUNDES 19/09/2019
		
		!@ ÁREA MÉDIA DOS PIXELS DE CADA USO DA MINIBACIA (KM^2)
		APIX = (ACEL(IC)*(PUSO(IC,IU)/100.))/NPXU(IC,IU) !@ (KM^2)
		APIX = APIX*100. !@ (Ha = hectares)


		!@ FATOR CONSTANTE DA MUSLE PARA CADA USO DA MINIBACIA
		FCTE = ALFsed(IB)*(APIX**BETsed(IB))*Kusle(IB,IU)*Cusle(IB,IU)*Pusle(IB,IU)*Rgros(IB,IU)*LSAcu(IC,IU)   ! HUGO FAGUNDES 04/10/17
		APIX = APIX/100. !@ (KM^2)


		IF (DSUP < 0.0) DSUP = 0.0 !@ CORREÇÃO DE DSUP NEGATIVO VINDO DA ROTINA "CELULA"

		!@ TAXA DE PICO DO ESCOAMENTO SUPERFICIAL PELA MÉTODO RACIONAL, CONSIDERANDO
        !@ A INTENSIDADE DA PRECIPITAÇÃO DENTRO DAS 24 HORAS (DIA) (M^3/S)
        IF (P(IC) > 0.0) QPIC = ((DSUP/P(IC))*(P(IC)/24.)*APIX)/36.0 !@ (M^3/S) !HOF 18.11.2020

		QPICaux(iaux) = QPICaux(iaux) + QPIC !@ acumula taxa de pico, apenas para gravação

		!@ PERDA DE SOLO PARA CADA USO DA MINIBACIA (TON)

		SLU(IC,IU) = SLU(IC,IU) + FCTE*(DSUP**BETsed(IB))*(QPIC**BETsed(IB)) !@ (TON)   HUGO FAGUNDES 04/10/17
        
        APIX = APIX*100. !@ (Ha = hectares)
        DSUPa(IC)=DSUPa(IC)+(ALFsed(IB)*(APIX**BETsed(IB))*(DSUP**BETsed(IB))*(QPIC**BETsed(IB)))/NU    !Calcula e armazena o Produto alfa(Q.qpico.A)^beta médio para cada minibacia HUGO FAGUNDES 23/09/2019
	ENDDO	! FIM DO LOOP DOS USOS


!    !@ MULTIPLICAR O SLU PELO SDR, ANTES DE FAZER A PROPAGAÇÃO NA REDE. O SDR AQUI É POR MINIBACIA, ENTÃO
!    !@ É CONSTANTE NOS BLOCOS. COM A OPÇÃO DE VARIAR POR BLOCO É SÓ MULTIPLICAR O SDRi DO BLODO i PELO SLUi DO
!    !@ BLOCO.
!    DO IU = 1,NU
!        !@ O SDR com o parametro de calibração não pode ser maior que 1!
!        SLU(IC,IU) = SLU(IC,IU)*min(Ksdr(IB,IU)*SDR(IC,NU+1),1.)
!    ENDDO


    FRAC = 0.0  !@ porcentagem de sedimento de cada bloco da minibacia
    if (sum(SLU(IC,:)) /= 0.0) then
        FRAC(:) = SLU(IC,:)/sum(SLU(IC,:))
    endif


	!@ PARÂMETRO DE RETARDO DO ESCOAMENTO SUPERFICIAL
	TKSed = AuxTKS(IB)*QSAUX(IC,1)! HUGO FAGUNDES 11/07/19
    IF (ICALIB==3) THEN
        TKS_SC(IC)=QSAUX(IC,1) !Hugo Fagundes 03/10/2019
    ENDIF
    
	!@ DESCARGA SÓLIDA DAS MINIBACIAS (TON/S) NO TEMPO IT
	!QSC(IC)   = SLC(IC)/TKSed
	QSU(IC,:) = SLU(IC,:)/TKSed	!@ 21/02/2011

	!@ ATUALIZA PERDA DE SOLO DAS MINIBACIAS (TON) NO TEMPO IT
	!SLX     = SLC(IC)   - QSC(IC)*DTP !@ verifica o que sobra de sedimentos na minibacia
	SLXU(:) = SLU(IC,:) - QSU(IC,:)*DTP !@ 21/02/2011


!@ **********************************************************************
!@ *** 21/02/2011 *******************************************************
!@ **********************************************************************
	DO IU=1,NU
		IF(SLXU(IU) < 0.0)THEN				!@ TKSed < DTP
			QSU(IC,IU) = SLU(IC,IU)/(DTP)
			SLU(IC,IU) = 0.0
		ELSE
			SLU(IC,IU) = SLXU(IU)
		ENDIF
	ENDDO
!@ **********************************************************************
!@ *** 21/02/2011 *******************************************************
!@ **********************************************************************



    !@ APORTE DE SEDIMENTOS DAS BLOCOS DAS MINIBACIAS PARA O RIO (TON)
	PSU(IC,:) = QSU(IC,:)*DTP	!@ 21/02/2011



	!@ APORTE DE SEDIMENTOS DAS MINIBACIAS PARA O RIO (TON)
	PSC(IC,1) = sum(QSU(IC,:))*DTP
	DO iUSO=1,NU-1
	    PSC(IC,2) = PSC(IC,2) + PSC(IC,1)*FRAC(iUSO)*Mareia(IB,iUSO)/100	!@ APORTE DE AREIA QUE CHEGA À DRENAGEM
	    PSC(IC,3) = PSC(IC,3) + PSC(IC,1)*FRAC(iUSO)*Msilte(IB,iUSO)/100	!@ APORTE DE SILTE QUE CHEGA À DRENAGEM
	    PSC(IC,4) = PSC(IC,4) + PSC(IC,1)*FRAC(iUSO)*Margila(IB,iUSO)/100	!@ APORTE DE ARGILA QUE CHEGA À DRENAGEM
    ENDDO   
    
    !PSC(IC,2) = PSC(IC,2)*DEFINIRCOEFICIENTE(IB)	!@ APORTE DE AREIA AJUSTADO HOF 01/04/2019
    
    !PSC(IC,3) = PSC(IC,3)*1.5
    !PSC(IC,4) = PSC(IC,4)*1.5

!@ VERIFICA O QUE OCORRE SE TUDO POR CONSIDERADO EM SUSPENSÃO
!PSC(IC,3) = PSC(IC,3) + PSC(IC,2)/2.
!PSC(IC,4) = PSC(IC,4) + PSC(IC,2)/2.
!PSC(IC,2) = 0.
!@ *************************************
	!@ ORDENANDO MINIBACIAS DAS SUBBACIAS DE INTERESSE
	PSCaux(iaux,1) = PSC(IC,1)	!@ ARMAZENA O APORTE DE SEDIMENTOS QUE CHEGA À DRENAGEM APENAS NAS MINIBACIAS DAS SUBBACIAS SELECIONADAS
	QSUPaux(iaux)   = QSUP(IC)	!@ ARMAZENA A VAZÃO SUPERFICIAL DA MINIBACIA QUE CHEGA À DRENAGEM APENAS NAS MINIBACIAS DAS SUBBACIAS SELECIONADAS
	QCELaux(iaux)   = QCEL2(IC)	!@ ARMAZENA A VAZÃO DA MINIBACIA QUE CHEGA À DRENAGEM APENAS NAS MINIBACIAS DAS SUBBACIAS SELECIONADAS 2013_02_06
    
	DO iUSO=1,NU-1
	    PSCaux(iaux,2) = PSCaux(iaux,2) + PSCaux(iaux,1)*FRAC(iUSO)*Mareia(IB,iUSO)/100     !@ ARMAZENA O APORTE DE AREIA QUE CHEGA À DRENAGEM APENAS NAS MINIBACIAS DAS SUBBACIAS SELECIONADAS
	    PSCaux(iaux,3) = PSCaux(iaux,3) + PSCaux(iaux,1)*FRAC(iUSO)*Msilte(IB,iUSO)/100     !@ ARMAZENA O APORTE DE SILTE QUE CHEGA À DRENAGEM APENAS NAS MINIBACIAS DAS SUBBACIAS SELECIONADAS
	    PSCaux(iaux,4) = PSCaux(iaux,4) + PSCaux(iaux,1)*FRAC(iUSO)*Margila(IB,iUSO)/100    !@ ARMAZENA O APORTE DE ARGILA QUE CHEGA À DRENAGEM APENAS NAS MINIBACIAS DAS SUBBACIAS SELECIONADAS
    ENDDO
!@ *************************************

ENDDO	!@ FIM DO LOOP DAS MINIBACIAS

!@ GRAVA CARGAS DE SEDIMENTO DAS MINIBACIAS QUE CHEGAM A REDE DE DRENAGEM
IF(ICALIB==0)THEN 
    WRITE(FILSED, '(4I10,<iaux>F15.3)') IT, IDIA, IMES, IANO, (PSCaux(IC,1),IC=1,iaux)  !@ SED_MINI.txt - total (TON/dia)
    WRITE(FILSEDa,'(4I10,<iaux>F15.3)') IT, IDIA, IMES, IANO, (PSCaux(IC,2),IC=1,iaux)  !@ SEDa_MINI.txt - silte (TON/dia)
    WRITE(FILSEDs,'(4I10,<iaux>F15.3)') IT, IDIA, IMES, IANO, (PSCaux(IC,3),IC=1,iaux)  !@ SEDs_MINI.txt - argila (TON/dia)
    WRITE(FILSEDc,'(4I10,<iaux>F15.3)') IT, IDIA, IMES, IANO, (PSCaux(IC,4),IC=1,iaux)  !@ SEDc_MINI.txt - areia (TON/dia)

    !@ GRAVA A TAXA DE PICO ACUMULADA NAS MINIBACIAS
    !WRITE(TAXPIC,'(4I10,<iaux>F15.3)')  IT, IDIA, IMES, IANO, (QPICaux(IC),IC=1,iaux)  !@ TAXpico_MINI.txt - acumulado (m3/s)

    !@ GRAVA VAZÃO SUPERFICIAL DAS MINIBACIAS QUE CHEGA A REDE DE DRENAGEM
    WRITE(FILQCE,'(4I10,<iaux>F15.3)') IT, IDIA, IMES, IANO, (QSUPaux(IC),IC=1,iaux) !@ QSUP_MINI.txt (m3/s)

    !@ GRAVA VAZÃO DAS MINIBACIAS QUE CHEGA A REDE DE DRENAGEM 2013_02_06
    WRITE(FILCEL,'(4I10,<iaux>F15.3)') IT, IDIA, IMES, IANO, (QCELaux(IC),IC=1,iaux) !@ QCEL_MINI.txt (m3/s) 2013_02_06
ENDIF

DEALLOCATE(PSCaux)
DEALLOCATE(DSUPa) !Hugo Fagundes 19/09/2019


RETURN

END