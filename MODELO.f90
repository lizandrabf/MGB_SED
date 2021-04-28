    !*********************************************************************************
    !
    !  SUBROUTINE MODELO controls the main time loop in MGB-IPH and call routines
	!					that realize catchment and river routing
    !
    !---------------------------------------------------------------------------------
    !  Discussion:
    ! 
    !    This routine has the main loop of the MGB-IPH, the time loop, from iT=1,nT
    !     where iT is the time interval and nT is the number of time steps.
	!
	!	 For each time interval date info is set, then rainfall and climate data are
	!	   loaded through subroutines LECHUVA and LECLIMA. Afterwards catchment flow
	!	   generation is done using subroutine CELULA. Finally, river routing is 
	!	   achieved by calling (i) REDE, for Muskingum-Cunge Method or
	!	   (ii) flood_inertial, for a 1D Hydrodynamic Model (w/ inertia and pressure)
	!	
	!	At the end discharge time series are stored in :
	!		QRB: calculated discharge time series in all subbasins
	!		QRG: calculated discharge time series in catchments pre-defined by user
	!		QR:  calculated discharge time series in subbasin outlet w observed data
	!	Those are recorded in files at the end when returns to SIMULA
	!
	!
	!	 * iT value should not be changed inside subroutines!
	!
	!    * AUX_MOD module from full hydrodynamic is deactivated (commented)!
	!    * Hidrodinamico2 subroutine from full hydrodynamic is deactivated (commented)!
    !
    !  Usage:
    !
    !    CALL MODELO
    !
    !    where
    !
    !    * no arguments are passed in this subroutine
    !
    !    uses modules and functions
    !
	!	 * function CALDAT	      		in	  caldat.f90
    !    * module     VARS_MAIN   		in    VARS_MAIN.f90
    !    * module     VARS_INERC  		in    VARSINERC.f90  !Quais? 
 	!	 * subroutine LECHUVA	  		in	  LECHUVA.f90
	!	 * subroutine LECLIMA	  		in	  LECLIMA.f90
	!	 * subroutine CELULA	  		in	  CELULA.f90
	!	 * subroutine REDE	      		in	  REDE.f90	
	!	 * subroutine flood_inercial	in	  flood_inercial.f90
	!
	!    Deactivated:	
	!	 * module	  AUX_MOD        	in	  HD_AUX_MOD.f90
	!	 * subroutine Hidrodinamico2	in	  Hidrodinamico2.f90
    !
    !	 opens
    !
    !    * Does not open files
    !
    !    reads
    !
    !    * Does not read files
    !
    !    creates
    !
    !    * Does not create files
    !    
    !
    !---------------------------------------------------------------------------------
    !  Licensing:
    !
    !    This code is distributed under the ...
    !
    !  Version/Modified: 
    !
    !    2014.25.11 - 25 November 2014 (By: Mino V. Sorribas)    
    !
    !  Authors:
    !
    !    Original fortran version by Walter Collischonn
    !    Present fortran version by:
    !    * Walter Collischonn
    !    * Rodrigo Cauduro Dias de Paiva
    !    * Diogo da Costa Buarque
    !    * Paulo Pontes Rógenes
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
    !   * Variables declarations and routines calls are all commented below.
    !
    !---------------------------------------------------------------------------------
	SUBROUTINE MODELO

	USE VARS_MAIN
	USE VARS_INERC
	USE SED_VARS !Hugo Fagundes 16/07/2019
	
	IMPLICIT NONE
    
	INTEGER K,KHID 					!indexes and counters
	INTEGER KC,JB,JC,KB,JCB,MWP, KSS, KSS2 	!... same
    INTEGER i, iUSO
	REAL(8) INmini(3), OUTmini(3), BALAreia, BALSilte, BALArgila


	! Initialize
	IT=0
    INmini=0.0
    OUTmini=0.0
    BALAreia=0.0
    BALSilte=0.0
    BALArgila=0.0
    
    WRITE(*,*)'SIMULANDO'
 
    ! Time Loop    
    DO WHILE (IT<NT)
		
		IT=IT+1
        if(mod(it,100)==0) write(*,*)100*iT/NT,'%'

		JDIA=IDINI+INT((IT+HORAINI-1)/(86400./DTP)) 	! Gets calendar day (big number)
		CALL CALDAT(JDIA,IMES,IDIA,IANO)

		DIAH(IT)=IDIA !stores day corresponding to iT time interval
		MESH(IT)=IMES !stores month corresponding to iT time interval
		ANOH(IT)=IANO !stores year corresponding to iT time interval
		HORAH(IT)=MOD(IT,24)+HORAINI 	!hour of day corresponding to iT time interval

		! Read rainfall data for all catchments in iT time-step
		ITCHUVA=IT
		CALL LECHUVA
	
		DO KC=1,NC
			PM(KC)=PM(KC)+P(KC) !Cumulative rainfall (used later for average)
		ENDDO		
	
		! Reads climate data
		TONTEM=TA
		CALL LECLIMA
	
		! Calls catchment routing/discharge for lateral inflow
		CALL CELULA		
	
        !*************************************************** Hugo Fagundes 04.07.19 ***************************************************!    
        IF(SED_INDEX==1)THEN
            CALL SED_BACIA
        ENDIF
		!*************************************************** Hugo Fagundes 04.07.19 ***************************************************!    
	
		! Saves detailed info for soil state variables in file NOSOLO.HIG.
		! Uses JC index for catchment and JB index for URH of interest.
		IF(ICALIB.EQ.0)THEN !only if not calibrating.
			JB=1
			JC=1
			!JC=57
			WRITE(FILSOL,75)IT,P(JC),W(JC,JB),SI(JC,JB),ET(JC,JB),CAF(JC,JB),QBAS(JC),QINT(JC),QSUP(JC)
			JB=2
			JC=2
			WRITE(FILSOL2,75)IT,P(JC),W(JC,JB),SI(JC,JB),ET(JC,JB),CAF(JC,JB),QBAS(JC),QINT(JC),QSUP(JC)			

!			write(971,66) (E0agua(iC),iC=1,nC)
!			write(972,66) (E0topo(iC),iC=1,nC)
!			write(973,66) (E0sup(iC),iC=1,nC)

		ENDIF

        ! Call main river routing routine using Muskingum-Cunge Method
		if(hdFLAG0==0)then
		    CALL REDE
		endif
		
		! Calculates lateral inflow ! necessary to run inertial model and to save files! 
        do ic=1,nc
		    QCEL2(IC)=QBAS(IC)+QINT(IC)+QSUP(IC) !sums surface, subsurface and 
            QITUDO(IC,IT)=QCEL2(IC) !save QCEL information at binary file !FMF and PETER 10/08/2016
		    QBASTUDO(IC,IT)=QBAS(IC) !save minibasin QBAS only information at binary file !FMF and PETER 10/08/2016
		enddo

        ! Calls river routing routine using Inertial Method
        if(hdFLAG0>0)then
            CALL flood_inercial
        endif

       !*************************************************** Hugo Fagundes 11.07.19 ***************************************************!    
        IF(SED_INDEX==1)THEN
            CALL SED_REDE
        ENDIF
		!*************************************************** Hugo Fagundes 11.07.19 ***************************************************!    

		! Stores calculated discharges in sub-basins with observed data - for calibration and assessment.
		DO K=1,NOBS
			KHID=IQOBS(K) 		! outlet catchment id
			QR(K,IT)=QJ2(KHID)  ! saves on QR(K,IT), that will be compared to QOBS in FOBJ routine.
		ENDDO
		
		! Stores discharges for file ouput when in simulation model
		IF(ICALIB.EQ.0)THEN 	! only if it is not calibration
			
			DO KB=1,NB				! store discharge in sub-basin
				JCB=IEXUT(KB) 		! outlet catchment of KB-th subbasin
				QRB(KB,IT)=QJ2(JCB)
			ENDDO
	
			DO K=1,NUMHIDG 				! store discharge in catchments defined in info.mgb
				KHID=IHIDG(K) 			! indexed catchment
				QRG(K,IT)=QJ2(KHID) 	! stored in QRG		
				
			    if(hdFLAG0>0)then !Stores the water level for inertial model 	           
                    HRG(K,IT)=Hfl(KHID) ! stored in HRG                  
                endif     
				
			ENDDO
            DO KSS2=1,AUXSS
                DO K=1,NUMSS(KSS2) 				! store sediment in catchments defined in info.mgb
				    KSS=ISS(KSS2,K) 			! indexed catchment

                    CAREIA(K,IT)=CSJ2(KSS,1)*(10.**(6.)) !@ (mg/l)
                    CSILTE(K,IT)=CSJ2(KSS,2)*(10.**(6.)) !@ (mg/l)
                    CARGIL(K,IT)=CSJ2(KSS,3)*(10.**(6.)) !@ (mg/l) 								
			    ENDDO
            ENDDO
!********************** !Hugo Fagundes 16/07/2019 ********************
               !!!     DO K=1,NUMSS(KSS2)    !NUMERO DE POSTOS RELACIONADOS AOS SS QUE FORAM INDICADOS PARA SEREM GRAVADOS
               !!!         KSS=ISS(KSS2,K)   !NUMERO DA MINI QUE CORRESPONDE AO POSTO DE SS A SER GRAVADO
               !!!
               !!!     DO K=1,NOBSS(KSS2)    !NUMERO DE POSTOS RELACIONADOS AOS SS QUE POSSUEM DADOS OBSERVADOS
               !!!         KSS=ISSOBS(K)     !NUMERO DA MINI QUE CORRESPONDE AO POSTO DE SS COM DADOS OBSERVADOS
               !!!
               !!!     DO K=1,NOBS           !NUMERO DE POSTOS COM DADOS DE VAZÃO OBSERVADOS
               !!!         KSS=IQOBS(K) 	 !NUMERO DA MINI QUE CORRESPONDE AO POSTO DE VAZÃO COM DADOS OBSERVADOS
               !!!
               !!!     DO K=1,NUMHIDG        !NUMERO DE POSTOS COM DADOS DE VAZÃO A SEREM GRAVADOS
               !!!         KSS=IHIDG(K)(K)   !NUMERO DA MINI QUE CORRESPONDE AO POSTO DE VAZÃO COM DADOS A SEREM GRAVADOS
            
            IF(SED_INDEX==1)THEN
                DO KSS2=1,AUXSS
                    DO K=1,NOBSS(KSS2) 
                         KSS=ISSOBS(KSS2,K) 
                         IF (FLAG_SEDVAR==1) THEN
                            CSSCAL(K,IT,KSS2)=(CSJ2(KSS,2)+CSJ2(KSS,3))*(10.**(6.))                              !CSS mg/L VAI SER COMPARADO A SSOBS(K,IT) NA ROTINA FOBJ   
                         ELSEIF(FLAG_SEDVAR==2) THEN
                            CSSCAL(K,IT,KSS2)=(CSJ2(KSS,2)+CSJ2(KSS,3))*(10.**(6.))*QJ2(KSS)*0.0864              !QSS ton/dia VAI SER COMPARADO A SSOBS(K,IT) NA ROTINA FOBJ 
                         ELSEIF(FLAG_SEDVAR==3) THEN
                            CSSCAL(K,IT,KSS2)=(CSJ2(KSS,1)+CSJ2(KSS,2)+CSJ2(KSS,3))*(10.**(6.))*QJ2(KSS)*0.0864  !QStot ton/dia VAI SER COMPARADO A SSOBS(K,IT) NA ROTINA FOBJ 
                         ENDIF
                    ENDDO
                    DO K=1,NUMSS(KSS2) 
                         KSS=ISS(KSS2,K) 
                         CAREIA(K,IT)=CSJ2(KSS,1)*(10.**(6.)) !@ (mg/l) 
                         CSILTE(K,IT)=CSJ2(KSS,2)*(10.**(6.)) !@ (mg/l) 
                         CARGIL(K,IT)=CSJ2(KSS,3)*(10.**(6.)) !@ (mg/l) 
                         IF (FLAG_SEDVAR==1) THEN
                            SEDstore(K,IT,KSS2)=(CSJ2(KSS,2)+CSJ2(KSS,3))*(10.**(6.))                              !CSS mg/L VAI SER COMPARADO A SSOBS(K,IT) NA ROTINA FOBJ   
                         ELSEIF(FLAG_SEDVAR==2) THEN
                            SEDstore(K,IT,KSS2)=(CSJ2(KSS,2)+CSJ2(KSS,3))*(10.**(6.))*QJ2(KSS)*0.0864              !QSS ton/dia VAI SER COMPARADO A SSOBS(K,IT) NA ROTINA FOBJ 
                         ELSEIF(FLAG_SEDVAR==3) THEN
                            SEDstore(K,IT,KSS2)=(CSJ2(KSS,1)+CSJ2(KSS,2)+CSJ2(KSS,3))*(10.**(6.))*QJ2(KSS)*0.0864  !QStot ton/dia VAI SER COMPARADO A SSOBS(K,IT) NA ROTINA FOBJ 
                         ENDIF
                    ENDDO     
                ENDDO
            ENDIF
            
            
            !************ DIOGO BUARQUE set/2012 ********************
			DO i = 1,nSEDmini
			    KHID=CONTIC(i) !CODIGO DA MINI-BACIA EM QUE SE DESEJA O HIDROGRAMA
				Qmini(i,IT)=QJ2(KHID) !QRG ARAMZENA OS HIDROGRAMAS NOS LOCAIS DESEJADOS
			ENDDO
            !************ DIOGO BUARQUE set/2012 ********************
            
			! Store discharge by water origin (i.e. surface, subsurface, groundwater) in a specified catchment
			MWP=1 						!catchment (this is manual)
			QB(IT)=QJ2(MWP)*PJB2(MWP)
			QBI(IT)=QB(IT)+QJ2(MWP)*PJI2(MWP)
			QBIS(IT)=QBI(IT)+QJ2(MWP)*PJS2(MWP)
            
		ELSE !Hugo Fagundes 16/07/2019 
!********************** !Hugo Fagundes 16/07/2019 ********************
            IF(SED_INDEX==1)THEN
                DO KSS2=1,AUXSS
                   DO K=1,NOBSS(KSS2) 
			            KSS=ISSOBS(KSS2,K) 
			            IF (FLAG_SEDVAR==1) THEN
                            CSSCAL(K,IT,KSS2)=(CSJ2(KSS,2)+CSJ2(KSS,3))*(10.**(6.))                              !CSS mg/L VAI SER COMPARADO A SSOBS(K,IT) NA ROTINA FOBJ   
                         ELSEIF(FLAG_SEDVAR==2) THEN
                            CSSCAL(K,IT,KSS2)=(CSJ2(KSS,2)+CSJ2(KSS,3))*(10.**(6.))*QJ2(KSS)*0.0864              !QSS ton/dia VAI SER COMPARADO A SSOBS(K,IT) NA ROTINA FOBJ 
                         ELSEIF(FLAG_SEDVAR==3) THEN
                            CSSCAL(K,IT,KSS2)=(CSJ2(KSS,1)+CSJ2(KSS,2)+CSJ2(KSS,3))*(10.**(6.))*QJ2(KSS)*0.0864  !QStot ton/dia VAI SER COMPARADO A SSOBS(K,IT) NA ROTINA FOBJ 
                         ENDIF
                    ENDDO     
                ENDDO
            ENDIF
!********************** !Hugo Fagundes 16/07/2019 ********************            
		ENDIF
        
	ENDDO !End time loop
	
!************ DIOGO BUARQUE ********************
IF(SED_INDEX==1)THEN
    OPEN(SEDSAI,FILE  = '.\output\SEDRIO_SAIDA.txt',STATUS='UNKNOWN')
    OPEN(SEDDEP,FILE  = '.\output\SEDRIO_DEPOSITO.txt',STATUS='UNKNOWN')
    OPEN(SEDERO,FILE  = '.\output\SEDRIO_EROSAO.txt',STATUS='UNKNOWN')
    OPEN(BALFLP,FILE  = '.\output\BAL_PLAN.txt',STATUS='UNKNOWN')  !@ DCB_HD_Sed - SALVA ENTRADAS E SAIDA NAS PLANICIES
        WRITE(SEDSAI,'(1A10,3A16)') 'IC', 'SAI_A', 'SAI_S', 'SAI_C'
        WRITE(SEDDEP,'(1A10,3A16)') 'IC', 'DEP_A', 'DEP_S', 'DEP_C'
        WRITE(SEDERO,'(1A10,3A16)') 'IC', 'ERO_A', 'ERO_S', 'ERO_C'
        WRITE(BALFLP,'(1A10,9A16)') 'IC', 'IN_FLPsand', 'IN_FLPsilt', 'IN_FLPclay', 'OUT_FLPsand', 'OUT_FLPsilt', 'OUT_FLPclay', 'QFL_in', 'QFL_out', 'DFL'      !@ DCB_HD_Sed
        DO i = 1,nSEDmini
            WRITE(SEDSAI,'(1I10,3ES16.7)'), CONTIC(i), iSAI(i,1), iSAI(i,2), iSAI(i,3)      !@ ton
            WRITE(SEDDEP,'(1I10,3ES16.7)'), CONTIC(i), DEPT(i,1), DEPT(i,2), DEPT(i,3)      !@ ton
            WRITE(SEDERO,'(1I10,3ES16.7)'), CONTIC(i), EROT(i,1), EROT(i,2), EROT(i,3)      !@ ton
            WRITE(BALFLP,'(1I10,9ES16.7)'), CONTIC(i), inFL(i,1), inFL(i,2), inFL(i,3), outFL(i,1), outFL(i,2), outFL(i,3), BalQFL(i,1), BalQFL(i,2), DFL(i,2)+DFL(i,3)  !@ DCB_HD_Sed (ton)
        ENDDO
    CLOSE(SEDSAI)
    CLOSE(SEDDEP)
    CLOSE(SEDERO)
    CLOSE(BALFLP)

    IF(ICALIB.EQ.0)THEN 	! only if it is not calibration
        INmini(1)  = sum(iENTRA(:,1)) !TON 
        INmini(2)  = sum(iENTRA(:,2)) !TON
        INmini(3)  = sum(iENTRA(:,3)) !TON
        !******************************************** Hugo Fagundes 13/08/2019 ***************************************************!
        DO IC=1,NC
            BALAreia= BALAreia + CSJ2(IC,1)*VolTREC2(IC)
            BALSilte= BALSilte + CSJ2(IC,2)*VolTREC2(IC)
            BALArgila= BALArgila + CSJ2(IC,3)*VolTREC2(IC)
        ENDDO
        !Soma os valores dos reservatórios + somatório do que saiu nos exutórios + o que ficou depositado nas planícies
    
        !somatório do que saiu nos exutórios + o que ficou depositado nas planícies + o que ainda está em suspensão nos trechos
        OUTmini(1) = OUTmini(1) + (sum(iSAI(:,1))) + BALAreia!TON  
        OUTmini(2) = OUTmini(2) + (sum(iSAI(:,2))) + ((Sum(inFL(:,2))) + (Sum(outFL(:,2)))) + BALSilte  !TON  
        OUTmini(3) = OUTmini(3) + (sum(iSAI(:,3))) + ((Sum(inFL(:,3))) + (Sum(outFL(:,3)))) + BALArgila!TON  

        write(*,*) 
        write(*,*) 
        write(*,*) 'VERFICACAO DE CONSERVACAO DE SEDIMENTOS' 
        write(*,*) 
        write(*,*) 'ENTRADA (Ton.)  = ', INmini(:)
        write(*,*) 'SAIDA   (Ton.)  = ', OUTmini(:)
        write(*,*) 'ERRO    (%)        = ', 100.*(INmini(1)-OUTmini(1))/INmini(1), 100.*(INmini(2)-OUTmini(2))/INmini(2), 100.*(INmini(3)-OUTmini(3))/INmini(3)
        write(*,*)
        !write(*,*)' BALAreia (tons), BALSilte (tons), BALArgila (tons)' !TON  
        !write(*,*) BALAreia, BALSilte, BALArgila                        !TON  
        pause
        !******************************************** Hugo Fagundes 13/08/2019 ***************************************************!
    ENDIF
ENDIF
	

75	FORMAT(I6,8F10.4)
66	FORMAT(<nC>F10.4)
67	FORMAT(<nC>F10.1)
68	FORMAT(<nC>F10.3)
69	FORMAT(<nC>F12.3)

     
	RETURN
	END