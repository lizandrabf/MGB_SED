	SUBROUTINE SED_Fob_Calib
	!ESTA SUBROTINA ANALISA A QUALIDADE DO AJUSTE ENTRE OS HIDROGRAMAS CALCULADOS E OBSERVADOS
	!COM BASE EM ALGUMAS FUNÇÕES OBJETIVO
	USE VARS_MAIN
    USE SED_VARS
	use vars_calib
    USE IFPORT
    
	IMPLICIT NONE
    
	INTEGER NANO,KB, KSS,K, IPOS, NPOSTOS
	REAL SR,ERRVT,SRL
	REAL POND(NBP,nf) 
	INTEGER INIAJU !DIA DO INICIO DA COMPARAÇÃO CALCULADO X OBSERVADO
	REAL SOMAOBS, SSOMAOBS !SOMATORIO DOS VALORES DE VAZÃO E DADOS DE RELATIVOS AOS SEDIMENTOS OBSERVADOS               !Hugo Fagundes
	REAL SOMACAL, SSOMACAL !SOMATORIO DOS VALORES DE VAZÃO E DADOS DE RELATIVOS AOS SEDIMENTOS CALCULADOS               !Hugo Fagundes
    REAL Resp(NOBSS(AUXCALIBSS),2)    ! MATRIZ PRA ARMAZENAR A SÉRIE DE DADOS DA CORRELAÇÃO ESPACIAL
	REAL SOMALOGS !SOMATORIO DOS VALORES DOS LOGARITMOS DE VAZÃO OBSERVADOS
	REAL SOMALCAL !SOMATORIO DOS VALORES DOS LOGARITMOS DE VAZÃO CALCULADOS
	INTEGER ITVAL, SITVAL !CONTA O NUMERO DE DIAS EM QUE EXISTE DADOS OBSERVADOS PARA COMPARAR
	REAL XMOBS,XLOGS, SSXOBS, CSSXCAL !VAZÃO MÉDIA OBSERVADA, MÉDIA DOS LOGARITMOS DAS VAZÕES OBSERVADAS,DADOS RELATIVOS AOS SED. OBSERVADOS, DADOS DE CSS CALCULADOS       !Hugo Fagundes
	REAL SQDES(NOBSS(AUXCALIBSS)), SSDES(NOBSS(AUXCALIBSS)),SSDPCAL, SSDPOBS !SOMATORIO DOS DESVIOS QUADRADOS ENTRE O VALOR OBSERVADO E CALCULADO, VARIANCIA DOS DADOS CALCULADOS E OBSERVADOS RELATIVOS AOS SED.    !Hugo Fagundes
	REAL SOMQ, SOMS,SSCOVAR !SOMATORIO DOS DESVIOS QUADRADOS ENTRE A VAZÃO OBSERVADA E A MÉDIA, COVARIANCIA DOS DADOS RELATIVOS AOS SED.        !Hugo Fagundes
	REAL SQLOG !IGUAL AO SQDES SÓ QUE PARA OS LOGARITMOS DAS VAZÕES
	REAL SOLOG !IGUAL AO SOMQ SÓ QUE PARA OS LOGARITMOS DAS VAZÕES
    integer i,j
    REAL Q10(2), Q50(2)
    REAL,ALLOCATABLE:: QSOROBS(:), QSORCALC(:) 		!observed discharges (without "missing data") for sorting and exceedance probaility curve

    
    SSCORRELtemp=0.0
    SSCORRELesp=0.0
    SSCORRELtudo=0.0
    R2=-100.0
	R2L=-100.0
	ERRV=1000.0
    R2SS=-100.0    
    RMSE=0.0
    ERRVS=0.0 
    SSKGE=-100.0
    KGEa=-100.0
    CoefB=-100
    DCPERM=1000
        
!************************************ HUGO FAGUNDES 21/09/17 - INICIO ******************************************   
    KSS=AUXCALIBSS !VAI FAZER APENAS PARA O DADO QUE SE ESTÁ CALIBRANDO
    iniaju=365 !HUGO FAGUNDES 03/10/2019
	SQDES=0.0
    SSDES=0.0       	      
    SSCORRELAC=0.0
    AUXOBS=0.0
    AUXCALC=0.0
    auxcont=0
    Resp=0.0
    ALLOCATE (QSOROBS(100000000))
    ALLOCATE (QSORCALC(100000000))
    QSOROBS=0.0
    QSORCALC=0.0
    VFO=0.0
    
    DO KB=1,NOBSS(KSS)           !PARA CADA POSTO
		IB=KB
        SOMS=0.0
        SSOMAOBS=0.0
		SSOMACAL=0.0
        SSXOBS=0.0
        CSSXCAL=0.0
        SSDPCAL=0.0
        SSDPOBS=0.0
        SSCOVAR=0.0
		SITVAL=0
        
        DO IT=iniaju,NT !PARA AVALIAR PERÍODO dd/mm/aaaa a dd/mm/aaaa
			IF(SSOBS(IB,IT,KSS).GE.0.0) THEN
                SSOMAOBS=SSOMAOBS+SSOBS(IB,IT,KSS)      !SOMATORIO DOS DADOS RELATIVOS AOS SEDIMENTOS OBSERVADOS        
				SSOMACAL=SSOMACAL+CSSCAL(IB,IT,KSS)           !SOMATORIO DOS DE SEDIMENTOS EM SUSPENSÃO CALCULADOS       
            	SITVAL=SITVAL+1                           !CONTADOR DE QUANTOS DADOS OBSERVADOS ENTRARAM NA ANALISE
                auxcont=auxcont+1
                AUXOBS(auxcont)=SSOBS(IB,IT,KSS)        !ARMAZENA OS DADOS OBSERVADOS PARA CALCULAR CORRELTUDO
                AUXCALC(auxcont)=CSSCAL(IB,IT,KSS)            !ARMAZENA OS DADOS CALCULADOS PARA CALCULAR CORRELTUDO    
                QSOROBS(auxcont)=SSOBS(IB,IT,KSS)    !ARMAZENA OS DADOS OBSERVADOS PARA CURVA DE PERMANENCIA
                QSORCALC(auxcont)=CSSCAL(IB,IT,KSS)    !ARMAZENA OS DADOS CALCULADOS PARA CURVA DE PERMANENCIA
            ENDIF                
        ENDDO
        !**********Curva de Permanencia**********
        Call SORTQQ (LOC(QSOROBS), auxcont, SRT$REAL4)            
        Call SORTQQ (LOC(QSORCALC), auxcont, SRT$REAL4)            
        !WRITE(*,*)  auxcont
        IPOS=INT(90*auxcont/100) 
        Q10(1)=QSOROBS(IPOS)  !Q10 observada 
        Q10(2)=QSORCALC(IPOS)  !Q10 calculada
        IPOS=INT(50*auxcont/100) 
        Q50(1)=QSOROBS(IPOS)  !Q50 observada 
        Q50(2)=QSORCALC(IPOS)  !Q50 calculada
        DCPERM(IB)=100*((Q50(2)-Q10(2))/(Q50(1)-Q10(1))-1)
        !**********Curva de Permanencia**********
        
        SSXOBS=SSOMAOBS/SITVAL                      !MEDIA DOS DADOS OBSERVADOS     
        CSSXCAL=SSOMACAL/SITVAL                     !MEDIA DOS DADOS CALCULADOS  
        !WRITE(*,*)'SITVAL,SSXOBS,CSSXCAL'
        !WRITE(*,*)SITVAL,SSXOBS,CSSXCAL 
        Resp(KB,1)=SSXOBS                           !ARMAZENA A MÉDIA DOS DADOS OBS. DO POSTO PARA CALCULAR CORREL ESPACIAL
        Resp(KB,2)=CSSXCAL                          !ARMAZENA A MÉDIA DOS DADOS CAL. DO POSTO PARA CALCULAR CORREL ESPACIAL

        DO IT=iniaju,NT !PARA AVALIAR PERÍODO dd/mm/aaaa a dd/mm/aaaa
            IF(SSOBS(IB,IT,KSS).GE.0.0) THEN
                SSDPOBS=SSDPOBS+(SSOBS(IB,IT,KSS)-SSXOBS)**2.0                    !VARIANCIA DOS DADOS OBSERVADOS
                SSDPCAL=SSDPCAL+(CSSCAL(IB,IT,KSS)-CSSXCAL)**2.0                        !VARIANCIA DOS DADOS CALCULADOS
                SSCOVAR=SSCOVAR+(SSOBS(IB,IT,KSS)-SSXOBS)*(CSSCAL(IB,IT,KSS)-CSSXCAL) !COVARIANCIA DOS DADOS OBSERVADOS E CALCULADOS
                SSDES(IB)=SSDES(IB)+(SSOBS(IB,IT,KSS)-CSSCAL(IB,IT,KSS))**2.0         !CALCULA O SOMATÓRIO DO NUMERADOR DA FORMULA DO COEFICIENTE DE NASH  
                SOMS=SOMS+(SSOBS(IB,IT,KSS)-SSXOBS)**2.0                          !CALCULA O SOMATÓRIO DO DENOMINADOR DA FORMULA DO COEFICIENTE DE NASH    
            ENDIF
        ENDDO
        !WRITE(*,*)'SSDPOBS,SSDPCAL,SSCOVAR'
        !WRITE(*,*)SSDPOBS,SSDPCAL,SSCOVAR 
        IF (SITVAL>0) THEN
            R2SS(IB)=1.-SSDES(IB)/SOMS                                        !CALCULA O NASH DOS DADOS DE SEDIMENTOS
            ERRVS(IB)=(SSOMACAL/SSOMAOBS-1.0)*100.                            !CALCULA O ERRO DOS DADOS DE SEDIMENTOS
            RMSE(IB)=(SSDES(IB)/SITVAL)**0.5                              !CALCULA A RAIZ DO ERRO QUADRATICO MÉDIO DOS SEDIMENTOS
            CoefB(IB)=(CSSXCAL-SSXOBS)/((SSDPOBS/SITVAL)**0.5)                !Beta: Calcula a relação entre a diferença relativa das medias das vazões, calculadas e observadas, normalizado pelo desvio
                                                                              !padrão das vazões observadas durante o período de investigação. Gupta et al. (2009) Wohling et al. (2013)
            !WRITE(*,*)'CoefB'
            !WRITE(*,*)CoefB(IB),(CSSXCAL-SSXOBS),SSDPOBS,SITVAL,(SSDPOBS/SITVAL)
        ENDIF
        
        SSCORRELAC(KB)=SSCOVAR/((SSDPCAL*SSDPOBS)**0.5)                       !COEFICIENTE DE CORRELAÇÃO LINEAR TEMPORAL PARA CADA POSTO
        !WRITE(*,*)'SSCORRELAC(KB)'
        !WRITE(*,*)SSCORRELAC(KB)
        SSCORRELtemp(KSS)=SSCORRELtemp(KSS) + SSCORRELAC(KB)/NOBSS(KSS)        !CALCULA A MÉDIA DAS CORRELAÇÕES TEMPORAIS 
        SSKGE(KB)=1-((1-DBLE(SSCORRELAC(KB)))**2 + (1-DBLE(CSSXCAL/SSXOBS))**2 + (1-DBLE((SSDPCAL/SSDPOBS)**0.5))**2)**0.5      !CALCULO DO COEFICIENTE KGE (GUPTA ET AL, 2009)
        KGEa(KB)=DBLE((SSDPCAL/SSDPOBS)**0.5)
        !pause
            
            !CALIBRA U
            !F1CORREL(KB)=(1-UCORRELAC)**2
            !F2DMEDIA(KB)=(1-(UXMCAL/UXMOBS))**2
            !F3DDESVPAD(KB)=(1-(UDPCAL/UDPOBS))**2

    ENDDO
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! CORREL ESPACIAL !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !      SSOMAOBS=0.0
		!SSOMACAL=0.0
  !      SSXOBS=0.0
  !      CSSXCAL=0.0
  !      SSDPCAL=0.0
  !      SSDPOBS=0.0
  !      SSCOVAR=0.0
  !      KB=0
  !      
  !      DO KB=1,NOBSS(KSS) ! A SÉRIE DE DADOS DA CORRELAÇÃO ESPACIAL É DO TAMANHO DO NUMERO DE POSTOS EXISTENTES
  !         SSOMAOBS=SSOMAOBS+Resp(KB,1)      !SOMATORIO DOS DADOS RELATIVOS AOS SEDIMENTOS OBSERVADOS        
  !         SSOMACAL=SSOMACAL+Resp(KB,2)      !SOMATORIO DOS DE SEDIMENTOS EM SUSPENSÃO CALCULADOS     
  !      ENDDO
  !      
  !      SSXOBS=SSOMAOBS/NOBSS(KSS)                      !MEDIA DOS DADOS OBSERVADOS     
  !      CSSXCAL=SSOMACAL/NOBSS(KSS)                     !MEDIA DOS DADOS CALCULADOS        
  !      DO KB=1,NOBSS(KSS)
  !         SSDPOBS=SSDPOBS+(Resp(KB,1)-SSXOBS)**2.0                    !VARIANCIA DOS DADOS OBSERVADOS
  !         SSDPCAL=SSDPCAL+(Resp(KB,2)-CSSXCAL)**2.0                  !VARIANCIA DOS DADOS CALCULADOS
  !         SSCOVAR=SSCOVAR+(Resp(KB,1)-SSXOBS)*(Resp(KB,2)-CSSXCAL) !COVARIANCIA DOS DADOS OBSERVADOS E CALCULADOS
  !      ENDDO
  !      
  !      SSCORRELesp(KSS)=SSCOVAR/((SSDPCAL*SSDPOBS)**0.5)                       !COEFICIENTE DE CORRELAÇÃO LINEAR ESPACIAL
  !      
  !      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! CORREL TUDO !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !      SSOMAOBS=0.0
		!SSOMACAL=0.0
  !      SSXOBS=0.0
  !      CSSXCAL=0.0
  !      SSDPCAL=0.0
  !      SSDPOBS=0.0
  !      SSCOVAR=0.0
  !      KB=0
  !      
  !      DO KB=1,auxcont                             ! A SÉRIE DE DADOS DA CORRELAÇÃO DE TUDO É DO TAMANHO DOS DADOS OBSERVADOS NOS POSTOS EXISTENTES
  !         SSOMAOBS=SSOMAOBS+AUXOBS(KB)             !SOMATORIO DOS DADOS RELATIVOS AOS SEDIMENTOS OBSERVADOS        
  !         SSOMACAL=SSOMACAL+AUXCALC(KB)            !SOMATORIO DOS DE SEDIMENTOS EM SUSPENSÃO CALCULADOS     
  !      ENDDO
  !      
  !      SSXOBS=SSOMAOBS/auxcont                      !MEDIA DOS DADOS OBSERVADOS     
  !      CSSXCAL=SSOMACAL/auxcont                     !MEDIA DOS DADOS CALCULADOS        
  !
  !      DO KB=1,auxcont
  !         SSDPOBS=SSDPOBS+(AUXOBS(KB)-SSXOBS)**2.0                    !VARIANCIA DOS DADOS OBSERVADOS
  !         SSDPCAL=SSDPCAL+(AUXCALC(KB)-CSSXCAL)**2.0                  !VARIANCIA DOS DADOS CALCULADOS
  !         SSCOVAR=SSCOVAR+(AUXOBS(KB)-SSXOBS)*(AUXCALC(KB)-CSSXCAL)   !COVARIANCIA DOS DADOS OBSERVADOS E CALCULADOS
  !      ENDDO
  !      
  !      SSCORRELtudo(KSS)=SSCOVAR/((SSDPCAL*SSDPOBS)**0.5)                       !COEFICIENTE DE CORRELAÇÃO LINEAR DE TUDO
  !
    !SR=0.0
	!SRL=0.0
    
	!ERRVT=0.0

    ! Considera funcao objetivo média ponderada dos postos que entram na calibracao
	IF (icalib==3) THEN
        !POND=0.0
        i=NOBSS(AUXCALIBSS)
	    !DO KB=1,i !USAR APENAS "NOBS" PARA CALIBRAR VAZÃO  HUGO FAGUNDES 27/09/17
            ! DO j=1,nf
                !  IF (calibFLAG(KB,j)>0.0) POND(KB,j)=calibFLAG(KB,j)/sum(calibFLAG(1:i,j)) 
            ! ENDDO
        !ENDDO
        
        DO KB=1,i    
            !VFO(1)=VFO(1)+(1.-SSCORRELAC(KB))
            VFO(1)=VFO(1)+(1-R2SS(KB))
            !VFO(1)=VFO(1)+(1-SSKGE(KB))
            !VFO(2)=VFO(2)+ABS(CoefB(KB))
            VFO(2)=VFO(2)+ABS(ERRVS(KB)) 
            !VFO(2)=VFO(2)+RMSE(KB) 
            !VFO(3)=VFO(3)+ABS(DCPERM(KB))
            !VFO(3)=VFO(3)+(KGEa(KB)-1)
            VFO(3)=VFO(3)+(1-SSKGE(KB))
        ENDDO    
    ENDIF
    
 
!	do kb=1,nobs
!!	write(*,*) kb,VFO(1)
!!	write(*,*) R2(kb)
!!	write(*,*) POND(kb,1)
!		IF (POND(KB,1)>0.0) VFO(1)=VFO(1)+(1.-R2(kb))*POND(kb,1)
!		IF (POND(KB,2)>0.0) VFO(2)=VFO(2)+(1.-R2L(kb))*POND(kb,2)
!		IF (POND(KB,3)>0.0) VFO(3)=VFO(3)+ABS(ERRV(kb))*POND(kb,3)
!	enddo
	

!	WRITE(*,*)(R2(IB),IB=1,NOBS)
!	WRITE(*,*)(R2L(IB),IB=1,NOBS)
!	WRITE(*,*)(ERRV(IB),IB=1,NOBS)
!************************************ HUGO FAGUNDES 27/09/17 - FIM ******************************************

	write(*,*) 'VFO(1),VFO(2),VFO(3)'
    write(*,*) VFO(1),VFO(2),VFO(3)
73  FORMAT(F10.3)
    
	RETURN
	END