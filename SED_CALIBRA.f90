	SUBROUTINE CALIBRA_SED
	! Subrotina Calibra organizada !!!! ADAPTADA POR HUGO FAGUNDES E, 22/09/17
	!-----------------------------------------------------------------------
	!
	! Descri��o das vari�veis locais:
	!
	! I,J = vari�veis auxiliares
	!
	!----------------------------------------------------------------------
	!
	use VARS_MAIN
	use VARS_CALIB
    use SED_VARS            
	IMPLICIT NONE
	! Declara��o de vari�veis:
	! Vari�veis locais:
	INTEGER I,L,J,KB,K,IW,KS,KSS2,IFUNC,IPAR,JF,IBOBO,IRMAX,CONT,j2,j3
	REAL XBOBO,RAN1
	INTEGER ICONG,ITESTE
	integer i1,i2,i3
	INTEGER NPCOMPLEX
	NPCOMPLEX=5
	!----------------------------------------------------------------------------
    
	!***************OTIMIZA��O***********************************************************************
	!************************************************************************************************
	OPEN(FILEVO,FILE='.\output\EVOLUTION.TXT',STATUS='UNKNOWN')
    OPEN(FILRESCALIB, FILE='.\output\RESULCALIB.TXT', STATUS='UNKNOWN') !ARQUIVO DE SA�DA COM OS RESULTADOS DA CALIBRA��O           !HUGO FAGUNDES 02/10/19


	!GERA A SEMENTE DO PROCESSO DE GERA��O DE NUMEROS ALEATORIOS
	CALL SEMENTE(ISEED) !GERA SEMENTE DE PROCESSO ALEAT�RIO	
	IBOBO=MOD(ISEED,100)
	DO IW=1,IBOBO
		XBOBO=RAN1(ISEED)
	ENDDO

	!GERA VALOR RELATIVO DO PARAMETRO E MULTIPLICA PELA FAIXA DE VALIDADE
	   
    PAR(1,1)=11.8   !Alfa da equa��o MUSLE (Williams, 1975) CHUTE INICIAL DOS PARAMETROS ALFA
    PAR(2,1)=0.56   !Beta da equa��o MUSLE (Williams, 1975) CHUTE INICIAL DOS PARAMETROS BETA
    PAR(3,1)=1.0    !Chute inicial para mudan�a proporcional do par�metro TKS
    DO I=1,NB-1
        PAR(3*I+1,1)=11.8
        PAR(3*I+2,1)=0.56
        PAR(3*I+3,1)=1.00
    ENDDO
    !PAR TER� DIMENSOES (NUMERO DE PARAMETROS X NUMERO DE INDIVIDUOS)
    DO I=2,NS
	    DO L=1,NPAR
            PPAR(L,I)=RAN1(ISEED)
            PAR(L,I)=PMIN(L)+PPAR(L,I)*(PMAX(L)-PMIN(L))
        ENDDO
    ENDDO
    xyz=1 !Vari�vel auxiliar para calibra��o

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!   HUGO FAGUNDES 22/09/2017 - INICIO    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    PARX=PAR(:,1)
    CALL CalibParam_SED    
    CALL CONDINIC 
    CALL SED_PARAM
    CALL SED_INICIAL
	CALL MODELO
    !!!CALL SED_Fob_Calib !FUN��O OBJETIVO PARA CALIBRA��O DO MODELO DE SEDIMENTOS
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!   HUGO FAGUNDES 22/09/2017 - FIM    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
	!-------------------------------------
	KOUNTF=0
	write(*,*) 'Gera populacao inicial:'
	!AVALIA A FUN��O EM CADA PONTO
	DO I=1,NS   !HOF 09/01/2019

		PARX=PAR(:,I)
		CALL CalibParam_SED
		!********************************************************************************************************
		CALL CONDINIC ! RP
        CALL SED_PARAM
        CALL SED_INICIAL
        !CALL MODELO
		CALL SED_CALIB_MODELO
		!********************************************************************************************************
		!!!CALL SED_Fob_Calib !FUN��O OBJETIVO
		!!!KOUNTF=KOUNTF+1 !CONTA QUANTAS VEZES A F.O. FOI AVALIADA (NUMERO DE SIMULA��ES)
		!!!!ARMAZENA O VALOR DAS FUN��ES OBJETIVO
		!!!DO JF=1,NF
		!!!	FO(I,JF)=VFO(JF)
		!!!ENDDO
		!!!write(*,'(A9,I5,A1,I3,A5,3F7.4)') 'Individuo',I,'/',NS,'F.O.=',(FO(I,JF),JF=1,NF)
  !!!      write(*,*) PARX(1),PARX(2),PARX(3)
	ENDDO	
	!-------------------------------------
	!-------------------------------------
	FMIN=999999999999999999.9
	ISHUFFLE=0
	IRMAX=10
	
	DO WHILE(IRMAX.GT.1.and.ISHUFFLE<iMaxGen) !ITERA ENQUANTO HOUVER MAIS DE 1 CATEGORIA HIERARQUICA
		
		! Itera enquanto ISHUFFLE<maximo de geracoes
		
		ISHUFFLE=ISHUFFLE+1
		IRUIM=0

		!INICIA HIERARQUIZA��O DAS SOLU��ES (RANKING)
		IPARET=1
		ISOMA=0
		IRMAX=1
3000	DO I=1,NS
			IF(IPARET(I).EQ.IRMAX)THEN
				DO J=1,NS
					IF(J.NE.I.AND.IPARET(J).EQ.IRMAX)THEN	
						DO K=1,NF
							IF(FO(I,K).GT.FO(J,K))THEN
								IDOMIN=1
								IDD=J
							ELSE
								IDOMIN=0
								GOTO 1000
							ENDIF
						ENDDO
						IF(IDOMIN.EQ.1)THEN !ESTE PONTO � DOMINADO
							GOTO 2000 
						ENDIF
1000				ENDIF
				ENDDO
2000			IF(IDOMIN.EQ.1)THEN
					!WRITE(2,*)I,X(I),' DOMINADO POR',IDD
					IPARET(I)=IPARET(I)+1
					ISOMA=1
				ELSE
					!WRITE(2,*)I,X(I),' NAO DOMINADO'
					IPARET(I)=IRMAX
				ENDIF
				IDOMIN=0
				IDD=0
			ENDIF
		ENDDO
		IF(ISOMA.EQ.1)THEN
			IRMAX=IRMAX+1
			ISOMA=0
			GOTO 3000
		ENDIF
		IF(IRMAX==1)EXIT !QUANDO TODOS OS PONTOS EST�O NO N�VEL 1, ENCERRA A CALIBRA��O
		!---------------------------------------------
		WRITE(FILEVO,*)
		WRITE(FILEVO,*)ISHUFFLE
		DO I=1,NS		   !GRAVA PONTOS NO EVOLUTION.TXT
		    WRITE(FILEVO,172)(PAR(L,I),L=1,NPAR),(FO(I,JF),JF=1,NF),IPARET(I)
		ENDDO
172		FORMAT(<NPAR>F10.3,<NF>F12.4,I8)
		!---------------------------------------------

		!QUANDO CHEGA AQUI EXISTEM NS PONTOS AVALIADOS, CADA UM COM UM RANKING DADO POR IPARET
		!IPARET(I)=1 INDICA OS MELHORES PONTOS E IPARET(I)=IRMAX INDICA OS PIORES PONTOS
		NPC=0
		DO I=1,NS
			IF(IPARET(I).EQ.IRMAX)THEN
				NPC=NPC+1 !CONTA O NUMERO DE PONTOS QUE TEM RANKING IGUAL AO PIOR DE TODOS (IRMAX)
				IRUIM(NPC)=I !GUARDA O LOCAL EM QUE EST� O PONTO RUIM
			ENDIF
		ENDDO

		!DEFINE PROBABILIDADE DE QUE CADA UM DOS NS PONTOS SEJA ESCOLHIDO PARA FAZER PARTE DE UM SIMPLEX
		SPARET=0.0
		
		RMAX=IRMAX
		DO I=1,NS
			SPARET=SPARET+IPARET(I)
		ENDDO
		DO I=1,NS
			PROB(I)=(RMAX-IPARET(I)+1.)/(NS*(RMAX+1.)-SPARET)
		ENDDO
		DO I=2,NS
			PROB(I)=PROB(I-1)+PROB(I) !CALCULA PROBABILIDADE ACUMULADA
		ENDDO
		IF(ABS(PROB(NS)-1.0).GT.0.001)STOP 'ALGO ERRADO PROBABILID'
		PROB(NS)=1.0 !ESTA LINHA CORRIGE EVENTUAIS ERROS DE ARREDONDAMENTO

		!VARI�VEL QUE CONTEM OS NPAR+1 PONTOS DO SIMPLEX 1...NPC, DEFINIDOS POR NPAR PAR�METROS
		!ALLOCATE (SPAR(NPC,NPAR+1,NPAR),FPLEX(NPC,NPAR+1,NF))
		ALLOCATE (SPAR(NPC,NPCOMPLEX,NPAR),FPLEX(NPC,NPCOMPLEX,NF))

		DO K=1,NPC !EVOLUI CADA UM DOS COMPLEXOS
			LRUIM=IRUIM(K)
			DO IPAR=1,NPAR
				SPAR(K,1,IPAR)=PAR(IPAR,LRUIM) !ARMAZENA PIOR PONTO
			ENDDO

			!DO J=2,NPAR+1 !ESCOLHE OS OUTROS NPAR PONTOS PARA CADA COMPLEXO
			DO J=2,NPCOMPLEX !ESCOLHE OS OUTROS NPAR PONTOS PARA CADA COMPLEXO
				XRAND=RAN1(ISEED) !GERA NUMERO ALEATORIO
				PINI=0.0
				DO I=1,NS
					IF(XRAND.GE.PINI.AND.XRAND.LE.PROB(I))THEN
						!I � O PONTO ESCOLHIDO
						DO IPAR=1,NPAR
							SPAR(K,J,IPAR)=PAR(IPAR,I)	!GUARDA VALORES DOS PAR�METROS
							DO IFUNC=1,NF
								FPLEX(K,J,IFUNC)=FO(I,IFUNC) !GUARDA VALOR DAS FUN��ES
							ENDDO
						ENDDO
					ENDIF
					PINI=PROB(I)
				ENDDO
			ENDDO
		ENDDO

107		FORMAT(I4,7F8.2)

		!CHEGANDO AQUI FORAM ESCOLHIDOS OS SIMPLEXES

		!----------------------------------------------------------			

		!QUANDO CHEGA AQUI, O PROGRAMA TEM UMA AMOSTRA DE NS PONTOS NO ESPA�O
		!DE DIMENSAO NPAR. A FUN�AO OBJETIVO FOI AVALIADA EM TODOS OS PONTOS
		!E A CADA PONTO FOI DADO UM RANKING MULTI-OBJETIVO
		!A AMOSTRA FOI SEPARADA EM NPC COMPLEXOS, CADA UM COM NS+1 PONTOS. 
		!A PARTIR DESTE PONTO INICIA O QUE OS INVENTORES DO SCE-UA CHAMAM DE 
		!COMPETITIVE COMPLEX EVOLUTION.
		
		DO K=1,NPC !EVOLUI CADA UM DOS SIMPLEXOS
				LRUIM=IRUIM(K)
				SOMAPAR=0.0
				IREJECT=0
				!DO J=2,NPAR+1 !CALCULA O CENTROIDE DOS MELHORES PONTOS
				DO J=2,NPCOMPLEX !CALCULA O CENTROIDE DOS MELHORES PONTOS
					DO L=1,NPAR
						SOMAPAR(L)=SOMAPAR(L)+SPAR(K,J,L)
					ENDDO
				ENDDO
				!SOMAPAR=SOMAPAR/(NPAR) !CALCULA O CENTROIDE DOS MELHORES PONTOS (MEDIA DAS COORDENADAS DO PONTO)
				SOMAPAR=SOMAPAR/(NPCOMPLEX-1) !CALCULA O CENTROIDE DOS MELHORES PONTOS (MEDIA DAS COORDENADAS DO PONTO)
			!write(*,*)' complexo:',k
			!write(*,*)' pior ponto: ',(spar(k,1,l),l=1,npar)
			!write(*,*)' centr�ide: ',(somapar(l),l=1,npar)
	
			DO L=1,NPAR
				DIFPAR=SOMAPAR(L)-SPAR(K,1,L) !DIST�NCIA DO PIOR PONTO AO CENTROIDE
				REFLEX(L)=SOMAPAR(L)+DIFPAR	!COORDENADAS DO PONTO DE REFLEX�O
				CONTRA(L)=SOMAPAR(L)-DIFPAR/2. !COORDENADAS DO PONTO DE CONTRA��O
			ENDDO
			!write(*,*)' reflex: ',(reflex(l),l=1,npar)
			!Verifica se o ponto gerado por reflex�o est� dentro da regi�o v�lida
			DO L=1,NPAR
				IF(REFLEX(L).LT.PMIN(L).OR.REFLEX(L).GE.PMAX(L))THEN
					write(*,*)' reflex rejeitado parametro',l
					!write(*,*)' contra: ',(contra(l),l=1,npar)
					IREJECT=1
					GOTO 4000
				ENDIF
			ENDDO


			PARX=REFLEX
            CALL CalibParam_SED
		    !********************************************************************************************************
		    CALL CONDINIC ! RP
            CALL SED_PARAM
            CALL SED_INICIAL
		    !CALL MODELO
            CALL SED_CALIB_MODELO
		    !********************************************************************************************************
		    CALL SED_Fob_Calib !FUN��O OBJETIVO
			KOUNTF=KOUNTF+1 !CONTA QUANTAS VEZES A F.O. FOI AVALIADA (NUMERO DE SIMULA��ES)
			!ARMAZENA O VALOR DAS FUN��ES OBJETIVO
			DO JF=1,NF
				FO(LRUIM,JF)=VFO(JF)
			ENDDO
			!VERIFICA SE O PONTO REFLEX � DOMINADO
			!SOMENTE ACEITA O PONTO SE REFLEX � N�O DOMINADO
			!DO J=2,NPAR+1 !VERIFICA SE O PONTO J DOMINA O NOVO PONTO
			DO J=2,NPCOMPLEX !VERIFICA SE O PONTO J DOMINA O NOVO PONTO
				IDOMIN=0
				DO IFUNC=1,NF
					IF(FPLEX(K,J,IFUNC).LT.FO(LRUIM,IFUNC))THEN
						IDOMIN=IDOMIN+1	!CONTA PARA QUANTAS FUN��ES OBJ. O NOVO PONTO � PIOR
					ENDIF
				ENDDO
				IF(IDOMIN.EQ.NF)THEN !SE O NOVO PONTO � PIOR PARA TODAS AS F. O.
					!ESTE PONTO J DOMINA O NOVO PONTO
					!REJEITA O NOVO PONTO
					WRITE(*,*)'PONTO REJEITADO POR',(FPLEX(K,J,IFUNC),IFUNC=1,NF)
					IREJECT=1
					GOTO 4000
				ENDIF
			ENDDO
4000		IF(IREJECT.EQ.0)THEN
				!SE CHEGOU AQUI O PONTO DE REFLEX�O FOI ACEITO
				WRITE(*,*)'PONTO ACEITO'
				DO IPAR=1,NPAR
					PAR(IPAR,LRUIM)=REFLEX(IPAR)
				ENDDO
	
			ELSE
				!SE CHEGOU AQUI O PONTO DE REFLEX�O FOI REJEITADO - USA CONTRA��O

				PARX=CONTRA
				CALL CalibParam_SED
		        !********************************************************************************************************
		        CALL CONDINIC ! RP
                CALL SED_PARAM
                CALL SED_INICIAL
		        !CALL MODELO
                CALL SED_CALIB_MODELO
		        !********************************************************************************************************
		        CALL SED_Fob_Calib !FUN��O OBJETIVO
				KOUNTF=KOUNTF+1 !CONTA QUANTAS VEZES A F.O. FOI AVALIADA (NUMERO DE SIMULA��ES)
				!ARMAZENA O VALOR DAS FUN��ES OBJETIVO
				DO JF=1,NF
					FO(LRUIM,JF)=VFO(JF)
				ENDDO
	
				DO IPAR=1,NPAR
					PAR(IPAR,LRUIM)=CONTRA(IPAR)
				ENDDO
			ENDIF

		ENDDO !FIM DO LOOP DOS SIMPLEXOS
	
		DEALLOCATE (SPAR,FPLEX)
	
		VMIN=99999999999990.0
		DO KS=1,NS
			DO JF=1,NF
				VMIN(JF)=MIN(VMIN(JF),FO(KS,JF))
			ENDDO
		ENDDO

		WRITE(*,79)ISHUFFLE,IRMAX,NPC,(VMIN(JF),JF=1,NF)
79		FORMAT(3I8,8F12.4)


    ENDDO !FIM DO LOOP DOS SHUFFLES (S� SAI QUANDO IRMAX=1)
	!!!!!************************************************************************************************
 !!!!   WRITE(FILRESCALIB,*)AUXCALIBSS !ESCREVE QUAL CONJUNTO DE DADOS ESTA SENDO USADO PARA CALIBRA��O
 !!!!   DO KS=1,NS
 !!!!       PARX=PAR(:,KS)
 !!!!       !ESCREVE OS PAR�METROS ALFA E BETA
 !!!!       !WRITE(FILRESCALIB,73)(PARX(K),K=1,NPAR)
 !!!!       CONT=0
 !!!!       CALL CalibParam_SED
	!!!!	!********************************************************************************************************
	!!!!	CALL CONDINIC ! RP
 !!!!       CALL SED_PARAM
 !!!!       CALL SED_INICIAL
	!!!!	!CALL MODELO
 !!!!       CALL SED_CALIB_MODELO
	!!!!	!********************************************************************************************************
 !!!!
 !!!!       CALL SED_Fob_Calib !FUN��O OBJETIVO PARA OS SEDIMENTOS PARA TODOS OS CONJUNTOS DE DADOS
 !!!!           
 !!!!       !ARMAZENA OS DADOS (PODERIA APENAS ESCREVER, SE N�O QUISER ARMAZENAR, COMENTAR AS LINHAS ABAIXO
 !!!!       RESULTSS(CONT+1,KS)=SSCORRELAC(AUXCALIBSS)
 !!!!       RESULTSS(CONT+2,KS)=KGEa(AUXCALIBSS)
 !!!!       RESULTSS(CONT+3,KS)=CoefB(AUXCALIBSS)
 !!!!       RESULTSS(CONT+4,KS)=R2SS(AUXCALIBSS)        !HOF 09/01/2019
 !!!!       RESULTSS(CONT+5,KS)=SSKGE(AUXCALIBSS)
 !!!!       RESULTSS(CONT+6,KS)=RMSE(AUXCALIBSS)
 !!!!       RESULTSS(CONT+7,KS)=ERRVS(AUXCALIBSS)
 !!!!       RESULTSS(CONT+8,KS)=DCPERM(AUXCALIBSS)
 !!!!
 !!!!       !ESCREVE OS PAR�METROS ALFA E BETA E OS RESULTADOS DAS FUN��ES OBJETIVO        
 !!!!       WRITE(FILRESCALIB,73)(PARX(K),K=1,NPAR),(RESULTSS(K,KS),K=1,8) !8 numero de fun��es objetivo  !(RESULTSS(K,KS),K=1,NF*AUXSS) HOF 09/01/2019
 !!!!   ENDDO  
    CLOSE (FILRESCALIB)
    73 FORMAT(F10.3)
    
	RETURN
	END