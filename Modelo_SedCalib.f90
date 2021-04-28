!	SUBROUTINE MODELO_SedCalib
!	!Esta subrotina comanda o loop do tempo do modelo hidrológico e chama as 
!	!rotinas de balanco e propagação nas células e de propagação na rede de drenagem
!	USE VARS_MAIN
!	USE VARS_INERC
!	USE SED_VARS
!	IMPLICIT NONE
!
!	INTEGER K,KHID !CONTADORES
!	INTEGER KC,JB,JC,KB,JCB,MWP,KSS, KSS2	!CONTADORES
!    integer count_daniel,dt_daniel
!	integer iTwrite
!	INTEGER i
!	REAL INmini(3), OUTmini(3)
!
!	IT=0
!    Dt_daniel=INT((NT-IT)/10)
!    itWrite = 0 !@ DCB ago/2012 (sem isso, o valor é grande no Fortran Intel)
!    WRITE(*,*)  !@ DCB set/2012
!    WRITE(*,*)'SIMULANDO'   !@ DCB set/2012
!    WRITE(*,*)  !@ DCB set/2012
!    DO WHILE (IT<NT)
!		IT=IT+1
!		if (icalib==0.and.itWrite<iT) then
!    	WRITE(*,*) '     Passo de tempo: ', IT   !@ DCB set/2012 ,maxval(QJ2),maxloc(QJ2),minval(QJ2),minloc(QJ2),QJ2(nC) !QJ2(6848)! RP
!!@ DCB ago/2012			write(7000,*) 'iT',iT
!!			write(*,*) 'iT',iT
!			itWrite=itWrite+1
!        endif
!        
!		if(it==count_daniel)then
!		   WRITE(*,701)'**'
!			701		FORMAT(A2,$)
!			count_daniel=Count_daniel+Dt_daniel
!        endif
!        if(it==NT)write(*,*)
!
!		JDIA=IDINI+INT((IT+HORAINI-1)/(86400./DTP)) !VERIFICA QUAL É O DIA DO CALENDÁRIO 
!		!SUBROTINA DE LEITURA E PREPARACAO DA CHUVA
!		
!!************ DIOGO BUARQUE ********************
!		CALL SED_BACIA
!		!@ Resultado é o aporte de cada fracao de sedimento para a rede
!		!@ PSC(IC,1:4) = 1 = total, 2 = areia, 3 = silte, 4 = argila
!!
!!************ DIOGO BUARQUE ********************
!		CALL SED_REDE
!!************ DIOGO BUARQUE ********************
!	
!		!ARMAZENA DADOS DE VAZÃO DAS CÉLULAS EM QUE EXISTE VAZÃO OBSERVADA
!		DO K=1,NOBS
!			KHID=IQOBS(K) !CÉLULA QUE CORRESPONDE AO POSTO
!			QR(K,IT)=QJ2(KHID) !QR(K,IT) VAI SER COMPARADO A QOBS(K,IT) NA ROTINA FOBJ
!        ENDDO
!
!		IF(ICALIB.EQ.0)THEN !SÓ ARMAZENA ESTES DADOS QUANDO NÃO ESTÁ CALIBRANDO
!!			DO KB=1,NB	!ARMAZENA VAZOES DAS SUB-BACIAS
!!				JCB=IEXUT(KB) !CELULA DO EXUTORIO DA SUB-BACIA KB
!!				QRB(KB,IT)=QJ2(JCB)
!!			ENDDO
!	
!			DO K=1,NUMHIDG !GUARDA DADOS PARA GRAVAR HIDROGRAMAS EM LOCAIS DEFINIDOS NO ARQUIVO PARHIG 
!				KHID=IHIDG(K) !IHIDG(K) É O NÚMERO DA CÉLULA EM QUE SE DESEJA O HIDROGRAMA
!				!QRG(K,IT)=QJ2(KHID) !QRG ARAMZENA OS HIDROGRAMAS NOS LOCAIS DESEJADOS
!!
!!
!!************ DIOGO BUARQUE ********************
!                CAREIA(K,IT)=CSJ2(KHID,1)*(10.**(6.)) !@ (mg/l)
!                CSILTE(K,IT)=CSJ2(KHID,2)*(10.**(6.)) !@ (mg/l)
!                CARGIL(K,IT)=CSJ2(KHID,3)*(10.**(6.)) !@ (mg/l)
!                
!                
!!                if (isNaN(CAREIA(K,IT)) .OR. isNaN(CSILTE(K,IT)) .OR. isNaN(CARGIL(K,IT))) then
!!                    write(*,*) 'IT, IC     = ', IT, IC
!!                    write(*,*) 'CSJ2       = ', CSJ2(KHID,1), CSJ2(KHID,2), CSJ2(KHID,3)
!!                    write(*,*) 'Ca, Cs, Cc = ', CAREIA(K,IT), CSILTE(K,IT), CARGIL(K,IT)
!!                    pause
!!                endif
!!************ DIOGO BUARQUE ********************
!!
!!
!            ENDDO
!!********************* HUGO FAGUNDES 18/08/17 ********************           
!            DO KSS2=1,AUXSS
!                DO K=1,NUMSS(KSS2) !NUMERO DE POSTOS RELACIONADOS AOS SS QUE FORAM INDICADOS PARA SEREM GRAVADOS
!			        KSS=ISS(KSS2,K) !NUMERO DA MINI QUE CORRESPONDE AO POSTO DE SS A SER GRAVADO
!			        CSSCAL(K,IT,KSS2)=(CSJ2(KSS,2)+CSJ2(KSS,3))*(10.**(6.))  !CSSCAL(K,IT,KSS2) EM mg/L VAI SER COMPARADO A SSOBS(K,IT) NA ROTINA FOBJ
!                ENDDO     
!            ENDDO     
!
!!
!!************ DIOGO BUARQUE set/2012 ********************
!			DO i = 1,nSEDmini
!			    KHID=CONTIC(i) !CODIGO DA MINI-BACIA EM QUE SE DESEJA O HIDROGRAMA
!				Qmini(i,IT)=QJ2(KHID) !QRG ARAMZENA OS HIDROGRAMAS NOS LOCAIS DESEJADOS
!			ENDDO
!!************ DIOGO BUARQUE set/2012 ********************			
!!
!!			
!	
!			!ARMAZENA VAZÕES SEGUNDO A ORIGEM PARA UMA CÉLULA - quando não está calibrando
!			MWP=1 !CÉLULA EM QUE SE DESEJAM OS RESULTADOS DE HIDROGRAMA SEPARADO POR ORIGEM
!			
!			MWP=2 !CÉLULA EM QUE SE DESEJAM OS RESULTADOS DE HIDROGRAMA SEPARADO POR ORIGEM
!			
!			QB(IT)=QJ2(MWP)*PJB2(MWP)
!			QBI(IT)=QB(IT)+QJ2(MWP)*PJI2(MWP)
!			QBIS(IT)=QBI(IT)+QJ2(MWP)*PJS2(MWP)
!
!
!        ELSE
!!********************* HUGO FAGUNDES 18/08/17 ********************           
!            DO KSS2=1,AUXSS
!                DO K=1,NUMSS(KSS2) !NUMERO DE POSTOS RELACIONADOS AOS SS QUE FORAM INDICADOS PARA SEREM GRAVADOS
!			        KSS=ISS(KSS2,K) !NUMERO DA MINI QUE CORRESPONDE AO POSTO DE SS A SER GRAVADO
!			        CSSCAL(K,IT,KSS2)=(CSJ2(KSS,2)+CSJ2(KSS,3))*(10.**(6.))  !CSSCAL(K,IT,KSS2) EM mg/L VAI SER COMPARADO A SSOBS(K,IT) NA ROTINA FOBJ
!                ENDDO     
!            ENDDO   
!     ENDIF
!
!!@ DCB ago/2012        write(316666,REC=iT)(QJ2(iC),iC=1,nC)
!!@ DCB ago/2012		write(974,67) (TWS(iC),iC=1,nC)
!!@ DCB ago/2012		write(988,68) (E0media(iC),iC=1,nC)
!!@ DCB ago/2012		write(989,68) (P(iC),iC=1,nC)
!
!!******************************************COMENTADO EM 19/11		
!!		QM2in=QM2
!!		if (IT==1) QM1in=QM2in
!!		do iC=1,nC
!!
!!			IF (hdFLAG(IC)==0) QM2in(iC)=QM2(iC)-(QSUP(IC)+QINT(IC)+QBAS(IC))
!!
!!			if (OD(iC)==1) then
!!				! Minibacia de cabeceira:
!!				DTWS(iC)=P(iC)+0.5*(-QJ1(iC)-QJ2(iC))*DTP/(ACEL(iC)*1000.0)-E0media(iC)
!!			else
!!				DTWS(iC)=P(iC)+0.5*(QM1in(iC)+QM2in(iC)-QJ1(iC)-QJ2(iC))*DTP/(ACEL(iC)*1000.0)-E0media(iC)
!!			endif
!!			! Tira valores absurdos em minibacias muito pequenas:
!!			if (ACEL(iC)<=10.0) DTWS(iC)=0.0
!!
!!		enddo
!!		QM1in=QM2in
!
!!@ DCB ago/2012		write(975,'(2F12.4)') Pbacia,Ebacia
!
!	
!	ENDDO !FIM DO LOOP DO TEMPO
!
!!************ DIOGO BUARQUE ********************
!
!75	FORMAT(I6,8F10.4)
!
!
!66	FORMAT(<nC>F10.4)
!67	FORMAT(<nC>F10.1)
!68	FORMAT(<nC>F10.3)
!69	FORMAT(<nC>F12.3)
!
!     
!	RETURN
!	END