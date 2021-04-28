  
	SUBROUTINE SED_CALIB_MODELO

	USE VARS_MAIN
	USE VARS_INERC
    USE SED_VARS 
	USE VARS_CALIB
    
    IMPLICIT NONE
	
	INTEGER k,KC,KB,KSS,KSS2 !indexes and counters

    ! Initialize
	IT=0    !3287 31/12/1998    !3652 31/12/1999 !5114 31/12/2003
        
    DO WHILE (IT<NT)
		
		IT=IT+1
		                 
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

        ! COMPUTES MUSLE SEDIMENT YIELD
        CALL SED_CALIB_BACIA
        
        ! ROUTING SEDIMENT LOAD
        CALL SED_CALIB_REDE

        IF(SED_INDEX==1)THEN
            DO KSS2=1,AUXSS
                DO K=1,NOBSS(KSS2) 
                    KSS=ISSOBS(KSS2,K) 
                    !CSSCAL(K,IT,KSS2)=((CSJ2(KSS,2)+CSJ2(KSS,3))*(10.**(6.)))                             !CSS mg/L VAI SER COMPARADO A SSOBS(K,IT) NA ROTINA FOBJ
                    CSSCAL(K,IT,KSS2)=(CSJ2(KSS,2)+CSJ2(KSS,3))*(10.**(6.))*QJ2(KSS)*0.0864              !QSS ton/dia VAI SER COMPARADO A SSOBS(K,IT) NA ROTINA FOBJ 
                    !CSSCAL(K,IT,KSS2)=(CSJ2(KSS,1)+CSJ2(KSS,2)+CSJ2(KSS,3))*(10.**(6.))*QJ2(KSS)*0.0864  !QStot ton/dia VAI SER COMPARADO A SSOBS(K,IT) NA ROTINA FOBJ 
                ENDDO     
            ENDDO
        ENDIF
      
	ENDDO !End time loop
     
	RETURN
	END
