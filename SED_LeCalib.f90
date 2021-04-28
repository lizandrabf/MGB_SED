	subroutine LeCalib_SED
	! Subrotina para leitura do arquivo ParCalibSED.txt
	!-----------------------------------------------------------------------
	!
	! Descrição das variáveis locais:
	!
	! I,J = variáveis auxiliares
	!
	!----------------------------------------------------------------------
	!
	use VARS_MAIN
	use VARS_CALIB
    use SED_VARS
	IMPLICIT NONE
	! Declaração de variáveis:
	! Variáveis locais:
	integer:: i,j,k,KB,j2
	character(60):: trashTex
	character(10):: trashTex2
	!----------------------------------------------------------------------------

    WRITE(*,*)'Le arquivo de parametros de calibracao:'

	OPEN(FILCALIB,FILE='.\input\ParCalibSED.txt',STATUS='OLD',ACTION='READ')

	READ(FILCALIB,75) trashTex
	WRITE(*,*)trashTex
	READ(FILCALIB,*) NS !NUMERO DE INDIVÍDUOS
	WRITE(*,*)NS
	Write(*,*)

	READ(FILCALIB,75) trashTex
	WRITE(*,*)trashTex
	READ(FILCALIB,*) NF !NÚMERO DE FUNÇÕES OBJETIVOS
	WRITE(*,*) NF
	Write(*,*)

	READ(FILCALIB,75) trashTex
	WRITE(*,*)trashTex
	READ(FILCALIB,*) iMaxGen
	WRITE(*,*)iMaxGen
	Write(*,*)
    
    READ(FILCALIB,75) trashTex
	WRITE(*,*)trashTex
	READ(FILCALIB,*) NPAR   !NÚMERO DE PARÂMETROS
	WRITE(*,*)NPAR      
    Write(*,*)
    
    READ(FILCALIB,75) trashTex
	!WRITE(*,*)trashTex
	ALLOCATE (PMIN(NPAR),PMAX(NPAR)) 
	DO i=1,NPAR
		READ(FILCALIB,*) PMIN(i),PMAX(i) !LIMITES DOS PARÂMETROS
		!WRITE(*,*) PMIN(i),PMAX(i)
	ENDDO
	Write(*,*)
	PAUSE
    
	READ(FILCALIB,75) trashTex
	WRITE(*,*)trashTex
    KB=NOBSS(AUXCALIBSS)
    ALLOCATE (CalibFLAG(NBP,nf))
    i=1
    READ(FILCALIB,*) (calibFLAG(i,j),j=1,NF)
    WRITE(*,*) calibFLAG(1,1),calibFLAG(1,2),calibFLAG(1,3)
    
	DO i=2,KB
        DO j=1,NF
            calibFLAG(i,j)=calibFLAG(1,j)
        ENDDO
	ENDDO
    PAUSE

    CLOSE (FILCALIB)

71	FORMAT(6I10)
72	FORMAT(5A10)
73	FORMAT(A10,7I6)
74	FORMAT(A20)
75	FORMAT(A60)
76	FORMAT(I10,F10.1)
77	FORMAT(A20)
78	FORMAT(A10,1I6)
!79	FORMAT(A10,2F6)

	RETURN
	END