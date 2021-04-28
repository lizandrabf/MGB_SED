
SUBROUTINE SED_CALIB_REDE

!@ **********************************
USE VARS_MAIN
USE SED_VARS
USE VARS_INERC
USE VARS_CALIB

IMPLICIT NONE

!@ **********************************
!@ Declaração de variáveis:
REAL HmJx, VmJx, RHmJx, SfTRECx, UatJx, Qlimx
REAL, ALLOCATABLE:: DEPaux2(:,:)
REAL QtFL, HFLSD

ALLOCATE (DEPaux2(NC,3))

!@ *****************************************************************************************
!@ EM CADA INTERVALO DE TEMPO O QUE ERA i+1 VIRA i
!@ -----------------------------------------------------------------------------------------
DO IC=1,NC
	CSM1(IC,:)      = CSM2(IC,:)
	CSJ1(IC,:)      = CSJ2(IC,:)
	CFL1(IC,:)      = CFL2(IC,:)
ENDDO

!@ INICIALIZANDO VARIÁVEIS
iSEDaux = 0     !@ Contador de minibacia com cálculo de sedimentos
CARGM   = 0.0   !@ Carga a montante do trecho
CSJ2    = 0.0
Qlimx   = 0.01
CFL2    = 0.0
DEPaux2 = 0.0 !Hugo Fagundes 19/08/2019

DO IC = 1,NC    !@ INICIO DO LOOP DAS MINIBACIAS
    IB=IBAC(IC)
    QtFL = 0.0  !@ Vazão de troca total na minibacia
    !@ CALCULA APENAS SUB-BACIAS SELECIONADAS
	IF(IB>SUBfim .OR. IB<SUBini)CYCLE 
    
    iSEDaux = iSEDaux + 1 !@ Contador de minibacia com cálculo de sedimentos
    
    VmJx = 0.0  ! Estou usando Meyer Peter e Muller que não usa essa variável, por isso ela pode ser zero
    UatJx = 0.0 ! Estou usando Meyer Peter e Muller que não usa essa variável, por isso ela pode ser zero
    
    !Atribui valores das variáveis usadas na propagação do modelo de sedimentos
    IF (IT==1) THEN
        VOLTREC1(IC)=0.0
        VFL1(IC)=0.0
        QM1(IC)=0.0
        QJ1(IC)=0.0
    ELSE
        VOLTREC1(IC)=VOLTREC2_SC(IC,IT-1)
        VFL1(IC)=VFL2_SC(IC,IT-1)
        
        !Sum of upstream IC flows
        Nentradas = MINIMONT(iC,1)
        IF(Nentradas==0)THEN
            QM1(IC)=0.0
        ELSE
            QM1(IC) = SUM(Q_SC(MINIMONT(IC,2:1+Nentradas),IT-1))    
        ENDIF
        
        QJ1(IC)=Q_SC(IC,IT-1)
    ENDIF
    
    QJ2(IC)=Q_SC(IC,IT)
    
    !Sum of upstream IC flows
    Nentradas = MINIMONT(iC,1)
    IF(Nentradas==0)THEN
        QM2(IC)=0.0
    ELSE
        QM2(IC) = SUM(Q_SC(MINIMONT(IC,2:1+Nentradas),IT))
    ENDIF
 
    VOLTREC2=VOLTREC2_SC(IC,IT)
    VFL2(IC)=VFL2_SC(IC,IT)
    SfTRECx=SfTRECx_SC(IC,IT)
    HmJx=HmJx_SC(IC,IT)
    HFLSD=HFLSD_SC(IC,IT)

    !@ Vazão de troca:
    QtFL = (VFL2(iC)-VFL1(iC))/DTP 
    CALL SED_PROPAG(HmJx, VmJx, SfTRECx, UatJx, Qlimx, QtFL, VFL1(IC), VFL2(IC), HFLSD, DEPaux2)

ENDDO   !@ FIM DO LOOP DAS MINIBACIAS

DEALLOCATE (DEPaux2)

RETURN
END SUBROUTINE