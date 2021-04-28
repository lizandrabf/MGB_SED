!@ *********************************************************************************************
!@ Criado por Diogo Buarque
!@ Data: Março de 2011
!@
!@ Atualizado: Abr 2011
!@
!@ >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MODULO SED_VARS >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!@
!@-------------------------------------------------------------------------------------
!@  VARIÁVEIS
!@ Kusle	= is the USLE soil erodibility factor [(0.013*ton*m2*hr)/(m3*ton*cm)]
!@ Culse	= is the USLE cover and management factor
!@ Pusle	= is the USLE support practice factor
!@ Rgros	= is the coarse fragment factor
!@
!@ Mareia	= is the percent sand content (0.05-2.00 mm diameter particles)
!@ Msilte	= is the percent silt content (0.002-0.05 mm diameter particles)
!@ Margila	= is the percent clay content (< 0.002 mm diameter particles)
!@ Morg		= is the percent of organic mater in the soil
!@ Mrocha	= is the percent rock in the first soil layer
!@
!@ LSAcu	= somatório do fator LS de cada uso em cada minibacia
!@ PSC		= aporte de sedimentos da minibacia para o rio (ton)
!@ PSU		= aporte de sedimentos dos blocos da minibacia para o rio (ton)
!@ SLC		= perda de solo total da minibacia (ton)
!@ SLU		= perda de solo total do bloco da minibacia (ton)
!@ QSC		= descarga sólida das minibacias (ton/s) no passo de tempo
!@ QSU      = descarga sólida dos blocos das minibacias (ton/s)
!@ QSAUX	= variável que armazena, para o escoamento superficial de cada minibacia, o tempo de retardo (s),
!@           o volume total (m3) e as lâminas (mm) para cada bloco
!@ NPXU		= número de pixels de cada uso em cada minibacia
!@ SDR      = valores do SDR de cada bloco de cada minibacia e de cada minibacia
!@
!@ DMa		= diâmetro característico para as partículas de areia
!@ DMs		= diâmetro característico para as partículas de silte
!@ DMg		= diâmetro característico para as partículas de argila
!@
!@ APIX     = área média dos pixels de uma minibacia (km2)
!@ FCTE     = fator constante da MUSLE para cada uso da minibacia
!@ TPIC     = taxa de pico do escoamento superficial (m3/s)
!@ SLX      = carga e sedimentos remanescente na minibacia no passo de tempo
!@ SLXU     = carga e sedimentos remanescente por bloco da minibacia no passo de tempo
!@
!@ RUGMAN   = rugosidade de Manning para MC
!@ VISC     = viscosidade cinemática da água (m2/s)
!@
!@ PSCaux   = aporte das frações de sedimentos que chegam à drenagem das minibacias das subbacias selecionadas
!@ DMP      = diâmetro nominal das partículas (m)
!@ WSP      = velocidade de queda das partículas (m/s)
!@ CTS      = capacidade de transporte do trecho (ton/m3)
!@ Fform    = fator de forma das partículas de sedimento
!@
!@ CONCTREC1= concentração de sedimentos no trecho no tempo 1 (ton/m3)
!@ CONCTREC2= concentração de sedimentos no trecho no tempo 2 (ton/m3)
!@ VolTREC1 = volume médio de água no trecho no tempo 1 (m3)
!@ VolTREC2 = volume médio de água no trecho no tempo 2 (m3)
!@
!@ DEPTREC  = carga acumulada de sedimento depositado no trecho (ton)
!@ EROSTREC = carga acumulada de sedimento erodido no trecho (ton)
!@ CARGM    = carga de sedimento saindo para o trecho de jusante(ton)
!@ CSdep    = quantidade de sedimentos depositada no trecho expressa em forma de concentração (ton/m3) !Hugo Fagundes 07/08/2017
!@
!@ HTR1     = PROFUNDIDADE NA SEÇÃO DE JUSANTE DO TRECHO NO TEMPO T
!@ HTR2     = PROFUNDIDADE NA SEÇÃO DE JUSANTE DO TRECHO NO TEMPO T+1
!@ VFL1     = VOLUME DE ÁGUA NA PLANÍCIENO TEMPO T
!@ VFL2     = VOLUME DE ÁGUA NA PLANÍCIENO TEMPO T+1
!@ CFL1     = CONCENTRAÇÃO NA PLANÍCIENO TEMPO T
!@ CFL2     = CONCENTRAÇÃO NA PLANÍCIENO TEMPO T+1
!@ DFL      = DEPÓSITO DE SILTE E ARGILA NA PLANÍCIE
!@-------------------------------------------------------------------------------------
!@
!@ *********************************************************************************************

MODULE SED_VARS

!@ JÁ LISTADAS
INTEGER,PARAMETER ::    FILPAR=300 !ARQUIVO DOS PARAMETROS DA MUSLE ASSOCIADOS AO USO
INTEGER,PARAMETER ::    FILTEX=301 !ARQUIVO DAS TEXTURAS DOS SOLOS ASSOCIADOS AO USO
INTEGER,PARAMETER ::    FILLSM=302 !ARQUIVO DOS LS ACUMULADOS DE CADA USO DE CADA MINIBACIA
INTEGER,PARAMETER ::    FILHRU=303 !ARQUIVO COM NUMERO DE PIXELS DE CADA USO DE CADA MINIBACIA
INTEGER,PARAMETER ::    FILSED=304 !ARQUIVO COM CARGA DE SEDIMENTOS DAS MINIBACIAS QUE CHEGA À REDE DE DRENAGEM
INTEGER,PARAMETER ::    FILSEDa=305 !ARQUIVO COM CARGA DE SEDIMENTOS DAS MINIBACIAS QUE CHEGA À REDE DE DRENAGEM
INTEGER,PARAMETER ::    FILSEDs=306 !ARQUIVO COM CARGA DE SEDIMENTOS DAS MINIBACIAS QUE CHEGA À REDE DE DRENAGEM
INTEGER,PARAMETER ::    FILSEDc=307 !ARQUIVO COM CARGA DE SEDIMENTOS DAS MINIBACIAS QUE CHEGA À REDE DE DRENAGEM
INTEGER,PARAMETER ::    FILQCE=308 !ARQUIVO COM A VAZÃO SUPERFICIAL DAS MINIBACIAS QUE CHEGA À REDE DE DRENAGEM
INTEGER,PARAMETER ::    FILCAR=309 !306 !ARQUIVO COM A CONCENTRAÇÃO DE AREIA DOS EXUTORIOS DESEJADOS(mg/L)
INTEGER,PARAMETER ::    FILCSI=310 !307 !ARQUIVO COM A CONCENTRAÇÃO DE SILTE DOS EXUTORIOS DESEJADOS(mg/L)
INTEGER,PARAMETER ::    FILCCL=311 !308 !ARQUIVO COM A CONCENTRAÇÃO DE ARGILA DOS EXUTORIOS DESEJADOS(mg/L)
INTEGER,PARAMETER ::    RIOCAR=312 !ARQUIVO COM A CONCENTRAÇÃO DE ARGILA DOS EXUTORIOS DAS MINIBACIAS(mg/L)
INTEGER,PARAMETER ::    RIOCSI=313 !ARQUIVO COM A CONCENTRAÇÃO DE ARGILA DOS EXUTORIOS DAS MINIBACIAS(mg/L)
INTEGER,PARAMETER ::    RIOCCL=314 !ARQUIVO COM A CONCENTRAÇÃO DE ARGILA DOS EXUTORIOS DAS MINIBACIAS(mg/L)
INTEGER,PARAMETER ::    FILSDR=315 !ARQUIVO COM A CONCENTRAÇÃO DE ARGILA DOS EXUTORIOS DAS MINIBACIAS(mg/L)
INTEGER,PARAMETER ::    SEDSAI=316 !ARQUIVOS COM A CARGA DE SEDIMENTOS NO TRECHO DA MINIBACIA (ton/dia)
INTEGER,PARAMETER ::    SEDDEP=317 !ARQUIVO COM A CARGA DE SEDIMENTOS DEPOSITADA NO TRECHO DA MINIBACIA (ton/dia)
INTEGER,PARAMETER ::    SEDERO=318 !ARQUIVO COM A CARGA DE SEDIMENTOS ERODIDA NO TRECHO DA MINIBACIA
INTEGER,PARAMETER ::    TAXPIC=319 !ARQUIVO COM O ACUMULADO DAS TAXAS DE PICO DOS BLOCOS DAS MINIBACIAS
INTEGER,PARAMETER ::    FILVAZ=320 !ARQUIVO DOS HIDROGRAMAS DAS MINI-BACIAS Hugo Fagundes 16/07/2019
INTEGER,PARAMETER ::    FILCEL=321 !317 !ARQUIVO COM A VAZÃO DAS MINIBACIAS QUE CHEGA À REDE DE DRENAGEM 2013_02_06
INTEGER,PARAMETER ::    BALFLP=322 !401 !@ DCB_HD_Sed
INTEGER,PARAMETER ::    FILSSCB=323
INTEGER,PARAMETER ::    FILSS=324 !ARQUIVO RELATIVO A DADOS DE SS OBSERVADOS                 !Hugo Fagundes
INTEGER,PARAMETER ::    FILSSB=325 !ARQUIVO RELATIVO A DADOS de CSS DE BASE DOS POSTOS        !Hugo Fagundes
INTEGER,PARAMETER ::    RIODEP=326 !ARQUIVO COM A QUANTIDADE DE SEDIMENTO DEPOSITADO EM CADA DIA EM CADA MINIBACIA (ton/dia)   
INTEGER,PARAMETER::     FILAJUS=327 !ARQUIVO DE AJUSTE PARA OS DADOS DE SEDIMENTOS       !HUGO FAGUNDES 19/08/17
INTEGER,PARAMETER::     FILRESCALIB=328 !ARQUIVO DE AJUSTE PARA OS DADOS DE SEDIMENTOS       !HUGO FAGUNDES 27/09/17
INTEGER,PARAMETER::     FILDSUP=329 !ARQUIVO PARA SALVAR DSUP MÉDIO DE CADA MINIBACIA PARA CADA DIA       !HUGO FAGUNDES 19/09/2019
INTEGER,PARAMETER::     FILCusle=330 !ARQUIVO PARA SALVAR DSUP MÉDIO DE CADA MINIBACIA PARA CADA DIA       !HUGO FAGUNDES 19/09/2019
INTEGER,PARAMETER::     FILKusle=331 !ARQUIVO PARA SALVAR DSUP MÉDIO DE CADA MINIBACIA PARA CADA DIA       !HUGO FAGUNDES 19/09/2019

CHARACTER               (10) Nuso(12) !NOMES DOS USOS DO SOLO (BLOCOS)
REAL                    DMa, DMs, DMg
REAL                    APIX, FCTE, QPIC, SLX
REAL                    VISC
REAL,ALLOCATABLE ::     Kusle(:,:),Cusle(:,:),Pusle(:,:),Rgros(:,:), Ksdr(:,:) !@ Ksdr(:) DCB 02-07-2011
REAL,ALLOCATABLE ::     Mareia(:,:),Msilte(:,:),Margila(:,:),Morg(:,:),Mrocha(:,:)
REAL,ALLOCATABLE ::     LSAcu(:,:), SDR(:,:) !@ SDR(:,:) DCB 02-07-2011
REAL(8),ALLOCATABLE ::  PSC(:,:),SLC(:),QSC(:),QSAUX(:,:)
REAL,ALLOCATABLE ::     PSU(:,:),QSU(:,:),SLXU(:)
REAL,ALLOCATABLE ::     SLU(:,:)
REAL,ALLOCATABLE ::     DSUPa(:)   !Armazena Dsup médio de cada minibacia para cada dia HUGO FAGUNDES 19/09/2019
REAL(4),ALLOCATABLE ::  PSCaux(:,:)
REAL,ALLOCATABLE ::     DMP(:), WSP(:), CTS(:,:), Fform(:)
REAL,ALLOCATABLE ::     CSM1(:,:), CSM2(:,:), CSJ1(:,:), CSJ2(:,:)
REAL,ALLOCATABLE ::     VolTREC1(:), VolTREC2(:)
REAL(8),ALLOCATABLE ::  DEPTREC(:,:), DEPT(:,:), EROSTREC(:,:), EROT(:,:), CARGM(:,:)
REAL(4), ALLOCATABLE :: CSdep(:,:) !Csdep - Hugo Fagundes em 07/08/17
REAL,ALLOCATABLE ::     QSJ2(:,:)
INTEGER,ALLOCATABLE ::  NPXU(:,:)
REAL,ALLOCATABLE ::     CARGAaux(:)
REAL,ALLOCATABLE ::     CAREIA(:,:)	!ARMAZENA CONCENTRACOES DE AREIA ONDE SE DESEJA GRAVAR
REAL,ALLOCATABLE ::     CSILTE(:,:)	!ARMAZENA CONCENTRACOES DE SILTE ONDE SE DESEJA GRAVAR
REAL,ALLOCATABLE ::     CARGIL(:,:)	!ARMAZENA CONCENTRACOES DE ARGILA ONDE SE DESEJA GRAVAR
REAL,ALLOCATABLE ::     CSSCAL(:,:,:),SEDstore(:,:,:)	!ARMAZENA CONCENTRACOES DE SEDIMENTOS EM SUSPENSÇAO QUE SERÃO COMPARADOS COM OS DADOS OBSERVADOS DE SSOB E ARMAZENA SEDIMENTOS EM LOCAIS DESEJADOS   !HUGO FAGUNDES 18/08/17
REAL,ALLOCATABLE ::     ALFsed(:),BETsed(:),AuxTKS(:)   !parametros calibraveis da MUSLE    !HUGO FAGUNDES 05/10/17
INTEGER                 iSEDaux, nSEDmini
INTEGER,ALLOCATABLE ::  CONTIC(:),CONTICSUB(:)
INTEGER                 AUXSS,AUXCALIBSS !NUMBER OF DATASETS RELATED TO SEDIMENT EVALUATION AND DEFINE WHAT VARIABLE WILL BE USED IN SEDIMENT MODEL CALIBRATION PROCESS.   1 TO SSC, 2 TO REFLECTANCE, 3 3 TO TURBIDITY, 4 TO TOTAL SOLIDS IN SUSPENSION
INTEGER                 flag_TC,Flag_SEDVAR ! FLAGS TO IDENTIFY WHAT TRANSPORT CAPACITY EQUATION AND VARAIBLE (SSC, QSS) WILL BE USED FOR SIMULATION AND COMPARISON WITH OBSERVED DATA
INTEGER,ALLOCATABLE::   NOBSS(:),ISSOBS(:,:) !NUMERO DE POSTOS RELATIVOS A SS, numero das células com dados observados
INTEGER,ALLOCATABLE::   NUMSS(:),ISS(:,:) !NUMERO DE POSTOS RELATIVOS A SS QUE SE DESEJA GRAVAR, NUMERO DAS CÉLULAS EM QUE SE DESEJA OS SS
CHARACTER (20)          ARQOBSSS(20)           !NOME DO ARQUIVO RELATIVO AOS SEDIMENTOS EM SUSPENSÃO
INTEGER                 SED_INDEX !Flag for sediment model      !Hugo Fagundes 28.06.19
INTEGER::               auxcont !variable to comput Slope Sed. duration curve in Fob_SedCalib
INTEGER::               xyz !Variável auxiliar para calibração
REAL,ALLOCATABLE::      SSCORRELesp(:), SSCORRELtudo(:), SSCORRELtemp(:) !CORREL ESPACIAL, DE TUDO E TEMPORAL                                       !Hugo Fagundes
REAL,ALLOCATABLE::      RMSE(:), CoefB(:), DCPERM(:), R2SS(:), ERRVS(:),SSCORRELAC(:),SSKGE(:),KGEa(:) !RMSE, Beta Coefficient, Slope Sed. duration curve, NASH SS, VOL. SS, CORREL TEMPORAL, KLING GUPTA, alfa Kling Gupta !Hugo Fagundes
REAL,ALLOCATABLE::      AUXOBS(:),AUXCALC(:), CSSbase(:,:) !ARMAZENAM OS VALORES MÉDIOS DE CADA POSTO PARA OS DADOS OBSERVADOS E CALCULADOS
REAL,ALLOCATABLE::      SSOBS(:,:,:) !DADOS RELATIVOS A SS  OBSERVADOS EM CADA MINI        !@HUGO FAGUNDES 18/08/17
INTEGER                 NPARsed, NSsed !Numero de parametros e numero de individuos calibração automática sedimentos
!REAL,ALLOCATABLE::      RESULTSS(:,:) !Results FO sediment calibration

!@ NÃO LISTADAS AINDA
REAL,ALLOCATABLE ::     FracS(:,:)

REAL,ALLOCATABLE ::     VFL1(:), VFL2(:)        !@ DCB_HD_Sed
REAL,ALLOCATABLE ::     CFL1(:,:), CFL2(:,:)    !@ DCB_HD_Sed
REAL,ALLOCATABLE ::     DFL(:,:)                !@ DCB_HD_Sed

!@ ------------------------------------------------------------------------
!@ BALANÇO DE SEDIMENTOS
REAL(8),ALLOCATABLE ::     iENTRA(:,:), iSAI(:,:)  !@ DCB_HD_Sed   ,iERRO(:,:)
REAL,ALLOCATABLE ::     inFL(:,:), outFL(:,:), BalQFL(:,:), VTOTAL(:,:)  !@ DCB_HD_Sed
!@ ------------------------------------------------------------------------




END MODULE