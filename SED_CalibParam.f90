subroutine CalibParam_SED
	! Subrotina para passar parametros do MOCOM-UA para MGB-SED
	!-----------------------------------------------------------------------
	!
	!
	USE VARS_MAIN
	USE VARS_CALIB
    USE SED_VARS
    
	IMPLICIT NONE
	! Declaração de variáveis:
	! Variáveis locais:
	integer:: k
    
    ALFsed(1)=PARX(1)          
    BETsed(1)=PARX(2)          
    AuxTKS(1)=PARX(3) 

    DO k=1,NB-1
        ALFsed(k+1)=PARX(3*k+1)       
        BETsed(k+1)=PARX(3*k+2)       
        AuxTKS(k+1)=PARX(3*k+3)    
    ENDDO

    !ALFsed(1:6)=PARX(1)          
    !BETsed(1:6)=PARX(2)          
    !AuxTKS(1:6)=PARX(3)          
    !ALFsed(7:8)=PARX(4)          
    !BETsed(7:8)=PARX(5) 
    !AuxTKS(7:8)=PARX(6)       
    !ALFsed(9:10)=PARX(7)         
    !BETsed(9:10)=PARX(8)      
    !AuxTKS(9:10)=PARX(9)       
    !ALFsed(11:14)=PARX(10)        
    !BETsed(11:14)=PARX(11)     
    !AuxTKS(11:14)=PARX(12)       
    !ALFsed(15:17)=PARX(13)        
    !BETsed(15:17)=PARX(14)
    !AuxTKS(15:17)=PARX(15)       
    
	RETURN
	END