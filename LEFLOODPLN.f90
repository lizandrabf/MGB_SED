    !---------------------------------------------------------------------------------
    !  Discussion:
    ! 
    !    This sub-routine reads the level x flooded area from PrePro (Inertial version).
    !
    !
    ! Usage:
    !
    ! *
    !
    ! uses modules, functions, and subroutines
    !
    ! * USE VARS_MAIN
    ! * USE VARS_INERC (only to Inertial version)
    !
    ! opens
    !
    ! * COTA_AREA.FLP
    !
    ! reads
    !
    ! * COTA_AREA.FLP
    !
    ! creates
    !
    ! * no files are created in this routine
    !
    !---------------------------------------------------------------------------------
    !  Licensing:
    !
    !   This program is free software: you can redistribute it and/or modify
    !   it under the terms of the GNU General Public License as published by
    !   the Free Software Foundation, either version 3 of the License, or
    !   (at your option) any later version.
    !
    !   This program is distributed in the hope that it will be useful,
    !   but WITHOUT ANY WARRANTY; without even the implied warranty of
    !   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    !   GNU General Public License for more details.
    !
    !   You should have received a copy of the GNU General Public License
    !   along with this program.  If not, see <http://www.gnu.org/licenses/>
    !
    !  Version/Modified: 
    !
    ! 2015.07.06 - 07 July 2015 (by Paulo Pontes)
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
    !
    !  Main Reference:
    !
    !    Walter Collischonn,
    !    Modelo de Grandes Bacias - Thesis
    !    Porto Alegre, 2001
    !    ISBN: XXXXXXXXXXX,
    !
    !---------------------------------------------------------------------------------
    ! Variables and Parameters:
    ! *Variables declarations and routines calls are all commented below.
    !---------------------------------------------------------------------------------
    ! End of header
    !---------------------------------------------------------------------------------
    
    SUBROUTINE LEFLOODPLN

    !  Variables and Parameters:
    USE VARS_MAIN
    USE VARS_INERC
    IMPLICIT NONE
    integer aux
    
    !Opening the level x flooded area file
    OPEN(41,FILE=INPUT_DIRECTORY // 'COTA_AREA.FLP',STATUS='OLD',action='READ')
    
    ! Reading data:
    READ(41,*)
    nPfl=0
    DO WHILE(.NOT.EOF(41))
	    READ(41,*) IC,ZFUNDOFL(IC),ZFL(NPFL(IC)+1,IC),AFL(NPFL(IC)+1,IC)
	    NPFL(iC)=NPFL(IC)+1
	ENDDO
    CLOSE(41)
   
    !ABRE ARQUIVO DE FACES
    if (1==0) then
    OPEN(42,FILE=INPUT_DIRECTORY // 'face.con',STATUS='OLD',action='READ')
    
    !Número da face / Minibacia 1 / Minibacia 2 / Ymin da fronteira (m) / DX entre as minibacias (m):
    aux=0
    DO WHILE(.NOT.EOF(42))
        aux=aux+1
	    READ(42,*) aux,nFACECAT1(aux),nFACEY1(aux),nFACECAT2(aux),nFACEY2(aux),nFACEDX(aux)
	ENDDO
    CLOSE(42)
    nFACE=aux
    endif

    RETURN
    END