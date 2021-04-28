	!---------------------------------------------------------------------------------
    !  Discussion:
    ! 
    !    This sub-routine generates initial conditions to MGB-IPH model (Inertial version).
    !
    !
    ! Usage:
    !
    ! *CALL HUNT or FINT
    !
    ! uses modules, functions, and subroutines
    !
    ! * USE VARS_MAIN
    ! * USE VARS_INERC (only to Inertial version)
    !
    ! opens
    !
    ! * no files are opened in this routine
    !
    ! reads
    !
    ! * no files are read in this routine
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

Subroutine flood_ini

!-------------------------------------------------------------------------------------
! Variables and Parameters:
use VARS_INERC
use VARS_MAIN
implicit none

integer ii
real*8 FINT,Vol2_ini,ALFA, BETA
!-------------------------------------------------------------------------------------
dtflood0=DBLE(DTP)

do iC=1,nC

        ALFA = (nMan(iC)*(DBLE(BRIO(iC))**(2.0/3.0)))**0.6/((DECL(iC)**0.5)**0.6) !Parameter ALFA: ALFA*(Eq. Manning)^BETA
        BETA = 0.6 !Parameter BETA: ALFA*(Eq. Manning)^BETA
        Q2fl(iC)=QREF(IC) !primeira vazão (base)
        Q2face=QREF(IC)      
        Hfl(iC)=ALFA*(Q2fl(iC)**BETA)/DBLE(BRIO(iC)) + 0.00001
        Yfl(iC)=ZTAB(1,iC)+Hfl(iC)
        Vol1(iC)=Hfl(iC)*DBLE(BRIO(iC))*(SRIO(iC)*1000.)
        Vol2(iC)=Vol1(iC)
        Area2(iC) = FINT(VTAB(:,iC),ATAB(:,iC),NPFL(IC)+2,Vol2(iC))
        Area2(iC) = mIN(Area2(iC),ACEL(iC))
        jtab(iC)=1  
enddo 

endsubroutine