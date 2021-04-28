      SUBROUTINE caldat(julian,mm,id,iyyy)
	  !subrotina que converte dia juliano em dia/mes/ano

!---------------------------------------------------------------------------------
! caldat.f90
!---------------------------------------------------------------------------------
! Discussion:
!
! This function determines the day/mounth/year of the Julian Calendar corresponding to the given Julian day
! Does the opposite of julday.f90
!
! Usage:
!
! CALL caldat(julian,mm,id,iyyy)
!
! uses modules, functions, and subroutines
!
! * no modules, functions, or subroutines are used in this funcation
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
! 2014.11.12 - 12 November 2014
! By: Fernando Mainardi Fan
!
! Authors:
!
! Original fortran version by Walter Collischonn
! Present fortran version by:
! * Walter Collischonn
! * Rodrigo Cauduro Dias de Paiva
! * Diogo da Costa Buarque
! * Paulo Pontes Rógenes
! * Mino Viana Sorribas
! * Fernando Mainardi Fan
! * Juan Martin Bravo
!
! Main References:
! COLLISCHONN, W. ; ALLASIA, D. G. ; SILVA, B. C. ; TUCCI, C. E. M. The MGB-IPH model for large-scale rainfall-runoff modelling. Hydrological Sciences Journal, v. 52, p. 878-895, 2007.
! COLLISCHONN, W., TUCCI, C. E. M. Simulação hidrológica de grandes bacias. Revista Brasileira de Recursos Hídricos, v. 6, n. 2, 2001.
! COLLISCHONN, W. Modelagem de Grandes Bacias - ph.d. Thesis. 2001
!
!---------------------------------------------------------------------------------
! Variables and Parameters:
! *Variables declarations and routines calls are all commented below.
!---------------------------------------------------------------------------------
! End of header
!---------------------------------------------------------------------------------
	
!Function:

	  implicit none
      INTEGER id,iyyy,julian,mm,IGREG
      PARAMETER (IGREG=2299161)
      INTEGER ja,jalpha,jb,jc,jd,je
      if(julian.ge.IGREG)then
        jalpha=int(((julian-1867216)-0.25)/36524.25)
        ja=julian+1+jalpha-int(0.25*jalpha)
      else
        ja=julian
      endif   
      jb=ja+1524
      jc=int(6680.+((jb-2439870)-122.1)/365.25)
      jd=365*jc+int(0.25*jc)
      je=int((jb-jd)/30.6001)
      id=jb-jd-int(30.6001*je)
      mm=je-1

      if(mm.gt.12)mm=mm-12
      iyyy=jc-4715
      if(mm.gt.2)iyyy=iyyy-1
      if(iyyy.le.0)iyyy=iyyy-1
      
      Return
      END
      
!End of function
