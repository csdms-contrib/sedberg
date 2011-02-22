////////////////////////////////////////////////////////////////////////////////////////
//	   SedBerg: simulates iceberg drift melt and sedimentation
//     Copyright (C) 2010  Ruth I Mugford
// 
//     This program is free software; you can redistribute it and/or modify
//     it under the terms of the GNU General Public License, version 2, as 
//	   published by the Free Software Foundation.
// 
//     This program is distributed in the hope that it will be useful,
//     but WITHOUT ANY WARRANTY; without even the implied warranty of
//     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//     GNU General Public License for more details.
// 
//     You should have received a copy of the GNU General Public License along
//     with this program; if not, write to the Free Software Foundation, Inc.,
//     51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
//     
//     Author:	R. I. Mugford, Scott Polar Research Institute, University of Cambridge, 
//     Lensfield Road, Cambridge, CB2 1ER, U.K.
//     R.I.Mugford@googlemail.com
//   
//	   Please cite the following manuscript in any publications about the use of or the
//     development of SedBerg:
//     Mugford, R. I., and J. A. Dowdeswell (2010), Modeling iceberg-rafted sedimentation
//     in high-latitude fjord environments, J. Geophys. Res., 115 (F03024), 
//     doi:10.1029/2009JF001564.  http://dx.doi.org/10.1029/2009JF001564
////////////////////////////////////////////////////////////////////////////////////////

#ifndef	IA
#define IA 16807
#endif

#ifndef	IM
#define IM 2147483647
#endif

#ifndef	AM
#define AM (1.0/IM)
#endif

#ifndef	IQ
#define IQ 127773
#endif

#ifndef	IR
#define IR 2836
#endif

#ifndef	NTAB
#define NTAB 32
#endif

#ifndef	NDIV
#define NDIV (1+(IM-1)/NTAB)
#endif

#ifndef	EPSRAN
#define EPSRAN 1.2e-7
#endif

#ifndef	RNMX 
#define RNMX (1.0-EPSRAN)
#endif

float ran(long *);
float ran1(long *);
float ran2(long *);
float ran3(long *);
float ran4(long *);
float ran5(long *);
float ran6(long *);
float ran7(long *);
float ran8(long *);
float ran9(long *);
float ran10(long *);
float ran11(long *);
float ran12(long *);
float f_ran_stayinfjord(long *);


float gausdev2(long *);
float gausdev4(long *);
float gausdev5(long *);



