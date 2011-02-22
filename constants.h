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

#ifndef CONSTANTS_H
#define CONSTANTS_H

double F_CORIOLIS;	/*coriolis parameter*/
FILE *out_runtime;
 
enum month_type {JAN, FEB, MAR, APR, MAY, JUN, JUL, AUG, SEP, OCT, NOV, DEC};

enum coord_type {X, Y, Z};

enum berg_idx_type {TOP, SIDES, BASE, FRONT_BACK, NO_BASAL};
			//		  0	   1	 2		3			 4
#endif

/*header file to go with berg1206.c */
#ifndef PI
#define PI 3.14159265
#endif

#ifndef g
#define g 9.81		//m/s^2
#endif


#ifndef OMEGA_0
#define OMEGA_0  0.2617 //rads per hour7.27E-5 //radians per second6.283185/*angular rotation of the Earth in radians per day*/
#endif

#ifndef RHO_ICE
#define RHO_ICE 850.0	  /* density of iceberg in kg/m3	*/
#endif
 
#ifndef TEMP_ICE
#define TEMP_ICE -4.0	  /* temp of ice in deg C	*/
#endif

#ifndef RHO_W
#define RHO_W 1027.0		/* density of sea water in kg/m3 */
#endif

#ifndef RHO_A
#define RHO_A 1.251		/* density of air in kg/m3 */	 
#endif

#ifndef RHO_S
#define RHO_S 	917.0	/* density of sea ice in kg/m3 */
#endif

#ifndef C_W
#define C_W 0.9		/* drag coefficient of sea water*/
#endif

#ifndef C_A
#define C_A 1.3	/* drag coefficient of air */
#endif

#ifndef C_S
#define C_S 0.9		/* drag coefficient of sea ice */
#endif

#ifndef C_surf
#define C_surf 0.002		/* drag coefficient of surface */
#endif

#ifndef K_AIR					  
#define K_AIR 0.0249		/* thermal conductivity of air W/(m deg C) at 283K from http://users.wpi.edu/~ierardi/FireTools/air_prop.html*/
#endif

#ifndef MU_AIR						
#define MU_AIR 1.460E-5		/* kinematic viscosity of air at 283K m^2/s from http://users.wpi.edu/~ierardi/FireTools/air_prop.html*/
#endif

#ifndef KAPPA_AIR			
#define KAPPA_AIR 2.160E-5		/* Thermal diffusivity of air at 283K m2/s from http://users.wpi.edu/~ierardi/FireTools/air_prop.html*/
#endif

#ifndef GAMMA_ICE			
#define GAMMA_ICE 3.34E5		/* Latent heat of melting of ice J/kg */
#endif

#ifndef MAX_BERG_SIZE			
#define MAX_BERG_SIZE 2000.		/* Maximum length/width of berg possible */
#endif

#ifndef PERC_BERG_WINTER			
#define PERC_BERG_WINTER 0.2		/* Percentage of annual production of bergs which stay in fjord...*/
#endif

#ifndef SPREAD_X			
#define SPREAD_X 4000.		/* Distance sed spread over when deposited (along fjord direction) */
#endif

#ifndef SPREAD_Y			
#define SPREAD_Y 2000.		/* Distance sed spread over when deposited (across fjord direction)  */
#endif

#ifndef RES_WATER_VEL_SUMMER			
#define RES_WATER_VEL_SUMMER -0.0432		/* Residual water velocity summer */
#endif

#ifndef RES_WATER_VEL_WINTER			
#define RES_WATER_VEL_WINTER -0.008		/* Residual water velocity winter */
#endif

#ifndef	P_calving
#define	P_calving 0.25 /* Probability of calving each timestep */ 
#endif





