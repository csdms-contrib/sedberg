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

#include <math.h>
#include <stdlib.h>
#include <stdio.h> 
#include "constants.h"

#ifndef NR_END
#define NR_END 1
#endif
#ifndef	FREE_ARG
#define FREE_ARG char*
#endif

typedef struct 
{
	double latitude; /*latitude of fjord in degrees */ 
	double h_front; /*height of ice front in metres  */  
	double w_front; /*width of ice front / fjord in metres*/
	double w_fjord; /*width of fjord in metres*/
	double v_front;  /* velocity of front in m/day    	  */
	double h_water; /*height of water   						*/
	double h_basal_sediment; /*thickness of basal sediment layer*/    
	double C_basal_i;  /*initial basal concentrations of each particle fraction*/    
	double C_englacial_i;
	double density_grain;
	double density_dep;
	double y_f; /* position of ice front (m)	   */
	double slope;  /*surface slope in degrees   */
	double bottomDepth;  /*height descended by glacier from edge of program to snout	 */
	double BrownsConst; /*constant in Brown's empirical calving equation ~ 27 per year*/
	double initial_snout_position;	/**/
	int calving_scenario; //SeasonalCalvingScenario:0=6timesinsummer,1=summer_only,2=sikussakbreakup(all_at_once),3=const_all_year
	double 	annual_berg_vol; //calving rate in km3/yr	 
	double 	mu; //mean of normal distribution for lognormal  iceberg size distribution	
	double 	sigma; //standard deviation of normal distribution for lognormal  iceberg size distribution	
	
}glacier;

typedef struct 
{
	long no_x_bins;
	long no_y_bins;
	long no_z_bins;
	int noLines;	/*number of lines in 'bathy.txt' - starts at line 0 */
	double bottom_x[100];	  	 
	double bottom_z[100]; 
	double max_depth;
	double max_length;
			 
}bathy;

typedef struct 
{
	double amplitude;
	double Z_var;
	double T_var;
	double delta_v;
	double eta_v;
				 
}water_vel;


typedef struct 
{
	double totnyears;
	double spinuptime;
	double dx;
	double dz; 
	int print_sediment;			
}options;

typedef struct 
{
	int	water_z_min;
	double **thickness;
				
}sediment;
	
glacier inputGlacier();
options inputOptions();
bathy inputBathy(glacier, options);
water_vel inputWater_vel();
sediment initialiseSediment(bathy, options, glacier);
void plumeThickness(options, bathy, glacier, int, int,double **,double **);
void gcThickness(options, bathy, glacier, int, int,double **,double **);
void inputSolarRadiation(double []);
void inputAirTemp(double []);
void inputSeaTempSummer(double *, bathy, options, sediment);
void inputSeaTempWinter(double *, bathy, options, sediment);
void inputSeaSalin(double *, bathy, options, sediment);
void calcMonth(double, int *);
void freezingpointwater(double,double *, double);


