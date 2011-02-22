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

#ifndef ICEBERG_H
#define  ICEBERG_H

#include "random.h"
#include "fjord.h"
#include "constants.h"

typedef struct 
{

	double tcalving_i;	/* time of calving t_c	*/
	double size_i[3];	/*0. length x	1. width y, 2. height z	*/
	double d_i;	/*4. depth of keel (calculated from volume and densities) */
	double position[2];	/*5.midpoint position  0. x direction along fjord 1. y direction across fjord   */
	double a_basal[4];	/*6. area of basal layer l_b	in x,y,z	*/
	double a_englacial[4];	/*7. area of englacial layer l_e  in x,y,z	 */
	double conc_sed[3];	/*8-12. average concentration of 5 sediment types in 3 melting indices 
							(0=top,1=sides+front+back,2=bottom	 */	
	double v_berg[2];	/*23. velocity of berg 0. x direction along fjord 1. y direction across fjord*/
	double v_mag;
	double R[3];	/*24. melting rate R 0=top;1=sides,front and back 2=base*/
	int basal_position;	/*25 where is basal layer on the iceberg? 1=side,0=top,2=bottom	3=front/back */
	double h_m; /*midpoint height of iceberg	*/
	int dirty_berg;
	double snout_i;/*initial position of snout*/
	double basal_t; /*thickness of baasal layer (useful when melting)*/
	double sed_top; //thickness of ice melted on surface of berg (to deposit sediment)
	double initialLength[3];
	double initialpos[2];
	double freeze_time;
	double stayinfjord;
	int left_fjord;
}iceberg;

void initialiseIceberg (iceberg *,int);	/* iceberg initialiseIceberg (iceberg *p_berg, int i)	 */
//double probCalving(glacier, double, double);		/*	double probCalvingBrown (glacier p_glacier, double dt)			 */
double dtCalving(glacier, double);
void icebergSize(iceberg *, glacier, int, long *);
void calving (double, iceberg *, int, glacier, options, long *, long *, long *, long *, long *, long *);	 /* iceberg calving (double time, iceberg *p_berg, int berg_no, glacier p_glacier) */
void bergBasalPos(iceberg *, int,glacier, long *);
void vectorMagnitude(double [], double *);
void waterVel(double [], double *, double []);
void airVel(double [], double *, long *, long *);

void waterVelTidal(double [], water_vel, double);

void perpLength_w(iceberg *, int, double [], double []);
void diffVect(double [],double [],double *,double,double);
void v_wDerivs(iceberg *, int, double [], water_vel, double);
void thicknessSeaIce(double *,double );
void perpLength_a(iceberg *, int, double *, double [])	;

void derivs(double,double [],double [], iceberg *, int, water_vel);
void bergMove(double [], iceberg *, int, water_vel, double);
void rk4(iceberg *, int, water_vel,double [],double [], double, double, double []);
void bounceOff(iceberg *,int,bathy, options, glacier, water_vel, double, long *, long *);
void meltWaterForcedConvection(double, double, double, double *);
void meltWaterVertConvection(double, double *);
void meltWaterWaveErosion(double,double *);
void meltAirSensHeat(double,double,double,double *);
void meltAirSolar(double, double *);
void calcMeanSeaTemp(double *,double *, iceberg *, int, options, bathy);
void calcMeanSeaSalin(double *,double *, iceberg *, int, options, bathy);

void meltBerg(iceberg *, int, glacier);
int stability(iceberg *, int);
void rollOver(iceberg *, int,glacier, long *);
int stability2(iceberg *, int);
void rollOver2(iceberg *, int, glacier, long *)	;
void depositionQuickWaterAir(sediment, bathy, glacier, options, iceberg *, int);
void depositionQuickWater(sediment, bathy, glacier, options, iceberg *, int);

double totbasal_berg, totbasal_glacier,toteng_berg, toteng_glacier;
double water_vel_res;

double tp, yp[4], hdid;

double p_v_w[2],p_v_a[2],p_v_w_mag,p_v_a_mag;

double p_v_w_tidal[2];
double v_w_total[2];
double mean_berg_volume, mean_berg_size, mean_berg_area;

double vol_sed_deposited;

//double warm_up_time=730. ; 

#endif

#ifndef BASALMELTWATER
#define BASALMELTWATER 	0.5	 /*factor decrease of debris rich ice melt rate vs englacial ice melt rate in water*/
#endif

#ifndef DEBRISMELTAIR
#define DEBRISMELTAIR 1.33	/*factor increase of ice with debris layer>0 and <2cm compared to clean ice in air*/
#endif
