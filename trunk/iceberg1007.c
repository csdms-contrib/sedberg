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

#include "iceberg.h"
#include "constants.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h> 


/*initialise info in iceberg structure to zero	*/
void initialiseIceberg (iceberg *p_berg, int i)
{
	int j;
	
	p_berg[i].tcalving_i = 0.0;	 				
	p_berg[i].d_i = 0.0;		
	for(j=0;j<2;j++)
	{
		p_berg[i].position[j] = 0.0;	
		p_berg[i].v_berg[j] = 0.0;	
		p_berg[i].initialpos[j] = 0.0;
	}
	p_berg[i].v_mag = 0.0;
	for(j=0;j<3;j++)
	{
		p_berg[i].size_i[j] = 0.0;	
		p_berg[i].R[j] = 0.0;		
		p_berg[i].conc_sed[j] = 0.0;
		p_berg[i].initialLength[j]=0.;
	}
	for(j=0;j<4;j++)
	{
		p_berg[i].a_basal[j] = 0.0;	
		p_berg[i].a_englacial[j] = 0.0;
	}
	p_berg[i].basal_position = NO_BASAL;
	p_berg[i].h_m=0.0;
	p_berg[i].dirty_berg =0;
	p_berg[i].snout_i=0.;
	p_berg[i].basal_t=0.;
	p_berg[i].sed_top=0.;
	p_berg[i].freeze_time=0.;
	p_berg[i].stayinfjord=0.;
	p_berg[i].left_fjord=0;

}

/*calculate calving timestep so that the probability of calving in each timestep is equal to P_calving (=0.25) by using Pelto and Warren's empirical calving relation*/ 
double dtCalving(glacier p_glacier, double annual_flux)
{
		/*dt in units of days*/
		double dt = 0.;
		double annualVol=annual_flux*1000000000.;//annual vol in m3 (annual_flux in km3per year)
  
		dt= (P_calving*mean_berg_volume*365.)/annualVol;

		return dt;			//dt in days

}

void icebergSize(iceberg *p_berg, glacier p_glacier, int berg_no_size, long *seed2)
{
	/*Calculates the berg size so that the size distribution
	of icebergs is reproduced using	inversion method of 
	random variable generation*/
	double random_size;
	double gausvar;

	gausvar= (double) gausdev2(seed2);
	random_size = (double) exp ( (double)(gausvar*p_glacier.sigma+p_glacier.mu));

	if (random_size<12.) 
	{
		random_size=12.;
	}
	/*size 0 is 1.62* size 1 and size 2 which are the same 
	the total volume is still the same as if all sides were size random size*/
	p_berg[berg_no_size].size_i[X]= random_size*1.3793;
	p_berg[berg_no_size].size_i[Y]=random_size*0.8514;
	p_berg[berg_no_size].size_i[Z]=random_size*0.8514;
}

//main iceberg initialisation routine - sets size, sediment content etc		   
void calving (double time, iceberg *p_berg, int berg_no_calving, glacier p_glacier, options p_options, long *seed2, long *seed, long *seed1, long *seed3, long *seed12, long *seed_stayinfjord)
{
    double random_basal = 0.0;	
	double random_y_pos = 0.0;
	double temp = 0.0;
	double nd_ntot=0.;
	double ran_freeze;
	double ran_no_stayinfjord;

	/*set all berg parameters to zero*/
	initialiseIceberg(p_berg, berg_no_calving);
	
	p_berg[berg_no_calving].tcalving_i=time;/*time of calving*/

	ran_freeze=(double) ran12(seed12);

	p_berg[berg_no_calving].freeze_time=ran_freeze;

	ran_no_stayinfjord=(double) f_ran_stayinfjord(seed_stayinfjord);

	p_berg[berg_no_calving].stayinfjord=ran_no_stayinfjord;

	icebergSize(p_berg, p_glacier, berg_no_calving, seed2);
	 
	p_berg[berg_no_calving].initialLength[X]=p_berg[berg_no_calving].size_i[X];
	p_berg[berg_no_calving].initialLength[Y]=p_berg[berg_no_calving].size_i[Y];
	p_berg[berg_no_calving].initialLength[Z]=p_berg[berg_no_calving].size_i[Z];
	
	//if the thickness of the berg is bigger than the ice front then adjust the thickness to the height of the front  
	if (p_berg[berg_no_calving].size_i[Z]>p_glacier.h_front)
	{
		p_berg[berg_no_calving].size_i[Z]=p_glacier.h_front;
		p_berg[berg_no_calving].dirty_berg = 1;
		p_berg[berg_no_calving].basal_t= p_glacier.h_basal_sediment;
		p_berg[berg_no_calving].basal_position = BASE;  //basal layer on bottom	
	}
	
	if 	(p_berg[berg_no_calving].size_i[X]>MAX_BERG_SIZE)
		p_berg[berg_no_calving].size_i[X]=MAX_BERG_SIZE;

	if(p_berg[berg_no_calving].size_i[Y]>MAX_BERG_SIZE)
		p_berg[berg_no_calving].size_i[Y]=MAX_BERG_SIZE;

/*if there is a basal sediment layer in the glacier*/
if (p_glacier.h_basal_sediment>0.)
{
	/*Find whether there is basal sediment in berg or not */											  
	if (p_berg[berg_no_calving].size_i[Z]<p_glacier.h_front)
	{
		nd_ntot=mean_berg_volume/mean_berg_area*p_glacier.h_front;
		random_basal = ran(seed);
		if(random_basal <= nd_ntot) 
		{
			p_berg[berg_no_calving].dirty_berg = 1;
			p_berg[berg_no_calving].basal_t= p_glacier.h_basal_sediment;
			//Basal position calculated as being on one of 4 possible sides after calving		
			bergBasalPos(p_berg,berg_no_calving,p_glacier, seed3);
		}
		else
		{
			p_berg[berg_no_calving].dirty_berg = 0;
			p_berg[berg_no_calving].basal_t= 0.;		
		}
	}  
}
/*if NO BASAL SED in icebergs*/
else
{
	p_berg[berg_no_calving].dirty_berg = 0;
	p_berg[berg_no_calving].basal_t= 0.;
	p_berg[berg_no_calving].basal_position = NO_BASAL;
  
	p_berg[berg_no_calving].a_englacial[TOP] = p_berg[berg_no_calving].size_i[X]*p_berg[berg_no_calving].size_i[Y];
	p_berg[berg_no_calving].a_englacial[SIDES] = 2.*p_berg[berg_no_calving].size_i[X]*p_berg[berg_no_calving].size_i[Z];
	p_berg[berg_no_calving].a_englacial[BASE] = p_berg[berg_no_calving].size_i[X]*p_berg[berg_no_calving].size_i[Y];
	p_berg[berg_no_calving].a_englacial[FRONT_BACK] = 2.*p_berg[berg_no_calving].size_i[Y]*p_berg[berg_no_calving].size_i[Z];
	p_berg[berg_no_calving].a_basal[TOP] = 0.;
	p_berg[berg_no_calving].a_basal[SIDES] = 0.;
	p_berg[berg_no_calving].a_basal[BASE] = 0.;
	p_berg[berg_no_calving].a_basal[FRONT_BACK] = 0.;
													
}
	p_berg[berg_no_calving].initialLength[X]=p_berg[berg_no_calving].size_i[X];
	p_berg[berg_no_calving].initialLength[Y]=p_berg[berg_no_calving].size_i[Y];
	p_berg[berg_no_calving].initialLength[Z]=p_berg[berg_no_calving].size_i[Z];
		
	/*average concentration of sediment on top of berg in kg/m */
	p_berg[berg_no_calving].conc_sed[TOP]=p_berg[berg_no_calving].a_basal[TOP]*p_glacier.C_basal_i*p_glacier.density_grain
		+p_berg[berg_no_calving].a_englacial[TOP]*p_glacier.C_englacial_i*p_glacier.density_grain;

	/*average concentration of sediment  on sides, front and back of berg in kg/m */
	p_berg[berg_no_calving].conc_sed[SIDES]=(p_berg[berg_no_calving].a_basal[SIDES]+p_berg[berg_no_calving].a_basal[FRONT_BACK])*p_glacier.C_basal_i*p_glacier.density_grain
		+(p_berg[berg_no_calving].a_englacial[SIDES]+p_berg[berg_no_calving].a_englacial[FRONT_BACK])*p_glacier.C_englacial_i*p_glacier.density_grain;

	/*average concentration of sediment  on bottom of berg in kg/m	*/
	p_berg[berg_no_calving].conc_sed[BASE]=p_berg[berg_no_calving].a_basal[BASE]*p_glacier.C_basal_i*p_glacier.density_grain
		+p_berg[berg_no_calving].a_englacial[BASE]*p_glacier.C_englacial_i*p_glacier.density_grain;

	if (time > p_options.spinuptime)
	{
		if (p_berg[berg_no_calving].basal_position == NO_BASAL)
		{
			temp=toteng_berg+(p_berg[berg_no_calving].size_i[X]*p_berg[berg_no_calving].size_i[Y]*p_berg[berg_no_calving].size_i[Z]);
			toteng_berg=temp;
		}

		else if (p_berg[berg_no_calving].basal_position == TOP)
		{
			temp= totbasal_berg+p_berg[berg_no_calving].a_basal[TOP]*p_berg[berg_no_calving].basal_t;
			totbasal_berg=temp;
			temp=toteng_berg+(p_berg[berg_no_calving].size_i[X]*p_berg[berg_no_calving].size_i[Y]*p_berg[berg_no_calving].size_i[Z]-p_berg[berg_no_calving].a_basal[TOP]*p_berg[berg_no_calving].basal_t);
			toteng_berg=temp;

		}
		else if (p_berg[berg_no_calving].basal_position == SIDES)
		{
			temp= totbasal_berg+p_berg[berg_no_calving].a_basal[SIDES]*p_berg[berg_no_calving].basal_t;	  	
			totbasal_berg=temp;	
			temp=toteng_berg+(p_berg[berg_no_calving].size_i[X]*p_berg[berg_no_calving].size_i[Y]*p_berg[berg_no_calving].size_i[Z]-p_berg[berg_no_calving].a_basal[SIDES]*p_berg[berg_no_calving].basal_t);
			toteng_berg=temp; 
		}																									   
		else if (p_berg[berg_no_calving].basal_position == BASE)
		{
			temp= totbasal_berg+p_berg[berg_no_calving].a_basal[BASE]*p_berg[berg_no_calving].basal_t;	
			totbasal_berg=temp;
			temp= toteng_berg+(p_berg[berg_no_calving].size_i[X]*p_berg[berg_no_calving].size_i[Y]*p_berg[berg_no_calving].size_i[Z]-p_berg[berg_no_calving].a_basal[BASE]*p_berg[berg_no_calving].basal_t);
			toteng_berg=temp;
		}
		
		else if (p_berg[berg_no_calving].basal_position == FRONT_BACK)
		{
			temp= totbasal_berg+p_berg[berg_no_calving].a_basal[FRONT_BACK]*p_berg[berg_no_calving].basal_t;	
			totbasal_berg=temp;
			temp= toteng_berg+(p_berg[berg_no_calving].size_i[X]*p_berg[berg_no_calving].size_i[Y]*p_berg[berg_no_calving].size_i[Z]-p_berg[berg_no_calving].a_basal[FRONT_BACK]*p_berg[berg_no_calving].basal_t);
			toteng_berg=temp;
		}
		
		temp=totbasal_glacier+(p_glacier.h_basal_sediment*p_berg[berg_no_calving].size_i[X]*p_berg[berg_no_calving].size_i[Y]*p_berg[berg_no_calving].size_i[Z]/p_glacier.h_front);
		totbasal_glacier= temp;

		temp= toteng_glacier+((p_glacier.h_front-p_glacier.h_basal_sediment)*p_berg[berg_no_calving].size_i[X]*p_berg[berg_no_calving].size_i[Y]*p_berg[berg_no_calving].size_i[Z]/p_glacier.h_front);
		toteng_glacier=	temp;
		
	}
	/*midpoint of iceberg in x-direction - iceberg starts depositing sediment one iceberg length in front of the snout*/
	p_berg[berg_no_calving].position[X] = p_glacier.y_f - (p_berg[berg_no_calving].size_i[X]/2.0);

	random_y_pos = ran1(seed1);
	/*y_position (parallel to glacier front).  Zero at left hand side of glacier when facing glacier*/ 
	p_berg[berg_no_calving].position[Y] = (p_glacier.w_fjord-p_glacier.w_front)/2.+random_y_pos*p_glacier.w_front;
			
	p_berg[berg_no_calving].initialpos[X] =	p_berg[berg_no_calving].position[X];
	p_berg[berg_no_calving].initialpos[Y] =	p_berg[berg_no_calving].position[Y];

	/*keel depth d	*/
	p_berg[berg_no_calving].d_i=0.833*p_berg[berg_no_calving].size_i[Z];	 
	
	/* initial velocity of berg u */
	p_berg[berg_no_calving].v_berg[X]=0.0; /*vel in km/hour*/
	p_berg[berg_no_calving].v_berg[Y]=0.0;

	p_berg[berg_no_calving].v_mag=0.0;
	p_berg[berg_no_calving].snout_i=p_glacier.y_f;

	p_berg[berg_no_calving].sed_top=0.;		
}	  


/*routine calculates position of basal layer after overturn*/
void bergBasalPos(iceberg *p_berg, int berg_no_basal, glacier p_glacier, long *seed3)
{
	double random_orientation = 0.0;
		
	random_orientation = ran3(seed3)*6.;
	
	/*if there is basal sediment in the iceberg*/
	if (p_berg[berg_no_basal].dirty_berg == 1)
	{
		/*basal layer on top*/
		if (random_orientation >= 0 && random_orientation <= 1)
		{
			p_berg[berg_no_basal].basal_position = TOP;	
					 
			p_berg[berg_no_basal].a_basal[TOP] = p_berg[berg_no_basal].size_i[X]*p_berg[berg_no_basal].size_i[Y];
			p_berg[berg_no_basal].a_basal[SIDES] = 2.*p_glacier.h_basal_sediment*p_berg[berg_no_basal].size_i[X];
			p_berg[berg_no_basal].a_basal[BASE] = 0.;
			p_berg[berg_no_basal].a_basal[FRONT_BACK] = 2.*p_glacier.h_basal_sediment*p_berg[berg_no_basal].size_i[Y];

			/*area of englacial layer */
			p_berg[berg_no_basal].a_englacial[TOP] = 0.;
			p_berg[berg_no_basal].a_englacial[SIDES] = 2.*p_berg[berg_no_basal].size_i[X]*p_berg[berg_no_basal].size_i[Z]-p_berg[berg_no_basal].a_basal[SIDES];
 			p_berg[berg_no_basal].a_englacial[BASE] = p_berg[berg_no_basal].size_i[X]*p_berg[berg_no_basal].size_i[Y];
			p_berg[berg_no_basal].a_englacial[FRONT_BACK] = 2.*p_berg[berg_no_basal].size_i[Y]*p_berg[berg_no_basal].size_i[Z]-p_berg[berg_no_basal].a_basal[FRONT_BACK];
			
		}	
		/*basal layer on either of 2 sides */
		else if (random_orientation > 1 && random_orientation <= 3)
		{
			p_berg[berg_no_basal].basal_position = SIDES;
		   
			/*area of basal layer*/			 
			p_berg[berg_no_basal].a_basal[TOP] = p_glacier.h_basal_sediment*p_berg[berg_no_basal].size_i[X];
			p_berg[berg_no_basal].a_basal[SIDES] = p_berg[berg_no_basal].size_i[X]*p_berg[berg_no_basal].size_i[Z];
			p_berg[berg_no_basal].a_basal[BASE] = p_glacier.h_basal_sediment*p_berg[berg_no_basal].size_i[X];
			p_berg[berg_no_basal].a_basal[FRONT_BACK] = 2.*p_glacier.h_basal_sediment*p_berg[berg_no_basal].size_i[Z];

			/*area of englacial layer */
			p_berg[berg_no_basal].a_englacial[TOP] = p_berg[berg_no_basal].size_i[X]*p_berg[berg_no_basal].size_i[Y]-p_berg[berg_no_basal].a_basal[TOP];
			p_berg[berg_no_basal].a_englacial[SIDES] = p_berg[berg_no_basal].size_i[X]*p_berg[berg_no_basal].size_i[Z];
 			p_berg[berg_no_basal].a_englacial[BASE] = p_berg[berg_no_basal].size_i[X]*p_berg[berg_no_basal].size_i[Y]-p_berg[berg_no_basal].a_basal[BASE];
			p_berg[berg_no_basal].a_englacial[FRONT_BACK] = 2.*p_berg[berg_no_basal].size_i[Y]*p_berg[berg_no_basal].size_i[Z]-p_berg[berg_no_basal].a_basal[FRONT_BACK];
		}	 

		/*basal layer on bottom*/
		else if (random_orientation > 3 && random_orientation <= 4)
		{
			p_berg[berg_no_basal].basal_position = BASE;	   

			/*area of basal layer*/					 
			p_berg[berg_no_basal].a_basal[TOP] = 0.;
			p_berg[berg_no_basal].a_basal[SIDES] = 2.*p_glacier.h_basal_sediment*p_berg[berg_no_basal].size_i[X];
			p_berg[berg_no_basal].a_basal[BASE] = p_berg[berg_no_basal].size_i[X]*p_berg[berg_no_basal].size_i[Y];
			p_berg[berg_no_basal].a_basal[FRONT_BACK] = 2.*p_glacier.h_basal_sediment*p_berg[berg_no_basal].size_i[Y];

			/*area of englacial layer */
			p_berg[berg_no_basal].a_englacial[TOP] = p_berg[berg_no_basal].size_i[X]*p_berg[berg_no_basal].size_i[Y]-p_berg[berg_no_basal].a_basal[TOP];
			p_berg[berg_no_basal].a_englacial[SIDES] = 2.*p_berg[berg_no_basal].size_i[X]*p_berg[berg_no_basal].size_i[Z]-p_berg[berg_no_basal].a_basal[SIDES];
 			p_berg[berg_no_basal].a_englacial[BASE] = 0.;
			p_berg[berg_no_basal].a_englacial[FRONT_BACK] = 2.*p_berg[berg_no_basal].size_i[Y]*p_berg[berg_no_basal].size_i[Z]-p_berg[berg_no_basal].a_basal[FRONT_BACK];
		}
			/*basal layer on FRONT/BACK*/

		else if (random_orientation > 4 && random_orientation <= 6)
		{
			p_berg[berg_no_basal].basal_position = FRONT_BACK;	   

			/*area of basal layer*/			 
			p_berg[berg_no_basal].a_basal[TOP] = p_glacier.h_basal_sediment*p_berg[berg_no_basal].size_i[Y];
			p_berg[berg_no_basal].a_basal[SIDES] = 	2.*p_glacier.h_basal_sediment*p_berg[berg_no_basal].size_i[Z];
			p_berg[berg_no_basal].a_basal[BASE] = p_glacier.h_basal_sediment*p_berg[berg_no_basal].size_i[Y];
			p_berg[berg_no_basal].a_basal[FRONT_BACK] = p_berg[berg_no_basal].size_i[Y]*p_berg[berg_no_basal].size_i[Z];

			/*area of englacial layer */
			p_berg[berg_no_basal].a_englacial[TOP] = p_berg[berg_no_basal].size_i[X]*p_berg[berg_no_basal].size_i[Y]-p_berg[berg_no_basal].a_basal[TOP];
			p_berg[berg_no_basal].a_englacial[SIDES] = 2.*p_berg[berg_no_basal].size_i[X]*p_berg[berg_no_basal].size_i[Z]-p_berg[berg_no_basal].a_basal[SIDES];
 			p_berg[berg_no_basal].a_englacial[BASE] = p_berg[berg_no_basal].size_i[X]*p_berg[berg_no_basal].size_i[Y]-p_berg[berg_no_basal].a_basal[BASE];
			p_berg[berg_no_basal].a_englacial[FRONT_BACK] = p_berg[berg_no_basal].size_i[Y]*p_berg[berg_no_basal].size_i[Z];			
		}
	}  /* end of case if there is basal sediment*/

	else if (p_berg[berg_no_basal].dirty_berg == 0)
	{
		p_berg[berg_no_basal].basal_position = NO_BASAL;

		p_berg[berg_no_basal].a_englacial[TOP] = p_berg[berg_no_basal].size_i[X]*p_berg[berg_no_basal].size_i[Y];
		p_berg[berg_no_basal].a_englacial[SIDES] = 2.*p_berg[berg_no_basal].size_i[X]*p_berg[berg_no_basal].size_i[Z];
		p_berg[berg_no_basal].a_englacial[BASE] = p_berg[berg_no_basal].size_i[X]*p_berg[berg_no_basal].size_i[Y];
		p_berg[berg_no_basal].a_englacial[FRONT_BACK] = 2.*p_berg[berg_no_basal].size_i[Y]*p_berg[berg_no_basal].size_i[Z];
		p_berg[berg_no_basal].a_basal[TOP] = 0.;
		p_berg[berg_no_basal].a_basal[SIDES] = 0.;
		p_berg[berg_no_basal].a_basal[BASE] = 0.;
		p_berg[berg_no_basal].a_basal[FRONT_BACK] = 0.;
	}

}

void vectorMagnitude(double vector[2], double *mag)
{
	*mag=sqrt(vector[X]*vector[X]+vector[Y]*vector[Y]);
}

void waterVel(double p_v_w[2], double *p_v_w_mag,double vel_random[2])
{
	double v_max1 = water_vel_res, v_max2 = 0.00;	/*maximum water speed km/hour (1m/s=3.6km/hr)  0.025*/
	
	p_v_w[X]=v_max1+v_max1*((vel_random[X]-0.5))/5.;
	p_v_w[Y]=v_max2+v_max1*((vel_random[Y]-0.5))/5.;
	vectorMagnitude(p_v_w,p_v_w_mag);
}

void airVel(double p_v_a[2], double *p_v_a_mag, long *seed4, long *seed5)
{
	double v_max1 = 15., v_max2 = 15.;	/*maximum wind speed km/hour (1m/s=3.6km/hr)*/
 	double U1, U2, b, sgn1, sgn2, find_log1, find_log2;
	//airvelocity varies with x max and min = vmax1, -vmax1 and y max and min = vmax2, -vmax2

	b=4.87;
	U1=ran4(seed4)-0.5;
	
	if (U1>0)
	{
		sgn1=1.;
	}
	else if (U1==0)
	{
		sgn1=0.;
	}
	else if (U1<0)
	{
		sgn1=-1.;
	}
	
	find_log1=1.-2.*fabs(U1);

	if (find_log1<=0.00000001)
	{
		find_log1=0.00000001;
		printf("Set find_log1 in calc of airvel[X] to avoid log of minus number\n");
	}
	
	U2=ran5(seed5)-0.5;

	if (U2>0)
	{
		sgn2=1.;
	}
	else if (U2==0)
	{
		sgn2=0.;
	}
	else if (U2<0)
	{
		sgn2=-1.;
	}

	find_log2=1.-2.*fabs(U2);

	if (find_log2<=0.00000001)
	{
		find_log2=0.00000001;
		printf("Set find_log2 in calc of airvel[Y] to avoid log of minus number\n");
	}

	p_v_a[X]=-b*sgn1*log(find_log1);
	p_v_a[Y]=-b*sgn2*log(find_log2);

	vectorMagnitude(p_v_a,p_v_a_mag);
}

void waterVelTidal(double p_v_w_tidal[2], water_vel p_water_vel,double t)
{
	 /*t is time in hours amplitude is in km/hour*/
	p_v_w_tidal[X]= p_water_vel.amplitude*sin(2.*PI*t/p_water_vel.T_var);
	p_v_w_tidal[Y]= p_water_vel.amplitude*cos(2.*PI*t/p_water_vel.T_var) ;

}


void perpLength_w(iceberg *p_berg, int berg_no_L_w, double perpL[2], double v_w_total[2])	 //length in km
{
	double theta;

	if (v_w_total[X]!=0.)
	{
		theta=atan2(v_w_total[Y],v_w_total[X]);
		if (theta<0.)
			theta=2.*PI+theta;
	}
	else
	theta=0.;

	/*length perpendicular to x-axis*/
	perpL[X]=fabs(p_berg[berg_no_L_w].size_i[Y]*cos(theta)+p_berg[berg_no_L_w].size_i[X]*sin(theta))/1000.;
	/*length perpendicular to y-axis*/
	perpL[Y]=fabs(p_berg[berg_no_L_w].size_i[X]*cos(theta)+p_berg[berg_no_L_w].size_i[Y]*sin(theta))/1000.;
} 

void v_wDerivs(iceberg *p_berg, int berg_no_vwderivs, double v_w_dt[2], water_vel p_water_vel, double t)
{
	v_w_dt[X]= p_water_vel.amplitude*2.*PI/p_water_vel.T_var*cos(2.*PI*t/p_water_vel.T_var);
	v_w_dt[Y]= -p_water_vel.amplitude*2.*PI/p_water_vel.T_var*sin(2.*PI*t/p_water_vel.T_var);
}

void thicknessSeaIce(double *thickness_si,double t)
{
	double yr_frac, days;
	days=t/24.;
	yr_frac=fmod(days,365.);
	if (yr_frac > 182. && yr_frac < 303.) *thickness_si=0.;
	else *thickness_si=0.001;	 /*thickness in km time in hours*/
}

void perpLength_a(iceberg *p_berg, int berg_no_perpL, double *perpL_a, double v_w_total[2])
{
	double theta;
	double lambda;
	
	if (v_w_total[X]!=0.)
	{
		theta=atan2(v_w_total[Y],v_w_total[X]);
		if (theta<0.)
			theta=2.*PI+theta;
	}
	else
		theta=0.;

	if (p_v_a[X]!=0.)
	{
	
		lambda=atan2(p_v_a[Y],p_v_a[X]);
		if (lambda<0.)
			lambda=2.*PI+lambda;
	}
	else
		lambda=0.;
	/*length perpendicular to air velocity*/
	*perpL_a=(fabs(p_berg[berg_no_perpL].size_i[X]*sin(fabs(lambda-theta)))+fabs(p_berg[berg_no_perpL].size_i[Y]*cos(fabs(lambda-theta))))/1000.;  //in km
} 
 
void diffVect(double vel[2], double diff[2], double *diff_mag, double y1, double y2)
{
	diff[X]=vel[X]-y1;
	diff[Y]=vel[Y]-y2;					   
	vectorMagnitude(diff,diff_mag);
}

//calculate derivatives of x,y,u and v
void derivs(double t, double y[4],double dydt[4], iceberg *p_berg, int berg_no_derivs, water_vel p_water_vel)
{	
	double p_v_w_diff[2], p_v_a_diff[2];
	double p_v_w_diff_mag, p_v_a_diff_mag;
	double perpL[2], v_w_dt[2], perpL_a;
	double berg_mass,Ek;
	double coriolis_x, coriolis_y, air_x,water_x, wave_x, press_x, air_y,water_y, wave_y, press_y; 
	double wave_amp;

	diffVect(v_w_total, p_v_w_diff, &p_v_w_diff_mag, y[2], y[3]);
	diffVect(p_v_a, p_v_a_diff, &p_v_a_diff_mag,  y[2], y[3]);
	
	/*Calculate length of iceberg perpendicular to x-axis [X] and y-axis[Y]*/
	perpLength_w(p_berg, berg_no_derivs, perpL, v_w_total);
	perpLength_a(p_berg, berg_no_derivs, &perpL_a,v_w_total);

	/*Calculate du_w/dt,dv_w/dt */
	v_wDerivs(p_berg, berg_no_derivs, v_w_dt, p_water_vel,t);

	berg_mass= p_berg[berg_no_derivs].size_i[X]*p_berg[berg_no_derivs].size_i[Y]*p_berg[berg_no_derivs].size_i[Z]*RHO_ICE;

	Ek=(p_berg[berg_no_derivs].d_i/1000.);
	if (p_berg[berg_no_derivs].d_i>90.)  Ek=0.09; //in km
	
	wave_amp=0.010125*p_v_a_mag*p_v_a_mag/12960000.; // in km^2

	/*[0] = X; [1] = Y; [2] = u (vel in X dir); [3] = v (vel in Y dir)*/

	dydt[0] = y[2];
	dydt[1] = y[3];

	coriolis_x= F_CORIOLIS*(y[3]-v_w_total[Y]);
	air_x = (0.5*RHO_A*1000000000.*C_A*perpL[X]*(p_berg[berg_no_derivs].size_i[Z]-p_berg[berg_no_derivs].d_i)/1000.+RHO_A*1000000000.*C_surf*p_berg[berg_no_derivs].size_i[X]*p_berg[berg_no_derivs].size_i[Y]/1000000.)*p_v_a_diff_mag*p_v_a_diff[X]/berg_mass;
	water_x = 0.5*(RHO_W*1000000000.*C_W*perpL[X]*p_berg[berg_no_derivs].d_i/1000.+RHO_W*1000000000.*C_surf*p_berg[berg_no_derivs].size_i[X]*p_berg[berg_no_derivs].size_i[Y]/1000000.)*p_v_w_diff_mag*p_v_w_diff[X]/berg_mass;
	if (p_v_a_mag>0.)
	wave_x = 0.25*RHO_W*1000000000.*(-g*12960.)*wave_amp*wave_amp*perpL_a*p_v_a[X]/(p_v_a_mag*berg_mass);	
	else
		wave_x=0.;
	press_x	= v_w_dt[X]+(1.5E-3*RHO_A*p_v_a_mag*p_v_a[X]/(Ek*RHO_W));	 //f x v_w in coriolis term

		
	dydt[2] = coriolis_x+air_x+water_x+wave_x+press_x;

	coriolis_y = -F_CORIOLIS*(y[2]-v_w_total[X]);
	air_y = (0.5*RHO_A*1000000000.*C_A*perpL[Y]*(p_berg[berg_no_derivs].size_i[Z]-p_berg[berg_no_derivs].d_i)/1000.+RHO_A*1000000000.*C_surf*p_berg[berg_no_derivs].size_i[X]*p_berg[berg_no_derivs].size_i[Y]/1000000.)*p_v_a_diff_mag*p_v_a_diff[Y]/berg_mass;
	water_y =  0.5*(RHO_W*1000000000.*C_W*perpL[Y]*p_berg[berg_no_derivs].d_i/1000.+RHO_W*1000000000.*C_surf*p_berg[berg_no_derivs].size_i[X]*p_berg[berg_no_derivs].size_i[Y]/1000000.)*p_v_w_diff_mag*p_v_w_diff[Y]/berg_mass;
	if (p_v_a_mag>0.)
		wave_y = 0.25*RHO_W*1000000000.*(-g*12960.)*wave_amp*wave_amp*perpL_a*p_v_a[Y]/(p_v_a_mag*berg_mass); 
	else
		wave_y=0.;
	press_y	=v_w_dt[Y]+(1.5E-3*RHO_A*p_v_a_mag*p_v_a[Y]/(Ek*RHO_W));	  //f x v_w in coriolis term

	dydt[3] = coriolis_y+air_y+water_y+wave_y+press_y; 	

}

//ODE Runge-Kutta solver
void rk4(iceberg *p_berg, int berg_no_rk4, water_vel p_water_vel, double y[4],double dydx[4], double x, double h, double yout[4])
{
	int i;
	double xh,hh,h6,dym[4],dyt[4],yt[4];

	hh=h*0.5;
	h6=h/6.0;		

	xh=x+hh;

	for (i=0;i<4;i++)
	{
		yt[i]=y[i]+hh*dydx[i];
	}
	derivs(xh,yt,dyt,p_berg,berg_no_rk4,p_water_vel);

	
	for (i=0;i<4;i++)
	{
		yt[i]=y[i]+hh*dyt[i];
	}
	derivs(xh,yt,dym,p_berg,berg_no_rk4,p_water_vel);

	for (i=0;i<4;i++)
	{
		yt[i]=y[i]+h*dym[i];
		dym[i] += dyt[i];
	}
	derivs(x+h,yt,dyt,p_berg,berg_no_rk4,p_water_vel);

	for (i=0;i<4;i++)
	{
		yout[i]=y[i]+h6*(dydx[i]+dyt[i]+2.0*dym[i]);
	}	
}

//routine to move iceberg
void bergMove(double ystart[4], iceberg *p_berg, int berg_no_move, water_vel p_water_vel, double h)
{	
	int i;
	double t;
	double y[4],yout[4],dy[4];

	t=tp;

	for (i=0;i<4;i++)
	{
		y[i]=ystart[i];	
	}

	derivs(t,y,dy,p_berg,berg_no_move,p_water_vel);

	rk4(p_berg,berg_no_move,p_water_vel,y,dy,t,h,yout);
	if ((double)(t+h) == t) 
	{
		printf("Error: Step size too small in routine rk4");
		exit(1);
	}

	for (i=0;i<4;i++)
	{
		y[i]=yout[i];
		yp[i]=y[i];
	}

}

//routine to simuulate deflection of iceberg off fjord walls and ice front
void bounceOff(iceberg *p_berg, int bounce_no,bathy fjord,options p_options, glacier p_glacier, water_vel p_water_vel, double t, long *seed10, long *seed11)
{
	double vel_initial[2], vel_final[2];	
	double v_initial_mag, v_final_mag, theta;
	double theta_final;
	double berg_mass, phi, lambda, beta;
	double r, r_y;
	double coeff_restitution=0.2, mom_inertia;
	double random_vel, random_deg, impulse;
	 
	vel_initial[X]=p_berg[bounce_no].v_berg[X];
	vel_initial[Y]=p_berg[bounce_no].v_berg[Y];
	berg_mass= p_berg[bounce_no].size_i[X]*p_berg[bounce_no].size_i[Y]*p_berg[bounce_no].size_i[Z]*RHO_ICE;

	vectorMagnitude(vel_initial, &v_initial_mag);

	mom_inertia=(1./12.)*berg_mass*(p_berg[bounce_no].size_i[X]*p_berg[bounce_no].size_i[X]+p_berg[bounce_no].size_i[Y]*p_berg[bounce_no].size_i[Y]);

	/*hits ice front*/
	if (floor((-p_berg[bounce_no].position[X]-p_berg[bounce_no].size_i[X]/2.)/p_options.dx)< ceil(-p_glacier.y_f/p_options.dx))
	{
		//for ice front
		theta=atan2(fabs(vel_initial[X]),fabs(vel_initial[Y]));  /*angle of berg moving into ice front */	

	
		if (((v_w_total[X]>0.) && (v_w_total[Y]>0.)) || ((v_w_total[X]<0.) && (v_w_total[Y]<0.)))
		{
			if(vel_initial[Y]<0.)
			{
				lambda=atan2(p_berg[bounce_no].size_i[X],p_berg[bounce_no].size_i[Y]);		//theta2
				phi=atan2(fabs(v_w_total[Y]),fabs(v_w_total[X])); /*angle between berg side and ice front (long axis parallel to water vel)*/
																 //phi2
			}
			else if(vel_initial[Y]>0.)
			{
				lambda=atan2(p_berg[bounce_no].size_i[Y],p_berg[bounce_no].size_i[X]);	//theta1
				phi=atan2(fabs(v_w_total[X]),fabs(v_w_total[Y])); /*angle between berg side and ice front (long axis parallel to water vel)*/
																 //phi1
			}
		}
		else if	(((v_w_total[X]>0.) && (v_w_total[Y]<0.)) || ((v_w_total[X]<0.) && (v_w_total[Y]>0.)))
		{
			if(vel_initial[Y]>0.)
			{
				lambda=atan2(p_berg[bounce_no].size_i[X],p_berg[bounce_no].size_i[Y]);	//theta2
				phi=atan2(fabs(v_w_total[Y]),fabs(v_w_total[X])); /*angle between berg side and ice front (long axis parallel to water vel)*/
																 //phi2
			}
			else if(vel_initial[Y]<0.)
			{
				lambda=atan2(p_berg[bounce_no].size_i[Y],p_berg[bounce_no].size_i[X]);	//theta1
				phi=atan2(fabs(v_w_total[X]),fabs(v_w_total[Y])); /*angle between berg side and ice front (long axis parallel to water vel)*/
																 //phi1
			}

		}

		r=0.5*sqrt(p_berg[bounce_no].size_i[X]*p_berg[bounce_no].size_i[X]+p_berg[bounce_no].size_i[Y]*p_berg[bounce_no].size_i[Y]);
		beta= lambda+phi;

		r_y=r*cos(beta);
		impulse	=vel_initial[X]*(1.+coeff_restitution)/((1./berg_mass)+r_y*r_y/mom_inertia);
		
		vel_final[X]= vel_initial[X]-impulse/berg_mass;

		vel_final[Y]= vel_initial[Y];

		theta_final=  atan2(fabs(vel_final[X]),fabs(vel_final[Y]));

		vectorMagnitude(vel_final, &v_final_mag);

		random_vel= (ran10(seed10)-0.5)/5.; //random no between -0.1 and 0.1

		v_final_mag=v_final_mag+random_vel*v_final_mag;	//final vel mag plus or minus 10%

		random_deg= (ran11(seed11)-0.5)*10.*PI/180.; //random no between -5 degrees (-5/180*pi) and +5 degrees (5/180*pi)
													  
		theta_final=theta_final+random_deg*theta_final;

		if (theta_final<0.)
		{
			theta_final=0.;
		}

		//calculate new velocity vectors from magnitude and direction of 'randomised' vectors, 
		//with direction depending on initial direction (so that berg 'bounces' of wall)
		if (vel_initial[X]<0.)
			p_berg[bounce_no].v_berg[X]= v_final_mag*sin(theta_final);
		else
			p_berg[bounce_no].v_berg[X]= -v_final_mag*sin(theta_final);
		 
		if (vel_initial[Y]<0.)
			p_berg[bounce_no].v_berg[Y]= -v_final_mag*cos(theta_final);
		else
			p_berg[bounce_no].v_berg[Y]= v_final_mag*cos(theta_final);


		p_berg[bounce_no].position[X]=p_glacier.y_f-p_berg[bounce_no].size_i[X]/2.;
	}
	/*if hits rhs fjord wall*/	
	else if (floor((p_berg[bounce_no].position[Y]-p_berg[bounce_no].size_i[X]/2.)/p_options.dx)<0. )
	{

			//for fjord walls
		theta=atan2(fabs(vel_initial[Y]),fabs(vel_initial[X]));  /*angle of berg moving into wall*/
	
		if (((v_w_total[X]>0.) && (v_w_total[Y]>0.)) || ((v_w_total[X]<0.) && (v_w_total[Y]<0.)))
		{
			if(vel_initial[X]>0.)
			{
				lambda=atan2(p_berg[bounce_no].size_i[X],p_berg[bounce_no].size_i[Y]);		//theta2
				phi=atan2(fabs(v_w_total[Y]),fabs(v_w_total[X])); /*angle between berg side and ice front (long axis parallel to water vel)*/
																 //phi2
			}
			else if(vel_initial[X]<0.)
			{
				lambda=atan2(p_berg[bounce_no].size_i[Y],p_berg[bounce_no].size_i[X]);	//theta1
				phi=atan2(fabs(v_w_total[X]),fabs(v_w_total[Y])); /*angle between berg side and ice front (long axis parallel to water vel)*/
																 //phi1
			}
		}
		else if	(((v_w_total[X]<0.) && (v_w_total[Y]<0.)) || ((v_w_total[X]<0.) && (v_w_total[Y]>0.)))
		{
			if(vel_initial[X]<0.)
			{
				lambda=atan2(p_berg[bounce_no].size_i[X],p_berg[bounce_no].size_i[Y]);	//theta2
				phi=atan2(fabs(v_w_total[Y]),fabs(v_w_total[X])); /*angle between berg side and ice front (long axis parallel to water vel)*/
																 //phi2
			}
			else if(vel_initial[X]>0.)
			{
				lambda=atan2(p_berg[bounce_no].size_i[Y],p_berg[bounce_no].size_i[X]);	//theta1
				phi=atan2(fabs(v_w_total[X]),fabs(v_w_total[Y])); /*angle between berg side and ice front (long axis parallel to water vel)*/
																 //phi1
			}
		}
	 
		r=0.5*sqrt(p_berg[bounce_no].size_i[X]*p_berg[bounce_no].size_i[X]+p_berg[bounce_no].size_i[Y]*p_berg[bounce_no].size_i[Y]);
		beta= lambda+phi;

		r_y=r*cos(beta);
		impulse	=vel_initial[Y]*(1.+coeff_restitution)/((1./berg_mass)+r_y*r_y/mom_inertia);
		
		vel_final[X]= vel_initial[X];

		vel_final[Y]= vel_initial[Y]-impulse/berg_mass;

		theta_final=  atan2(fabs(vel_final[Y]),fabs(vel_final[X]));

		vectorMagnitude(vel_final, &v_final_mag);

		random_vel= (ran10(seed10)-0.5)/5.; //random no between -0.1 and 0.1

		v_final_mag=v_final_mag+random_vel*v_final_mag;	//final vel mag plus or minus 10%

		random_deg= (ran11(seed11)-0.5)*10.*PI/180.; //random no between -5 degrees (-5/180*pi) and +5 degrees (5/180*pi)
													  
		theta_final=theta_final+random_deg*theta_final;

		if (theta_final<0.)
		{
			theta_final=0.;
		}

		//calculate new velocity vectors from magnitude and direction of 'randomised' vectors, 
		//with direction depending on initial direction (so that berg 'bounces' of wall)
		if (vel_initial[X]<0.)
			p_berg[bounce_no].v_berg[X]= -v_final_mag*cos(theta_final);
		else
			p_berg[bounce_no].v_berg[X]= v_final_mag*cos(theta_final);
		 
		if (vel_initial[Y]<0.)
			p_berg[bounce_no].v_berg[Y]= v_final_mag*sin(theta_final);
		else
			p_berg[bounce_no].v_berg[Y]= -v_final_mag*sin(theta_final);


		p_berg[bounce_no].position[Y]=0.+p_berg[bounce_no].size_i[X]/2.;
	}
	/*if hits lhs / east fjord wall*/
	else if ((int)ceil((p_berg[bounce_no].position[Y]+p_berg[bounce_no].size_i[X]/2.)/p_options.dx)>fjord.no_y_bins )
	{
	
		//for fjord walls
		theta=atan2(fabs(vel_initial[Y]),fabs(vel_initial[X]));  /*angle of berg moving into wall*/

		if (((v_w_total[X]>0.) && (v_w_total[Y]>0.)) || ((v_w_total[X]<0.) && (v_w_total[Y]<0.)))
		{
			if(vel_initial[X]<0.)
			{
				lambda=atan2(p_berg[bounce_no].size_i[X],p_berg[bounce_no].size_i[Y]);		//theta2
				phi=atan2(fabs(v_w_total[Y]),fabs(v_w_total[X])); /*angle between berg side and ice front (long axis parallel to water vel)*/
																 //phi2
			}
			else if(vel_initial[X]>0.)
			{
				lambda=atan2(p_berg[bounce_no].size_i[Y],p_berg[bounce_no].size_i[X]);	//theta1
				phi=atan2(fabs(v_w_total[X]),fabs(v_w_total[Y])); /*angle between berg side and ice front (long axis parallel to water vel)*/
																 //phi1
			}
		}
		else if	(((v_w_total[X]<0.) && (v_w_total[Y]<0.)) || ((v_w_total[X]<0.) && (v_w_total[Y]>0.)))
		{
			if(vel_initial[X]>0.)
			{
				lambda=atan2(p_berg[bounce_no].size_i[X],p_berg[bounce_no].size_i[Y]);	//theta2
				phi=atan2(fabs(v_w_total[Y]),fabs(v_w_total[X])); /*angle between berg side and ice front (long axis parallel to water vel)*/
																 //phi2
			}
			else if(vel_initial[X]<0.)
			{
				lambda=atan2(p_berg[bounce_no].size_i[Y],p_berg[bounce_no].size_i[X]);	//theta1
				phi=atan2(fabs(v_w_total[X]),fabs(v_w_total[Y])); /*angle between berg side and ice front (long axis parallel to water vel)*/
																 //phi1
			}
		}

		r=0.5*sqrt(p_berg[bounce_no].size_i[X]*p_berg[bounce_no].size_i[X]+p_berg[bounce_no].size_i[Y]*p_berg[bounce_no].size_i[Y]);
		beta= lambda+phi;

		r_y=r*cos(beta);
		impulse	=vel_initial[Y]*(1.+coeff_restitution)/((1./berg_mass)+r_y*r_y/mom_inertia);
		
		//calculate velocity after collision
		vel_final[X]= vel_initial[X];
		vel_final[Y]= vel_initial[Y]-impulse/berg_mass;

		theta_final=  atan2(fabs(vel_final[Y]),fabs(vel_final[X]));

		vectorMagnitude(vel_final, &v_final_mag);

		//add random noise to velocity after collision
		random_vel= (ran10(seed10)-0.5)/5.; //random no between -0.1 and 0.1

		v_final_mag=v_final_mag+random_vel*v_final_mag;	//final vel mag plus or minus 10%

		random_deg= (ran11(seed11)-0.5)*10.*PI/180.; //random no between -5 degrees (-5/180*pi) and +5 degrees (5/180*pi)
													  
		theta_final=theta_final+random_deg*theta_final;

		//check that theta isn't less than zero
		if (theta_final<0.)
		{
			theta_final=0.;
		}

		//calculate new velocity vectors from magnitude and direction of 'randomised' vectors, 
		//with direction depending on initial direction (so that berg 'bounces' of wall)
		if (vel_initial[X]<0.)
			p_berg[bounce_no].v_berg[X]= -v_final_mag*cos(theta_final);
		else
			p_berg[bounce_no].v_berg[X]= v_final_mag*cos(theta_final);
		 
		if (vel_initial[Y]<0.)
			p_berg[bounce_no].v_berg[Y]= v_final_mag*sin(theta_final);
		else
			p_berg[bounce_no].v_berg[Y]= -v_final_mag*sin(theta_final);

		p_berg[bounce_no].position[Y]=p_glacier.w_fjord-p_berg[bounce_no].size_i[X]/2.;
	}
	 
}

/*Melt rate due to forced convection (Weeks-Campbell equation).  Produces melt rate in m/day*/
/*this melt rate is too large!!*/
void meltWaterForcedConvection(double vel_diff, double tempWater, double L_a, double *meltForcedC_w)
{
	/*meltForcedC_w units of m/day
	vel_diff in km/hr
	L_a in m*/
	*meltForcedC_w = (6.74E-6*pow((vel_diff*1000./(60.*60.)),0.8)*(tempWater-TEMP_ICE)/pow(L_a,0.2))*60.*60.*24.;
 
}  
 
void meltWaterVertConvection(double tempWater,double *meltVert_w)
{
	/*meltVert_wunits of m/day*/
	*meltVert_w = 7.62E-3*tempWater+1.29E-3*tempWater*tempWater;
}

void meltWaterWaveErosion(double windSpeed,double *meltWave_w)
{
	/*units of m/day*/
	double seaState;

	seaState=pow(((windSpeed*1000./(60.*60.))/0.836),(2./3.));
		
	*meltWave_w = seaState*0.5;	
}
void meltAirSensHeat(double windSpeed, double longAxis, double tempAir, double *meltSens_a)
{	
	/*meltSens_a units of m/day*/
	/*windSpeed in km/hr, longAxis in m, tempAir in deg C*/
	/*MU_AIRin m^2/s
	KAPPA_AIR in m2/s
	K_AIR in W/(m deg C)
	GAMMA_ICE J/kg
	RHO_ICE kg/m3*/
	double Nusselt, Reynolds, Prandtl, heatFlowRate;

	if (tempAir>0.)
	{
		Prandtl= MU_AIR/ KAPPA_AIR;		/*no units*/
		Reynolds=windSpeed*1000./(60.*60.)*longAxis/MU_AIR;	/*no units*/
		Nusselt=0.0296*pow(Reynolds,0.8)*pow(Prandtl, 1./3.);
		
		heatFlowRate=Nusselt*K_AIR*tempAir/longAxis; 
		
		*meltSens_a = heatFlowRate*60.*60.*24./(RHO_ICE*GAMMA_ICE);
	}
	else
		*meltSens_a = 0.;
}
void meltAirSolar(double solarFlux, double *meltSolar_a)
{
	/*meltSolar_a units of m/day*/
	*meltSolar_a = solarFlux*60.*60.*24./(RHO_ICE*GAMMA_ICE);
}

void calcMeanSeaTemp(double *avTempWater, double *seaTemp, iceberg *p_berg, int berg_no_mseat, options p_options, bathy fjord)
{
	int max_for_av_temp = 0, l, max_depth;
	double total_temp = 0.0;

	max_depth=  fjord.no_z_bins-1;
	max_for_av_temp = (int)(ceil(p_berg[berg_no_mseat].d_i/p_options.dz));	  //max depth of berg in bins 
						
	*avTempWater = 0.0;
	if (max_for_av_temp > max_depth) 
	{
		printf("There is no temperature registered for this depth of iceberg keel\n");
		max_for_av_temp = max_depth;
	}

	total_temp = 0.0;
	
	//seaTemp index from 0 to fjord.no_z_bins-1 
	for (l=0;l<=max_for_av_temp;l++)
	{				
		total_temp = total_temp + seaTemp[l];
	}
	//average = total temperature/total number of measurements
	*avTempWater = total_temp/(max_for_av_temp+1);

}
  
void calcMeanSeaSalin(double *avSalinWater, double *seaSalin, iceberg *p_berg, int berg_no_mseat, options p_options, bathy fjord)
{
	int max_for_av_Salin = 0, l, max_depth;
	double total_Salin = 0.0;

	max_depth=  fjord.no_z_bins-1;
	max_for_av_Salin = (int)(ceil(p_berg[berg_no_mseat].d_i/p_options.dz));	  //max depth of berg in bins 
						
	*avSalinWater = 0.0;
	if (max_for_av_Salin > max_depth) 
	{
		printf("There is no temperature registered for this depth of iceberg keel\n");
		max_for_av_Salin = max_depth;
	}

	total_Salin = 0.0;
	
	//seaTemp index from 0 to fjord.no_z_bins-1 
	for (l=0;l<=max_for_av_Salin;l++)
	{				
		total_Salin = total_Salin + seaSalin[l];
	}
	//average = total temperature/total number of measurements
	*avSalinWater = total_Salin/(max_for_av_Salin+1);

}

void meltBerg(iceberg *p_berg, int berg_no_melt, glacier p_glacier)
{
	int i;
	double expcoeff_a=1.11, expcoeff_b=0.09, air_debris_melt;
	double sed_surf1, sed_surf2;  //boundaries between the criteria for enhanced melt by surface debris cover in metres of surface ice melt
		
	/*'if' loop to check whether berg size greater than 0 after melting (above)
		to calculate new sediment concs*/
	if ((p_berg[berg_no_melt].size_i[X]>10.) && (p_berg[berg_no_melt].size_i[Y]>10.) && (p_berg[berg_no_melt].size_i[Z]>10.))
	{	 
		/*If there is a BASAL layer BEFORE MELTING*/		/*If there is a BASAL layer BEFORE MELTING*/

		if(p_berg[berg_no_melt].dirty_berg==1)
		{	 	
			/*Basal layer is on the top of the iceberg*/
			if (p_berg[berg_no_melt].basal_position==TOP)
			{
				sed_surf1=0.02/(p_glacier.C_basal_i); //convert 0.02 m of sediment into m of ice equivalent
				sed_surf2=0.3/(p_glacier.C_basal_i); //convert 0.3 m of sediment into m of ice equivalent
	
				if ((p_berg[berg_no_melt].sed_top>=0.) && (p_berg[berg_no_melt].sed_top<=sed_surf1))
				{
					air_debris_melt=DEBRISMELTAIR*p_berg[berg_no_melt].R[TOP]*(hdid/24.); //melt rate of ice adjusted for sediment
					/*with thin debris layer/dirty ice melt rate increased by factor 1.33 (BASALMELTAIR) */
					p_berg[berg_no_melt].sed_top=p_berg[berg_no_melt].sed_top+air_debris_melt;
				}
				else if((p_berg[berg_no_melt].sed_top>sed_surf1) && (p_berg[berg_no_melt].sed_top<=sed_surf2))
				{
					air_debris_melt=expcoeff_a*exp(-expcoeff_b*p_berg[berg_no_melt].sed_top)*p_berg[berg_no_melt].R[TOP]*(hdid/24.); //melt rate adjusted for sediment
						p_berg[berg_no_melt].sed_top=p_berg[berg_no_melt].sed_top+air_debris_melt;
				}	
				else if(p_berg[berg_no_melt].sed_top>sed_surf2)
				{	
					air_debris_melt=0.;
					p_berg[berg_no_melt].sed_top=sed_surf2;
				}

				/*MELTING:  Change y,x and z size of iceberg*/
				p_berg[berg_no_melt].size_i[X]=p_berg[berg_no_melt].size_i[X]-2.*p_berg[berg_no_melt].R[SIDES]*hdid/24.;
				p_berg[berg_no_melt].size_i[Y]=p_berg[berg_no_melt].size_i[Y]-2.*p_berg[berg_no_melt].R[SIDES]*hdid/24.;
				p_berg[berg_no_melt].size_i[Z]=p_berg[berg_no_melt].size_i[Z]-(air_debris_melt+(p_berg[berg_no_melt].R[BASE]*hdid/24.));
				/*in units of m^3*/

				/*Calculate the new thickness of the basal layer*/
				p_berg[berg_no_melt].basal_t = p_berg[berg_no_melt].basal_t-air_debris_melt; 
				/*if the basal layer has melted away*/
				if (p_berg[berg_no_melt].basal_t<0.)
				{
					p_berg[berg_no_melt].dirty_berg=0;
					p_berg[berg_no_melt].basal_t = 0.;
					p_berg[berg_no_melt].basal_position=NO_BASAL;

					for (i=0;i<4;i++)
					{
						p_berg[berg_no_melt].a_basal[i] = 0.;
					}  	
					p_berg[berg_no_melt].a_englacial[TOP] = p_berg[berg_no_melt].size_i[X]*p_berg[berg_no_melt].size_i[Y];
					p_berg[berg_no_melt].a_englacial[SIDES] = 2.*p_berg[berg_no_melt].size_i[X]*p_berg[berg_no_melt].size_i[Z];
 					p_berg[berg_no_melt].a_englacial[BASE] = p_berg[berg_no_melt].size_i[X]*p_berg[berg_no_melt].size_i[Y];
					p_berg[berg_no_melt].a_englacial[FRONT_BACK] = 2.*p_berg[berg_no_melt].size_i[Y]*p_berg[berg_no_melt].size_i[Z];
				}
				else
				{
					p_berg[berg_no_melt].a_basal[TOP] = p_berg[berg_no_melt].size_i[X]*p_berg[berg_no_melt].size_i[Y];
					p_berg[berg_no_melt].a_basal[SIDES] = 2.*p_berg[berg_no_melt].basal_t*p_berg[berg_no_melt].size_i[X];
					p_berg[berg_no_melt].a_basal[BASE] = 0.;   	
					p_berg[berg_no_melt].a_basal[FRONT_BACK] = 2.*p_berg[berg_no_melt].basal_t*p_berg[berg_no_melt].size_i[Y];
						/*area of englacial layer */
					p_berg[berg_no_melt].a_englacial[TOP] = 0.;
					p_berg[berg_no_melt].a_englacial[SIDES] = 2.*p_berg[berg_no_melt].size_i[X]*p_berg[berg_no_melt].size_i[Z]-p_berg[berg_no_melt].a_basal[SIDES];
 					p_berg[berg_no_melt].a_englacial[BASE] = p_berg[berg_no_melt].size_i[X]*p_berg[berg_no_melt].size_i[Y]-p_berg[berg_no_melt].a_basal[BASE];
					p_berg[berg_no_melt].a_englacial[FRONT_BACK] = 2.*p_berg[berg_no_melt].size_i[Y]*p_berg[berg_no_melt].size_i[Z]-p_berg[berg_no_melt].a_basal[FRONT_BACK];
				}									  
			}

			/*basal layer on either of 2 sides */
			if (p_berg[berg_no_melt].basal_position==SIDES)
			{
				sed_surf1=0.02/((p_glacier.C_basal_i*p_berg[berg_no_melt].a_basal[SIDES]/(p_berg[berg_no_melt].a_basal[SIDES]+p_berg[berg_no_melt].a_englacial[SIDES]))+(p_glacier.C_englacial_i*p_berg[berg_no_melt].a_englacial[SIDES]/(p_berg[berg_no_melt].a_basal[SIDES]+p_berg[berg_no_melt].a_englacial[SIDES]))); //convert 0.02 m of sediment into m of ice equivalent
				sed_surf2=0.3/((p_glacier.C_basal_i*p_berg[berg_no_melt].a_basal[SIDES]/(p_berg[berg_no_melt].a_basal[SIDES]+p_berg[berg_no_melt].a_englacial[SIDES]))+(p_glacier.C_englacial_i*p_berg[berg_no_melt].a_englacial[SIDES]/(p_berg[berg_no_melt].a_basal[SIDES]+p_berg[berg_no_melt].a_englacial[SIDES]))); //convert 0.3 m of sediment into m of ice equivalent
	
				if ((p_berg[berg_no_melt].sed_top>=0.) && (p_berg[berg_no_melt].sed_top<=sed_surf1))
				{
					air_debris_melt=DEBRISMELTAIR*p_berg[berg_no_melt].R[TOP]*(hdid/24.); //melt rate adjusted for sediment
					/*with thin debris layer/dirty ice melt rate increased by factor 1.33 (BASALMELTAIR) */
					p_berg[berg_no_melt].sed_top=p_berg[berg_no_melt].sed_top+air_debris_melt;
				}
				else if((p_berg[berg_no_melt].sed_top>sed_surf1) && (p_berg[berg_no_melt].sed_top<=sed_surf2))
				{
					air_debris_melt=expcoeff_a*exp(-expcoeff_b*p_berg[berg_no_melt].sed_top)*p_berg[berg_no_melt].R[TOP]*(hdid/24.); //melt rate adjusted for sediment
						p_berg[berg_no_melt].sed_top=p_berg[berg_no_melt].sed_top+air_debris_melt;
				}	
				else if(p_berg[berg_no_melt].sed_top>sed_surf2)
				{	
					air_debris_melt=0.;
					p_berg[berg_no_melt].sed_top=sed_surf2;
				}
				/*If the berg has a basal and englacial layer*/
					/*MELTING:  Change y,x and z size of iceberg*/
				p_berg[berg_no_melt].size_i[X]=p_berg[berg_no_melt].size_i[X]-2.*p_berg[berg_no_melt].R[SIDES]*hdid/24.;
				p_berg[berg_no_melt].size_i[Y]=p_berg[berg_no_melt].size_i[Y]-(1.+BASALMELTWATER)*p_berg[berg_no_melt].R[SIDES]*hdid/24.;
				p_berg[berg_no_melt].size_i[Z]=p_berg[berg_no_melt].size_i[Z]-(air_debris_melt+(p_berg[berg_no_melt].R[BASE]*hdid/24.));
				/*in units of m^3*/

				/*Calculate the new thickness of the basal layer*/
				p_berg[berg_no_melt].basal_t = p_berg[berg_no_melt].basal_t-BASALMELTWATER*p_berg[berg_no_melt].R[SIDES]*hdid/24.; 
				/*if the basal layer has melted away*/
				if (p_berg[berg_no_melt].basal_t<0.)
				{
					p_berg[berg_no_melt].dirty_berg=0;
					p_berg[berg_no_melt].basal_t = 0.;
					p_berg[berg_no_melt].basal_position=NO_BASAL;
					for (i=0;i<4;i++)
					{
						p_berg[berg_no_melt].a_basal[i] = 0.;
					}
					p_berg[berg_no_melt].a_englacial[TOP] = p_berg[berg_no_melt].size_i[X]*p_berg[berg_no_melt].size_i[Y];
					p_berg[berg_no_melt].a_englacial[SIDES] = 2.*p_berg[berg_no_melt].size_i[X]*p_berg[berg_no_melt].size_i[Z];
 					p_berg[berg_no_melt].a_englacial[BASE] = p_berg[berg_no_melt].size_i[X]*p_berg[berg_no_melt].size_i[Y];
					p_berg[berg_no_melt].a_englacial[FRONT_BACK] = 2.*p_berg[berg_no_melt].size_i[Y]*p_berg[berg_no_melt].size_i[Z];
				}
				else
				{
					p_berg[berg_no_melt].a_basal[TOP] = p_berg[berg_no_melt].basal_t*p_berg[berg_no_melt].size_i[X];
					p_berg[berg_no_melt].a_basal[SIDES] = p_berg[berg_no_melt].size_i[X]*p_berg[berg_no_melt].size_i[Z];
					p_berg[berg_no_melt].a_basal[BASE] = p_berg[berg_no_melt].basal_t*p_berg[berg_no_melt].size_i[X];
					p_berg[berg_no_melt].a_basal[FRONT_BACK] = 2.*p_berg[berg_no_melt].basal_t*p_berg[berg_no_melt].size_i[Z];

					p_berg[berg_no_melt].a_englacial[TOP] = p_berg[berg_no_melt].size_i[X]*p_berg[berg_no_melt].size_i[Y]-p_berg[berg_no_melt].a_basal[TOP];
					p_berg[berg_no_melt].a_englacial[SIDES] = 2.*p_berg[berg_no_melt].size_i[X]*p_berg[berg_no_melt].size_i[Z]-p_berg[berg_no_melt].a_basal[SIDES];
 					p_berg[berg_no_melt].a_englacial[BASE] = p_berg[berg_no_melt].size_i[X]*p_berg[berg_no_melt].size_i[Y]-p_berg[berg_no_melt].a_basal[BASE];
					p_berg[berg_no_melt].a_englacial[FRONT_BACK] = 2.*p_berg[berg_no_melt].size_i[Y]*p_berg[berg_no_melt].size_i[Z]-p_berg[berg_no_melt].a_basal[FRONT_BACK];
				}
			}
 
			/*Basal layer is on the bottom of the iceberg*/
			if (p_berg[berg_no_melt].basal_position==BASE)
			{	  
				sed_surf1=0.02/(p_glacier.C_englacial_i); //convert 0.02 m of sediment into m of ice equivalent
				sed_surf2=0.3/(p_glacier.C_englacial_i); //convert 0.3 m of sediment into m of ice equivalent
	
				if ((p_berg[berg_no_melt].sed_top>=0.) && (p_berg[berg_no_melt].sed_top<=sed_surf1))
				{
					air_debris_melt=DEBRISMELTAIR*p_berg[berg_no_melt].R[TOP]*(hdid/24.); //melt rate adjusted for sediment
					/*with thin debris layer/dirty ice melt rate increased by factor 1.33 (BASALMELTAIR) */
					p_berg[berg_no_melt].sed_top=p_berg[berg_no_melt].sed_top+air_debris_melt;
				}
				else if((p_berg[berg_no_melt].sed_top>sed_surf1) && (p_berg[berg_no_melt].sed_top<=sed_surf2))
				{
					air_debris_melt=expcoeff_a*exp(-expcoeff_b*p_berg[berg_no_melt].sed_top)*p_berg[berg_no_melt].R[TOP]*(hdid/24.); //melt rate adjusted for sediment
					p_berg[berg_no_melt].sed_top=p_berg[berg_no_melt].sed_top+air_debris_melt;
				}	
				else if(p_berg[berg_no_melt].sed_top>sed_surf2)
				{	
					air_debris_melt=0.;
					p_berg[berg_no_melt].sed_top=sed_surf2;
				}

				/*If the berg has a basal and englacial layer*/
				/*MELTING:  Change y,x and z size of iceberg*/
				p_berg[berg_no_melt].size_i[X]=p_berg[berg_no_melt].size_i[X]-2.*p_berg[berg_no_melt].R[SIDES]*hdid/24.;
				p_berg[berg_no_melt].size_i[Y]=p_berg[berg_no_melt].size_i[Y]-2.*p_berg[berg_no_melt].R[SIDES]*hdid/24.;
				p_berg[berg_no_melt].size_i[Z]=p_berg[berg_no_melt].size_i[Z]-(air_debris_melt+(p_berg[berg_no_melt].R[BASE]*hdid/24.));
				/*in units of m^3*/
			
				/*Calculate the new thickness of the basal layer*/
				p_berg[berg_no_melt].basal_t = p_berg[berg_no_melt].basal_t-BASALMELTWATER*p_berg[berg_no_melt].R[BASE]*hdid/24.; 
				/*if the basal layer has melted away*/
				if (p_berg[berg_no_melt].basal_t<0.)
				{	
					p_berg[berg_no_melt].dirty_berg=0;
					p_berg[berg_no_melt].basal_t = 0.;
					p_berg[berg_no_melt].basal_position=NO_BASAL;

					for (i=0;i<4;i++)
					{
						p_berg[berg_no_melt].a_basal[i] = 0.;
					}
					p_berg[berg_no_melt].a_englacial[TOP] = p_berg[berg_no_melt].size_i[X]*p_berg[berg_no_melt].size_i[Y];
					p_berg[berg_no_melt].a_englacial[SIDES] = 2.*p_berg[berg_no_melt].size_i[X]*p_berg[berg_no_melt].size_i[Z];
 					p_berg[berg_no_melt].a_englacial[BASE] = p_berg[berg_no_melt].size_i[X]*p_berg[berg_no_melt].size_i[Y];
					p_berg[berg_no_melt].a_englacial[FRONT_BACK] = 2.*p_berg[berg_no_melt].size_i[Y]*p_berg[berg_no_melt].size_i[Z];
				}

				else
				{
					p_berg[berg_no_melt].a_basal[TOP] = 0.;
					p_berg[berg_no_melt].a_basal[SIDES] = 2.*p_berg[berg_no_melt].basal_t*p_berg[berg_no_melt].size_i[X];
					p_berg[berg_no_melt].a_basal[BASE] = p_berg[berg_no_melt].size_i[X]*p_berg[berg_no_melt].size_i[Y];
					p_berg[berg_no_melt].a_basal[FRONT_BACK] = 2.*p_berg[berg_no_melt].basal_t*p_berg[berg_no_melt].size_i[Y];
					/*area of englacial layer */
					p_berg[berg_no_melt].a_englacial[TOP] = p_berg[berg_no_melt].size_i[X]*p_berg[berg_no_melt].size_i[Y];
					p_berg[berg_no_melt].a_englacial[SIDES] = 2.*p_berg[berg_no_melt].size_i[X]*p_berg[berg_no_melt].size_i[Z]-p_berg[berg_no_melt].a_basal[SIDES];
 					p_berg[berg_no_melt].a_englacial[BASE] = 0.;
					p_berg[berg_no_melt].a_englacial[FRONT_BACK] = 2.*p_berg[berg_no_melt].size_i[Y]*p_berg[berg_no_melt].size_i[Z]-p_berg[berg_no_melt].a_basal[FRONT_BACK];
				}						
			}

			/*Basal layer is on the front/back of the iceberg*/
			if (p_berg[berg_no_melt].basal_position==FRONT_BACK)
			{	  
				sed_surf1=0.02/((p_glacier.C_basal_i*p_berg[berg_no_melt].a_basal[FRONT_BACK]/(p_berg[berg_no_melt].a_basal[FRONT_BACK]+p_berg[berg_no_melt].a_englacial[FRONT_BACK]))+(p_glacier.C_englacial_i*p_berg[berg_no_melt].a_englacial[FRONT_BACK]/(p_berg[berg_no_melt].a_basal[FRONT_BACK]+p_berg[berg_no_melt].a_englacial[FRONT_BACK]))); //convert 0.02 m of sediment into m of ice equivalent
				sed_surf2=0.3/((p_glacier.C_basal_i*p_berg[berg_no_melt].a_basal[FRONT_BACK]/(p_berg[berg_no_melt].a_basal[FRONT_BACK]+p_berg[berg_no_melt].a_englacial[FRONT_BACK]))+(p_glacier.C_englacial_i*p_berg[berg_no_melt].a_englacial[FRONT_BACK]/(p_berg[berg_no_melt].a_basal[FRONT_BACK]+p_berg[berg_no_melt].a_englacial[FRONT_BACK]))); //convert 0.3 m of sediment into m of ice equivalent
	
				if ((p_berg[berg_no_melt].sed_top>=0.) && (p_berg[berg_no_melt].sed_top<=sed_surf1))
				{
					air_debris_melt=DEBRISMELTAIR*p_berg[berg_no_melt].R[TOP]*(hdid/24.); //melt rate adjusted for sediment
					/*with thin debris layer/dirty ice melt rate increased by factor 1.33 (BASALMELTAIR) */
					p_berg[berg_no_melt].sed_top=p_berg[berg_no_melt].sed_top+air_debris_melt;
				}
				else if((p_berg[berg_no_melt].sed_top>sed_surf1) && (p_berg[berg_no_melt].sed_top<=sed_surf2))
				{
					air_debris_melt=expcoeff_a*exp(-expcoeff_b*p_berg[berg_no_melt].sed_top)*p_berg[berg_no_melt].R[TOP]*(hdid/24.); //melt rate adjusted for sediment
					p_berg[berg_no_melt].sed_top=p_berg[berg_no_melt].sed_top+air_debris_melt;
				}	
				else if(p_berg[berg_no_melt].sed_top>sed_surf2)
				{	
					air_debris_melt=0.;
					p_berg[berg_no_melt].sed_top=sed_surf2;
				}

				/*If the berg has a basal and englacial layer*/
				/*MELTING:  Change y,x and z size of iceberg*/
				p_berg[berg_no_melt].size_i[X]=p_berg[berg_no_melt].size_i[X]-(1.+BASALMELTWATER)*p_berg[berg_no_melt].R[SIDES]*hdid/24.;
				p_berg[berg_no_melt].size_i[Y]=p_berg[berg_no_melt].size_i[Y]-2.*p_berg[berg_no_melt].R[SIDES]*hdid/24.;
				p_berg[berg_no_melt].size_i[Z]=p_berg[berg_no_melt].size_i[Z]-(air_debris_melt+(p_berg[berg_no_melt].R[BASE]*hdid/24.));
				/*in units of m^3*/
			
				/*Calculate the new thickness of the basal layer*/
				p_berg[berg_no_melt].basal_t = p_berg[berg_no_melt].basal_t-BASALMELTWATER*p_berg[berg_no_melt].R[SIDES]*hdid/24.; 
				/*if the basal layer has melted away*/
				if (p_berg[berg_no_melt].basal_t<0.)
				{
					p_berg[berg_no_melt].dirty_berg=0;
					p_berg[berg_no_melt].basal_t = 0.;
					p_berg[berg_no_melt].basal_position=NO_BASAL;

					for (i=0;i<4;i++)
					{
						p_berg[berg_no_melt].a_basal[i] = 0.;
					}
					p_berg[berg_no_melt].a_englacial[TOP] = p_berg[berg_no_melt].size_i[X]*p_berg[berg_no_melt].size_i[Y];
					p_berg[berg_no_melt].a_englacial[SIDES] = 2.*p_berg[berg_no_melt].size_i[X]*p_berg[berg_no_melt].size_i[Z];
 					p_berg[berg_no_melt].a_englacial[BASE] = p_berg[berg_no_melt].size_i[X]*p_berg[berg_no_melt].size_i[Y];
					p_berg[berg_no_melt].a_englacial[FRONT_BACK] = 2.*p_berg[berg_no_melt].size_i[Y]*p_berg[berg_no_melt].size_i[Z];
				}

				else
				{
					p_berg[berg_no_melt].a_basal[TOP] = p_berg[berg_no_melt].basal_t*p_berg[berg_no_melt].size_i[Y];
					p_berg[berg_no_melt].a_basal[SIDES] = 2.*p_berg[berg_no_melt].basal_t*p_berg[berg_no_melt].size_i[Z];
					p_berg[berg_no_melt].a_basal[BASE] = p_berg[berg_no_melt].basal_t*p_berg[berg_no_melt].size_i[Y];
					p_berg[berg_no_melt].a_basal[FRONT_BACK] = p_berg[berg_no_melt].size_i[Y]*p_berg[berg_no_melt].size_i[Z];
						/*area of englacial layer */
					p_berg[berg_no_melt].a_englacial[TOP] = p_berg[berg_no_melt].size_i[X]*p_berg[berg_no_melt].size_i[Y]-p_berg[berg_no_melt].a_basal[TOP];
					p_berg[berg_no_melt].a_englacial[SIDES] = 2.*p_berg[berg_no_melt].size_i[X]*p_berg[berg_no_melt].size_i[Z]-p_berg[berg_no_melt].a_basal[SIDES];
 					p_berg[berg_no_melt].a_englacial[BASE] = p_berg[berg_no_melt].size_i[X]*p_berg[berg_no_melt].size_i[Y]-p_berg[berg_no_melt].a_basal[BASE];
					p_berg[berg_no_melt].a_englacial[FRONT_BACK] = p_berg[berg_no_melt].size_i[Y]*p_berg[berg_no_melt].size_i[Z];

				}					
			}
		}	 /*end of case if there is basal sediment present*/

		/* if there is ONLY ENGLACIAL sediment before melting */
		else if (p_berg[berg_no_melt].dirty_berg==0)
		{	
			sed_surf1=0.02/(p_glacier.C_englacial_i); //convert 0.02 m of sediment into m of ice equivalent
			sed_surf2=0.3/(p_glacier.C_englacial_i); //convert 0.3 m of sediment into m of ice equivalent

			if ((p_berg[berg_no_melt].sed_top>=0.) && (p_berg[berg_no_melt].sed_top<=sed_surf1))
			{
				air_debris_melt=DEBRISMELTAIR*p_berg[berg_no_melt].R[TOP]*(hdid/24.); //melt rate adjusted for sediment
				/*with thin debris layer/dirty ice melt rate increased by factor 1.33 (BASALMELTAIR) */
				p_berg[berg_no_melt].sed_top=p_berg[berg_no_melt].sed_top+air_debris_melt;
			}
			else if((p_berg[berg_no_melt].sed_top>sed_surf1) && (p_berg[berg_no_melt].sed_top<=sed_surf2))
			{
				air_debris_melt=expcoeff_a*exp(-expcoeff_b*p_berg[berg_no_melt].sed_top)*p_berg[berg_no_melt].R[TOP]*(hdid/24.); //melt rate adjusted for sediment
				p_berg[berg_no_melt].sed_top=p_berg[berg_no_melt].sed_top+air_debris_melt;
			}	
			else if(p_berg[berg_no_melt].sed_top>sed_surf2)
			{	
				air_debris_melt=0.;
				p_berg[berg_no_melt].sed_top=sed_surf2;
			}
			/*MELTING:  Change y,x and z size of iceberg*/
			p_berg[berg_no_melt].size_i[X]=p_berg[berg_no_melt].size_i[X]-2.*p_berg[berg_no_melt].R[SIDES]*hdid/24.;
			p_berg[berg_no_melt].size_i[Y]=p_berg[berg_no_melt].size_i[Y]-2.*p_berg[berg_no_melt].R[SIDES]*hdid/24.;
			p_berg[berg_no_melt].size_i[Z]=p_berg[berg_no_melt].size_i[Z]-(air_debris_melt+(p_berg[berg_no_melt].R[BASE]*hdid/24.));
			/*in units of m^3*/
				
			p_berg[berg_no_melt].a_englacial[TOP] = p_berg[berg_no_melt].size_i[X]*p_berg[berg_no_melt].size_i[Y];
			p_berg[berg_no_melt].a_englacial[SIDES] = 2.*p_berg[berg_no_melt].size_i[X]*p_berg[berg_no_melt].size_i[Z];
 			p_berg[berg_no_melt].a_englacial[BASE] = p_berg[berg_no_melt].size_i[X]*p_berg[berg_no_melt].size_i[Y];
			p_berg[berg_no_melt].a_englacial[FRONT_BACK] = 2.*p_berg[berg_no_melt].size_i[Y]*p_berg[berg_no_melt].size_i[Z];
		}
	
		/*average concentration of sediment on top of berg in kg/m */
		p_berg[berg_no_melt].conc_sed[TOP]=p_berg[berg_no_melt].a_basal[TOP]*p_glacier.C_basal_i*p_glacier.density_grain
			+p_berg[berg_no_melt].a_englacial[TOP]*p_glacier.C_englacial_i*p_glacier.density_grain;

		/*average concentration of sediment  on sides, front and back of berg in kg/m */
		p_berg[berg_no_melt].conc_sed[SIDES]=(p_berg[berg_no_melt].a_basal[SIDES]+p_berg[berg_no_melt].a_basal[FRONT_BACK])*p_glacier.C_basal_i*p_glacier.density_grain
			+(p_berg[berg_no_melt].a_englacial[SIDES]+p_berg[berg_no_melt].a_englacial[FRONT_BACK])*p_glacier.C_englacial_i*p_glacier.density_grain;

		/*average concentration of sediment  on bottom of berg in kg/m	*/
		p_berg[berg_no_melt].conc_sed[BASE]=p_berg[berg_no_melt].a_basal[BASE]*p_glacier.C_basal_i*p_glacier.density_grain
			+p_berg[berg_no_melt].a_englacial[BASE]*p_glacier.C_englacial_i*p_glacier.density_grain;

		/*keel depth d	*/
		p_berg[berg_no_melt].d_i=0.833*p_berg[berg_no_melt].size_i[Z];	
	
	} /*End of 'if' loop to check whether berg size greater than 10 before
		calculating new sediment concs*/

}

void depositionQuickWater(sediment p_sediment, bathy fjord, glacier p_glacier, options p_options, iceberg *p_berg, int berg_no_dep)
{
	long min_y_dep = 0, min_x_dep = 0, max_x_dep=0,max_y_dep=0, temp, cellx, celly;
	double  area_dep=0., cellThickness=0., temp_thickness, diff_x, diff_y;
	double dep_spread_dist_x=SPREAD_X;  //distance over which sed spread by currents (m)
	double dep_spread_dist_y=SPREAD_Y; 

	min_x_dep = (int)(floor(-p_berg[berg_no_dep].position[X]/p_options.dx));
	max_x_dep = (int)(floor((-p_berg[berg_no_dep].position[X]+p_berg[berg_no_dep].v_berg[X]*hdid)/p_options.dx));
	min_y_dep = (int)(floor(p_berg[berg_no_dep].position[Y]/p_options.dx));
	max_y_dep = (int)(floor((p_berg[berg_no_dep].position[Y]+p_berg[berg_no_dep].v_berg[Y]*hdid)/p_options.dx));

	/* x and y distance over which deposition occurs in m	*/
	if (min_x_dep>max_x_dep)
	{
		temp=max_x_dep;
		max_x_dep=min_x_dep;
		min_x_dep=temp;
	}

	if (min_y_dep>max_y_dep)
	{
		temp=max_y_dep;
		max_y_dep=min_y_dep;
		min_y_dep=temp;
	}

	min_x_dep = (int)(floor((min_x_dep*p_options.dx-p_berg[berg_no_dep].size_i[X]/2.-dep_spread_dist_x)/p_options.dx));
	max_x_dep = (int)(floor((max_x_dep*p_options.dx+p_berg[berg_no_dep].size_i[X]/2.+dep_spread_dist_x)/p_options.dx));
	min_y_dep = (int)(floor((min_y_dep*p_options.dx-p_berg[berg_no_dep].size_i[X]/2.-dep_spread_dist_y)/p_options.dx));
	max_y_dep = (int)(floor((max_y_dep*p_options.dx+p_berg[berg_no_dep].size_i[X]/2.+dep_spread_dist_y)/p_options.dx));

	if (min_x_dep< (int)floor(-p_glacier.y_f/p_options.dx)) 
		min_x_dep=(int)(floor(-p_glacier.y_f/p_options.dx));
	
	if (max_y_dep>=fjord.no_y_bins) 
		max_y_dep=fjord.no_y_bins-1;
	
	if (min_y_dep<0) 
		min_y_dep=0;
	
	diff_x =  (double)((max_x_dep-min_x_dep+1)*p_options.dx); /*length in m*/
	diff_y =  (double)((max_y_dep-min_y_dep+1)*p_options.dx);

	if (diff_x<=0. || diff_y<=0.)
		printf("oops\n");

	if (max_x_dep>=fjord.no_x_bins) 
		max_x_dep=fjord.no_x_bins-1;

	area_dep = diff_x*diff_y;
	
	//vol dep in m^3
	vol_sed_deposited=vol_sed_deposited+((p_berg[berg_no_dep].R[SIDES]*(hdid/24.)*p_berg[berg_no_dep].conc_sed[SIDES]+p_berg[berg_no_dep].R[BASE]*(hdid/24.)*p_berg[berg_no_dep].conc_sed[BASE])/p_glacier.density_grain);
	//in cm
	cellThickness= (p_berg[berg_no_dep].R[SIDES]*(hdid/24.)*p_berg[berg_no_dep].conc_sed[SIDES]+p_berg[berg_no_dep].R[BASE]*(hdid/24.)*p_berg[berg_no_dep].conc_sed[BASE])/(p_glacier.density_dep*area_dep)*100.;
		
	for (cellx=min_x_dep;cellx<=max_x_dep;cellx++)
	{
		for	(celly=min_y_dep;celly<=max_y_dep;celly++)
		{
			temp_thickness=p_sediment.thickness[cellx][celly];		
			p_sediment.thickness[cellx][celly]=temp_thickness+cellThickness;

		}
	}

}

void depositionQuickWaterAir(sediment p_sediment, bathy fjord, glacier p_glacier, options p_options, iceberg *p_berg, int berg_no_dep)
{
	long min_y_dep = 0, min_x_dep = 0, max_x_dep=0,max_y_dep=0, temp, cellx, celly;
	double  area_dep=0., cellThickness=0., temp_thickness, diff_x, diff_y;
	double dep_spread_dist_x=SPREAD_X;  //distance over which sed spread by currents (m)
	double dep_spread_dist_y=SPREAD_Y; 

	min_x_dep = (int)(floor(-p_berg[berg_no_dep].position[X]/p_options.dx));
	max_x_dep = (int)(ceil((-p_berg[berg_no_dep].position[X]+p_berg[berg_no_dep].v_berg[X]*hdid)/p_options.dx));
	min_y_dep = (int)(floor(p_berg[berg_no_dep].position[Y]/p_options.dx));
	max_y_dep = (int)(ceil((p_berg[berg_no_dep].position[Y]+p_berg[berg_no_dep].v_berg[Y]*hdid)/p_options.dx));

	/* x and y distance over which deposition occurs in m	*/
	if (min_x_dep>max_x_dep)
	{
		temp=max_x_dep;
		max_x_dep=min_x_dep;
		min_x_dep=temp;
	}

	if (min_y_dep>max_y_dep)
	{
		temp=max_y_dep;
		max_y_dep=min_y_dep;
		min_y_dep=temp;
	}

	min_x_dep = (int)(floor((min_x_dep*p_options.dx-p_berg[berg_no_dep].size_i[X]/2.-dep_spread_dist_x)/p_options.dx));
	max_x_dep = (int)(floor((max_x_dep*p_options.dx+p_berg[berg_no_dep].size_i[X]/2.+dep_spread_dist_x)/p_options.dx));
	min_y_dep = (int)(floor((min_y_dep*p_options.dx-p_berg[berg_no_dep].size_i[X]/2.-dep_spread_dist_y)/p_options.dx));
	max_y_dep = (int)(floor((max_y_dep*p_options.dx+p_berg[berg_no_dep].size_i[X]/2.+dep_spread_dist_y)/p_options.dx));

	if (min_x_dep< (int)ceil(-p_glacier.y_f/p_options.dx)) 
		min_x_dep=(int)(floor(-p_glacier.y_f/p_options.dx));

	if (max_y_dep>=fjord.no_y_bins) 
		max_y_dep=fjord.no_y_bins-1;

	if (min_y_dep<0) 
		min_y_dep=0;
	
	diff_x =  (double)((max_x_dep-min_x_dep+1)*p_options.dx); /*length in m*/
	diff_y =  (double)((max_y_dep-min_y_dep+1)*p_options.dx);
	
	if (diff_x<=0. || diff_y<=0.)
		printf("oops\n");

	if (max_x_dep>=fjord.no_x_bins) 
		max_x_dep=fjord.no_x_bins-1;

	area_dep = diff_x*diff_y;

	//vol dep in m^3
	vol_sed_deposited=vol_sed_deposited+((p_berg[berg_no_dep].sed_top*p_berg[berg_no_dep].conc_sed[TOP]+p_berg[berg_no_dep].R[SIDES]*(hdid/24.)*p_berg[berg_no_dep].conc_sed[SIDES]+p_berg[berg_no_dep].R[BASE]*(hdid/24.)*p_berg[berg_no_dep].conc_sed[BASE])/p_glacier.density_grain);	
	//in cm
	cellThickness= (p_berg[berg_no_dep].sed_top*p_berg[berg_no_dep].conc_sed[TOP]+p_berg[berg_no_dep].R[SIDES]*(hdid/24.)*p_berg[berg_no_dep].conc_sed[SIDES]+p_berg[berg_no_dep].R[BASE]*(hdid/24.)*p_berg[berg_no_dep].conc_sed[BASE])/(p_glacier.density_dep*area_dep)*100.;
		
	for (cellx=min_x_dep;cellx<=max_x_dep;cellx++)
	{
		for	(celly=min_y_dep;celly<=max_y_dep;celly++)
		{
			temp_thickness=p_sediment.thickness[cellx][celly];		
			p_sediment.thickness[cellx][celly]=temp_thickness+cellThickness;
		}
	}
	p_berg[berg_no_dep].sed_top=0.;
}

int stability(iceberg *p_berg, int berg_no_stab)
{
	double lhs,rhs;

	lhs=p_berg[berg_no_stab].size_i[X]/p_berg[berg_no_stab].size_i[Z];

	rhs=0.7;

	if (lhs>rhs) /*stable*/
	{
		return 0;
	}
	else /*unstable*/
		return 1;
}

void rollOver(iceberg *p_berg, int berg_no_roll, glacier p_glacier, long *seed9)
{
	double temp_size, temp_basal, temp_englacial;

	if (p_berg[berg_no_roll].size_i[Z]>=p_berg[berg_no_roll].size_i[Y])
	{	
		temp_size=p_berg[berg_no_roll].size_i[X];
		p_berg[berg_no_roll].size_i[X]=p_berg[berg_no_roll].size_i[Z];
		p_berg[berg_no_roll].size_i[Z]=temp_size;

		if ((p_berg[berg_no_roll].basal_position==TOP) || (p_berg[berg_no_roll].basal_position==SIDES) || (p_berg[berg_no_roll].basal_position==BASE) ||(p_berg[berg_no_roll].basal_position==NO_BASAL))
		{
			temp_basal=p_berg[berg_no_roll].a_basal[FRONT_BACK];
			p_berg[berg_no_roll].a_basal[FRONT_BACK]=p_berg[berg_no_roll].a_basal[TOP]+p_berg[berg_no_roll].a_basal[BASE];
			p_berg[berg_no_roll].a_basal[TOP]=temp_basal/2.;
			p_berg[berg_no_roll].a_basal[BASE]=temp_basal/2.;

			temp_englacial=p_berg[berg_no_roll].a_englacial[FRONT_BACK];
			p_berg[berg_no_roll].a_englacial[FRONT_BACK]=p_berg[berg_no_roll].a_englacial[TOP]+p_berg[berg_no_roll].a_englacial[BASE];
			p_berg[berg_no_roll].a_englacial[TOP]=temp_englacial/2.;
			p_berg[berg_no_roll].a_englacial[BASE]=temp_englacial/2.;
			
			if (p_berg[berg_no_roll].basal_position==TOP)
				p_berg[berg_no_roll].basal_position=FRONT_BACK;

			else if (p_berg[berg_no_roll].basal_position==BASE)
				p_berg[berg_no_roll].basal_position=FRONT_BACK;
				
		}

		else if (p_berg[berg_no_roll].basal_position==FRONT_BACK)
		{
			if(ran9(seed9)<0.5)
			{
				temp_basal=p_berg[berg_no_roll].a_basal[FRONT_BACK];
				p_berg[berg_no_roll].a_basal[FRONT_BACK]=p_berg[berg_no_roll].a_basal[TOP]+p_berg[berg_no_roll].a_basal[BASE];
				p_berg[berg_no_roll].a_basal[TOP]=temp_basal;
				p_berg[berg_no_roll].a_basal[BASE]=0.;

				temp_englacial=p_berg[berg_no_roll].a_englacial[FRONT_BACK];
				p_berg[berg_no_roll].a_englacial[FRONT_BACK]=p_berg[berg_no_roll].a_englacial[TOP]+p_berg[berg_no_roll].a_englacial[BASE];
				p_berg[berg_no_roll].a_englacial[TOP]=0.;
				p_berg[berg_no_roll].a_englacial[BASE]=temp_englacial;
					
				p_berg[berg_no_roll].basal_position=TOP;
			}
			else if(ran9(seed9)>0.5)
			{
				temp_basal=p_berg[berg_no_roll].a_basal[FRONT_BACK];
				p_berg[berg_no_roll].a_basal[FRONT_BACK]=p_berg[berg_no_roll].a_basal[TOP]+p_berg[berg_no_roll].a_basal[BASE];
				p_berg[berg_no_roll].a_basal[TOP]=0.;
				p_berg[berg_no_roll].a_basal[BASE]=temp_basal;

				temp_englacial=p_berg[berg_no_roll].a_englacial[FRONT_BACK];
				p_berg[berg_no_roll].a_englacial[FRONT_BACK]=p_berg[berg_no_roll].a_englacial[TOP]+p_berg[berg_no_roll].a_englacial[BASE];
				p_berg[berg_no_roll].a_englacial[TOP]=temp_englacial;
				p_berg[berg_no_roll].a_englacial[BASE]=0.;

				p_berg[berg_no_roll].basal_position=BASE;
			}
		}				
	}
	else
	{	
		temp_size=p_berg[berg_no_roll].size_i[X];
		p_berg[berg_no_roll].size_i[X]=p_berg[berg_no_roll].size_i[Z];
		p_berg[berg_no_roll].size_i[Z]=p_berg[berg_no_roll].size_i[Y];
		p_berg[berg_no_roll].size_i[Y]=temp_size;

		if ((p_berg[berg_no_roll].basal_position==TOP) || (p_berg[berg_no_roll].basal_position==SIDES) || (p_berg[berg_no_roll].basal_position==BASE)||(p_berg[berg_no_roll].basal_position==NO_BASAL))
		{
		
			temp_basal=p_berg[berg_no_roll].a_basal[FRONT_BACK];
			p_berg[berg_no_roll].a_basal[FRONT_BACK]=p_berg[berg_no_roll].a_basal[SIDES];
			p_berg[berg_no_roll].a_basal[SIDES]=p_berg[berg_no_roll].a_basal[TOP]+p_berg[berg_no_roll].a_basal[BASE];
			p_berg[berg_no_roll].a_basal[TOP]=temp_basal/2.;
			p_berg[berg_no_roll].a_basal[BASE]=temp_basal/2.;
			
			temp_englacial=p_berg[berg_no_roll].a_englacial[FRONT_BACK];
			p_berg[berg_no_roll].a_englacial[FRONT_BACK]=p_berg[berg_no_roll].a_englacial[SIDES];
			p_berg[berg_no_roll].a_englacial[SIDES]=p_berg[berg_no_roll].a_englacial[TOP]+p_berg[berg_no_roll].a_englacial[BASE];
			p_berg[berg_no_roll].a_englacial[TOP]=temp_englacial/2.;
			p_berg[berg_no_roll].a_englacial[BASE]=temp_englacial/2.;

			if (p_berg[berg_no_roll].basal_position==TOP)
			 p_berg[berg_no_roll].basal_position=SIDES;
			else if (p_berg[berg_no_roll].basal_position==SIDES)
				p_berg[berg_no_roll].basal_position=FRONT_BACK;
			else if (p_berg[berg_no_roll].basal_position==BASE)
				 p_berg[berg_no_roll].basal_position=SIDES;

		}
		else if (p_berg[berg_no_roll].basal_position==FRONT_BACK)
		{
			if(ran9(seed9)<0.5)
			{
				temp_basal=p_berg[berg_no_roll].a_basal[FRONT_BACK];
				p_berg[berg_no_roll].a_basal[FRONT_BACK]=p_berg[berg_no_roll].a_basal[SIDES];
				p_berg[berg_no_roll].a_basal[SIDES]=p_berg[berg_no_roll].a_basal[TOP]+p_berg[berg_no_roll].a_basal[BASE];
				p_berg[berg_no_roll].a_basal[TOP]=temp_basal;
				p_berg[berg_no_roll].a_basal[BASE]=0.;

				temp_englacial=p_berg[berg_no_roll].a_englacial[FRONT_BACK];
				p_berg[berg_no_roll].a_englacial[FRONT_BACK]=p_berg[berg_no_roll].a_englacial[SIDES];
				p_berg[berg_no_roll].a_englacial[SIDES]=p_berg[berg_no_roll].a_englacial[TOP]+p_berg[berg_no_roll].a_englacial[BASE];
				p_berg[berg_no_roll].a_englacial[TOP]=0.;
				p_berg[berg_no_roll].a_englacial[BASE]=temp_englacial;
				p_berg[berg_no_roll].basal_position=TOP;
			}
			if(ran9(seed9)>0.5)
			{
				temp_basal=p_berg[berg_no_roll].a_basal[FRONT_BACK];
				p_berg[berg_no_roll].a_basal[FRONT_BACK]=p_berg[berg_no_roll].a_basal[SIDES];
				p_berg[berg_no_roll].a_basal[SIDES]=p_berg[berg_no_roll].a_basal[TOP]+p_berg[berg_no_roll].a_basal[BASE];
				p_berg[berg_no_roll].a_basal[TOP]=0.;
				p_berg[berg_no_roll].a_basal[BASE]=temp_basal;

				temp_englacial=p_berg[berg_no_roll].a_englacial[FRONT_BACK];
				p_berg[berg_no_roll].a_englacial[FRONT_BACK]=p_berg[berg_no_roll].a_englacial[SIDES];
				p_berg[berg_no_roll].a_englacial[SIDES]=p_berg[berg_no_roll].a_englacial[TOP]+p_berg[berg_no_roll].a_englacial[BASE];
				p_berg[berg_no_roll].a_englacial[TOP]=temp_englacial;
				p_berg[berg_no_roll].a_englacial[BASE]=0.;

				p_berg[berg_no_roll].basal_position=BASE;
			}
		}
	}
			 
	/*average concentration of sediment on top of berg in kg/m */
	p_berg[berg_no_roll].conc_sed[TOP]=p_berg[berg_no_roll].a_basal[TOP]*p_glacier.C_basal_i*p_glacier.density_grain
		+p_berg[berg_no_roll].a_englacial[TOP]*p_glacier.C_englacial_i*p_glacier.density_grain;

	/*average concentration of sediment  on sides, front and back of berg in kg/m */
	p_berg[berg_no_roll].conc_sed[SIDES]=(p_berg[berg_no_roll].a_basal[SIDES]+p_berg[berg_no_roll].a_basal[FRONT_BACK])*p_glacier.C_basal_i*p_glacier.density_grain
		+(p_berg[berg_no_roll].a_englacial[SIDES]+p_berg[berg_no_roll].a_englacial[FRONT_BACK])*p_glacier.C_englacial_i*p_glacier.density_grain;

	/*average concentration of sediment  on bottom of berg in kg/m	*/
	p_berg[berg_no_roll].conc_sed[BASE]=p_berg[berg_no_roll].a_basal[BASE]*p_glacier.C_basal_i*p_glacier.density_grain
		+p_berg[berg_no_roll].a_englacial[BASE]*p_glacier.C_englacial_i*p_glacier.density_grain;
}

int stability2(iceberg *p_berg, int berg_no_stab)
{ 
	double lhs, rhs;
	lhs=p_berg[berg_no_stab].size_i[Y]/p_berg[berg_no_stab].size_i[Z];

	rhs=0.7;
	if (lhs>rhs)
	{
		return 0;
		/*iceberg stable*/
	}
	else  /*unstable*/
		return 1;

}

void rollOver2(iceberg *p_berg, int berg_no_roll, glacier p_glacier, long *seed9)
{
	double temp_size, temp_basal, temp_englacial;

	if (p_berg[berg_no_roll].size_i[X]>=p_berg[berg_no_roll].size_i[Z])
	{	
		temp_size=p_berg[berg_no_roll].size_i[Y];
		p_berg[berg_no_roll].size_i[Y]=p_berg[berg_no_roll].size_i[Z];
		p_berg[berg_no_roll].size_i[Z]=temp_size;

		if ((p_berg[berg_no_roll].basal_position==TOP) || (p_berg[berg_no_roll].basal_position==FRONT_BACK) || (p_berg[berg_no_roll].basal_position==BASE)||(p_berg[berg_no_roll].basal_position==NO_BASAL))
		{
			temp_basal=p_berg[berg_no_roll].a_basal[SIDES];
			p_berg[berg_no_roll].a_basal[SIDES]=p_berg[berg_no_roll].a_basal[TOP]+p_berg[berg_no_roll].a_basal[BASE];
			p_berg[berg_no_roll].a_basal[TOP]=temp_basal/2.;
			p_berg[berg_no_roll].a_basal[BASE]=temp_basal/2.;

			temp_englacial=p_berg[berg_no_roll].a_englacial[SIDES];
			p_berg[berg_no_roll].a_englacial[SIDES]=p_berg[berg_no_roll].a_englacial[TOP]+p_berg[berg_no_roll].a_englacial[BASE];
			p_berg[berg_no_roll].a_englacial[TOP]=temp_englacial/2.;
			p_berg[berg_no_roll].a_englacial[BASE]=temp_englacial/2.;
			
			if (p_berg[berg_no_roll].basal_position==TOP)
				p_berg[berg_no_roll].basal_position=SIDES;

			else if (p_berg[berg_no_roll].basal_position==BASE)
				p_berg[berg_no_roll].basal_position=SIDES;
				
		}

		else if (p_berg[berg_no_roll].basal_position==SIDES)
		{
			if(ran9(seed9)<0.5)
			{
				temp_basal=p_berg[berg_no_roll].a_basal[SIDES];
				p_berg[berg_no_roll].a_basal[SIDES]=p_berg[berg_no_roll].a_basal[TOP]+p_berg[berg_no_roll].a_basal[BASE];
				p_berg[berg_no_roll].a_basal[TOP]=temp_basal;
				p_berg[berg_no_roll].a_basal[BASE]=0.;

				temp_englacial=p_berg[berg_no_roll].a_englacial[SIDES];
				p_berg[berg_no_roll].a_englacial[SIDES]=p_berg[berg_no_roll].a_englacial[TOP]+p_berg[berg_no_roll].a_englacial[BASE];
				p_berg[berg_no_roll].a_englacial[TOP]=0.;
				p_berg[berg_no_roll].a_englacial[BASE]=temp_englacial;
				p_berg[berg_no_roll].basal_position=TOP;
			}
			else if(ran9(seed9)>0.5)
			{
				temp_basal=p_berg[berg_no_roll].a_basal[SIDES];
				p_berg[berg_no_roll].a_basal[SIDES]=p_berg[berg_no_roll].a_basal[TOP]+p_berg[berg_no_roll].a_basal[BASE];
				p_berg[berg_no_roll].a_basal[TOP]=0.;
				p_berg[berg_no_roll].a_basal[BASE]=temp_basal;

				temp_englacial=p_berg[berg_no_roll].a_englacial[SIDES];
				p_berg[berg_no_roll].a_englacial[SIDES]=p_berg[berg_no_roll].a_englacial[TOP]+p_berg[berg_no_roll].a_englacial[BASE];
				p_berg[berg_no_roll].a_englacial[TOP]=temp_englacial;
				p_berg[berg_no_roll].a_englacial[BASE]=0.;

				p_berg[berg_no_roll].basal_position=BASE;
			}
		}				
	}

	else
	{	
		temp_size=p_berg[berg_no_roll].size_i[Y];
		p_berg[berg_no_roll].size_i[Y]=p_berg[berg_no_roll].size_i[X];
		p_berg[berg_no_roll].size_i[X]=p_berg[berg_no_roll].size_i[Z];
		p_berg[berg_no_roll].size_i[Z]=temp_size;

		if ((p_berg[berg_no_roll].basal_position==TOP) || (p_berg[berg_no_roll].basal_position==FRONT_BACK) || (p_berg[berg_no_roll].basal_position==BASE)||(p_berg[berg_no_roll].basal_position==NO_BASAL))
		{
			temp_basal=p_berg[berg_no_roll].a_basal[SIDES];
			p_berg[berg_no_roll].a_basal[SIDES]=p_berg[berg_no_roll].a_basal[FRONT_BACK];
			p_berg[berg_no_roll].a_basal[FRONT_BACK]=p_berg[berg_no_roll].a_basal[TOP]+p_berg[berg_no_roll].a_basal[BASE];
			p_berg[berg_no_roll].a_basal[TOP]=temp_basal/2.;
			p_berg[berg_no_roll].a_basal[BASE]=temp_basal/2.;
			
			temp_englacial=p_berg[berg_no_roll].a_englacial[SIDES];
			p_berg[berg_no_roll].a_englacial[SIDES]=p_berg[berg_no_roll].a_englacial[FRONT_BACK];
			p_berg[berg_no_roll].a_englacial[FRONT_BACK]=p_berg[berg_no_roll].a_englacial[TOP]+p_berg[berg_no_roll].a_englacial[BASE];
			p_berg[berg_no_roll].a_englacial[TOP]=temp_englacial/2.;
			p_berg[berg_no_roll].a_englacial[BASE]=temp_englacial/2.;

			if (p_berg[berg_no_roll].basal_position==TOP)
			 p_berg[berg_no_roll].basal_position=FRONT_BACK;
			else if (p_berg[berg_no_roll].basal_position==FRONT_BACK)
				p_berg[berg_no_roll].basal_position=SIDES;
			else if (p_berg[berg_no_roll].basal_position==BASE)
				 p_berg[berg_no_roll].basal_position=FRONT_BACK;

		}
		else if (p_berg[berg_no_roll].basal_position==SIDES)
		{
			if(ran9(seed9)<0.5)
			{
				temp_basal=p_berg[berg_no_roll].a_basal[SIDES];
				p_berg[berg_no_roll].a_basal[SIDES]=p_berg[berg_no_roll].a_basal[FRONT_BACK];
				p_berg[berg_no_roll].a_basal[FRONT_BACK]=p_berg[berg_no_roll].a_basal[TOP]+p_berg[berg_no_roll].a_basal[BASE];
				p_berg[berg_no_roll].a_basal[TOP]=temp_basal;
				p_berg[berg_no_roll].a_basal[BASE]=0.;

				temp_englacial=p_berg[berg_no_roll].a_englacial[SIDES];
				p_berg[berg_no_roll].a_englacial[SIDES]=p_berg[berg_no_roll].a_englacial[FRONT_BACK];
				p_berg[berg_no_roll].a_englacial[FRONT_BACK]=p_berg[berg_no_roll].a_englacial[TOP]+p_berg[berg_no_roll].a_englacial[BASE];
				p_berg[berg_no_roll].a_englacial[TOP]=0.;
				p_berg[berg_no_roll].a_englacial[BASE]=temp_englacial;

				p_berg[berg_no_roll].basal_position=TOP;
			}
			if(ran9(seed9)>0.5)
			{
				temp_basal=p_berg[berg_no_roll].a_basal[SIDES];
				p_berg[berg_no_roll].a_basal[SIDES]=p_berg[berg_no_roll].a_basal[FRONT_BACK];
				p_berg[berg_no_roll].a_basal[FRONT_BACK]=p_berg[berg_no_roll].a_basal[TOP]+p_berg[berg_no_roll].a_basal[BASE];
				p_berg[berg_no_roll].a_basal[TOP]=0.;
				p_berg[berg_no_roll].a_basal[BASE]=temp_basal;

				temp_englacial=p_berg[berg_no_roll].a_englacial[SIDES];
				p_berg[berg_no_roll].a_englacial[SIDES]=p_berg[berg_no_roll].a_englacial[FRONT_BACK];
				p_berg[berg_no_roll].a_englacial[FRONT_BACK]=p_berg[berg_no_roll].a_englacial[TOP]+p_berg[berg_no_roll].a_englacial[BASE];
				p_berg[berg_no_roll].a_englacial[TOP]=temp_englacial;
				p_berg[berg_no_roll].a_englacial[BASE]=0.;

				p_berg[berg_no_roll].basal_position=BASE;
			}
		}
	}
		 
	/*average concentration of sediment on top of berg in kg/m */
	p_berg[berg_no_roll].conc_sed[TOP]=p_berg[berg_no_roll].a_basal[TOP]*p_glacier.C_basal_i*p_glacier.density_grain
		+p_berg[berg_no_roll].a_englacial[TOP]*p_glacier.C_englacial_i*p_glacier.density_grain;

	/*average concentration of sediment  on sides, front and back of berg in kg/m */
	p_berg[berg_no_roll].conc_sed[SIDES]=(p_berg[berg_no_roll].a_basal[SIDES]+p_berg[berg_no_roll].a_basal[FRONT_BACK])*p_glacier.C_basal_i*p_glacier.density_grain
		+(p_berg[berg_no_roll].a_englacial[SIDES]+p_berg[berg_no_roll].a_englacial[FRONT_BACK])*p_glacier.C_englacial_i*p_glacier.density_grain;

	/*average concentration of sediment  on bottom of berg in kg/m	*/
	p_berg[berg_no_roll].conc_sed[BASE]=p_berg[berg_no_roll].a_basal[BASE]*p_glacier.C_basal_i*p_glacier.density_grain
		+p_berg[berg_no_roll].a_englacial[BASE]*p_glacier.C_englacial_i*p_glacier.density_grain;
}

#undef BASALMELTWATER
#undef DEBRISMELTAIR



