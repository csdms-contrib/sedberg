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

#include "fjord.h"		

glacier inputGlacier()
{
	FILE *in_glacier;
	glacier p_glacier;
	char temp[120];
	/*GLACIER CHARACTERISTICS file */

	in_glacier = fopen("Input/glacier.txt", "r"); 

	if (in_glacier == NULL)
	{
		printf("Error: The file 'glacier.txt' could not be opened\n");
		exit(1);
	}
	else  
	{
		fscanf(in_glacier, "%lf	%s	%lf	%s	%lf	%s	%lf	%s	%lf	%s	%lf	%s	%lf	%s	%lf	%s	%lf	%s	%lf	%s	%lf	%s	%lf	%s	%lf	%s	%lf	%s	%d	%s	%lf	%s	%lf	%s	%lf",
			&p_glacier.latitude,temp,&p_glacier.h_front,temp,&p_glacier.w_front,temp,&p_glacier.w_fjord,temp,&p_glacier.v_front,
			temp,&p_glacier.h_water,temp,&p_glacier.h_basal_sediment,temp,
			&p_glacier.C_basal_i,temp,&p_glacier.C_englacial_i,temp,&p_glacier.y_f, temp, &p_glacier.slope,temp,
			&p_glacier.density_grain, temp, &p_glacier.density_dep, temp, &p_glacier.BrownsConst,temp, 
			&p_glacier.calving_scenario,temp, &p_glacier.annual_berg_vol,temp, &p_glacier.mu,temp, &p_glacier.sigma);
	}
	fclose(in_glacier);	  /*Close "glacier.txt"*/
	p_glacier.initial_snout_position =  p_glacier.y_f;
	
	return p_glacier;

}

/**********************************************************************************************/

options inputOptions()
{
	options p_options;
	FILE *inputoptions;
	char temp[120];

	/*open OPTIONS file	*/
	   
	inputoptions = fopen ("Input/options.txt","r"); 

	if (inputoptions == NULL)
	{   		
   	
		printf("Error: The file 'options.txt' could not be opened\n");
		exit(1);
	}
  
	else  
	{	

		fscanf(inputoptions, "%lf	%s	%lf	%s	%lf	%s	%lf	%s	%d", 
			&p_options.totnyears,temp, &p_options.spinuptime,temp, &p_options.dx,temp, &p_options.dz,temp, &p_options.print_sediment);
			
	}
	/*Close "options.txt" */
printf("p_options.totnyears=%lf\n",p_options.totnyears);
	fclose(inputoptions);

	return p_options;
}

/**********************************************************************************************/

bathy inputBathy(glacier p_glacier, options p_options) 
{
	FILE *inputbathy;
	bathy fjord;
	int checkSnoutPosition = 0;
	int o;

	/*BATHYMETRY file */

 	inputbathy = fopen("Input/bathy.txt", "r"); 

	if (inputbathy == NULL)
	{
		printf("Error: The file 'bathy.txt' could not be opened\n");
		exit(1);
	}
  
	else  
	{
		fjord.noLines=0;
		fscanf(inputbathy, "%lf	%lf", &fjord.bottom_x[fjord.noLines],&fjord.bottom_z[fjord.noLines]);
	
		while(!feof(inputbathy))
		{	
			fjord.noLines++; 
			
			if (fjord.noLines>100)  
			{
				printf("Error: The file 'bathy.txt' has too many lines of data -\n");
				printf("you must increase the size of array 'noLines' in program\n");
				exit(1);
			}

			fscanf(inputbathy, "%lf	%lf", &fjord.bottom_x[fjord.noLines],&fjord.bottom_z[fjord.noLines]);
	//		printf ("bottom_x= %lf\nbottom_z= %lf\nno_lines=%d\n", fjord.bottom_x[fjord.noLines],fjord.bottom_z[fjord.noLines],fjord.noLines);
//#		printf("Number of bathy read: %d\nx=%lf	z=%lf\n", fjord.noLines, fjord.bottom_x[fjord.noLines],fjord.bottom_z[fjord.noLines]);

		
		}
	//	printf("Number of bathy read: %d\nx=%lf	z=%lf\n", fjord.noLines, fjord.bottom_x[fjord.noLines],fjord.bottom_z[fjord.noLines]);
	}
	
	fclose(inputbathy);	/*Close "Bathy.txt"	*/


/*DEFINE number of x and z bins	   */
		 
	fjord.max_depth=0.;
	fjord.max_length=0.;

	for (o=0;o<fjord.noLines;o++)
	{
	//	printf("bottom depth %lf\n", fjord.bottom_z[o]);
		if (fjord.bottom_z[o]>fjord.max_depth)
		{
			fjord.max_depth= fjord.bottom_z[o];
		}
		if(-fjord.bottom_x[o]>fjord.max_length)
		{
		   fjord.max_length= fjord.bottom_x[o];
		}
/*		printf(	"p_glacier.y_f = %lf\n", p_glacier.y_f);		  */
		if (p_glacier.y_f==fjord.bottom_x[o])
		{
		   checkSnoutPosition = 1;
		   p_glacier.bottomDepth = fjord.bottom_z[o];
  //printf("stored bottom depth %lf\n", p_glacier.bottomDepth );	 
		}
		if ((o==fjord.noLines-1) && (checkSnoutPosition==0))
		{
			printf("Error: there must be a depth value (z) in the file bathy.txt\n");
			printf("which corresponds to the initial position (x) of the glacier snout\n");
			exit(1);
		}
	}	
		
 //	printf ("max depth= %lf\n max length= %lf\np_glacier.y_f= %lf\np_glacier.h_front= %lf\np_glacier.bottomDepth= %lf\np_options.dz= %lf\n",
//		fjord.max_depth, fjord.max_length, p_glacier.y_f,p_glacier.h_front,p_glacier.bottomDepth,p_options.dz);	   

	fjord.no_x_bins = -(int)(ceil(fjord.max_length/p_options.dx)); 
	fjord.no_y_bins =	(int)(ceil(p_glacier.w_fjord/p_options.dx));
//	fjord.no_z_bins = (int)(ceil((fjord.max_depth+(p_glacier.y_f*tan(p_glacier.slope*2.*PI/360.))+(p_glacier.h_front-p_glacier.bottomDepth))/p_options.dz));
	fjord.no_z_bins = (int)(ceil((p_glacier.h_water-p_glacier.bottomDepth+fjord.max_depth)/p_options.dz));


//	fprintf (stdout,"no_y_bins= %ld\nno_x_bins= %ld\n no_z_bins=%ld\nfjord.max_depth=%lf\np_glacier.w_fjord=%lf\n", fjord.no_y_bins, fjord.no_x_bins,fjord.no_z_bins,fjord.max_depth,p_glacier.w_fjord);
//	fflush(stdout);
	fprintf (out_runtime,"no_y_bins= %ld\nno_x_bins= %ld\n no_z_bins=%ld\nfjord.max_depth=%lf\np_glacier.w_fjord=%lf\n", fjord.no_y_bins, fjord.no_x_bins,fjord.no_z_bins,fjord.max_depth,p_glacier.w_fjord);
	fflush(out_runtime);
	
	return fjord;
}

/**********************************************************************************************/
water_vel inputWater_vel()
{
	FILE *inputWater_vel;
	water_vel p_water_vel;
	char temp[80];

		
	inputWater_vel = fopen ("Input/waterVelConst.txt","r"); 

	if (inputWater_vel == NULL)
	{   		
   	
		printf("Error: The file 'waterVelConst.txt' could not be opened\n");
		exit(1);
	}
  
	else  
	{	

		fscanf(inputWater_vel, "%lf	%s	%lf	%s	%lf	%s	%lf	%s	%lf", 
			&p_water_vel.amplitude,temp,&p_water_vel.Z_var,temp,&p_water_vel.T_var,temp,&p_water_vel.delta_v,temp,&p_water_vel.eta_v);
			
	}


	fclose(inputWater_vel);

//	fprintf(stdout,"p_water_vel.eta_v=%lf\n",p_water_vel.eta_v);
//	fflush(stdout);
	
	fprintf(out_runtime,"p_water_vel.eta_v=%lf\n",p_water_vel.eta_v);
	fflush(out_runtime);

	return p_water_vel;
}

sediment initialiseSediment(bathy fjord, options p_options, glacier p_glacier)
{
	sediment p_sediment;
	//int row, r,c,o,l,m,d;
//	int x_min,x_max,r,c,m,o,l;
	int r,c;
	/*SEDIMENT*/
	p_sediment.thickness= (double **) malloc(fjord.no_x_bins * sizeof(double *) );
 
	for (r = 0; r < fjord.no_x_bins; r++) 
	{
		p_sediment.thickness[r]= malloc (fjord.no_y_bins * sizeof(double));
	}

/*INITIALISE SEDIMENT ARRAY*/
	for (r=0;r<fjord.no_x_bins;r++)
	{	
		for (c=0;c<fjord.no_y_bins;c++)
		{ 
   		  	p_sediment.thickness[r][c]=0.;	
		}
	}
					
	p_sediment.water_z_min = 0;

	return p_sediment;
}		  	


void inputSolarRadiation(double solarRad[12])
{
	FILE *in_solar;

	/*solar radiation file - monthly mean climatology values in W/m^2 */
	in_solar = fopen("Input/solrad.txt", "r"); 

	if (in_solar == NULL)
	{
		printf("Error: The file 'solrad.txt' could not be opened\n");
		exit(1);
	}
	else  
	{
		fscanf(in_solar, "%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf",
			&solarRad[JAN],&solarRad[FEB],&solarRad[MAR],&solarRad[APR],&solarRad[MAY],&solarRad[JUN],
			&solarRad[JUL],&solarRad[AUG],&solarRad[SEP],&solarRad[OCT],&solarRad[NOV],&solarRad[DEC]);
	}
	fclose(in_solar);	 

	fprintf(out_runtime,"solar radiation (Dec)=%lf\n",solarRad[DEC]);
	fflush(out_runtime);
}

void inputAirTemp(double airTemp[12])
{
	FILE *in_airtemp;

	/*air temperature file - monthly mean climatology values in deg C */

	in_airtemp = fopen("Input/airtemp.txt", "r"); 

	if (in_airtemp == NULL)
	{
		printf("Error: The file 'airtemp.txt' could not be opened\n");
		exit(1);
	}
	else  
	{
		fscanf(in_airtemp, "%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf",
			&airTemp[JAN],&airTemp[FEB],&airTemp[MAR],&airTemp[APR],&airTemp[MAY],&airTemp[JUN],
			&airTemp[JUL],&airTemp[AUG],&airTemp[SEP],&airTemp[OCT],&airTemp[NOV],&airTemp[DEC]);

	}
	fclose(in_airtemp);
	
	fprintf(out_runtime, "air temp (Dec)=%lf\n",airTemp[DEC]);
	fflush(out_runtime);

}
/*Reads in data from seatemp_winter.txt then interpolates data to produce 
seaTempSummer - an array of water temp for each grid cell*/
void inputSeaTempWinter(double *seaTemp, bathy fjord, options p_options, sediment p_sediment)
{
	FILE *water_temp_f;
	int l,m,o, noLines_t, z_min_t = 0,z_max_t = 0;
	double water_temp_z[100];	
	double water_temp_t[100];
	double dist_temp = 0.;

//	WATER TEMPERATURE file

	water_temp_f = fopen("Input/seatemp_winter.txt", "r"); 
	if (water_temp_f == NULL)
	{
		printf("Error: The file 'seatemp_winter.txt' could not be opened\n");
		exit(1);
	}
	else  
	{
		noLines_t=0;  

		fscanf(water_temp_f, "%lf	%lf", &water_temp_z[noLines_t],&water_temp_t[noLines_t]);
		while(!feof(water_temp_f))
		{	
			noLines_t++; 	
			if (noLines_t>100)  
			{
				printf("Error: The file 'seatemp_winter.txt' has too many lines of data -\nyou must increase the size of array 'noLines_t' in program\n");
				exit(1);
			}
			fscanf(water_temp_f, "%lf	%lf", &water_temp_z[noLines_t],&water_temp_t[noLines_t]);

		}
	}
	fclose(water_temp_f);	 

	   fprintf(out_runtime,"winter sea temp 1st depth=%lf	temp=%lf\n",water_temp_z[0],water_temp_t[0]);
	   fflush(out_runtime);

	   fprintf(out_runtime,"winter sea temp final depth=%lf	temp=%lf\n",water_temp_z[noLines_t],water_temp_t[noLines_t]);
	   fflush(out_runtime);

/*initialise water temp array		 */	
	for (m=0;m<fjord.no_z_bins;m++)
		seaTemp[m] = 0.0;
		 	
/*INTERPOLATE WATER TEMPERATURE	   */  //	 set a value for water temperature in each vertical bin	

	for (o=0;o<noLines_t;o++)
		{		  
			z_min_t = (int)(ceil(water_temp_z[o]/p_options.dz));
			z_max_t = (int)(ceil(water_temp_z[o+1]/p_options.dz)); //water_temp_z[noLines_t]/dz	

			if (z_max_t>=fjord.no_z_bins)
			{
				printf("The water temperature data is too deep for this fjord\n");
				exit(1);
			}

			dist_temp = z_max_t-z_min_t;

			for (l=z_min_t;l<z_max_t;l++)
			{	
			   //interpolate sea temperature in water_temp so that there is a value at each grid cell
				 seaTemp[l] = water_temp_t[o]+((l-z_min_t)*(water_temp_t[o+1]-water_temp_t[o])/dist_temp);
			}	 
		}

	//Inserts value of water velocity at maximum depth in input file to every depth greater than this in fjord

	z_max_t = (int)(ceil(water_temp_z[noLines_t]/p_options.dz));

	for (l=z_max_t;l<fjord.no_z_bins;l++)
	{
		seaTemp[l]=water_temp_t[noLines_t];	
	}	


}

 /*Reads in data from seatemp_summer.txt then interpolates data to produce 
seaTempSummer - an array of water temp for each grid cell*/
void inputSeaTempSummer(double *seaTemp, bathy fjord, options p_options, sediment p_sediment)
{
	FILE *water_temp_f;
	int l,m,o, noLines_t, z_min_t = 0,z_max_t = 0;
	double water_temp_z[100];	
	double water_temp_t[100];
	double dist_temp = 0.;

//	WATER TEMPERATURE file

	water_temp_f = fopen("Input/seatemp_summer.txt", "r"); 

	if (water_temp_f == NULL)
	{
		printf("Error: The file 'seatemp_summer.txt' could not be opened\n");
		exit(1);
	}
	else  
	{
		noLines_t=0;  

		fscanf(water_temp_f, "%lf	%lf", &water_temp_z[noLines_t],&water_temp_t[noLines_t]);
		while(!feof(water_temp_f))
		{	
			noLines_t++; 	
			if (noLines_t>100)  
			{
				printf("Error: The file 'seatemp_summer.txt' has too many lines of data -\nyou must increase the size of array 'noLines_t' in program\n");
				exit(1);
			}
			fscanf(water_temp_f, "%lf	%lf", &water_temp_z[noLines_t],&water_temp_t[noLines_t]);
		}
	}
	fclose(water_temp_f);	 

	fprintf(out_runtime,"summer sea temp 1st depth=%lf	temp=%lf\n",water_temp_z[0],water_temp_t[0]);
   	fflush(out_runtime);
	
	fprintf(out_runtime,"summer sea temp final depth=%lf	temp=%lf\n", water_temp_z[noLines_t], water_temp_t[noLines_t]);
	fflush(out_runtime); 
 
/*initialise water temp array		 */	
	for (m=0;m<fjord.no_z_bins;m++)
		seaTemp[m] = 0.0;
		
/*INTERPOLATE WATER TEMPERATURE	   */  //	 set a value for water temperature in each vertical bin	

	for (o=0;o<noLines_t;o++)
	{	
		z_min_t = (int)(ceil(water_temp_z[o]/p_options.dz));
		z_max_t = (int)(ceil(water_temp_z[o+1]/p_options.dz)); //water_temp_z[noLines_t]/dz	

		if (z_max_t>=fjord.no_z_bins)
		{
			printf("The water temperature data is too deep for this fjord\n");
			exit(1);
		}

		dist_temp = z_max_t-z_min_t;

		for (l=z_min_t;l<z_max_t;l++)
		{	
		   //interpolate sea temperature in water_temp so that there is a value at each grid cell
			 seaTemp[l] = water_temp_t[o]+((l-z_min_t)*(water_temp_t[o+1]-water_temp_t[o])/dist_temp);
		}	 
	}

	//Inserts value of water velocity at maximum depth in input file to every depth greater than this in fjord

	z_max_t = (int)(ceil(water_temp_z[noLines_t]/p_options.dz));

	for (l=z_max_t;l<fjord.no_z_bins;l++)
	{
		seaTemp[l]=water_temp_t[noLines_t];	
	}	
}


void inputSeaSalin(double *seaSalin, bathy fjord, options p_options, sediment p_sediment)
{
	FILE *water_salin_f;
	int l,m,o, noLines_t, z_min_t = 0,z_max_t = 0;
	double water_salin_z[100];	
	double water_salin_t[100];
	double dist_salin = 0.;

//	WATER SALINITY file

	water_salin_f = fopen("Input/salin.txt", "r"); 

	if (water_salin_f == NULL)
	{
		printf("Error: The file 'salin.txt' could not be opened\n");
		exit(1);
	}
	else  
	{
		noLines_t=0;  
		fscanf(water_salin_f, "%lf	%lf", &water_salin_z[noLines_t],&water_salin_t[noLines_t]);
		while(!feof(water_salin_f))
		{	
			noLines_t++; 	
			if (noLines_t>100)  
			{
				printf("Error: The file 'water_salin.txt' has too many lines of data -\nyou must increase the size of array 'noLines_t' in program\n");
				exit(1);
			}
			fscanf(water_salin_f, "%lf	%lf", &water_salin_z[noLines_t],&water_salin_t[noLines_t]);
		}		
	}
	fclose(water_salin_f);	 

   fprintf(out_runtime,"salin 1st depth=%lf	salin=%lf\n",water_salin_z[0], water_salin_t[0]);
   fflush(out_runtime);

   fprintf(out_runtime,"salin final depth=%lf	salin=%lf\n",water_salin_z[noLines_t], water_salin_t[noLines_t]);
   fflush(out_runtime);
	   
/*initialise water salin array		 */	
	for (m=0;m<fjord.no_z_bins;m++)
		seaSalin[m] = 0.0;		
   	
/*INTERPOLATE WATER SALINITY	   */  //	 set a value for water temperature in each vertical bin	

	for (o=0;o<noLines_t;o++)
		{	  
			z_min_t = (int)(ceil(water_salin_z[o]/p_options.dz));
			z_max_t = (int)(ceil(water_salin_z[o+1]/p_options.dz)); //water_salin_z[noLines_t]/dz	

			if (z_max_t>=fjord.no_z_bins)
			{
				printf("The water salinity data is too deep for this fjord\n");
				exit(1);
			}

			dist_salin = z_max_t-z_min_t;

			for (l=z_min_t;l<z_max_t;l++)
			{	
			   //interpolate sea salinerature in water_salin so that there is a value at each grid cell
				 seaSalin[l] = water_salin_t[o]+((l-z_min_t)*(water_salin_t[o+1]-water_salin_t[o])/dist_salin);
			}	 
		}
	//Inserts value of water velocity at maximum depth in input file to every depth greater than this in fjord
	z_max_t = (int)(ceil(water_salin_z[noLines_t]/p_options.dz));

	for (l=z_max_t;l<fjord.no_z_bins;l++)
	{
		seaSalin[l]=water_salin_t[noLines_t];	
	}	

}



void calcMonth(double time, int *month)
{
	if (fmod(time,365.)>=0. && fmod(time,365.)<=31.)
			*month=JAN;
	else if (fmod(time,365.)>31. && fmod(time,365.)<=59.)
			*month=FEB;

	else if (fmod(time,365.)>59. && fmod(time,365.)<=90.)
			*month=MAR;

	else if (fmod(time,365.)>90. && fmod(time,365.)<=120.)
			*month=APR;

	else if (fmod(time,365.)>120. && fmod(time,365.)<=151.)
			*month=MAY;

	else if (fmod(time,365.)>151. && fmod(time,365.)<=181.)
			*month=JUN;

	else if (fmod(time,365.)>181. && fmod(time,365.)<=212.)
			*month=JUL;

	else if (fmod(time,365.)>212. && fmod(time,365.)<=243.)
			*month=AUG;

	else if (fmod(time,365.)>243. && fmod(time,365.)<=273.)
			*month=SEP;
								 
	else if (fmod(time,365.)>273. && fmod(time,365.)<=304.)
			*month=OCT;

	else if (fmod(time,365.)>304. && fmod(time,365.)<=334.)
			*month=NOV;

	else if (fmod(time,365.)>334. && fmod(time,365.)<=365.)
			*month=DEC;

}

void freezingpointwater(double salin,double *freezingpt, double press)
{
	double a0 = -0.0575, a1 = 1.710523e-3, a2 = -2.154996e-4, b  = -7.53e-4;
	 //  salin = salinity    [psu      (PSS-78)]
	//   press = pressure    [db] / m
		
	*freezingpt = (a0*salin + a1*salin*sqrt(salin) + a2*salin*salin + b*press) / 1.00024;
}

