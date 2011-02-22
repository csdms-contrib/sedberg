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
#include <string.h>
#include "iceberg.h"
#include "constants.h"

int main()
{

	double endTime = 0.;	//maximum time in days
	double n_years;
    double fracyear;
	double winterdt=0.;   //winter timestep in hours
 	double velchange_dt=0.083;   //timestep that wind and water vel changed in hours

	double fracvelcalc;
	double n_velcalc;
	double frac_calving;
	double n_calving;
	double calvingdt=0.;//calving timestep in days
	double calvingdt_winter= 0.;//calving timestep in days
	double calvingdt_summer=0.;//calving timestep in days

	double time = 0.; /*time  in days*/
	double oldtime = 0.;
	double time_i=0.;
	int loop_bergs;
	int berg_no = -1; /*Berg array index - zero for the first berg*/
	int total_no_bergs = 0;
 	int n_bergs = 10000;
	double tempvol = 0.,totvol=0., theoreticalvol=0.,theoreticalnobergs=0.;
	double random_calving = 0.0;  /*initialise random number			 */

	iceberg *p_berg; 
	glacier p_glacier;
	bathy fjord;
	options p_options;
	sediment p_sediment;
	water_vel p_water_vel;

	int i=0;
	int r=0;
	int c=0;
	int del_berg_loop;

	double timefrac;
	double n_time;

	int no_bad_bergs = 0;

	FILE *out_sed;
	FILE *out_sed2;
	FILE *out_sed_avalongfjord;
	FILE *out_sed_avcrossfjord;

	double ystart[4];	
	double h1=.03;	// initial time step in hours
	long seed,seed1,seed2, seed3,seed4, seed5, seed6, seed7, seed8, seed9, seed10, seed11, seed12, seed_stayinfjord;

	double vel_random[2];
	double meltForcedC_w;
	double meltVert_w;
	double meltSens_a;
	double meltSolar_a;
	int max_depth_berg;
	double solarRad[12], airTemp[12], *seaTemp,*seaTempSummer,*seaTempWinter,avTempWater, *seaSalin;
	int month, loop_temp;
	double p_v_w_diff_m[2],p_v_w_diff_mag_m;
	double *tot_crossfjord, *tot_alongfjord;
	int min_av, max_av;	 
	double width_av;
	double tempgtfreez;
	double freezingpt, avSalinWater;
	long no_bergs_left_fjord=0;
	double residence_time_sum=0.;
	double residence_time=0.;
	double prob_berg_stays_overwinter=0.;
	double no_bergs_end_summer=0.;
	int calc_prob=0;
	double volume_sed_out_fjord=0.;
	int endwarmup=0;
	double meansize=0.;
	double av_dep_rate=0.;

	vol_sed_deposited=0.;

	totbasal_berg=0., totbasal_glacier=0.,toteng_berg=0., toteng_glacier=0.;

	for (i=0;i<4;i++)
	{
		ystart[i]=0.0;	
		yp[i]=0.;
	}

/*open runtime output file*/	
  	out_runtime = fopen("Output/logfile.txt","w");  
	if(out_runtime==NULL) 
	{
	   printf("Error: can't open file logfile.txt\n");
 	   /* DON'T PASS A NULL POINTER TO fclose !! */
		exit(1);
	}


/* end of variable definitions; start of prog*/	
/* INITIALISE VARIABLES	 */

/*
	Random number generator seeds
*/
	seed = -604;
	seed1 = -692;
	seed2 = -386;
	seed3 = -491;
	seed4 = -811;
	seed5 = -260;
	seed6 = -137;
	seed7 = -595;
	seed8 = -379;
	seed9 = -181;
	seed10 = -555;
	seed11 = -204;
	seed12 = -376;
	seed_stayinfjord = -456;
										  
	fprintf(stdout,"SedBerg, Copyright (C) 2010 R. I. Mugford\nSedBerg comes with ABSOLUTELY NO WARRANTY;\nThis is free software, and you are welcome\nto redistribute it under certain conditions\n");
	fflush(stdout);
 										  
	fprintf(out_runtime,"SedBerg, Copyright (C) 2010 R. I. Mugford\nSedBerg comes with ABSOLUTELY NO WARRANTY;\nThis is free software, and you are welcome\nto redistribute it under certain conditions\n");
	fflush(out_runtime);
		
/*INPUT FILES */	 
	p_glacier = inputGlacier();	 
	p_options = inputOptions();	
	fjord = inputBathy(p_glacier,p_options);
	p_sediment = initialiseSediment(fjord, p_options, p_glacier);
	p_water_vel = inputWater_vel();
						
	endTime = p_options.spinuptime+p_options.totnyears*365.;	//maximum time in days

	seaTemp = (double *) malloc( fjord.no_z_bins* sizeof(double) );
  	seaTempSummer= (double *) malloc( fjord.no_z_bins* sizeof(double) );
	seaTempWinter= (double *) malloc( fjord.no_z_bins* sizeof(double) );
	seaSalin= (double *) malloc( fjord.no_z_bins* sizeof(double) );

	tot_crossfjord  = (double *) malloc( fjord.no_y_bins * sizeof(double) );
	if(tot_crossfjord==NULL)
	{
		printf("Error: Couldn't allocate memory for tot_crossfjord\n");
		exit(1);
	}

	tot_alongfjord  = (double *) malloc( fjord.no_x_bins * sizeof(double) );
	if(tot_alongfjord==NULL)
	{
		printf("Error: Couldn't allocate memory for tot_alongfjord\n");
		exit(1);
	}

	for (i=0;i<fjord.no_y_bins;i++)
	{
		tot_crossfjord[i]=0.;
	}

    for (i=0;i<fjord.no_x_bins;i++)
	{
		tot_alongfjord[i]=0.;
	}

	inputSolarRadiation(solarRad);						
   	inputAirTemp(airTemp);
	inputSeaTempSummer(seaTempSummer, fjord, p_options, p_sediment);
	inputSeaTempWinter(seaTempWinter, fjord, p_options, p_sediment);

	inputSeaSalin(seaSalin, fjord, p_options, p_sediment);

/*ASSIGN POINTERS */
		
/*BERG*/

	p_berg = (iceberg*) malloc(n_bergs*sizeof(iceberg));

	if(p_berg == NULL)
    {
	    printf("out of memory\n");
        exit(1);
    }
  
	F_CORIOLIS=2.*OMEGA_0*sin(PI*p_glacier.latitude/180.); 
	 
	fprintf(stdout,"no years sed deposited over (%.0lf day spin up added to this) =	%lf\n",p_options.spinuptime,p_options.totnyears);
	fprintf(out_runtime,"no years sed deposited over (%.0lf day spin up added to this) =	%lf\n",p_options.spinuptime,p_options.totnyears);

//	fprintf(stdout,"F_CORIOLIS %lf\n",F_CORIOLIS);
	fprintf(out_runtime,"F_CORIOLIS %lf\n",F_CORIOLIS);

	mean_berg_size = exp(p_glacier.mu+0.5*p_glacier.sigma*p_glacier.sigma);
	mean_berg_area = exp(2.*p_glacier.mu+2.*p_glacier.sigma*p_glacier.sigma);
	mean_berg_volume= exp(3.*p_glacier.mu+4.5*p_glacier.sigma*p_glacier.sigma);
	
//	fprintf(stdout,"begin h=	%lf hours\n",h1);	 										  
	fprintf(stdout,"mean berg size=	%lf\n",mean_berg_size );
	fflush(stdout);
	
	fprintf(out_runtime,"begin h=	%lf hours\n",h1);	 										  
	fprintf(out_runtime,"mean berg size=	%lf\n",mean_berg_size );
	fflush(out_runtime);
	
	theoreticalvol = p_glacier.annual_berg_vol;
	theoreticalnobergs= (p_glacier.annual_berg_vol*1000000000.)/mean_berg_volume;

	time=time_i;
	tp=time_i;
	
	vel_random[0]=ran6(&seed6);
	vel_random[1]=ran7(&seed7);
	waterVel(p_v_w, &p_v_w_mag,vel_random);
	airVel(p_v_a, &p_v_a_mag, &seed4, &seed5);

	if (p_glacier.calving_scenario==0) //seasonal calving rate variation scenario 6 times higher in summer than in winter
	{
		calvingdt_summer = dtCalving(p_glacier, (9./4.)*p_glacier.annual_berg_vol);
		calvingdt_winter = dtCalving(p_glacier, (3./8.)*p_glacier.annual_berg_vol);
		winterdt=0.3;   //winter timestep in hours
	}
	else if (p_glacier.calving_scenario==1)  //seasonal calving rate variation scenario summer calving only
	{
		calvingdt_winter=0.;
		calvingdt_summer = dtCalving(p_glacier, p_glacier.annual_berg_vol*3.);
		winterdt=0.5;   //winter timestep in hours
	}
	else if (p_glacier.calving_scenario==2)  //seasonal calving rate variation scenario sikussak breakup scenario: all at once
	{
		calvingdt_summer = dtCalving(p_glacier, (p_glacier.annual_berg_vol*10.*365./28.));
		calvingdt_winter=calvingdt_summer;
		winterdt=0.5;

	}
	else if (p_glacier.calving_scenario==3)  //seasonal calving rate variation scenario constant all year round
	{
		calvingdt_summer = dtCalving(p_glacier, p_glacier.annual_berg_vol);
		calvingdt_winter = calvingdt_summer;
		double winterdt=0.1; 
	}
	
	if (p_glacier.calving_scenario==0)
	{
		fprintf(stdout,"Seasonal calving rate variation scenario:\ncalving rate 6 times higher in summer than in winter\n");	 										  
		fflush(stdout);
		
		fprintf(out_runtime,"Seasonal calving rate variation scenario:\ncalving rate 6 times higher in summer than in winter\n",calvingdt_summer,calvingdt_winter);	 										  
		fflush(out_runtime);
	}
	else if (p_glacier.calving_scenario==1)
	{
		fprintf(stdout,"Seasonal calving rate variation scenario:\ncalving only in summer\n");	 										  
		fflush(stdout);
		
		fprintf(out_runtime,"Seasonal calving rate variation scenario:\ncalving only in summer\n",calvingdt_summer,calvingdt_winter);	 										  
		fflush(out_runtime);
	}
		if (p_glacier.calving_scenario==2)
	{
		fprintf(stdout,"Seasonal calving rate variation scenario:\nsikussak breakup scenario:\n10 years of icebergs calved in 28 days from 17th August to 13th September\n");	 										  
		fflush(stdout);
		
		fprintf(out_runtime,"Seasonal calving rate variation scenario:\nsikussak breakup scenario:\n10 years of icebergs calved in 28 days from 17th August to 13th September\n",calvingdt_summer,calvingdt_winter);	 										  
		fflush(out_runtime);
	}
		if (p_glacier.calving_scenario==3)
	{
		fprintf(stdout,"Seasonal calving rate variation scenario:\ncalving rate constant all year round\n");	 										  
		fflush(stdout);
		
		fprintf(out_runtime,"Seasonal calving rate variation scenario:\ncalving rate constant all year round\n",calvingdt_summer,calvingdt_winter);	 										  
		fflush(out_runtime);
	}

	
	fprintf(out_runtime,"calvingdt_summer=	%lf days\ncalvingdt_winter=	%lf days\n",calvingdt_summer,calvingdt_winter);	 										  
	fflush(out_runtime);


	fprintf(out_runtime,"P_calving=	%lf\n",P_calving);	 										  
	fflush(out_runtime);


	fprintf(out_runtime,"p_glacier.annual_berg_vol=	%lf	km3\n",p_glacier.annual_berg_vol);	 										  
	fflush(out_runtime);

	
	fprintf(out_runtime,"Wind and water velocity change timestep=	%lf	minutes\n\n",velchange_dt*60.);	 										  
	fflush(out_runtime);

	/* beginning of main TIME loop for model */
	while (time<endTime)
	{	 
		waterVelTidal(p_v_w_tidal, p_water_vel,time*24.);
		v_w_total[0]=p_v_w_tidal[0]+p_v_w[0];
		v_w_total[1]=p_v_w_tidal[1]+p_v_w[1];

		calcMonth(time, &month);
	
		fracyear= modf((time-time_i)/365., &n_years);

		if ((month>=JUN) && (month<=SEP))
		{
			for (loop_temp=0;loop_temp<fjord.no_z_bins;loop_temp++)
				seaTemp[loop_temp]=seaTempSummer[loop_temp]; 
			
			calvingdt=calvingdt_summer;			
			if (p_glacier.calving_scenario!=2)
			{
				hdid=h1;
			}
			else
			{
				h1=calvingdt*24.;
				hdid=h1;
				if ((time<(p_options.spinuptime+229.)) && (time>(p_options.spinuptime+257.))) //calving occurs over 4 weeks for sikussak breakup scenario- calving timestep zero outside this time 
				{
					calvingdt=0.;
				}
			}
			water_vel_res=RES_WATER_VEL_SUMMER;

		}
		else
		{  
			for (loop_temp=0;loop_temp<fjord.no_z_bins;loop_temp++)
			
			seaTemp[loop_temp]=seaTempWinter[loop_temp];
			calvingdt=calvingdt_winter;	
						
			if (p_glacier.calving_scenario!=2)
			{
				hdid=winterdt;
			}
			else
			{
				h1=calvingdt*24.;
				hdid=h1;
				if ((time<(p_options.spinuptime+229.)) && (time>(p_options.spinuptime+257.))) //calving occurs over 4 weeks for sikussak breakup scenario- calving timestep zero outside this time 
				{
					calvingdt=0.;
				}
			}
			water_vel_res=RES_WATER_VEL_WINTER;
		 }

			/*calculate calving probability : this is number of icebergs/day  */
			/* The iceberg calving rate is calculated based on empirical calving relation from Pelto and Warren - 
			it is a probability between 0 and 1 per timestep*/
			if (calvingdt>0.) //no calving if calvingdt<=0
			{
				frac_calving= modf((time-time_i)/calvingdt, &n_calving);
		
				if(((time-time_i)>0.) && ((time-time_i)>=(n_calving*calvingdt)) && ((oldtime-time_i)<(n_calving*calvingdt)))
				{
				
					/*Set random number for whether calving even happens to between 0 and 1		 */
					random_calving = ran8(&seed8);
			
					/*  create new iceberg	  */
					if (random_calving <= P_calving)
					{
						berg_no++; 	
						if(berg_no>=n_bergs)
						{
							printf("The number of bergs exceeds the size of the berg array\n");
							exit(1);
						}	
						calving (time, p_berg, berg_no, p_glacier,p_options, &seed2, &seed, &seed1, &seed3, &seed12, &seed_stayinfjord);	
						tempvol = totvol + p_berg[berg_no].size_i[0]*p_berg[berg_no].size_i[1]*p_berg[berg_no].size_i[2];
						totvol=tempvol;
						total_no_bergs++;
						meansize=meansize+((p_berg[berg_no].size_i[0]/1.3793)+(p_berg[berg_no].size_i[1]/0.8514)+(p_berg[berg_no].size_i[2]/0.8514))/(3.);
					
				
					} //end of if loop for new berg 
				}   /*end of calving if statement: */
			}				  
	
		/*new wind and current calculated every velchange_dt hours */
		fracvelcalc= modf((time-time_i)/((velchange_dt)/24.), &n_velcalc);

		if(((time-time_i)>0.) && ((time-time_i)>=n_velcalc*((velchange_dt)/24.)) && ((oldtime-time_i)<n_velcalc*((velchange_dt)/24.)))
		{
			vel_random[0]=ran6(&seed6);
			vel_random[1]=ran7(&seed7);
			waterVel(p_v_w, &p_v_w_mag,vel_random);
			airVel(p_v_a, &p_v_a_mag, &seed4, &seed5);
		 }
		
		/*loop melting and moving iceberg	*/
		no_bad_bergs=0;
		if (berg_no>=0)
		{			
			tp=time*24.;	/*tp in hours, time in days*/
			
			for (loop_bergs=0;loop_bergs<=berg_no;loop_bergs++)
			{	
				/*check that berg 'moving' - either summer or one of bergs which doesn't stay over winter*/
				if ( ((month>=JUN) && (month<=SEP)) && (time<=((n_years*365.)+272.+p_berg[loop_bergs].freeze_time)) )
				{
					
					ystart[0]=p_berg[loop_bergs].position[0]/1000.;
					ystart[1]=p_berg[loop_bergs].position[1]/1000.;				
					ystart[2]=p_berg[loop_bergs].v_berg[0];
					ystart[3]=p_berg[loop_bergs].v_berg[1];
					
					bergMove(ystart, p_berg, loop_bergs, p_water_vel, hdid);
				
					p_berg[loop_bergs].position[0]=yp[0]*1000.;
					p_berg[loop_bergs].position[1]=yp[1]*1000.;				
					p_berg[loop_bergs].v_berg[0]=yp[2];
					p_berg[loop_bergs].v_berg[1]=yp[3];

					bounceOff(p_berg, loop_bergs, fjord, p_options, p_glacier, p_water_vel,time*24., &seed10, &seed11);

					if (p_berg[loop_bergs].left_fjord==0)
						{
						calcMeanSeaTemp(&avTempWater, seaTemp, p_berg, loop_bergs, p_options, fjord);
						calcMeanSeaSalin(&avSalinWater, seaSalin, p_berg, loop_bergs, p_options, fjord);
						freezingpointwater(avSalinWater,&freezingpt, p_berg[loop_bergs].d_i/2.)	;
						tempgtfreez= avTempWater-freezingpt;
						
						meltAirSensHeat(p_v_a_mag, p_berg[loop_bergs].size_i[0], airTemp[month], &meltSens_a);
						meltAirSolar(solarRad[month], &meltSolar_a);			
							
						// Calculate Melt rate for base, sides and top in m/day 
						 
						diffVect(v_w_total, p_v_w_diff_m, &p_v_w_diff_mag_m, p_berg[loop_bergs].v_berg[0],p_berg[loop_bergs].v_berg[1] );

						meltWaterForcedConvection(p_v_w_diff_mag_m, avTempWater, p_berg[loop_bergs].size_i[0], &meltForcedC_w);
						meltWaterVertConvection(tempgtfreez, &meltVert_w);
						p_berg[loop_bergs].R[TOP]=meltSens_a+meltSolar_a;	 //top
						p_berg[loop_bergs].R[SIDES]= meltForcedC_w+meltVert_w;  //sides (and front and back)

						max_depth_berg = (int)(ceil(p_berg[loop_bergs].d_i/p_options.dz));
				   		meltWaterForcedConvection(p_v_w_diff_mag_m, seaTemp[max_depth_berg], p_berg[loop_bergs].size_i[0], &meltForcedC_w);
					
						p_berg[loop_bergs].R[BASE]=meltForcedC_w;	//base
	 								  
						if (stability(p_berg,loop_bergs)==1)
						{	
							rollOver(p_berg, loop_bergs, p_glacier, &seed9);
							if (time>p_options.spinuptime)
							{  
								depositionQuickWaterAir(p_sediment, fjord, p_glacier, p_options, p_berg,loop_bergs);					
							}	
						}				
						else if (stability2(p_berg,loop_bergs)==1)
						{					
							rollOver2(p_berg, loop_bergs, p_glacier, &seed9);
							if (time>p_options.spinuptime)
							{  
								depositionQuickWaterAir(p_sediment, fjord, p_glacier, p_options, p_berg,loop_bergs);
							}	
						}			
						else
						{	 
							if (time>p_options.spinuptime)
							{
								depositionQuickWater(p_sediment, fjord, p_glacier, p_options, p_berg,loop_bergs);					   	
							}
						}
						meltBerg(p_berg,loop_bergs, p_glacier);	 	
					} 	/*END IF LOOP: IS BERG IN FJORD (ONLY MELT IF IT IS)*/
				}

				//if berg frozen in one place - during winter and at a randomly selected time during the last day of 'summer' (end sept) for bergs frozen in winter
				else
				{
					p_berg[loop_bergs].v_berg[0]=0.;
					p_berg[loop_bergs].v_berg[1]=0.;

					if (p_berg[loop_bergs].left_fjord==0)
					{
						calcMeanSeaTemp(&avTempWater, seaTemp, p_berg, loop_bergs, p_options, fjord);
						calcMeanSeaSalin(&avSalinWater, seaSalin, p_berg, loop_bergs, p_options, fjord);
						freezingpointwater(avSalinWater,&freezingpt, p_berg[loop_bergs].d_i/2.)	;
						tempgtfreez= avTempWater-freezingpt;
						//if the berg has a size greater than 0			
						
						meltAirSensHeat(p_v_a_mag, p_berg[loop_bergs].size_i[0], airTemp[month], &meltSens_a);
						meltAirSolar(solarRad[month], &meltSolar_a);
					
						// Calculate Melt rate for base, sides and top in m/day 
	
						diffVect(v_w_total, p_v_w_diff_m, &p_v_w_diff_mag_m, p_berg[loop_bergs].v_berg[0],p_berg[loop_bergs].v_berg[1] );

						meltWaterForcedConvection(p_v_w_diff_mag_m, avTempWater, p_berg[loop_bergs].size_i[0], &meltForcedC_w);
						meltWaterVertConvection(tempgtfreez, &meltVert_w);
						p_berg[loop_bergs].R[TOP]=meltSens_a+meltSolar_a;	 //top
						p_berg[loop_bergs].R[SIDES]= meltForcedC_w+meltVert_w;  //sides (and front and back)

						max_depth_berg = (int)(ceil(p_berg[loop_bergs].d_i/p_options.dz));
				   		meltWaterForcedConvection(p_v_w_diff_mag_m, seaTemp[max_depth_berg], p_berg[loop_bergs].size_i[0], &meltForcedC_w);
				
						p_berg[loop_bergs].R[BASE]=meltForcedC_w;	//base
	 									 
						if (time>p_options.spinuptime)
						{
							depositionQuickWater(p_sediment, fjord, p_glacier, p_options, p_berg,loop_bergs);	
						}			
						meltBerg(p_berg,loop_bergs, p_glacier);
					}	/*END IF LOOP: IS BERG IN FJORD (ONLY MELT IF IT IS)*/

				}	/*end if berg frozen in fjord on last day of summer (melting only)*/
			
			
				/*If the berg has EXITED THE RECORDED DEPOSITION MODEL AREA*/
				if (((int)(ceil((-p_berg[loop_bergs].position[0]-p_berg[loop_bergs].size_i[0]-SPREAD_X)/p_options.dx))>(int)(ceil(fjord.no_x_bins))) && (p_berg[loop_bergs].left_fjord==0))
				{
					if (time>p_options.spinuptime)
					{
						volume_sed_out_fjord=volume_sed_out_fjord+p_berg[loop_bergs].size_i[0]*p_berg[loop_bergs].size_i[1]*p_berg[loop_bergs].size_i[2]*p_glacier.C_englacial_i;
					}
					p_berg[loop_bergs].left_fjord=1;
				}

					//If the berg has MELTED COMPLETELY
				if (((p_berg[loop_bergs].size_i[0])<=10.) || ((p_berg[loop_bergs].size_i[1])<=10.) ||((p_berg[loop_bergs].size_i[2])<=10.))
				{
					if (berg_no>=0)
					{
						//deleted berg from array p_berg
			
						if (time>p_options.spinuptime)
						{
							volume_sed_out_fjord=volume_sed_out_fjord+p_berg[loop_bergs].size_i[0]*p_berg[loop_bergs].size_i[1]*p_berg[loop_bergs].size_i[2]*p_glacier.C_englacial_i;
						}

						for (del_berg_loop=loop_bergs;del_berg_loop<berg_no;del_berg_loop++)
						{
							p_berg[del_berg_loop]=p_berg[del_berg_loop+1];
						} 
						initialiseIceberg (p_berg, berg_no);
						berg_no = berg_no-1;
						loop_bergs=loop_bergs-1;
						
						no_bad_bergs++;
					}
				}	  /*end of if loop when MELTED COMPLETELY*/	 
				
				/*If the berg has EXITED THE MODEL AREA*/
				else if ((int)(ceil((-p_berg[loop_bergs].position[0]-p_berg[loop_bergs].size_i[0])/p_options.dx))>(int)(ceil(1.25*fjord.no_x_bins)))
				{
					if (berg_no>=0)
					{
					
						if (((time-time_i)-p_berg[loop_bergs].tcalving_i)<230.)
						{
							residence_time_sum=residence_time_sum+((time-time_i)-p_berg[loop_bergs].tcalving_i);
							no_bergs_left_fjord++;
						}
						//delete berg from array p_berg
			
						for (del_berg_loop=loop_bergs;del_berg_loop<berg_no;del_berg_loop++)
						{
							p_berg[del_berg_loop]=p_berg[del_berg_loop+1];
						} 
						initialiseIceberg (p_berg, berg_no);
						berg_no = berg_no-1;
						loop_bergs=loop_bergs-1;
						no_bad_bergs++;
					}
				} //end of 'if' the berg has EXITED THE MODEL AREA		*/

			}	  /*END	FOR LOOP THROUGH ALL BERGS*/
		 }//end of if clause to check there is an iceberg in the system	
		
		 
		 timefrac= modf((time-time_i)/((endTime-time_i)/10.), &n_time);

		/*PRINT FRACTION OF PROGRAM RUN SO FAR*/
		if(((time-time_i)>(n_time*((endTime-time_i)/10.))) && ((oldtime-time_i)<=(n_time*((endTime-time_i)/10.))))
		{
			fprintf(stdout,"fraction of time passed = %0.lf%%, time=	%lf\n",n_time*10.,time);
			fprintf(stdout,"Number of bergs in fjord=%d\n",berg_no+1);
			fprintf(stdout,"total no bergs =	%d\n",total_no_bergs);
			fflush(stdout);

			fprintf(out_runtime,"fraction of time passed = %0.lf%%, time=	%lf\n",n_time*10.,time);
			fprintf(out_runtime,"Number of bergs in fjord=%d\n",berg_no+1);
			fprintf(out_runtime,"total no bergs =	%d\n",total_no_bergs);
			fflush(out_runtime);		
		}
		//time in days (h1 in hours)
		oldtime=time;
		//INCREMENT TIME
		time+=(hdid/24.);

		if ((time>p_options.spinuptime) && (endwarmup==0))
		{

			fprintf(stdout,"No bergs in fjord after warm up=	%d\n", berg_no+1);
			fprintf(out_runtime,"No bergs in fjord after warm up=	%d\n", berg_no+1);
			endwarmup=1;
		}
	
	}//end of main TIME for loop 


	fprintf(stdout,"No bergs in fjord at end of program=	%d\n", berg_no+1);
	fprintf(out_runtime,"No bergs in fjord at end of program=	%d\n", berg_no+1) ;
	
	residence_time=(residence_time_sum/no_bergs_left_fjord)/1.25;
	
	fprintf(stdout,"average residence_time=%lf\n\n",residence_time);
	fprintf(stdout,"tot vol basal ice bergs in %lf years=	%e\ntot vol basal ice glacier in %lf years =	%e\ntot vol eng ice berg in %lf years=	%e\ntot vol eng ice glacier in %lf years=	%e\n",
		p_options.totnyears,totbasal_berg,p_options.totnyears,totbasal_glacier,p_options.totnyears,toteng_berg,p_options.totnyears,toteng_glacier);
	fprintf(stdout,"tot vol sed from glacier in %lf years=	%e\n",p_options.totnyears,(totbasal_glacier+toteng_glacier)*p_glacier.C_englacial_i);

	fprintf(stdout,"Total vol sediment produced in bergs in %lf years=	%e\n",(time-time_i)/365., totvol*p_glacier.C_englacial_i);

	fprintf(stdout,"Volume sediment deposited in fjord in %lf years=	%e\nVolume sediment floated out of fjord in bergs in %lf years=	%e\nTotal volume sediment 'out' in %lf years =	%e\n",
		p_options.totnyears,vol_sed_deposited,p_options.totnyears, volume_sed_out_fjord,p_options.totnyears, volume_sed_out_fjord+vol_sed_deposited);

	fprintf(stdout,"final time = %lf days\n",time);
	fprintf(stdout,"total no bergs = %d\n",total_no_bergs);
	fflush(stdout);

	fprintf(out_runtime,"average residence_time=%lf\n\n",residence_time);
	fprintf(out_runtime,"tot vol basal ice bergs in %lf years=	%e\ntot vol basal ice glacier in %lf years =	%e\ntot vol eng ice berg in %lf years=	%e\ntot vol eng ice glacier in %lf years=	%e\n",
		p_options.totnyears,totbasal_berg,p_options.totnyears,totbasal_glacier,p_options.totnyears,toteng_berg,p_options.totnyears,toteng_glacier);
	fprintf(out_runtime,"tot vol sed from glacier in %lf years=	%e\n",p_options.totnyears,(totbasal_glacier+toteng_glacier)*p_glacier.C_englacial_i);

	fprintf(out_runtime,"Total vol sediment produced in bergs in %lf years=	%e\n",(time-time_i)/365.,totvol*p_glacier.C_englacial_i);

	fprintf(out_runtime,"Volume sediment deposited in fjord in %lf years=	%e\nVolume sediment floated out of fjord in bergs in %lf years=	%e\nTotal volume sediment 'out' in %lf years =	%e\n",
		p_options.totnyears, vol_sed_deposited,p_options.totnyears, volume_sed_out_fjord,p_options.totnyears, volume_sed_out_fjord+vol_sed_deposited);

	fprintf(out_runtime,"final time = %lf days\n",time);
	fprintf(out_runtime,"total no bergs = %d\n",total_no_bergs);
	fflush(out_runtime);

	min_av=(int)floor((fjord.no_y_bins-1)/2.-5.);
	max_av=(int)ceil((fjord.no_y_bins-1)/2.+5.);	 
	width_av= max_av-min_av+1.;

	fprintf(stdout,"min_av=%d\nmax_av=%d\nwidth_av=%lf\n",min_av, max_av, width_av);
	fflush(stdout);

 	fprintf(out_runtime,"min_av=%d\nmax_av=%d\nwidth_av=%lf\n",min_av, max_av, width_av);
	fflush(out_runtime);
	/*OUTPUT FILES */
	  
	out_sed = fopen("Output/sedthickness.txt","w");  
	if(out_sed==NULL) 
	{
	   printf("Error: can't open file sedthickness.txt\n");
    /* DON'T PASS A NULL POINTER TO fclose !! */
		exit(1);
	}	   

	out_sed2 = fopen("Output/sedthickness_small.txt","w");  
	if(out_sed2==NULL) 
	{
	   printf("Error: can't open file sedthickness_small.txt\n");
    /* DON'T PASS A NULL POINTER TO fclose !! */
		exit(1);
	}

	out_sed_avalongfjord = fopen("Output/avalongfjord.txt","w");  
	if(out_sed_avalongfjord==NULL) 
	{
	   printf("Error: can't open file avalongfjord.txt\n");
    /* DON'T PASS A NULL POINTER TO fclose !! */
		exit(1);
	}	   

  	out_sed_avcrossfjord = fopen("Output/avcrossfjord.txt","w");  
	if(out_sed_avcrossfjord==NULL) 
	{
	   printf("Error: can't open file avcrossfjord.txt\n");
    /* DON'T PASS A NULL POINTER TO fclose !! */
		exit(1);
	}

	/* writing the sediment array to file  	*/
	fprintf(stdout,"Writing sediment array... please wait\n");
	fflush(stdout);

	fprintf(out_runtime,"Writing sediment array... please wait\n");
	fflush(out_runtime);
	
	if (p_options.totnyears>0.)
	{
		for (r=0;r<fjord.no_x_bins;r++)
		{
			for (c=0;c<fjord.no_y_bins;c++)
	 		{ 
				tot_crossfjord[c]=tot_crossfjord[c]+p_sediment.thickness[r][c]/(p_options.totnyears*fjord.no_x_bins);			
				fprintf(out_sed, "%lf	",p_sediment.thickness[r][c]/p_options.totnyears);
			}
			for (c=min_av;c<=max_av;c++)
	 		{ 
				tot_alongfjord[r]=tot_alongfjord[r]+p_sediment.thickness[r][c]/(p_options.totnyears*width_av);	
			}

			fprintf(out_sed_avalongfjord, "%lf	%lf\n",r*p_options.dx,tot_alongfjord[r]);
			fprintf(out_sed, "\n");
		}	
		
		fprintf(stdout,"Writing smaller sediment array... please wait\n");
		fprintf(out_runtime,"Writing smaller sediment array... please wait\n");
		
		for (r=0;r<fjord.no_x_bins;r=r+10)
		{
			for (c=0;c<fjord.no_y_bins;c=c+10)
	 		{ 	
				fprintf(out_sed2, "%lf	",p_sediment.thickness[r][c]/p_options.totnyears);
			}
			fprintf(out_sed2, "\n");
		}
			
		for (c=0;c<fjord.no_y_bins;c++)
		{ 		
			fprintf(out_sed_avcrossfjord, "%lf	%lf\n",c*p_options.dx,tot_crossfjord[c]);
			av_dep_rate=av_dep_rate+(tot_crossfjord[c]/fjord.no_y_bins);
		}


	}

	fprintf(stdout,"Average deposition rate over whole fjord	%lf	cm/yr\n\n",av_dep_rate);
	fprintf(out_runtime,"Average deposition rate over whole fjord	%lf	cm/yr\n\n",av_dep_rate);

	fprintf(stdout,"Total volume actually calved per year in %lf years=	%lfkm^3\nTheoretical volume from calving rate per year=	%lfkm^3\n",
		p_options.totnyears,toteng_berg/(p_options.totnyears*1000000000.),theoreticalvol);
	fprintf(stdout,"Total no bergs actually calved per year in %lf years=	%lf\nTheoretical no bergs from calving rate per year=	%lf\n",
		(time-time_i)/365.,(float)(total_no_bergs*365./(time-time_i)), theoreticalnobergs);
	fprintf(stdout,"Summer Calving Six times winter calving\n");	
	fflush(stdout);

	fprintf(out_runtime,"Total volume actually calved per year in %lf years=	%lfkm^3\nTheoretical volume from calving rate per year=	%lfkm^3\n",
		p_options.totnyears,toteng_berg/(p_options.totnyears*1000000000.),theoreticalvol);
	fprintf(out_runtime,"Total no bergs actually calved per year in %lf years=	%lf\nTheoretical no bergs from calving rate per year=	%lf\n",
		(time-time_i)/365.,(float)(total_no_bergs*365./(time-time_i)), theoreticalnobergs);
	fprintf(out_runtime,"Summer Calving Six times winter calving\n");
	fflush(out_runtime);

	meansize= (double) meansize/total_no_bergs;

	fprintf(stdout,"meansize=	%lf\n",meansize);
	fflush(stdout);


	fprintf(out_runtime,"meansize=	%lf\n",meansize);
	fflush(out_runtime);

	fclose(out_sed);
	fclose(out_sed2);
	fclose(out_sed_avalongfjord);
	fclose(out_sed_avcrossfjord);
	fclose(out_runtime);

	free(p_berg);
	
  for (r = 0; r < fjord.no_x_bins; r++) 
  {
    free(p_sediment.thickness[r]);
  }
  /* now free the array of pointers */
  free(p_sediment.thickness);

	free(seaTemp);
	free(seaTempSummer);
	free(seaTempWinter);
	free(seaSalin);			
	free(tot_crossfjord);
	free(tot_alongfjord);

    return 0;
}


#undef PI
#undef g
#undef OMEGA_0
#undef RHO_ICE
#undef TEMP_ICE
#undef RHO_W
#undef RHO_A
#undef RHO_S
#undef C_W
#undef C_A
#undef C_S
#undef C_surf
#undef K_AIR
#undef MU_AIR
#undef KAPPA_AIR
#undef GAMMA_ICE
#undef MAX_BERG_SIZE
#undef PERC_BERG_WINTER
#undef SPREAD_X
#undef SPREAD_Y			
#undef RES_WATER_VEL_SUMMER			
#undef RES_WATER_VEL_WINTER			
#undef	P_calving
#undef BASALMELTWATER
#undef DEBRISMELTAIR
#undef NR_END
#undef FREE_ARG
#undef	IA
#undef	IM
#undef	AM
#undef	IQ
#undef	IR
#undef	NTAB
#undef	NDIV
#undef	EPSRAN
#undef	RNMX