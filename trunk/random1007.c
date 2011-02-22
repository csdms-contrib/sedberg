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
#include "random.h"

 /*"Minimal" random number generator of Park and Miller with Bays-Durham shuffle and added safeguards.  
Returns random deviate between 0.0 and 1.0.  From Numerical Recipes for C*/
float ran(long *itemp)
{
	int j;
	long k;
	static long	iy=0;
	static long iv[NTAB];
	float temp;

	if (*itemp <= 0 || !iy)	  /*Initialise*/
	{
		if (-(*itemp) < 1) *itemp=1;	/*Be sure to prevent itemp=0*/
		else *itemp = -(*itemp);
		for (j=NTAB+7;j>=0;j--)		   /*Load the shuffle table (after 8 warm ups)*/
		{
			k=(*itemp)/IQ;
			*itemp=IA*(*itemp-k*IQ)-IR*k;
			if (*itemp < 0) *itemp += IM;
			if (j < NTAB) iv[j] = *itemp;
		}
		iy=iv[0];
	}
	k=(*itemp)/IQ;					/*Start here when not initialising*/
	*itemp=IA*(*itemp-k*IQ)-IR*k;		 /*Compute itemp=(IA*itemp) % IM withough overflows by Schrage's method*/
	if (*itemp < 0) *itemp += IM;
	j=iy/NDIV;						 /*Will be in the range 0..NTAB-1*/
	iy = iv[j];						/*Output previously stored value and refil the shuffle table*/
	iv[j] = *itemp;
	temp=(float)AM*iy;
	if (temp > RNMX) return (float) RNMX;	   /*Because users don't expect endpoint values*/
//	if ((temp=AM*iy) > RNMX) return RNMX;
	else return temp;

}

float ran1(long *itemp1)
{
	int j;
	long k;
	static long	iy1=0;
	static long iv1[NTAB];
	float temp;

	if (*itemp1 <= 0 || !iy1)	  /*Initialise*/
	{
		if (-(*itemp1) < 1) *itemp1=1;	/*Be sure to prevent itemp1=0*/
		else *itemp1 = -(*itemp1);
		for (j=NTAB+7;j>=0;j--)		   /*Load the shuffle table (after 8 warm ups)*/
		{
			k=(*itemp1)/IQ;
			*itemp1=IA*(*itemp1-k*IQ)-IR*k;
			if (*itemp1 < 0) *itemp1 += IM;
			if (j < NTAB) iv1[j] = *itemp1;
		}
		iy1=iv1[0];
	}
	k=(*itemp1)/IQ;					/*Start here when not initialising*/
	*itemp1=IA*(*itemp1-k*IQ)-IR*k;		 /*Compute itemp1=(IA*itemp1) % IM withough overflows by Schrage's method*/
	if (*itemp1 < 0) *itemp1 += IM;
	j=iy1/NDIV;						 /*Will be in the range 0..NTAB-1*/
	iy1 = iv1[j];						/*Output previously stored value and refil the shuffle table*/
	iv1[j] = *itemp1;
	temp=(float)AM*iy1;
	if (temp > RNMX) return (float) RNMX;	   /*Because users don't expect endpoint values*/
//	if ((temp=AM*iy1) > RNMX) return RNMX;
	else return temp;

}

float ran2(long *itemp2)
{
	int j;
	long k;
	static long	iy2=0;
	static long iv2[NTAB];
	float temp;

	if (*itemp2 <= 0 || !iy2)	  /*Initialise*/
	{
		if (-(*itemp2) < 1) *itemp2=1;	/*Be sure to prevent itemp2=0*/
		else *itemp2 = -(*itemp2);
		for (j=NTAB+7;j>=0;j--)		   /*Load the shuffle table (after 8 warm ups)*/
		{
			k=(*itemp2)/IQ;
			*itemp2=IA*(*itemp2-k*IQ)-IR*k;
			if (*itemp2 < 0) *itemp2 += IM;
			if (j < NTAB) iv2[j] = *itemp2;
		}
		iy2=iv2[0];
	}
	k=(*itemp2)/IQ;					/*Start here when not initialising*/
	*itemp2=IA*(*itemp2-k*IQ)-IR*k;		 /*Compute itemp2=(IA*itemp2) % IM withough overflows by Schrage's method*/
	if (*itemp2 < 0) *itemp2 += IM;
	j=iy2/NDIV;						 /*Will be in the range 0..NTAB-1*/
	iy2 = iv2[j];						/*Output previously stored value and refil the shuffle table*/
	iv2[j] = *itemp2;
	temp=(float)AM*iy2;
	if (temp > RNMX) return (float) RNMX;	   /*Because users don't expect endpoint values*/
//	if ((temp=AM*iy2) > RNMX) return RNMX;
	else return temp;

}

float ran3(long *itemp3)
{
	int j;
	long k;
	static long	iy3=0;
	static long iv3[NTAB];
	float temp;

	if (*itemp3 <= 0 || !iy3)	  /*Initialise*/
	{
		if (-(*itemp3) < 1) *itemp3=1;	/*Be sure to prevent itemp3=0*/
		else *itemp3 = -(*itemp3);
		for (j=NTAB+7;j>=0;j--)		   /*Load the shuffle table (after 8 warm ups)*/
		{
			k=(*itemp3)/IQ;
			*itemp3=IA*(*itemp3-k*IQ)-IR*k;
			if (*itemp3 < 0) *itemp3 += IM;
			if (j < NTAB) iv3[j] = *itemp3;
		}
		iy3=iv3[0];
	}
	k=(*itemp3)/IQ;					/*Start here when not initialising*/
	*itemp3=IA*(*itemp3-k*IQ)-IR*k;		 /*Compute itemp3=(IA*itemp3) % IM withough overflows by Schrage's method*/
	if (*itemp3 < 0) *itemp3 += IM;
	j=iy3/NDIV;						 /*Will be in the range 0..NTAB-1*/
	iy3 = iv3[j];						/*Output previously stored value and refil the shuffle table*/
	iv3[j] = *itemp3;
	temp=(float)AM*iy3;
	if (temp > RNMX) return (float) RNMX;	   /*Because users don't expect endpoint values*/
//	if ((temp=AM*iy3) > RNMX) return RNMX;
	else return temp;

}

float ran4(long *itemp4)
{
	int j;
	long k;
	static long	iy4=0;
	static long iv4[NTAB];
	float temp;

	if (*itemp4 <= 0 || !iy4)	  /*Initialise*/
	{
		if (-(*itemp4) < 1) *itemp4=1;	/*Be sure to prevent itemp4=0*/
		else *itemp4 = -(*itemp4);
		for (j=NTAB+7;j>=0;j--)		   /*Load the shuffle table (after 8 warm ups)*/
		{
			k=(*itemp4)/IQ;
			*itemp4=IA*(*itemp4-k*IQ)-IR*k;
			if (*itemp4 < 0) *itemp4 += IM;
			if (j < NTAB) iv4[j] = *itemp4;
		}
		iy4=iv4[0];
	}
	k=(*itemp4)/IQ;					/*Start here when not initialising*/
	*itemp4=IA*(*itemp4-k*IQ)-IR*k;		 /*Compute itemp4=(IA*itemp4) % IM withough overflows by Schrage's method*/
	if (*itemp4 < 0) *itemp4 += IM;
	j=iy4/NDIV;						 /*Will be in the range 0..NTAB-1*/
	iy4 = iv4[j];						/*Output previously stored value and refil the shuffle table*/
	iv4[j] = *itemp4;
	temp=(float)AM*iy4;
	if (temp > RNMX) return (float) RNMX;	   /*Because users don't expect endpoint values*/
//	if ((temp=AM*iy4) > RNMX) return RNMX;
	else return temp;

}

float ran5(long *itemp5)
{
	int j;
	long k;
	static long	iy5=0;
	static long iv5[NTAB];
	float temp;

	if (*itemp5 <= 0 || !iy5)	  /*Initialise*/
	{
		if (-(*itemp5) < 1) *itemp5=1;	/*Be sure to prevent itemp5=0*/
		else *itemp5 = -(*itemp5);
		for (j=NTAB+7;j>=0;j--)		   /*Load the shuffle table (after 8 warm ups)*/
		{
			k=(*itemp5)/IQ;
			*itemp5=IA*(*itemp5-k*IQ)-IR*k;
			if (*itemp5 < 0) *itemp5 += IM;
			if (j < NTAB) iv5[j] = *itemp5;
		}
		iy5=iv5[0];
	}
	k=(*itemp5)/IQ;					/*Start here when not initialising*/
	*itemp5=IA*(*itemp5-k*IQ)-IR*k;		 /*Compute itemp5=(IA*itemp5) % IM withough overflows by Schrage's method*/
	if (*itemp5 < 0) *itemp5 += IM;
	j=iy5/NDIV;						 /*Will be in the range 0..NTAB-1*/
	iy5 = iv5[j];						/*Output previously stored value and refil the shuffle table*/
	iv5[j] = *itemp5;
	temp=(float)AM*iy5;
	if (temp > RNMX) return (float) RNMX;	   /*Because users don't expect endpoint values*/
//	if ((temp=AM*iy5) > RNMX) return RNMX;
	else return temp;

}

float ran6(long *itemp6)
{
	int j;
	long k;
	static long	iy6=0;
	static long iv6[NTAB];
	float temp;

	if (*itemp6 <= 0 || !iy6)	  /*Initialise*/
	{
		if (-(*itemp6) < 1) *itemp6=1;	/*Be sure to prevent itemp6=0*/
		else *itemp6 = -(*itemp6);
		for (j=NTAB+7;j>=0;j--)		   /*Load the shuffle table (after 8 warm ups)*/
		{
			k=(*itemp6)/IQ;
			*itemp6=IA*(*itemp6-k*IQ)-IR*k;
			if (*itemp6 < 0) *itemp6 += IM;
			if (j < NTAB) iv6[j] = *itemp6;
		}
		iy6=iv6[0];
	}
	k=(*itemp6)/IQ;					/*Start here when not initialising*/
	*itemp6=IA*(*itemp6-k*IQ)-IR*k;		 /*Compute itemp6=(IA*itemp6) % IM withough overflows by Schrage's method*/
	if (*itemp6 < 0) *itemp6 += IM;
	j=iy6/NDIV;						 /*Will be in the range 0..NTAB-1*/
	iy6 = iv6[j];						/*Output previously stored value and refil the shuffle table*/
	iv6[j] = *itemp6;
	temp=(float)AM*iy6;
	if (temp > RNMX) return (float) RNMX;	   /*Because users don't expect endpoint values*/
//	if ((temp=AM*iy6) > RNMX) return RNMX;
	else return temp;

}

float ran7(long *itemp7)
{
	int j;
	long k;
	static long	iy7=0;
	static long iv7[NTAB];
	float temp;

	if (*itemp7 <= 0 || !iy7)	  /*Initialise*/
	{
		if (-(*itemp7) < 1) *itemp7=1;	/*Be sure to prevent itemp7=0*/
		else *itemp7 = -(*itemp7);
		for (j=NTAB+7;j>=0;j--)		   /*Load the shuffle table (after 8 warm ups)*/
		{
			k=(*itemp7)/IQ;
			*itemp7=IA*(*itemp7-k*IQ)-IR*k;
			if (*itemp7 < 0) *itemp7 += IM;
			if (j < NTAB) iv7[j] = *itemp7;
		}
		iy7=iv7[0];
	}
	k=(*itemp7)/IQ;					/*Start here when not initialising*/
	*itemp7=IA*(*itemp7-k*IQ)-IR*k;		 /*Compute itemp7=(IA*itemp7) % IM withough overflows by Schrage's method*/
	if (*itemp7 < 0) *itemp7 += IM;
	j=iy7/NDIV;						 /*Will be in the range 0..NTAB-1*/
	iy7 = iv7[j];						/*Output previously stored value and refil the shuffle table*/
	iv7[j] = *itemp7;
	temp=(float)AM*iy7;
	if (temp > RNMX) return (float) RNMX;	   /*Because users don't expect endpoint values*/
//	if ((temp=AM*iy7) > RNMX) return RNMX;
	else return temp;

}

float ran8(long *itemp8)
{
	int j;
	long k;
	static long	iy8=0;
	static long iv8[NTAB];
	float temp;

	if (*itemp8 <= 0 || !iy8)	  /*Initialise*/
	{
		if (-(*itemp8) < 1) *itemp8=1;	/*Be sure to prevent itemp8=0*/
		else *itemp8 = -(*itemp8);
		for (j=NTAB+7;j>=0;j--)		   /*Load the shuffle table (after 8 warm ups)*/
		{
			k=(*itemp8)/IQ;
			*itemp8=IA*(*itemp8-k*IQ)-IR*k;
			if (*itemp8 < 0) *itemp8 += IM;
			if (j < NTAB) iv8[j] = *itemp8;
		}
		iy8=iv8[0];
	}
	k=(*itemp8)/IQ;					/*Start here when not initialising*/
	*itemp8=IA*(*itemp8-k*IQ)-IR*k;		 /*Compute itemp8=(IA*itemp8) % IM withough overflows by Schrage's method*/
	if (*itemp8 < 0) *itemp8 += IM;
	j=iy8/NDIV;						 /*Will be in the range 0..NTAB-1*/
	iy8 = iv8[j];						/*Output previously stored value and refil the shuffle table*/
	iv8[j] = *itemp8;
	temp=(float)AM*iy8;
	if (temp > RNMX) return (float) RNMX;	   /*Because users don't expect endpoint values*/
//	if ((temp=AM*iy8) > RNMX) return RNMX;
	else return temp;

}
float ran9(long *itemp9)
{
	int j;
	long k;
	static long	iy9=0;
	static long iv9[NTAB];
	float temp;

	if (*itemp9 <= 0 || !iy9)	  /*Initialise*/
	{
		if (-(*itemp9) < 1) *itemp9=1;	/*Be sure to prevent itemp9=0*/
		else *itemp9 = -(*itemp9);
		for (j=NTAB+7;j>=0;j--)		   /*Load the shuffle table (after 8 warm ups)*/
		{
			k=(*itemp9)/IQ;
			*itemp9=IA*(*itemp9-k*IQ)-IR*k;
			if (*itemp9 < 0) *itemp9 += IM;
			if (j < NTAB) iv9[j] = *itemp9;
		}
		iy9=iv9[0];
	}
	k=(*itemp9)/IQ;					/*Start here when not initialising*/
	*itemp9=IA*(*itemp9-k*IQ)-IR*k;		 /*Compute itemp9=(IA*itemp9) % IM withough overflows by Schrage's method*/
	if (*itemp9 < 0) *itemp9 += IM;
	j=iy9/NDIV;						 /*Will be in the range 0..NTAB-1*/
	iy9 = iv9[j];						/*Output previously stored value and refil the shuffle table*/
	iv9[j] = *itemp9;
	temp=(float)AM*iy9;
	if (temp > RNMX) return (float) RNMX;	   /*Because users don't expect endpoint values*/
//	if ((temp=AM*iy9) > RNMX) return RNMX;
	else return temp;

}
float ran10(long *itemp10)
{
	int j;
	long k;
	static long	iy10=0;
	static long iv10[NTAB];
	float temp;

	if (*itemp10 <= 0 || !iy10)	  /*Initialise*/
	{
		if (-(*itemp10) < 1) *itemp10=1;	/*Be sure to prevent itemp9=0*/
		else *itemp10 = -(*itemp10);
		for (j=NTAB+7;j>=0;j--)		   /*Load the shuffle table (after 8 warm ups)*/
		{
			k=(*itemp10)/IQ;
			*itemp10=IA*(*itemp10-k*IQ)-IR*k;
			if (*itemp10 < 0) *itemp10 += IM;
			if (j < NTAB) iv10[j] = *itemp10;
		}
		iy10=iv10[0];
	}
	k=(*itemp10)/IQ;					/*Start here when not initialising*/
	*itemp10=IA*(*itemp10-k*IQ)-IR*k;		 /*Compute itemp9=(IA*itemp9) % IM withough overflows by Schrage's method*/
	if (*itemp10 < 0) *itemp10 += IM;
	j=iy10/NDIV;						 /*Will be in the range 0..NTAB-1*/
	iy10 = iv10[j];						/*Output previously stored value and refil the shuffle table*/
	iv10[j] = *itemp10;
	temp=(float)AM*iy10;
	if (temp > RNMX) return (float) RNMX;	   /*Because users don't expect endpoint values*/
//	if ((temp=AM*iy9) > RNMX) return RNMX;
	else return temp;

}

float ran11(long *itemp11)
{
	int j;
	long k;
	static long	iy11=0;
	static long iv11[NTAB];
	float temp;

	if (*itemp11 <= 0 || !iy11)	  /*Initialise*/
	{
		if (-(*itemp11) < 1) *itemp11=1;	/*Be sure to prevent itemp9=0*/
		else *itemp11 = -(*itemp11);
		for (j=NTAB+7;j>=0;j--)		   /*Load the shuffle table (after 8 warm ups)*/
		{
			k=(*itemp11)/IQ;
			*itemp11=IA*(*itemp11-k*IQ)-IR*k;
			if (*itemp11 < 0) *itemp11 += IM;
			if (j < NTAB) iv11[j] = *itemp11;
		}
		iy11=iv11[0];
	}
	k=(*itemp11)/IQ;					/*Start here when not initialising*/
	*itemp11=IA*(*itemp11-k*IQ)-IR*k;		 /*Compute itemp9=(IA*itemp9) % IM withough overflows by Schrage's method*/
	if (*itemp11 < 0) *itemp11 += IM;
	j=iy11/NDIV;						 /*Will be in the range 0..NTAB-1*/
	iy11 = iv11[j];						/*Output previously stored value and refil the shuffle table*/
	iv11[j] = *itemp11;
	temp=(float)AM*iy11;
	if (temp > RNMX) return (float) RNMX;	   /*Because users don't expect endpoint values*/
//	if ((temp=AM*iy9) > RNMX) return RNMX;
	else return temp;

}

float ran12(long *itemp12)
{
	int j;
	long k;
	static long	iy12=0;
	static long iv12[NTAB];
	float temp;

	if (*itemp12 <= 0 || !iy12)	  /*Initialise*/
	{
		if (-(*itemp12) < 1) *itemp12=1;	/*Be sure to prevent itemp9=0*/
		else *itemp12 = -(*itemp12);
		for (j=NTAB+7;j>=0;j--)		   /*Load the shuffle table (after 8 warm ups)*/
		{
			k=(*itemp12)/IQ;
			*itemp12=IA*(*itemp12-k*IQ)-IR*k;
			if (*itemp12 < 0) *itemp12 += IM;
			if (j < NTAB) iv12[j] = *itemp12;
		}
		iy12=iv12[0];
	}
	k=(*itemp12)/IQ;					/*Start here when not initialising*/
	*itemp12=IA*(*itemp12-k*IQ)-IR*k;		 /*Compute itemp9=(IA*itemp9) % IM withough overflows by Schrage's method*/
	if (*itemp12 < 0) *itemp12 += IM;
	j=iy12/NDIV;						 /*Will be in the range 0..NTAB-1*/
	iy12 = iv12[j];						/*Output previously stored value and refil the shuffle table*/
	iv12[j] = *itemp12;
	temp=(float)AM*iy12;
	if (temp > RNMX) return (float) RNMX;	   /*Because users don't expect endpoint values*/
//	if ((temp=AM*iy9) > RNMX) return RNMX;
	else return temp;

}


float f_ran_stayinfjord(long *itemp_stayinfjord)
{
	int j;
	long k;
	static long	iy_stayinfjord=0;
	static long iv_stayinfjord[NTAB];
	float temp;

	if (*itemp_stayinfjord <= 0 || !iy_stayinfjord)	  /*Initialise*/
	{
		if (-(*itemp_stayinfjord) < 1) *itemp_stayinfjord=1;	/*Be sure to prevent itemp9=0*/
		else *itemp_stayinfjord = -(*itemp_stayinfjord);
		for (j=NTAB+7;j>=0;j--)		   /*Load the shuffle table (after 8 warm ups)*/
		{
			k=(*itemp_stayinfjord)/IQ;
			*itemp_stayinfjord=IA*(*itemp_stayinfjord-k*IQ)-IR*k;
			if (*itemp_stayinfjord < 0) *itemp_stayinfjord += IM;
			if (j < NTAB) iv_stayinfjord[j] = *itemp_stayinfjord;
		}
		iy_stayinfjord=iv_stayinfjord[0];
	}
	k=(*itemp_stayinfjord)/IQ;					/*Start here when not initialising*/
	*itemp_stayinfjord=IA*(*itemp_stayinfjord-k*IQ)-IR*k;		 /*Compute itemp9=(IA*itemp9) % IM withough overflows by Schrage's method*/
	if (*itemp_stayinfjord < 0) *itemp_stayinfjord += IM;
	j=iy_stayinfjord/NDIV;						 /*Will be in the range 0..NTAB-1*/
	iy_stayinfjord = iv_stayinfjord[j];						/*Output previously stored value and refil the shuffle table*/
	iv_stayinfjord[j] = *itemp_stayinfjord;
	temp=(float)AM*iy_stayinfjord;
	if (temp > RNMX) return (float) RNMX;	   /*Because users don't expect endpoint values*/
//	if ((temp=AM*iy9) > RNMX) return RNMX;
	else return temp;

}







/*Returns a normally distributed deviate with zero mean and unit variance using ran2(itemp) as the source of uniform deviates*/
float gausdev2(long *itemp2)
{
	float ran2(long *itemp2);
	static int iset = 0;
	static float gset2;									   
	float fac,rsq,v1,v2;

	if (*itemp2 < 0) iset=0;	 /*Reinitialise*/
	if (iset == 0)				  /*We don't have an extra deviate handy, so*/
	{
		do 
		{
			v1=(float)(2.0*ran2(itemp2)-1.0);				/*pick two uniform numbers in the square*/
			v2=(float)(2.0*ran2(itemp2)-1.0);				/* extending from -1 to +1 in each direction, */
		//	v1=2.0*ran1(itemp2)-1.0;				/*pick two uniform numbers in the square*/
		//	v2=2.0*ran1(itemp2)-1.0;	
			rsq=v1*v1+v2*v2;					/*see if they are in the unit circle*/
		} while (rsq >= 1.0 || rsq == 0.0);		/*and if they are not, try again*/
		fac=(float) sqrt(-2.0*log(rsq)/rsq);				/*Now make the Box-Muller transformation to get two normal deviates.  */
	//	fac=sqrt(-2.0*log(rsq)/rsq);
		gset2=v1*fac;							/*Return one and save the other for next time*/
		iset=1;									/*Set flag*/
		return ((float) v2*fac);
	} 
	else 						 
	{							 /*We have an extra deviate handy*/ 
		iset=0;					 /* so unset the flag	*/
		return ((float) gset2);			 /*  and return it. */

	}
}
/*Returns a normally distributed deviate with zero mean and unit variance using ran1(itemp2) as the source of uniform deviates*/
float gausdev4(long *itemp4)
{
	float ran4(long *itemp4);
	static int iset = 0;
	static float gset4;									   
	float fac,rsq,v1,v2;

	if (*itemp4 < 0) iset=0;	 /*Reinitialise*/
	if (iset == 0)				  /*We don't have an extra deviate handy, so*/
	{
		do 
		{
			v1=(float)(2.0*ran4(itemp4)-1.0);				/*pick two uniform numbers in the square*/
			v2=(float)(2.0*ran4(itemp4)-1.0);				/* extending from -1 to +1 in each direction, */
		//	v1=2.0*ran1(itemp4)-1.0;				/*pick two uniform numbers in the square*/
		//	v2=2.0*ran1(itemp4)-1.0;	
			rsq=v1*v1+v2*v2;					/*see if they are in the unit circle*/
		} while (rsq >= 1.0 || rsq == 0.0);		/*and if they are not, try again*/
		fac=(float) sqrt(-2.0*log(rsq)/rsq);				/*Now make the Box-Muller transformation to get two normal deviates.  */
	//	fac=sqrt(-2.0*log(rsq)/rsq);
		gset4=v1*fac;							/*Return one and save the other for next time*/
		iset=1;									/*Set flag*/
		return ((float) v2*fac);
	} 
	else 						 
	{							 /*We have an extra deviate handy*/ 
		iset=0;					 /* so unset the flag	*/
		return ((float) gset4);			 /*  and return it. */

	}
}

/*Returns a normally distributed deviate with zero mean and unit variance using ran1(itemp4) as the source of uniform deviates*/
float gausdev5(long *itemp5)
{
	float ran5(long *itemp5);
	static int iset = 0;
	static float gset5;									   
	float fac,rsq,v1,v2;

	if (*itemp5 < 0) iset=0;	 /*Reinitialise*/
	if (iset == 0)				  /*We don't have an extra deviate handy, so*/
	{
		do 
		{
			v1=(float)(2.0*ran5(itemp5)-1.0);				/*pick two uniform numbers in the square*/
			v2=(float)(2.0*ran5(itemp5)-1.0);				/* extending from -1 to +1 in each direction, */
		//	v1=2.0*ran1(itemp5)-1.0;				/*pick two uniform numbers in the square*/
		//	v2=2.0*ran1(itemp5)-1.0;	
			rsq=v1*v1+v2*v2;					/*see if they are in the unit circle*/
		} while (rsq >= 1.0 || rsq == 0.0);		/*and if they are not, try again*/
		fac=(float) sqrt(-2.0*log(rsq)/rsq);				/*Now make the Box-Muller transformation to get two normal deviates.  */
	//	fac=sqrt(-2.0*log(rsq)/rsq);
		gset5=v1*fac;							/*Return one and save the other for next time*/
		iset=1;									/*Set flag*/
		return ((float) v2*fac);
	} 
	else 						 
	{							 /*We have an extra deviate handy*/ 
		iset=0;					 /* so unset the flag	*/
		return ((float) gset5);			 /*  and return it. */

	}
}



#undef IA
#undef IM
#undef AM
#undef IQ
#undef IR
#undef NTAB
#undef NDIV
#undef EPSRAN
#undef RNMX

