
/*--------------------------------------------- 
	First Attempt to :								
 	RADIX-2 FFT, DECIMATION IN TIME 	
	v 0.1	Emant, 2014					
---------------------------------------------*/

#include <stdlib.h>
#include <math.h>

#include "fft_c.h"

unsigned int reverse_bit(unsigned int c, unsigned int depth)
{
	unsigned int i=0,o=0;
	o = o | (c & 0x01);
	for(i = 0; i < depth - 1; i++)
	{
		o <<= 1 ; c >>= 1;
		o = o | (c & 0x01);	
	}
	return o;
}


int fft_c(	float *in_r  ,			/* Float Input  Array (Real Part)*/  
			float *in_i  ,			/* Float Input  Array (Imag Part)*/  
			float *out_r , 			/* Float Output Array (Real Part)*/ 
			float *out_i , 			/* Float Output Array (Imag Part)*/ 
			unsigned char level )	/* N Bins = 2^level 			 */
{
	float *twid_r, *twid_i, tw;
	unsigned int i,j,k,l,m,tw_max,step, bins;
	
	tw_max = 1<<(level - 1);
	bins = tw_max*2;
	
	/* Allocating Memory For Tweedles. 
	I'm Using Optimized Butterfly Structure 
	(See "R.G. Lyons Understanding DSP"  §4.6 )
	Only half factors needed ! */
	
	twid_r = (float*)calloc(bins,sizeof(float));
	twid_i = twid_r + tw_max;
	
	
	tw = (2.0*M_PI)/bins;
	
	/*Butterfly Level 1 , Reverse Bit Index, Calc. Tweedles*/
	
	for (i = 0; i < tw_max; i++)
	{
		unsigned int k, i1 , i2;
		k = 2*i;
		
		/* Pre - Calculating Twiddles*/
		twid_r[i] = cos(tw*i);
		twid_i[i] = -sin(tw*i);
		
		/* Bit Inversion for Input Array Index */
		i1 = reverse_bit(k,level);
		i2 = reverse_bit(k+1,level);
		
		/* Butterfly Level 1 (Twiddle Factors are +1, -1) */
		out_r[k]   = in_r[i1] + in_r[i2];
		out_r[k+1] = in_r[i1] - in_r[i2];
		out_i[k]   = in_i[i1] + in_i[i2];
		out_i[k+1] = in_i[i1] - in_i[i2];
		#ifdef DEBUG
				printf("tw[%i] : (%f,%f)\ni1=%i k=%i\ni2=%i k+1=%i\n",i,twid_r[i],twid_i[i],i1,k,i2,k+1);
				printf("-- out[%i] = (%f,%f) , out[%i] = (%f,%f) \n",k,out_r[k],
							out_i[k],k+1,out_r[k+1],out_i[k+1]);
		#endif
	}
	
	/* Remaining Butterfy Layers... */
	step = tw_max;						/* Step for traveling Tweedles Array */
	l = 2;								/* Elements to be processed in a Butterfly Level */ 
	#ifdef DEBUG
		printf("tw_max=%i\nbins=%i\n\n",tw_max,bins);
	#endif
	for(i=1; i < level ; i++)			/* Outer Loop, Thru Layers of Butterflies*/
	{
		float a_r, a_i , b_r ,b_i;
		step >>= 1; 					/* Divide Tweedles Step by 2*/			
		#ifdef DEBUG
			printf("i=%i step = %i\n",i,step);
		#endif
		for(j=0; j< bins ; j += l*2 )	/* Center Loop, Thru Butterflies in a Level */	
		{
			#ifdef DEBUG
				printf("	j=%i\n",j);
			#endif
			m = 0; 						/* Twiddle Index = m += step while < bins/2*/
			for(k=j ; k < j+l ; k++)	/* Inner Loop, Thru Element in a Butterfly */
			{
				/* FFT CORE: Optimized Structure - Only One Complex Mult. is needed ! */
				a_r = out_r[k+l]*twid_r[m] - out_i[k+l]*twid_i[m];
				a_i = out_r[k+l]*twid_i[m] + out_i[k+l]*twid_r[m];
				b_r = out_r[k];
				b_i = out_i[k];
				out_r[k] += a_r;
				out_i[k] += a_i;
				out_r[k+l] = b_r - a_r;
				out_i[k+l] = b_i - a_i;
				#ifdef DEBUG
					printf("		k   = %i\n",k);
					printf("		k+l = %i\n",k+l);
					printf("-- out[%i] = (%f,%f) , out[%i] = (%f,%f) \n",k,out_r[k],
						out_i[k],k+l,out_r[k+l],out_i[k+l]);
				#endif
				m += step;				/* Increment Tweedle Index */
			}		
		}
		l <<=1;						/* Double Number of Elements in a Butterfly Layer */
	}
	return 0;	
}

int ifft_c (float *in_r  ,			/* Float Input  Array (Real Part)*/  
			float *in_i  ,			/* Float Input  Array (Imag Part)*/  
			float *out_r , 			/* Float Output Array (Real Part)*/ 
			float *out_i , 			/* Float Output Array (Imag Part)*/ 
			unsigned char level )	/* N Bins = 2^level 			 */
{
	int i;
	int k = 1<<level;
	for(i = 0; i < k; ++i)
	{
		in_i[i] = -in_i[i]/k;
		in_r[i] /= k;
	}		
	fft_c(in_r,in_i,out_r,out_i,level);
	return 0;
}

