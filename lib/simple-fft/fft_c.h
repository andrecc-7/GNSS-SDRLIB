
/*--------------------------------------------- 
	First Attempt to :								
 	RADIX-2 FFT, DECIMATION IN TIME 	
	v 0.1	Emant, 2014					
---------------------------------------------*/


int fft_c(	float *in_r  ,			/* Float Input  Array (Real Part)*/  
			float *in_i  ,			/* Float Input  Array (Imag Part)*/  
			float *out_r , 			/* Float Output Array (Real Part)*/ 
			float *out_i , 			/* Float Output Array (Imag Part)*/ 
			unsigned char level );	/* N Bins = 2^level 			 */
			
int ifft_c (float *in_r  ,			/* Float Input  Array (Real Part)*/  
			float *in_i  ,			/* Float Input  Array (Imag Part)*/  
			float *out_r , 			/* Float Output Array (Real Part)*/ 
			float *out_i , 			/* Float Output Array (Imag Part)*/ 
			unsigned char level );	/* N Bins = 2^level 			 */
