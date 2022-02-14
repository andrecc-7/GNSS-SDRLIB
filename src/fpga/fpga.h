#ifndef FPGA_H
#define FPGA_H

#include "xaxidma.h"
#include "xparameters.h"
#include "xil_printf.h"
#include "xtime_l.h"

/*
 * short int : 2 bytes
 * int : 4 bytes
 * long int : 4 bytes
 * long long int : 8 bytes
 * */

//#define CHECK_STEPS

/**********************************/
/*         ALLOCATIONS            */
/**********************************/

/* SIZES */
#define NSAMPLES	16368//524288//131072//10000000//16368 	/* SAMPLES IN 1 CODE PERIOD */
#define NCHIPS		1023	                                /* CHIPS   "              " */
#define SAMP_CHIP 	16
#define CORRN 1
#define CORRP (2*CORRN + 1)
#define CODE_COPIES 3

/* ADDRESSES */
#define ALIGN(X) (((X) + sizeof(int)) & 0xFFFFFFE0)

#define BASE_ADDR_0 0x100000
// 255 MB
#define BASE_ADDR_1	0x10000000 // IN FILE 1 : SIGNAL									/* 16 MB */
#define BASE_ADDR_2	0x11000000 // IN FILE 2 : CODE										/* 16 MB */
// to check results
#define CORRECT_II	0x12000000 // IN FILE 3 : In Phase									/* 16 MB */
#define CORRECT_QQ	0x13000000 // IN FILE 4 : Quadratue									/* 16 MB */
#define CORRECT_CE	0x14000000 // IN FILE 5 : Early										/* 16 MB */
#define CORRECT_CP	0x15000000 // IN FILE 6 : Prompt									/* 16 MB */
#define CORRECT_CL	0x16000000 // IN FILE 7 : Late										/* 16 MB */

#define GPSGIOVEA01 0x1A000000 // IN FILE 8 : GPS_and_GIOVE_A-NN-fs16_3676-if4_1304.bin	/* 16 MB */
#define GPSGIOVEA02 0x1B000000 // IN FILE 9 : GPS_fs16_3676E6-if4_1304E6_8bit_Real.bin	/* 16 MB */

#define BASE_ADDR_3	0x17000000 // ALLOCATIONS											/* 32 MB */
#define BASE_ADDR_4 0x19100000 // AUX ALLOCATIONS										/*  1 MB */
#define BASE_ADDR_5 0x19200000 // RESULT												/*  1 MB */

#define CC 8

#define DATA   		BASE_ADDR_1//GPSGIOVEA02 //BASE_ADDR_1
#define CODE   		BASE_ADDR_2

#define DATAI  		BASE_ADDR_3
#define DATAQ    	ALIGN(DATAI   +sizeof(short)*NSAMPLES*CC)
#define C_PROMPT 	ALIGN(DATAQ   +sizeof(short)*NSAMPLES*CC)
#define C_EARLY  	ALIGN(C_PROMPT+sizeof(short)*(SAMP_CHIP+NSAMPLES)*CC)
#define C_LATE   	ALIGN(C_EARLY +sizeof(short)*(SAMP_CHIP+NSAMPLES)*CC)
#define CORRI	 	ALIGN(C_LATE  +sizeof(short)*(SAMP_CHIP+NSAMPLES)*CC)
#define CORRQ	 	ALIGN(CORRI   +sizeof(short)*CORRP*CC)

/**********************************/
/*           MACROS               */
/**********************************/

#define SHORT(X) ((short) ( !(X & 0x80) ?  X : X | 0xFF00 ))
#define TRY(error_code, str) {if(error_code) {xil_printf(str, (error_code)); exit(-1);}}
#define TRYALLOC(pt, str) {if((pt) == NULL) {xil_printf(str); exit(-1);}}

/**********************************/
/*         DECLARATIONS           */
/**********************************/

void transaction_tracking(int n, char* data, int dtype_, double phase_, double phase_step_,
		short int* code, int clen, double coff_, double crate_, double csperiod_, int smax, int ns,
		double* CI, double* CQ);

/* VECTORS */
extern volatile char *data;
extern volatile short *code;
extern volatile short *dataI;
extern volatile short *dataQ;
extern volatile short *code_e;
extern volatile short *code_p;
extern volatile short *code_l;
extern volatile double *corrI;
extern volatile double *corrQ;

#endif
