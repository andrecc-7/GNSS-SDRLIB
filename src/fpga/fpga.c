#include "fpga.h"

/**********************************/
/*         ALLOCATIONS            */
/**********************************/

/* Define vectors */
volatile char *data = DATA;
volatile short *code = CODE;
volatile short *dataI = DATAI;
volatile short *dataQ = DATAQ;
volatile short *code_e = C_EARLY;
volatile short *code_p = C_PROMPT;
volatile short *code_l = C_LATE;
								   //  EARLY  PROMPT  LATE
volatile double *corrI = CORRI;    //  {0,    0,      0};
volatile double *corrQ = CORRQ;	   //  {0,    0,      0};

volatile int* aux = BASE_ADDR_4;
volatile int* res = BASE_ADDR_5;

/**********************************/
/*         DECLARATIONS           */
/**********************************/

void transaction_rescode(XAxiDma *dma, short int* code, int clen, int n, double coff_, double csperiod_, int smax, int ns);
void transaction_mixcarr(XAxiDma *dma, char* data, int n, int dtype_, double phase_, double phase_step_, int ns);

/**********************************/
/*        TRANSACTIONS            */
/**********************************/

/* Initialize DMA */
int init_dma(XAxiDma *dma, u32 Device)
{
	XAxiDma_Config *cfg;

	TRYALLOC(
		cfg = XAxiDma_LookupConfig(Device),
		"")
	TRY(
		XAxiDma_CfgInitialize(dma, cfg),
		"ERROR: DMA: Init")

    XAxiDma_IntrDisable(dma, XAXIDMA_IRQ_ALL_MASK, XAXIDMA_DEVICE_TO_DMA);
    XAxiDma_IntrDisable(dma, XAXIDMA_IRQ_ALL_MASK, XAXIDMA_DMA_TO_DEVICE);

    return 0;
}

#define O_TYP_B 8
void transaction_tracking(int n, char* data, int dtype_, double phase_, double phase_step_,
		short int* code, int clen, double coff_, double crate_, double csperiod_, int smax, int ns,
		double* CI, double* CQ)
{
	int status, size;

	XAxiDma dma0, dma1;
	XAxiDma_Config *cfg0, *cfg1;

	init_dma(&dma0, XPAR_AXI_DMA_0_DEVICE_ID);
	init_dma(&dma1, XPAR_AXI_DMA_1_DEVICE_ID);

Xil_DCacheDisable(); // --------------------------------------------------

    // ** OUTPUT **

	xil_printf("Requesting %d results from each correlator\n", 2*ns+1);

    TRY(
    	XAxiDma_SimpleTransfer(&dma0, CI, (2*ns+1)*O_TYP_B, XAXIDMA_DEVICE_TO_DMA),
		"DMA 0: RX: error: %d\n")
	TRY(
		XAxiDma_SimpleTransfer(&dma1, CQ, (2*ns+1)*O_TYP_B, XAXIDMA_DEVICE_TO_DMA),
		"DMA 1: RX: error: %d\n")

	// ** INPUT **

	xil_printf("TRANSACTION: RESAMPLE CODE\n");
	transaction_rescode(&dma1, code, clen, n, coff_ - csperiod_*smax, csperiod_, smax, ns);
	xil_printf("TRANSACTION: MIX CARRIER\n");
    transaction_mixcarr(&dma0, data, n, dtype_, phase_, phase_step_, ns);

	// ** RESULTS **

	xil_printf("Wait for resuslts...\n");
	while(XAxiDma_Busy(&dma0, XAXIDMA_DEVICE_TO_DMA) == TRUE);
	while(XAxiDma_Busy(&dma1, XAXIDMA_DEVICE_TO_DMA) == TRUE);

Xil_DCacheEnable(); // --------------------------------------------------

}

#define RS_I_BUS_b 32
#define RS_I_BUS_B 4
#define RS_O_BUS_b 64
#define RS_O_BUS_B 8
#define RS_I_TYP_B 2
#define RS_O_TYP_B 2

void transaction_rescode(XAxiDma *dma, short int* code, int clen, int n, double coff_, double csperiod_, int smax, int ns)
{
	int status, size;

	// ** ARGUMENTS **

	int last = RS_O_TYP_B*(n/RS_O_BUS_B) - 1;
	long long int coff = (long long int) (coff_ * pow(2,48)); // Q16.48
	long long int csperiod = (long long int) (csperiod_ * pow(2,48)); // Q16.48
	short int codes = CODE_COPIES;

	aux[0] = last;
	aux[1] = ((int*) &coff)[0];
	aux[2] = ((int*) &coff)[1];
	aux[3] = ((int*) &csperiod)[0];
	aux[4] = ((int*) &csperiod)[1];
	((char*) aux)[18 /*4*4+2*/] = (char) smax; // code displacement
	((char*) aux)[19 /*4*4+3*/] = (char) codes; // number of code copies

	size = 5*RS_I_BUS_B;

	TRY(
		XAxiDma_SimpleTransfer(dma, aux, size, XAXIDMA_DMA_TO_DEVICE),
		"ERROR TX DMA: VAR: last: %d\n")

	while(XAxiDma_Busy(dma, XAXIDMA_DMA_TO_DEVICE));

	// ** INPUT **

	size = RS_I_TYP_B*clen;

	TRY(
		XAxiDma_SimpleTransfer(dma, code, size, XAXIDMA_DMA_TO_DEVICE),
		"ERROR TX DMA: VAR: signal: %d\n")

	while(XAxiDma_Busy(dma, XAXIDMA_DMA_TO_DEVICE));
}

#define MC_I_BUS_b 32
#define MC_I_BUS_B 4
#define MC_O_BUS_b 64
#define MC_O_BUS_B 8
#define MC_I_TYP_B 1
#define MC_O_TYP_B 2

void transaction_mixcarr(XAxiDma *dma, char* data, int n, int dtype_, double phase_, double phase_step_, int ns)
{
	int status, size, i;

	// ** ARGUMENTS **

	int last = MC_I_TYP_B*(n/MC_I_BUS_B) - 1;
    int data_type = dtype_;
    long long int phase = (long long int) (phase_ * pow(2,44)); // Q20.44 (64b);
    long long int phase_step = (long long int) (phase_step_ * pow(2,44)); // Q20.44 (64b)

    aux[0] = last;
    aux[1] = data_type;
    aux[2] = ((int*) &phase)[0];
    aux[3] = ((int*) &phase)[1];
    aux[4] = ((int*) &phase_step)[0];
    aux[5] = ((int*) &phase_step)[1];
    aux[6] = CODE_COPIES;

    size = 7*MC_I_BUS_B; // Last (32) + Data Type (32) + Phase (32 + 32) + Phase Step (32 + 32)

	TRY(
		XAxiDma_SimpleTransfer(dma, aux, size, XAXIDMA_DMA_TO_DEVICE),
		"DMA: TX: last: error: %d\n")

	while(XAxiDma_Busy(dma, XAXIDMA_DMA_TO_DEVICE) == TRUE); /* Wait */

	// ** INPUT **

	size = MC_I_TYP_B*n;

	for(i = 0; i < 2*ns+1; i++) {
		TRY(
			XAxiDma_SimpleTransfer(dma, data, size, XAXIDMA_DMA_TO_DEVICE),
			"DMA: TX: signal: error: %d\n")

		while(XAxiDma_Busy(dma, XAXIDMA_DMA_TO_DEVICE) == TRUE); /* Wait */
	}
}
