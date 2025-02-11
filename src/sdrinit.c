/*------------------------------------------------------------------------------
* sdrinit.c : SDR initialize/cleanup functions
*
* Copyright (C) 2014 Taro Suzuki <gnsssdrlib@gmail.com>
*-----------------------------------------------------------------------------*/
#include "sdr.h"

/* sdr structs */
sdrini_t sdrini={0};
sdrstat_t sdrstat={0};
sdrch_t sdrch[MAXSAT]={{0}};
sdrspec_t sdrspec={0};
sdrout_t sdrout={0};

/* read ini file -------------------------------------------------------------*/
static int GetFileAttributes(const char *file)
{
#if defined(VITIS)
#else
    FILE *fp;
    if (!(fp=fopen(file,"r"))) return -1;
    fclose(fp);
#endif
    return 0;
}

static void GetPrivateProfileString(const char *sec, const char *key,
    const char *def, char *str, int len, const char *file)
{
    FILE *fp;
    char buff[1024],*p,*q;
    int enter=0;
#if defined(VITIS)

#else
    strncpy(str,def,len-1); str[len-1]='\0';

    if (!(fp=fopen(file,"r"))) {
        fprintf(stderr,"ini file open error [%s]\n",file);
        return;
    }
    while (fgets(buff,sizeof(buff),fp)) {
        if ((p=strchr(buff,';'))) *p='\0';
        if ((p=strchr(buff,'['))&&(q=strchr(p+1,']'))) {
            *q='\0';
            enter=!strcmp(p+1,sec);
        }
        else if (enter&&(p=strchr(buff,'='))) {
            *p='\0';
            for (q=p-1;q>=buff&&(*q==' '||*q=='\t');) *q--='\0';
            if (strcmp(buff,key)) continue;
            for (q=p+1+strlen(p+1)-1;q>=p+1&&(*q=='\r'||*q=='\n');) *q--='\0';
            strncpy(str,p+1,len-1); str[len-1]='\0';
            break;
        }
    }
    fclose(fp);
#endif
}

/* functions used in read ini file ---------------------------------------------
* note : these functions are only used in CLI application
*-----------------------------------------------------------------------------*/
int splitint(char *src, char *dlm, int *out, int n)
{
    int i;
    char *str;
    for (i=0;i<n;i++) {
        if ((str=strtok((i==0)?src:NULL,dlm))==NULL) return -1;
        out[i]=atoi(str);
    }
    return 0;
}
int splitdouble(char *src, char *dlm, double *out, int n)
{
    int i;
    char *str;
    for (i=0;i<n;i++) {
        if ((str=strtok((i==0)?src:NULL,dlm))==NULL) return -1;
        out[i]=atof(str);
    }
    return 0;
}
int readiniint(char *file, char *sec, char *key)
{
    char str[256];
    GetPrivateProfileString(sec,key,"",str,256,file);
    return atoi(str);
}
int readiniints(char *file, char *sec, char *key, int *out, int n)
{
    char str[256];
    GetPrivateProfileString(sec,key,"",str,256,file);
    return splitint(str,",",out,n);
}
double readinidouble(char *file, char *sec, char *key)
{
    char str[256];
    GetPrivateProfileString(sec,key,"",str,256,file);
    return atof(str);
}
int readinidoubles(char *file, char *sec, char *key, double *out, int n)
{
    char str[256];
    GetPrivateProfileString(sec,key,"",str,256,file);
    return splitdouble(str,",",out,n);
}
void readinistr(char *file, char *sec, char *key, char *out)
{
    GetPrivateProfileString(sec,key,"",out,256,file);
}
/* read ini file ---------------------------------------------------------------
* read ini file and set value to sdrini struct
* args   : sdrini_t *ini    I/0 sdrini struct
* return : int                  0:okay -1:error
* note : this function is only used in CLI application
*-----------------------------------------------------------------------------*/
extern int readinifile(sdrini_t *ini)
{
    int i;

    ini->fend	= FEND_FILE;
    //strcpy(ini->file1, "../../../../data/GPS_and_GIOVE_A-NN-fs16_3676-if4_1304.bin");
    strcpy(ini->file1, "../../../../data/gn3sv3_l1/gn3sv3_l1.bin");
    ini->useif1	= ON;
	// strcpy(ini->file2, "");
	ini->useif2	= OFF;
    ini->f_cf[0]= 1575.42e6;
    ini->f_sf[0]= 16.3676e6;
    ini->f_if[0]= 4.092e6;//4.1304e6;
    ini->dtype[0]=1;
    ini->f_cf[1]= 0.0;
    ini->f_sf[1]= 0.0;
    ini->f_if[1]= 0.0;
    ini->dtype[1]=0;

    /* RTL-SDR only */
    ini->rtlsdrppmerr=30;

    /* tracking parameter setting */
    ini->trkcorrn=3;
    ini->trkcorrd=3;
    ini->trkcorrp=6;
    ini->trkdllb[0]=5.0;
    ini->trkpllb[0]=30.0;
    ini->trkfllb[0]=200.0;
    ini->trkdllb[1]=1.0;
    ini->trkpllb[1]=10.0;
    ini->trkfllb[1]=50.0;
    
    /* channel setting */
    ini->nch=1;

	ini->prn[0]   = 1;//3;
	//ini->prn[1]   = 4;//15;
	//ini->prn[2]   = 6;//16;
	//ini->prn[3]   = 20;//18;
    for(int i = 0; i < 1; i++) {
		ini->sys[i]   = 1;
		ini->ctype[i] = 1;
		ini->ftype[i] = 1;
	}

    /* output setting */
    ini->outms   =100;
    ini->rinex   =1;
    ini->rtcm    =0;
    ini->lex     =0;
    ini->sbas    =0;
    ini->log     =0;
    strcpy(ini->rinexpath, "./");
    ini->rtcmport=9999;
    ini->lexport =9998;
    ini->sbasport=9997;

    /* sdr channel setting */
    for (i=0;i<ini->nch;i++) {
        if (ini->ctype[i]==CTYPE_L1CA ||
            ini->ctype[i]==CTYPE_G1 ||
            ini->ctype[i]==CTYPE_E1B ||
            ini->ctype[i]==CTYPE_B1I
            ) {
            ini->nchL1++;
        }
        if (ini->ctype[i]==CTYPE_LEXS) {
            ini->nchL6++;
        }
    }
    return 0;
}
/* check initial value ---------------------------------------------------------
* checking value in sdrini struct
* args   : sdrini_t *ini    I   sdrini struct
* return : int                  0:okay -1:error
*-----------------------------------------------------------------------------*/
extern int chk_initvalue(sdrini_t *ini)
{
    int ret;

    /* checking frequency input */
    if ((ini->f_sf[0]<=0||ini->f_sf[0]>100e6) ||
        (ini->f_if[0]<0 ||ini->f_if[0]>100e6)) {
            SDRPRINTF("error: wrong freq. input sf1: %.0f if1: %.0f\n",
                ini->f_sf[0],ini->f_if[0]);
            return -1;
    }

    /* checking frequency input */
    if(ini->useif2||ini->fend==FEND_STEREO) {
        if ((ini->f_sf[1]<=0||ini->f_sf[1]>100e6) ||
            (ini->f_if[1]<0 ||ini->f_if[1]>100e6)) {
                SDRPRINTF("error: wrong freq. input sf2: %.0f if2: %.0f\n",
                    ini->f_sf[1],ini->f_if[1]);
                return -1;
        }
    }

    /* checking filepath */
    if (ini->fend==FEND_FILE   ||ini->fend==FEND_FSTEREO||
        ini->fend==FEND_FGN3SV2||ini->fend==FEND_FGN3SV2||
        ini->fend==FEND_FRTLSDR||ini->fend==FEND_FBLADERF) {
        if (ini->useif1&&((ret=GetFileAttributes(ini->file1))<0)){
            SDRPRINTF("error: file1 doesn't exist: %s\n",ini->file1);
            return -1;
        }
        if (ini->useif2&&((ret=GetFileAttributes(ini->file2))<0)){
            SDRPRINTF("error: file2 doesn't exist: %s\n",ini->file2);
            return -1;
        }
        if ((!ini->useif1)&&(!ini->useif2)) {
            SDRPRINTF("error: file1 or file2 are not selected\n");
            return -1;
        }
    }

    /* checking rinex directory */
    if (ini->rinex) {
        if ((ret=GetFileAttributes(ini->rinexpath))<0) {
            SDRPRINTF("error: rinex output directory doesn't exist: %s\n",
                ini->rinexpath);
            return -1;
        }
    }

    return 0;
}
/* initialize acquisition struct -----------------------------------------------
* set value to acquisition struct
* args   : int sys          I   system type (SYS_GPS...)
*          int ctype        I   code type (CTYPE_L1CA...)
*          int prn          I   PRN
*          sdracq_t *acq    I/0 acquisition struct
* return : none
*-----------------------------------------------------------------------------*/
extern void initacqstruct(int sys, int ctype, int prn, sdracq_t *acq)
{
    if (ctype==CTYPE_L1CA) acq->intg=ACQINTG_L1CA;
    if (ctype==CTYPE_G1)   acq->intg=ACQINTG_G1;
    if (ctype==CTYPE_E1B)  acq->intg=ACQINTG_E1B;
    if (ctype==CTYPE_B1I)  acq->intg=ACQINTG_B1I;
    if (ctype==CTYPE_L1SAIF||ctype==CTYPE_L1SBAS) acq->intg=ACQINTG_SBAS;

    acq->hband=ACQHBAND;
    acq->step=ACQSTEP;
    acq->nfreq=2*(ACQHBAND/ACQSTEP)+1;
}
/* initialize tracking parameter struct ----------------------------------------
* set value to tracking parameter struct
* args   : sdrtrk_t *trk    I/0 tracking struct
* return : none
*-----------------------------------------------------------------------------*/
extern void inittrkprmstruct(sdrtrk_t *trk)
{
    /* set tracking parameter */
    trk->prm1.dllb=sdrini.trkdllb[0];
    trk->prm1.pllb=sdrini.trkpllb[0];
    trk->prm1.fllb=sdrini.trkfllb[0];
    trk->prm2.dllb=sdrini.trkdllb[1];
    trk->prm2.pllb=sdrini.trkpllb[1];
    trk->prm2.fllb=sdrini.trkfllb[1];

    /* calculation loop filter parameters (before nav frame synchronization) */
    trk->prm1.dllw2=(trk->prm1.dllb/0.53)*(trk->prm1.dllb/0.53);
    trk->prm1.dllaw=1.414*(trk->prm1.dllb/0.53);
    trk->prm1.pllw2=(trk->prm1.pllb/0.53)*(trk->prm1.pllb/0.53);
    trk->prm1.pllaw=1.414*(trk->prm1.pllb/0.53);
    trk->prm1.fllw =trk->prm1.fllb/0.25;
    
    /* calculation loop filter parameters (after nav frame synchronization) */
    trk->prm2.dllw2=(trk->prm2.dllb/0.53)*(trk->prm2.dllb/0.53);
    trk->prm2.dllaw=1.414*(trk->prm2.dllb/0.53);
    trk->prm2.pllw2=(trk->prm2.pllb/0.53)*(trk->prm2.pllb/0.53);
    trk->prm2.pllaw=1.414*(trk->prm2.pllb/0.53);
    trk->prm2.fllw =trk->prm2.fllb/0.25;
}
/* initialize tracking struct --------------------------------------------------
* set value to tracking struct
* args   : int    sat       I   satellite number
*          int    ctype     I   code type (CTYPE_L1CA...)
*          double ctime     I   code period (s)
*          sdrtrk_t *trk    I/0 tracking struct
* return : int                  0:okay -1:error
*-----------------------------------------------------------------------------*/
extern int inittrkstruct(int sat, int ctype, double ctime, sdrtrk_t *trk)
{
    int i,prn;
    int ctimems=(int)(ctime*1000);
    int sys=satsys(sat,&prn);

    /* set tracking parameter */
    inittrkprmstruct(trk);

    /* correlation point */
    trk->corrn=sdrini.trkcorrn;
    trk->corrp=(int *)malloc(sizeof(int)*trk->corrn);
    for (i=0;i<trk->corrn;i++) {
        trk->corrp[i]=sdrini.trkcorrd*(i+1);
        if (trk->corrp[i]==sdrini.trkcorrp){
            trk->ne=2*(i+1)-1; /* Early */
            trk->nl=2*(i+1);   /* Late */
        }
    }
    /* correlation point for plot */
    (trk->corrx=(double *)calloc(2*trk->corrn+1,sizeof(double)));
    for (i=1;i<=trk->corrn;i++) {
        trk->corrx[2*i-1]=-sdrini.trkcorrd*i;
        trk->corrx[2*i  ]= sdrini.trkcorrd*i;
    }

    trk->II     =(double*)calloc(1+2*trk->corrn,sizeof(double));
    trk->QQ     =(double*)calloc(1+2*trk->corrn,sizeof(double));
    trk->oldI   =(double*)calloc(1+2*trk->corrn,sizeof(double));
    trk->oldQ   =(double*)calloc(1+2*trk->corrn,sizeof(double));
    trk->sumI   =(double*)calloc(1+2*trk->corrn,sizeof(double));
    trk->sumQ   =(double*)calloc(1+2*trk->corrn,sizeof(double));
    trk->oldsumI=(double*)calloc(1+2*trk->corrn,sizeof(double));
    trk->oldsumQ=(double*)calloc(1+2*trk->corrn,sizeof(double));

    if (ctype==CTYPE_L1CA)   trk->loop=LOOP_L1CA;
    if (ctype==CTYPE_G1)     trk->loop=LOOP_G1;
    if (ctype==CTYPE_E1B)    trk->loop=LOOP_E1B;
    if (ctype==CTYPE_L1SAIF) trk->loop=LOOP_SBAS;
    if (ctype==CTYPE_L1SBAS) trk->loop=LOOP_SBAS;
    if (ctype==CTYPE_B1I&&prn>5 ) trk->loop=LOOP_B1I;
    if (ctype==CTYPE_B1I&&prn<=5) trk->loop=LOOP_B1IG;
    
    /* for LEX */
    if (sys==SYS_QZS&&ctype==CTYPE_L1CA&&sdrini.nchL6) trk->loop=LOOP_LEX;
    
    /* loop interval (ms) */
    trk->loopms=trk->loop*ctimems;

    if (!trk->II||!trk->QQ||!trk->oldI||!trk->oldQ||!trk->sumI||!trk->sumQ||
        !trk->oldsumI||!trk->oldsumQ) {
        SDRPRINTF("error: inittrkstruct memory allocation\n");
        return -1;
    }
    return 0;
}
/* initialize navigation struct ------------------------------------------------
* set value to navigation struct
* args   : int sys          I   system type (SYS_GPS...)
*          int ctype        I   code type (CTYPE_L1CA...)
*          int prn          I   PRN (or SV) number
*          sdrnav_t *nav    I/0 navigation struct
* return : int                  0:okay -1:error
*-----------------------------------------------------------------------------*/
extern int initnavstruct(int sys, int ctype, int prn, sdrnav_t *nav)
{
    int i;
    int pre_l1ca[8]= { 1,-1,-1,-1, 1,-1, 1, 1}; /* L1CA preamble*/
    int pre_e1b[10]= { 1,-1, 1,-1,-1, 1, 1, 1, 1, 1}; /* E1B preamble */
    int pre_g1[30]=  {-1,-1,-1,-1,-1, 1, 1, 1,-1,-1,
                       1,-1,-1,-1, 1,-1, 1,-1, 1, 1,
                       1, 1,-1, 1, 1,-1, 1,-1,-1, 1}; /* G1 preamble */
    int pre_b1i[11]= {-1,-1,-1, 1, 1, 1,-1, 1, 1,-1, 1}; /* B1I preamble */
    int pre_sbs[24]= { 1,-1, 1,-1, 1, 1,-1,-1,-1, 1,
                       1,-1,-1, 1,-1, 1,-1,-1, 1, 1,
                       1 -1,-1, 1}; /* SBAS L1/QZS L1SAIF preamble */

    int poly[2]={V27POLYA,V27POLYB};

    nav->ctype=ctype;
    nav->sdreph.ctype=ctype;
    nav->sdreph.prn=prn;
    nav->sdreph.eph.iodc=-1;

    /* GPS/QZS L1CA */
    if (ctype==CTYPE_L1CA) {
        nav->rate=NAVRATE_L1CA;
        nav->flen=NAVFLEN_L1CA;
        nav->addflen=NAVADDFLEN_L1CA;
        nav->prelen=NAVPRELEN_L1CA;
        nav->sdreph.cntth=NAVEPHCNT_L1CA;
        nav->update=(int)(nav->flen*nav->rate);
        memcpy(nav->prebits,pre_l1ca,sizeof(int)*nav->prelen);

        /* overlay code (all 1) */
        nav->ocode=(short *)calloc(nav->rate,sizeof(short));
        for (i=0;i<nav->rate;i++) nav->ocode[i]=1;

	/* GEO (D2 NAV) */
	} else {
		nav->rate=NAVRATE_B1IG;
		nav->flen=NAVFLEN_B1IG;
		nav->addflen=NAVADDFLEN_B1IG;
		nav->prelen=NAVPRELEN_B1IG;
		nav->sdreph.cntth=NAVEPHCNT_B1IG;
		nav->update=(int)(nav->flen*nav->rate);
		memcpy(nav->prebits,pre_b1i,sizeof(int)*nav->prelen);

		/* overlay code (all 1) */
		nav->ocode=(short *)calloc(nav->rate,sizeof(short));
		for (i=0;i<nav->rate;i++) nav->ocode[i]=1;
	}

    if (!(nav->bitsync= (int *)calloc(nav->rate,sizeof(int))) || 
        !(nav->fbits=   (int *)calloc(nav->flen+nav->addflen,sizeof(int))) ||
        !(nav->fbitsdec=(int *)calloc(nav->flen+nav->addflen,sizeof(int)))) {
            SDRPRINTF("error: initnavstruct memory alocation\n");
            return -1;
    }
    return 0;
}
/* initialize sdr channel struct -----------------------------------------------
* set value to sdr channel struct
* args   : int    chno      I   channel number (1,2,...)
*          int    sys       I   system type (SYS_***)
*          int    prn       I   PRN number
*          int    ctype     I   code type (CTYPE_***)
*          int    dtype     I   data type (DTYPEI or DTYPEIQ)
*          int    ftype     I   front end type (FTYPE1 or FTYPE2)
*          double f_cf      I   center (carrier) frequency (Hz)
*          double f_sf      I   sampling frequency (Hz)
*          double f_if      I   intermidiate frequency (Hz)
*          sdrch_t *sdr     I/0 sdr channel struct
* return : int                  0:okay -1:error
*-----------------------------------------------------------------------------*/
extern int initsdrch(int chno, int sys, int prn, int ctype, int dtype,
                     int ftype, double f_cf, double f_sf, double f_if,
                     sdrch_t *sdr)
{
    int i;
    short *rcode;

    sdr->no=chno;
    sdr->sys=sys;
    sdr->prn=prn;
    sdr->sat=satno(sys,prn);
    sdr->ctype=ctype;
    sdr->dtype=dtype;
    sdr->ftype=ftype;
    sdr->f_sf=f_sf;
    sdr->f_if=f_if;
    sdr->ti=1/f_sf;

    sdr->flagacq=OFF;
    
    /* code generation */
    sdr->code=gencode(prn,ctype,&sdr->clen,&sdr->crate);
    sdr->ci=sdr->ti*sdr->crate;
    sdr->ctime=sdr->clen/sdr->crate;
    sdr->nsamp=(int)(f_sf*sdr->ctime);
    sdr->nsampchip=(int)(sdr->nsamp/sdr->clen);
    satno2id(sdr->sat,sdr->satstr);

    /* set carrier frequency */
    if (ctype==CTYPE_G1) {
        sprintf(sdr->satstr,"R%d",prn); /* frequency number instead of PRN */
        sdr->f_cf=FREQ1_GLO+DFRQ1_GLO*prn; /* FDMA */
        sdr->foffset=DFRQ1_GLO*prn; /* FDMA */
    } else if (sdrini.fend==FEND_FRTLSDR) {
        sdr->foffset=f_cf*sdrini.rtlsdrppmerr*1e-6;
    } else {
        sdr->f_cf=f_cf; /* carrier frequency */
        sdr->foffset=0.0; /* frequency offset */
    }

    /* for BeiDou B1I */
    if (ctype==CTYPE_B1I) sdr->nsampchip*=2; /* for BOC code */

    /* acqisition struct */
    initacqstruct(sys,ctype,prn,&sdr->acq);
    sdr->acq.nfft=2*sdr->nsamp;//calcfftnum(2*sdr->nsamp,0);

    /* memory allocation */
    if (!(sdr->acq.freq=(double*)malloc(sizeof(double)*sdr->acq.nfreq))) {
        SDRPRINTF("error: initsdrch memory alocation\n"); return -1;
    }

    /* doppler search frequency */
    for (i=0;i<sdr->acq.nfreq;i++)
        sdr->acq.freq[i]=sdr->f_if+((i-(sdr->acq.nfreq-1)/2)*sdr->acq.step)
                            +sdr->foffset;

    /* tracking struct */
    if (inittrkstruct(sdr->sat,ctype,sdr->ctime,&sdr->trk)<0) return -1;

    /* navigation struct */
    if (initnavstruct(sys,ctype,prn,&sdr->nav)<0) {
        return -1;
    }
    /* memory allocation */
    if (!(rcode=(short *)sdrmalloc(sizeof(short)*sdr->acq.nfft)) || 
        !(sdr->xcode=cpxmalloc(sdr->acq.nfft))) {
            SDRPRINTF("error: initsdrch memory alocation\n"); return -1;
    }
    /* other code generation */
    for (i=0;i<sdr->acq.nfft;i++) rcode[i]=0; /* zero padding */
    rescode(sdr->code,sdr->clen,0,0,sdr->ci,sdr->nsamp,rcode); /* resampling */
    cpxcpx(rcode,NULL,1.0,sdr->acq.nfft,sdr->xcode); /* FFT for acquisition */
    cpxfft(sdr->xcode,sdr->acq.nfft);

    sdrfree(rcode);
    return 0;
}
/* free sdr channel struct -----------------------------------------------------
* free memory in sdr channel struct
* args   : sdrch_t *sdr     I/0 sdr channel struct
* return : none 
*-----------------------------------------------------------------------------*/
extern void freesdrch(sdrch_t *sdr)
{
    free(sdr->code);
    cpxfree(sdr->xcode);
    free(sdr->nav.fbits);
    free(sdr->nav.fbitsdec);
    free(sdr->nav.bitsync);
    free(sdr->trk.II);
    free(sdr->trk.QQ);
    free(sdr->trk.oldI);
    free(sdr->trk.oldQ);
    free(sdr->trk.sumI);
    free(sdr->trk.sumQ);
    free(sdr->trk.oldsumI);
    free(sdr->trk.oldsumQ);
    free(sdr->trk.corrp);
    free(sdr->acq.freq);

    if (sdr->nav.fec!=NULL)
        // delete_viterbi27_port(sdr->nav.fec);

    if (sdr->nav.ocode!=NULL)
        free(sdr->nav.ocode);
}
