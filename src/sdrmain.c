/*------------------------------------------------------------------------------
* sdrmain.c : SDR main functions
*
* Copyright (C) 2014 Taro Suzuki <gnsssdrlib@gmail.com>
*-----------------------------------------------------------------------------*/
#include "sdr.h"

/* global variables ----------------------------------------------------------*/
#ifdef GUI
GCHandle hform;
#endif

/* thread handle and mutex */
thread_t hmainthread;
thread_t hsyncthread;
thread_t hspecthread;
thread_t hkeythread;
mlock_t hbuffmtx;
mlock_t hreadmtx;
mlock_t hfftmtx;
mlock_t hpltmtx;
mlock_t hobsmtx;
mlock_t hlexmtx;
event_t hlexeve;

/* sdr structs */
sdrini_t sdrini={0};
sdrstat_t sdrstat={0};
sdrch_t sdrch[MAXSAT]={{0}};
sdrspec_t sdrspec={0};
sdrout_t sdrout={0};

/* initsdrgui ------------------------------------------------------------------
* initialize sdr gui application  
* args   : maindlg^ form       I   main dialog class
*          sdrini_t* sdrinigui I   sdr init struct
* return : none
* note : This function is only used in GUI application 
*-----------------------------------------------------------------------------*/
#ifdef GUI
extern void initsdrgui(maindlg^ form, sdrini_t* sdrinigui)
{
    /* initialization global structs */
    memset(&sdrini,0,sizeof(sdrini));
    memset(&sdrstat,0,sizeof(sdrstat));
    memset(&sdrch,0,sizeof(sdrch));

    hform=GCHandle::Alloc(form);
    memcpy(&sdrini,sdrinigui,sizeof(sdrini_t)); /* copy setting from GUI */
}
#else
/* keyboard thread -------------------------------------------------------------
* keyboard thread for program termination  
* args   : void   *arg      I   not used
* return : none
* note : this thread is only created in CLI application
*-----------------------------------------------------------------------------*/
#ifdef WIN32
extern void keythread(void * arg) 
#else
extern void *keythread(void * arg) 
#endif
{
    do {
        switch(getchar()) {
        case 'q':
        case 'Q':
            sdrstat.stopflag=1;
            break;
        default:
            SDRPRINTF("press 'q' to exit...\n");
            break;
        }
    } while (!sdrstat.stopflag);

    return THRETVAL;
}
/* main function ---------------------------------------------------------------
* main entry point in CLI application  
* args   : none
* return : none
* note : This function is only used in CLI application 
*-----------------------------------------------------------------------------*/
int main(int argc, char **argv)
{
    /* read ini file */
    if (readinifile(&sdrini)<0) {
        return -1; 
    }
    //cratethread(hkeythread,keythread,NULL);

    startsdr();

    return 0;
}
#endif
/* sdr start -------------------------------------------------------------------
* start sdr function  
* args   : void   *arg      I   not used
* return : none
* note : This function is called as thread in GUI application and is called as
*        function in CLI application
*-----------------------------------------------------------------------------*/
#ifdef GUI
extern void startsdr(void *arg) /* call as thread */
#else
extern void startsdr(void) /* call as function */
#endif
{
    //sdrini.nch = 1;

    int i;
    SDRPRINTF("GNSS-SDRLIB start!\n");

    /* check initial value */
    if (chk_initvalue(&sdrini)<0) {
        SDRPRINTF("error: chk_initvalue\n");
        quitsdr(&sdrini,1);
        return;
    }

    /* receiver initialization */
    if (rcvinit(&sdrini)<0) {
        SDRPRINTF("error: rcvinit\n");
        quitsdr(&sdrini,1);
        return;
    }
    /* initialize sdr channel struct */
    for (i=0;i<sdrini.nch;i++) {
        if (initsdrch(i+1,sdrini.sys[i],sdrini.prn[i],sdrini.ctype[i],
            sdrini.dtype[sdrini.ftype[i]-1],sdrini.ftype[i],
            sdrini.f_cf[sdrini.ftype[i]-1],sdrini.f_sf[sdrini.ftype[i]-1],
            sdrini.f_if[sdrini.ftype[i]-1],&sdrch[i])<0) {
            
            SDRPRINTF("error: initsdrch\n");
            quitsdr(&sdrini,2);
            return;
        }
    }

    /* mutexes and events */
    openhandles();

    /* start grabber */
    if (rcvgrabstart(&sdrini)<0) {
        quitsdr(&sdrini,4);
        return;
    }

    //cratethread(sdrch[0].hsdr,sdrthread,&sdrch[0]);

    sdrthread(NULL);

    //while(!sdrstat.stopflag) {
    //    if(rcvgrabdata(&sdrini) < 0)
    //        sdrstat.stopflag = ON;
    //}

    /* wait thereds */
    //waitthread(hsyncthread);
    //for (i=0;i<sdrini.nch;i++) 
    //    waitthread(sdrch[i].hsdr);

    /* sdr termination */
    quitsdr(&sdrini,0);

    SDRPRINTF("GNSS-SDRLIB is finished!\n");
}
/* sdr termination -------------------------------------------------------------
* sdr termination process  
* args   : sdrini_t *ini    I   sdr init struct
* args   : int    stop      I   stop position in function 0: run all
* return : none
*-----------------------------------------------------------------------------*/
extern void quitsdr(sdrini_t *ini, int stop)
{
    int i;

    if (stop==1) return;

    /* sdr termination */
    rcvquit(ini);
    if (stop==2) return;

    /* free memory */
    for (i=0;i<ini->nch;i++) freesdrch(&sdrch[i]);
    if (stop==3) return;

    /* mutexes and events */
    closehandles();
    if (stop==4) return;
}
/* sdr channel thread ----------------------------------------------------------
* sdr channel thread for signal acquisition and tracking  
* args   : void   *arg      I   sdr channel struct
* return : none
* note : This thread handles the acquisition and tracking of one of the signals. 
*        The thread is created at startsdr function.
*-----------------------------------------------------------------------------*/
extern void *sdrthread(void *arg)
{
    sdrch_t *sdr;
    uint64_t *cnt = malloc(sizeof(uint64_t)*sdrini.nch), *loopcnt = malloc(sizeof(uint64_t)*sdrini.nch);
    double *acqpower=NULL;
    FILE** fp = (FILE**) malloc(sizeof(FILE*)*sdrini.nch);
    char fname[100];
    
    int k;
    int ontrk = 0;
    int request = ontrk;

    /* create tracking log file */
    for(k = 0; k < sdrini.nch; k++)
    {
        cnt[k] = 0; loopcnt[k] = 0;
        sdr = &sdrch[k];

        if (sdrini.log) {
            sprintf(fname,"log%s.csv",sdr->satstr);
            if((fp[k]=createlog(fname,&sdr->trk))==NULL) {
                SDRPRINTF("error: invailed log file: %s\n",fname);
                return THRETVAL;
            }
        }
    }

    //sleepms(sdr->no*500);
    //SDRPRINTF("**** %s sdr thread %d start! ****\n",sdr->satstr,sdr->no);

    while (!sdrstat.stopflag) {       

        if(request == ontrk)
            file_pushtomembuf();

        request = 0;
        ontrk = 0;

        for(k = 0, sdr = &sdrch[0], syncthread(NULL);
                k < sdrini.nch; sdr = &sdrch[++k])
        {
            /* acquisition */
            
            if (!sdr->flagacq) {               
                /* memory allocation */
                if (acqpower!=NULL) free(acqpower);
                acqpower = (double*) calloc(sizeof(double),sdr->nsamp*sdr->acq.nfreq);

                /* fft correlation */
                sdr->trk.buffloc=sdraccuisition(sdr,acqpower);
            }

            /* tracking */

            if (sdr->flagacq) {
                ontrk ++;

                sdr->trk.bufflocnow = sdrtracking(sdr,sdr->trk.buffloc,cnt[k]);
                
                if(sdr->trk.bufflocnow<=sdr->trk.buffloc)
                    request ++;

                if (sdr->flagtrk) {

                    /* correlation output accumulation */
                    cumsumcorr(&sdr->trk,sdr->nav.ocode[sdr->nav.ocodei]);

                    sdr->trk.flagloopfilter=0;
                    if (!sdr->nav.flagsync) {
                        pll(sdr,&sdr->trk.prm1,sdr->ctime); /* PLL */
                        dll(sdr,&sdr->trk.prm1,sdr->ctime); /* DLL */
                        sdr->trk.flagloopfilter=1;
                    }
                    else if (sdr->nav.swloop) {
                        pll(sdr,&sdr->trk.prm2,(double)sdr->trk.loopms/1000);
                        dll(sdr,&sdr->trk.prm2,(double)sdr->trk.loopms/1000);
                        sdr->trk.flagloopfilter=2;

                        mlock(hobsmtx);

                        /* calculate observation data */
                        if (loopcnt[k]%(SNSMOOTHMS/sdr->trk.loopms)==0)
                            setobsdata(sdr,sdr->trk.buffloc,cnt[k],&sdr->trk,1);
                        else
                            setobsdata(sdr,sdr->trk.buffloc,cnt[k],&sdr->trk,0);

                        unmlock(hobsmtx);

                        /* LEX thread */
                        if (sdrini.nchL6!=0&&sdr->no==sdrini.nch+1&&loopcnt[k]>250) 
                            setevent(hlexeve);

                        loopcnt[k]++;
                    }

                    if (cnt[k]%(1000*10)==0)
                        SDRPRINTF("process %d sec...\n",(int)cnt[k]/(1000));

                    /* write tracking log */
                    if (sdrini.log) writelog(fp[k],&sdr->trk,&sdr->nav);

                    if (sdr->trk.flagloopfilter) clearcumsumcorr(&sdr->trk);
                    cnt[k]++;
                    sdr->trk.buffloc+=sdr->currnsamp;
                }
            }
        }
    }
    
    //if (sdrini.nchL6!=0&&sdr->no==sdrini.nch+1) 
    //    setevent(hlexeve);
    
    /* plot termination */
    // quitpltstruct(&pltacq,&plttrk);

    /* cloase tracking log file */
    if (sdrini.log) {
        for(k = 0; k < sdrini.nch; k ++)
            closelog(fp[k]);
    }

    if (sdr->flagacq) {
        SDRPRINTF("SDR channel %s thread finished! Delay=%d [ms]\n",
            sdr->satstr,(int)(sdr->trk.bufflocnow-sdr->trk.buffloc)/sdr->nsamp);
    } else {
        SDRPRINTF("SDR channel %s thread finished!\n",sdr->satstr);
    }

    return THRETVAL;
}
