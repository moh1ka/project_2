/*************************************************************************************
			       DEPARTMENT OF ELECTRICAL AND ELECTRONIC ENGINEERING
					   		     IMPERIAL COLLEGE LONDON 

 				      EE 3.19: Real Time Digital Signal Processing
					       Dr Paul Mitcheson and Daniel Harvey

				        		 PROJECT: Frame Processing

 				            ********* ENHANCE. C **********
							 Shell for speech enhancement 

  		Demonstrates overlap-add frame processing (interrupt driven) on the DSK. 

 *************************************************************************************
 				             By Danny Harvey: 21 July 2006
							 Updated for use on CCS v4 Sept 2010
 ************************************************************************************/

/*
 *	You should modify the code so that a speech enhancement project is built 
 *  on top of this template.
 */
/**************************** Pre-processor statements ******************************/
//  library required when using calloc
#include <stdlib.h>
//  Included so program can make use of DSP/BIOS configuration tool.  
#include "dsp_bios_cfg.h"

/* The file dsk6713.h must be included in every program that uses the BSL.  This 
   example also includes dsk6713_aic23.h because it uses the 
   AIC23 codec module (audio interface). */
#include "dsk6713.h"
#include "dsk6713_aic23.h"

// math library (trig functions)
#include <math.h>

/* Some functions to help with Complex algebra and FFT. */
#include "cmplx.h"      
#include "fft_functions.h"  

// Some functions to help with writing/reading the audio ports when using interrupts.
#include <helper_functions_ISR.h>

#define WINCONST 0.85185			/* 0.46/0.54 for Hamming window */
#define FSAMP 8000.0		/* sample frequency, ensure this matches Config for AIC */
#define FFTLEN 256					/* fft length = frame length 256/8000 = 32 ms*/
#define NFREQ (1+FFTLEN/2)			/* number of frequency bins from a real FFT */
#define OVERSAMP 4					/* oversampling ratio (2 or 4) */  
#define FRAMEINC (FFTLEN/OVERSAMP)	/* Frame increment */
#define CIRCBUF (FFTLEN+FRAMEINC)	/* length of I/O buffers */

#define OUTGAIN 16000.0				/* Output gain for DAC */
#define INGAIN  (1.0/16000.0)		/* Input gain for ADC  */
// PI defined here for use in your code 
#define PI 3.141592653589793
#define TFRAME FRAMEINC/FSAMP       /* time between calculation of each frame */
/******************************* Enhancement Switches*********************/
//#define enhance_1 0;
//#define enhance_2 1;
//#define enhance_3 0;
/******************************* Global declarations ********************************/

/* Audio port configuration settings: these values set registers in the AIC23 audio 
   interface to configure it. See TI doc SLWS106D 3-3 to 3-10 for more info. */
DSK6713_AIC23_Config Config = { \
			 /**********************************************************************/
			 /*   REGISTER	            FUNCTION			      SETTINGS         */ 
			 /**********************************************************************/\
    0x0017,  /* 0 LEFTINVOL  Left line input channel volume  0dB                   */\
    0x0017,  /* 1 RIGHTINVOL Right line input channel volume 0dB                   */\
    0x01f9,  /* 2 LEFTHPVOL  Left channel headphone volume   0dB                   */\
    0x01f9,  /* 3 RIGHTHPVOL Right channel headphone volume  0dB                   */\
    0x0011,  /* 4 ANAPATH    Analog audio path control       DAC on, Mic boost 20dB*/\
    0x0000,  /* 5 DIGPATH    Digital audio path control      All Filters off       */\
    0x0000,  /* 6 DPOWERDOWN Power down control              All Hardware on       */\
    0x0043,  /* 7 DIGIF      Digital audio interface format  16 bit                */\
    0x008d,  /* 8 SAMPLERATE Sample rate control        8 KHZ-ensure matches FSAMP */\
    0x0001   /* 9 DIGACT     Digital interface activation    On                    */\
			 /**********************************************************************/
};

// Codec handle:- a variable used to identify audio interface  
DSK6713_AIC23_CodecHandle H_Codec;

float *inbuffer, *outbuffer;   		/* Input/output circular buffers */
float *inframe, *outframe;          /* Input and output frames *///and intermediate frame for processing
float *inwin, *outwin;              /* Input and output windows */
float ingain, outgain;				/* ADC and DAC gains */ 
float cpufrac; 						/* Fraction of CPU time used */
volatile int io_ptr=0;              /* Input/ouput pointer for circular buffers */
volatile int frame_ptr=0;           /* Frame pointer */
float *mag;
float *mag_min;

float *p;
float *x;

float *M1;
float *M2;
float *M3;
float *M4;
float *G;
float *G_old;
float *G_older;

float lambda = 0.1;
float alpha = 4;
float tau = 0.02;
float lambda_6 = 0.1;
float alpha_6 = 4;

int counter = 0;          
  
int process_on = 1;//0 = no processing
int enhance_1 = 0;
int enhance = 0;
int enhance_6 = 0;
/*********************************************************************************************/
float min (float a, float b){
	if(b<a){
		a=b;
	}
	return a;
}

float max (float a, float b){
	if(b>a){
		a=b;
	}
	return a;
}
 /******************************* Function prototypes *******************************/
void init_hardware(void);    	/* Initialize codec */ 
void init_HWI(void);            /* Initialize hardware interrupts */
void ISR_AIC(void);             /* Interrupt service routine for codec */
void process_frame(void);       /* Frame processing routine */
float min (float a, float b);     
float max (float a, float b);        
void lpf_1 (float *x, float *pbuf);
/********************************** Main routine ************************************/
void main()
{      

  	int k; // used in various for loops
/*  Initialize and zero fill arrays */  

	inbuffer	= (float *) calloc(CIRCBUF, sizeof(float));	/* Input array */
    outbuffer	= (float *) calloc(CIRCBUF, sizeof(float));	/* Output array */
	inframe		= (float *) calloc(FFTLEN, sizeof(float));	/* Array for processing*/
    outframe	= (float *) calloc(FFTLEN, sizeof(float));	/* Array for processing*/
    inwin		= (float *) calloc(FFTLEN, sizeof(float));	/* Input window */
    outwin		= (float *) calloc(FFTLEN, sizeof(float));	/* Output window */
	mag	        = (float *) calloc(FFTLEN, sizeof(float)); /* magnitude spectrum*/
	mag_min	    = (float *) calloc(FFTLEN, sizeof(float)); /* old mag spectrum*/
	
	p = (float *) calloc(FFTLEN, sizeof(float));
	x = (float *) calloc(FFTLEN, sizeof(float)); /* magnitude spectrum*/
	
	//Min noise buffers
	M1	        = (float *) calloc(FFTLEN, sizeof(float)); /* magnitude spectrum*/
	M2	        = (float *) calloc(FFTLEN, sizeof(float)); /* magnitude spectrum*/
	M3	        = (float *) calloc(FFTLEN, sizeof(float)); /* magnitude spectrum*/
	M4	        = (float *) calloc(FFTLEN, sizeof(float)); /* magnitude spectrum*/
	G	        = (float *) calloc(FFTLEN, sizeof(float)); /* magnitude spectrum*/
	G_old	        = (float *) calloc(FFTLEN, sizeof(float)); /* magnitude spectrum*/
	G_older	        = (float *) calloc(FFTLEN, sizeof(float)); /* magnitude spectrum*/
	/* initialize board and the audio port */
  	init_hardware();
  
  	/* initialize hardware interrupts */
  	init_HWI();    
  
/* initialize algorithm constants */  
                       
  	for (k=0;k<FFTLEN;k++)
	{                           
	inwin[k] = sqrt((1.0-WINCONST*cos(PI*(2*k+1)/FFTLEN))/OVERSAMP);
	outwin[k] = inwin[k]; 
	} 
  	ingain=INGAIN;
  	outgain=OUTGAIN;        

 							
  	/* main loop, wait for interrupt */  
  	while(1) 	process_frame();
}
    
/********************************** init_hardware() *********************************/  
void init_hardware()
{
    // Initialize the board support library, must be called first 
    DSK6713_init();
    
    // Start the AIC23 codec using the settings defined above in config 
    H_Codec = DSK6713_AIC23_openCodec(0, &Config);

	/* Function below sets the number of bits in word used by MSBSP (serial port) for 
	receives from AIC23 (audio port). We are using a 32 bit packet containing two 
	16 bit numbers hence 32BIT is set for  receive */
	MCBSP_FSETS(RCR1, RWDLEN1, 32BIT);	

	/* Configures interrupt to activate on each consecutive available 32 bits 
	from Audio port hence an interrupt is generated for each L & R sample pair */	
	MCBSP_FSETS(SPCR1, RINTM, FRM);

	/* These commands do the same thing as above but applied to data transfers to the 
	audio port */
	MCBSP_FSETS(XCR1, XWDLEN1, 32BIT);	
	MCBSP_FSETS(SPCR1, XINTM, FRM);	
	

}
/********************************** init_HWI() **************************************/ 
void init_HWI(void)
{
	IRQ_globalDisable();			// Globally disables interrupts
	IRQ_nmiEnable();				// Enables the NMI interrupt (used by the debugger)
	IRQ_map(IRQ_EVT_RINT1,4);		// Maps an event to a physical interrupt
	IRQ_enable(IRQ_EVT_RINT1);		// Enables the event
	IRQ_globalEnable();				// Globally enables interrupts

}
        
/******************************** process_frame() ***********************************/  
void process_frame(void)
{
	int k, m, i; 
	int io_ptr0;  
	float snr = 0;
	complex *c;

	/* work out fraction of available CPU time used by algorithm */    
	cpufrac = ((float) (io_ptr & (FRAMEINC - 1)))/FRAMEINC;  
		
	/* wait until io_ptr is at the start of the current frame */ 	
	while((io_ptr/FRAMEINC) != frame_ptr); 
	
	/* then increment the framecount (wrapping if required) */ 
	if (++frame_ptr >= (CIRCBUF/FRAMEINC)) frame_ptr=0;
 	
 	/* save a pointer to the position in the I/O buffers (inbuffer/outbuffer) where the 
 	data should be read (inbuffer) and saved (outbuffer) for the purpose of processing */
 	io_ptr0=frame_ptr * FRAMEINC;
	
	/* copy input data from inbuffer into inframe (starting from the pointer position) */ 
	 
	m=io_ptr0;
    for (k=0;k<FFTLEN;k++)
	{                           
		inframe[k] = inbuffer[m] * inwin[k]; 
		if (++m >= CIRCBUF) m=0; /* wrap if required */
	} 
	
	/************************* DO PROCESSING OF FRAME  HERE **************************/
	
	
	/* please add your code, at the moment the code simply copies the input to the 
	ouptut with no processing */	 
	
	c = (complex *) calloc(FFTLEN, sizeof(complex)); /* Array for processing*/
	
	//copy inframe into complex valued array c	
	for(i=0; i<FFTLEN; i++)
	{
		c[i].r = inframe[i];
	}
	
	fft(FFTLEN, c);
	
	if(process_on == 1)
	{
	
		counter++;
		//calculate magnitude
		for(i=0; i<FFTLEN; i++)
		{
			mag[i] = cabs(c[i]);
		}
		
		lpf_1(mag,p); //p stores low pass filtered input (mag)
	
		if(enhance_1 == 1)
		{	
			memcpy (x, p, FFTLEN*sizeof(float));
	
		}
		else
		{
			memcpy (x, mag, FFTLEN*sizeof(float));
		}
		
		memcpy (M1, x, FFTLEN*sizeof(float));
		
		for (i=0;i<FFTLEN;i++)
		{
			M1[i] = min(M1[i], x[i]);
		}
			
		if (counter == 79)
		{
			counter = 0;
			memcpy (M4, M3, FFTLEN*sizeof(float));
			memcpy (M3, M2, FFTLEN*sizeof(float));
			memcpy (M2, M1, FFTLEN*sizeof(float));
			memcpy (M1, M4, FFTLEN*sizeof(float));
		
		}
		
		for( i=0; i<FFTLEN;i++)
		{
			mag_min[i] =  min(M1[i], M2[i]); //finding min across M1-M4
			mag_min[i] =  min(M3[i], mag_min[i]);//finding min across M1-M4
			mag_min[i] = ( alpha * (min(M4[i], mag_min[i])) );//finding min across M1-M4
			switch (enhance)
			{
				case 0:
				G[i] = max( lambda,  ( 1 - (mag_min[i]/mag[i]) ) );
				if(enhance_8 == 1)
				{
					if( (mag_min[i]/mag[i]) > th)
					{
						G[i] = min(G[i],G_old[i]);
						G[i] = min(G[i],G_older[i]);
					}
					G_old[i] = G[i];
					G_older[i] = G_old[i];
				}
				break;
				case 1:
				G[i] = max( lambda * (mag_min[i]/mag[i]) ,  ( 1 - (mag_min[i]/mag[i]) ) );
				break ;
				case 2:
				G[i] = max( lambda * (p[i]/mag[i]) ,  ( 1 - (mag_min[i]/mag[i]) ) );
				break;
				case 3:
				G[i] = max( lambda * (mag_min[i]/p[i]) ,  ( 1 - (mag_min[i]/p[i]) ) );
				break;
				case 4:
				G[i] = max( lambda ,  ( 1 - (mag_min[i]/p[i]) ) );
				break;
				case 6:
				snr = (mag[i]/mag_min[i]);
				//G[i] = max( lambda_6 ,  ( 1 - alpha_6* pow(mag_min[i]/mag[i],2)) );
				G[i] = max( lambda_6 ,  ( 1 - (alpha_6*(1/snr)*mag_min[i]) ) );
				break;
				default:
				G[i] = max( lambda,  ( 1 - (mag_min[i]/mag[i]) ) );
			}
			//free(p);
			c[i] = rmul( G[i], c[i]);
			//free(G);
		}
	}
	
	ifft(FFTLEN, c);
	
	//copy into outframe
	for(i=0; i<FFTLEN; i++)
	{
		outframe[i] = c[i].r;
	}
	
	free(c);

	
	
	/********************************************************************************/
	
    /* multiply outframe by output window and overlap-add into output buffer */  
                           
	m=io_ptr0;
    
    for (k=0;k<(FFTLEN-FRAMEINC);k++) 
	{    										/* this loop adds into outbuffer */                       
	  	outbuffer[m] = outbuffer[m]+outframe[k]*outwin[k];   
		if (++m >= CIRCBUF) m=0; /* wrap if required */
	}         
    for (;k<FFTLEN;k++) 
	{                           
		outbuffer[m] = outframe[k]*outwin[k];   /* this loop over-writes outbuffer */        
	    m++;
	}	                                   
}        
/*************************** INTERRUPT SERVICE ROUTINE  *****************************/

// Map this to the appropriate interrupt in the CDB file
   
void ISR_AIC(void)
{       
	short sample;
	/* Read and write the ADC and DAC using inbuffer and outbuffer */
	
	sample = mono_read_16Bit();
	inbuffer[io_ptr] = ((float)sample)*ingain;
		/* write new output data */
	mono_write_16Bit((int)(outbuffer[io_ptr]*outgain)); 
	
	/* update io_ptr and check for buffer wraparound */    
	
	if (++io_ptr >= CIRCBUF) io_ptr=0;
}

/************************************************************************************/
void lpf_1 (float *x, float *pbuf)
{	
	int j;
	double lpf_k = exp( -(TFRAME/tau) );
	
	for  (j=0; j<FFTLEN;j++)
	{
		pbuf[j] = ( (1-lpf_k) * x[j] ) + ( lpf_k * pbuf[j] );
	}
}
/*******************************************************************************************/
