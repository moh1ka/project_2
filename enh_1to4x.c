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
#define NFREQ (1+FFTLEN/2)			/* number of frequency bins from a real FFT */	//129
#define OVERSAMP 4					/* oversampling ratio (2 or 4) */  
#define FRAMEINC (FFTLEN/OVERSAMP)	/* Frame increment */		//64
#define CIRCBUF (FFTLEN+FRAMEINC)	/* length of I/O buffers */		//320

#define OUTGAIN 16000.0				/* Output gain for DAC */
#define INGAIN  (1.0/16000.0)		/* Input gain for ADC  */
// PI defined here for use in your code 
#define PI 3.141592653589793
#define TFRAME FRAMEINC/FSAMP       /* time between calculation of each frame */  //8ms

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
float *inframe, *outframe;          /* Input and output frames */
float *inwin, *outwin;              /* Input and output windows */
float ingain, outgain;				/* ADC and DAC gains */ 
float cpufrac; 						/* Fraction of CPU time used */
volatile int io_ptr=0;              /* Input/ouput pointer for circular buffers */
volatile int frame_ptr=0;           /* Frame pointer */
float *intermediate;
float *mag;
float *copy_x;
float *mag_min;
float *M1; 
float *M2; 
float *M3; 
float *M4; 
float *G; 
float *P;

float lambda = 0.1; 
float alpha = 2.0;  // alpha will probably be around 2 to 4
int counter = 0; 
float min(float a, float b){
	if(b<a)a=b;
	return a;
}
float max(float a, float b){
	if(b>a)a=b;
	return a;
}

float tau1 = 0.02;
float tau2 = 0.08;
int enhance_1 = 1;	//0 off; 1 on
int enhance_3 = 0;	//0 off; 1 on
int enhance_4 = 0;	//0 off; 1 on
int enhance_6 = 0;
int G_index = 1;	//1 to 4

float alpha_6 = 4.0; 
 /******************************* Function prototypes *******************************/
void init_hardware(void);    	/* Initialize codec */ 
void init_HWI(void);            /* Initialize hardware interrupts */
void ISR_AIC(void);             /* Interrupt service routine for codec */
void process_frame(void);       /* Frame processing routine */
  
void lpf_1(float *data, float t);     
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
	intermediate = (float *) calloc(FFTLEN, sizeof(float)); /* Array for processing*/
	
	mag	     	= (float *) calloc(FFTLEN, sizeof(float)); /* Array for processing*/
	mag_min	    = (float *) calloc(FFTLEN, sizeof(float)); 
	copy_x	= (float *) calloc(FFTLEN, sizeof(float));
	P = (float *) calloc(FFTLEN, sizeof(float));
		//Min noise buffers
	M1	        = (float *) calloc(FFTLEN, sizeof(float)); /* magnitude spectrum*/
	M2	        = (float *) calloc(FFTLEN, sizeof(float)); /* magnitude spectrum*/
	M3	        = (float *) calloc(FFTLEN, sizeof(float)); /* magnitude spectrum*/
	M4	        = (float *) calloc(FFTLEN, sizeof(float)); /* magnitude spectrum*/
	G	        = (float *) calloc(FFTLEN, sizeof(float)); /* magnitude spectrum*/
	
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
	complex *c;
	c = (complex *) calloc(FFTLEN, sizeof(complex));
	/* work out fraction of available CPU time used by algorithm */    
	cpufrac = ((float) (io_ptr & (FRAMEINC - 1)))/FRAMEINC;  
	
	/* wait until io_ptr is at the start of the current frame */ 	
	while((io_ptr/FRAMEINC) != frame_ptr); 
	
	/* then increment the framecount (wrapping if required) */ 
	if (++frame_ptr >= (CIRCBUF/FRAMEINC)) frame_ptr=0;		//if >= 5, ->0.
 	
 	/* save a pointer to the position in the I/O buffers (inbuffer/outbuffer) where the 
 	data should be read (inbuffer) and saved (outbuffer) for the purpose of processing */
 	io_ptr0 = frame_ptr * FRAMEINC;
	
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
							      	
	/*inframe processing*/
	for(i=0; i<FFTLEN; i++){
		c[i].r = inframe[i];
	}
	fft(FFTLEN,c);
	
	/*deal with noise here*/
	counter++;
	
	for(i=0; i<FFTLEN; i++){
		mag[i] = cabs(c[i]);
		copy_x[i] = cabs(c[i]);
	}
	
	if(enhance_1 == 1)lpf_1(mag,tau1);
	if(enhance_4 == 1)lpf_1(mag,tau2);

	memcpy (M1, mag, FFTLEN*sizeof(float));
	
	for(i=0; i<FFTLEN; i++){
		M1[i] = min(M1[i],mag[i]);
	}
	if(counter == 78){
		counter = 0;
		memcpy (M4, M3, FFTLEN*sizeof(float));
		memcpy (M3, M2, FFTLEN*sizeof(float));
		memcpy (M2, M1, FFTLEN*sizeof(float));
		memcpy (M1, mag, FFTLEN*sizeof(float));
	}
	
	for(i=0; i<FFTLEN; i++){
		mag_min[i] =  min(M1[i], M2[i]);
		mag_min[i] =  min(M3[i], mag_min[i]);
		mag_min[i] = ( alpha * (min(M4[i], mag_min[i])) );
		
		if(enhance_3 == 1)lpf_1(mag_min, tau1);
		if(enhance_4 == 1){
			switch(G_index){
				case 1:
				G[i] = max( lambda*(mag_min[i]/copy_x[i]), ( 1 - (mag_min[i]/copy_x[i])));
				break;
				case 2:
				G[i] = max( lambda*(mag[i]/copy_x[i]),  ( 1 - (mag_min[i]/copy_x[i]) ) );
				break;
				case 3:
				G[i] = max( lambda*(mag_min[i]/mag[i]),  ( 1 - (mag_min[i]/mag[i]) ) );
				break;
				case 4:
				G[i] = max( lambda,  ( 1 - (mag_min[i]/mag[i]) ) );
				break;
			}
		}
		if(enhance_6 == 1){
			for(i=0;i<20;i++){
			mag_min[i]= alpha_6* (mag_min[i]/alpha);
			}
			G[i] = max( lambda,  ( 1 - (mag_min[i]/mag[i])*(mag_min[i]/mag[i])) );
		}
		else{
			G[i] = max( lambda, ( 1 - (mag_min[i]/copy_x[i])));
		}

		c[i] = rmul( G[i], c[i]);
	}
	
	ifft(FFTLEN,c);
	for(i=0; i<FFTLEN; i++){
		outframe[i] = c[i].r;
	}
	free(c);	
    for (k=0;k<FFTLEN;k++)
	{                           
		outframe[k] = inframe[k];/* copy input straight into output */ 
	} 
	
	/********************************************************************************/
	
    /* multiply outframe by output window and overlap-add into output buffer */  
                           
	m=io_ptr0;
    
    for (k=0;k<(FFTLEN-FRAMEINC);k++) //k<(256-64)
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
void ISR_AIC(void){       
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

void lpf_1(float *data, float t){
	
	int i;
	double k = exp(-(TFRAME/t));
	 
	for(i=0; i<FFTLEN; i++){
		P[i] = (1 - k)*data[i] + k*P[i];
	}
	data [i] = P[i];
}
