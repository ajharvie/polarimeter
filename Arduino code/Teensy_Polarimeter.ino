/* 
Software for Teensy-powered polarimeter. Performs sliding DFT on periodic intensities 
from two OPT101s, with reference for laser intensity on a light-frequency converter
(TSL235R-LF). After "initialisation" step, where sDFT bins are filled, peaks are located 
to find the shared fundamental frequency of the periodic  signals. The sDFT is then 
continuously updated at the peak location only, and a value of the phase difference 
between the two signals is calculated at the peak and output over serial connection.

Uses single-precision floats to take advantage of Teensy 3.6's FPU.

Sliding DFT libary adapted slightly from the implementation by Bronson Philippa 
https://github.com/bronsonp/SlidingDFT

ADC library for Teensy by Pedro Villanueva, pedvide https://github.com/pedvide/ADC 
*/


#include "sliding_dft_narrow.h"
#include <math.h>
#include <ADC.h>
#define pi 3.14159265358979

//Choose sample rate in microseconds (600)
#define looptime 600

//Choose update rate (number of samples between each output)
#define up_rate 64 

// Number of initialisation samples
#define number_loops 1025

//LED pin
#define ledPin 13

//Pin for light-frequency converter input
#define ltfPin 23

//Pins for OPT101s 
#define pin1 A1
#define pin2 A2

//Smoothing parameter for L-F converter.
float alpha = 0.01;

//Smoothing parameter for output filter (multiplier = 1 for approx 6s time constant)
float lambda = 1*0.0128;

float min_lambda = 1-lambda;

//Variables///////////////////////////
int i;
int up_index = 0;
int peak_location = 0; 
int initialiser0[number_loops];
int initialiser1[number_loops];
int initialiser_refs[number_loops];
unsigned long current_time = 0;
unsigned long old_time = 0;
unsigned long gap;
float oldref = 0.0;
float ref_value; 
float phase1;      
float phase2;
float lag;
float re1;
float re2;
float im1;
float im2;
bool wait = true;
bool first = true;
float checkmag;
float maxmag = 0.0;
volatile int count = 0;
int lastcount = 0;
int grabcount = 0;
int avstart = 0;
float degconv;
float peakLeft;
float peakCentre;
float peakRight;
float gapLeft;
float gapCentre;
float gapRight;
float delta;
int stabilised = 0;
float filtered;


// wow look at all these g l o b a l s //
/////////////////////////////////////////

//ADC object
ADC *adc = new ADC();

// Intitialise 2 DFTs in parallel. 
static SlidingDFT<float, 1024> dft1;
static SlidingDFT<float, 1024> dft2;

//ISR for Light-Frequency converter
void irq1(){
  ++count;
}

//Function to fill up sDFT
void initialise(){
  
  digitalWrite(ledPin, HIGH); //LED on for setup phase
  
  // For "initialisation" phase - loop fills arrays
  for (i=0; i < (number_loops); ++i){
    wait = true;

    //timing loop. Don't use delay() as it doesn't always play nicely with interrupts
    while(wait == true){
      current_time = micros();
      gap = current_time - old_time;
      if (gap > looptime){
        old_time = current_time;      
        wait = false;
        break;
      }
    }
    
    //Get values from ADCs and l-f converter
    ADC::Sync_result result = adc->analogSyncRead(pin1, pin2);
    grabcount = count;
    count = 0;
    
    initialiser_refs[i] = grabcount;
    initialiser0[i] = result.result_adc0;
    initialiser1[i] = result.result_adc1;
  }

  //apply corrections and update sDFTs
  for (i=0; i < (number_loops); ++i){
    float voltage0 = (float)initialiser0[i] / (float)initialiser_refs[i];               
    float voltage1 = (float)initialiser1[i] / (float)initialiser_refs[i];
    dft1.update(voltage0);
    dft2.update(voltage1);      
  }
  
  //Find max. Starts at bin 2 to avoid spurious signal at 0 frequency              
  for (i=2; i<number_loops/2; ++i){ 
    std::complex<float> DC_bin = dft1.dft[i];
    re1 = DC_bin.real();
    im1 = DC_bin.imag();
    checkmag = sqrt((re1*re1)+(im1*im1));
    
    if (checkmag > maxmag){
      peak_location = i;
      maxmag = checkmag;
    }
  }

  //interpolating for true freq - 3 bins around peak
  std::complex<float> leftBin = dft1.dft[peak_location-1];
  std::complex<float> centreBin = dft1.dft[peak_location];
  std::complex<float> rightBin = dft1.dft[peak_location+1];
  //magnitudes
  peakLeft = abs(leftBin);
  peakCentre = abs(centreBin);
  peakRight = abs(rightBin);
  //using Candan's method
  delta = 2*((peakRight-peakLeft)/(2*peakCentre+peakLeft+peakRight));
  
  digitalWrite(ledPin, LOW); //LED off after initialisation step
}

void setup() {
  pinMode(ledPin, OUTPUT);
  pinMode(ltfPin, INPUT);
  attachInterrupt(digitalPinToInterrupt(ltfPin), irq1, RISING);
  adc->setResolution(12,0);
  adc->setResolution(12,1);
  adc->setSamplingSpeed(ADC_SAMPLING_SPEED::VERY_LOW_SPEED, 0);
  adc->setSamplingSpeed(ADC_SAMPLING_SPEED::VERY_LOW_SPEED, 1);
  Serial.begin(9600); 
  delay(50);

  //run function that fills sDFTs and finds peak
  initialise();
  
}

void loop() {

  //"Output" phase
     
  wait=true;
  
  //timing circuitry    
  while(wait==true){     
    current_time = micros();
    gap = current_time - old_time;
    if (gap > looptime){
      old_time = current_time;
      wait = false;
      break;
    }
  }
  
  //Get values from ADCs and l-f converter
  ADC::Sync_result result = adc->analogSyncRead(pin1, pin2);
  grabcount = count;
  count = 0;

  //Exponential smoothing on light-frequency converter reading
  if (first==true){
    ref_value = (float)grabcount;
  }
  else{
    ref_value = ((float)grabcount * alpha) + ((1.0 - alpha)*oldref);
  }
  
  int value0 = result.result_adc0; 
  int value1 = result.result_adc1; 
      
   //Apply correction
  float voltage0 = (float)value0 / ref_value;             
  float voltage1 = (float)value1 / ref_value;  

  //Update sDFT around peak location only
  dft1.updatepeak(voltage0, peak_location);
  dft2.updatepeak(voltage1, peak_location);
  
  //Extracting phase data at peak location
  if(delta > 0.0){ //max to rhs of peak
    gapCentre = arg(dft2.dft[peak_location])-arg(dft1.dft[peak_location]);
    gapRight = arg(dft2.dft[peak_location + 1])-arg(dft1.dft[peak_location + 1]);
    
    if (gapCentre < 0.0){
      gapCentre = gapCentre + 2*pi;
    }
    if (gapRight < 0.0){
      gapRight = gapRight + 2*pi;
    }
    
    lag = gapCentre + delta*(gapRight - gapCentre); //linear interpolation of phase

 }
  else{ //to LHS of peak
    gapCentre = arg(dft2.dft[peak_location])-arg(dft1.dft[peak_location]);
    gapLeft = arg(dft2.dft[peak_location - 1])-arg(dft1.dft[peak_location - 1]);

    if (gapCentre < 0.0){
      gapCentre = gapCentre + 2*pi;
    }
    if (gapLeft < 0.0){
      gapLeft = gapLeft + 2*pi;
    }
  
    lag = gapCentre - delta*(gapLeft - gapCentre);

 }

  oldref = ref_value;
  first = false;

  lag = lag*(180.0/pi); //convert to degrees
  
  if (lag < 0){ 
    lag = lag + 360.0; 
  }

  //don't output first 1024 samples
  if (stabilised < 1024){
    up_index = 0;
    ++stabilised;
    filtered = lag;
  }

  
  ++up_index;
  if (up_index == up_rate){
    
    filtered = filtered*min_lambda + lambda*lag;
    Serial.println(filtered,5);

    up_index = 0;
  }
   
}
