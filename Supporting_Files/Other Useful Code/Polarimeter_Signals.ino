#include <ADC.h>

#define ledPin 13
#define ltfPin 23

int current_time = 0;
int old_time = 0;
int gap = 0;
int looptime = 600;
bool wait = true;
bool first = true;
volatile int count = 0;
int grabcount;
float corr;
float lastcorr;
float alpha = 0.01;


ADC *adc = new ADC();

void irq1(){
  ++count;
}

void setup() {
  attachInterrupt(digitalPinToInterrupt(ltfPin), irq1, RISING); 
  adc->setResolution(12,0);
  adc->setResolution(12,1);
  //sampling speed
  adc->setSamplingSpeed(ADC_SAMPLING_SPEED::VERY_LOW_SPEED, 0);
  adc->setSamplingSpeed(ADC_SAMPLING_SPEED::VERY_LOW_SPEED, 1);
  pinMode(ledPin, OUTPUT); 
  // serial port for monitoring
  Serial.begin(9600);
}

void loop() {

  while(wait == true){
    current_time = micros();
    gap = current_time - old_time;
    if (gap > looptime){
      old_time = current_time;      
      wait = false;
      break;
    }
  }
  
  ADC::Sync_result result = adc->analogSyncRead(A1, A2);
  grabcount = count;
  count = 0;
  int value0 = result.result_adc0; 
  int value1 = result.result_adc1; 
  if(first == true){
    corr = (float)grabcount;
    lastcorr = corr;    
  }
  else{
    corr = ((float)grabcount * alpha) + ((1.0 - alpha)*lastcorr);
  }
  
  float voltage0 = value0/corr;
  float voltage1 = value1/corr; 

  Serial.print(voltage0);
  Serial.print(" ");
  Serial.println(voltage1);

  wait = true; 
  first = false; 
  
}
