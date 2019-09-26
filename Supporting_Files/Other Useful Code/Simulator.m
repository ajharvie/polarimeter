%generate signals
timestep = 0.0006;
freqspace = (0:1023)*(1/(timestep*1024));
sigFreq = 16.504;

noise = 0.05; %amount of noise. Half this value to find approx. percentage

phase = pi/3;

time = timestep*linspace(0,999999,1000000);
timeav = 64*time;

maxpos = (timestep*1024)*sigFreq;
maxbin = round(maxpos);
delta = maxpos - maxbin;



%digitised
signal1 = round(2048*(0.9*sin(2*pi*sigFreq*time) + noise*0.9*randn(1,1000000)/200)+ 2060);
signal2 = round(2048*(0.9*sin(2*pi*sigFreq*time + phase) + noise*0.9*randn(1,1000000)/200) + 2060);

%not digitised
% signal1 = 2048*(0.9*sin(2*pi*sigFreq*time) + 0*0.9*randn(1,300000)/200)+ 2060;
% signal2 = 2048*(0.9*sin(2*pi*sigFreq*time + phase) + 0*0.9*randn(1,300000)/200) + 2060;

%with extra component
% signal1 = round(2048*(0.9*sin(2*pi*sigFreq*time) + 0.09*sin(4*pi*sigFreq*time) + noise*0.9*randn(1,1000000)/200)+ 2060);
% signal2 = round(2048*(0.9*sin(2*pi*sigFreq*time + phase) + 0.09*sin(4*pi*sigFreq*time + 2*phase) + noise*0.9*randn(1,1000000)/200) + 2060);

%plot signals
subplot(3,1,1)
plot(time(1:1001), signal1(1000:2000),time(1:1001), signal2(1000:2000))


%sdft stuff
gaps = zeros(1, length(signal1));
%extremely horrible way of initialising complex array
a = zeros(1,1024);
twiddle = complex(a,0);
maxind1 = 10; %just initialising it

windowsize = 1024;
windowed1 = zeros(windowsize,1);
windowed2 = zeros(windowsize,1);

%circular buffers for time domain values
live_samples1 = zeros(windowsize,1);
live_samples2 = zeros(windowsize,1); 

dft1 = zeros(windowsize,1);
dft2 = zeros(windowsize,1);
sample_index = 1;
r = 0.99999999999999; %damping factor
rton = power(r,windowsize);

%compute twiddle factors
for k = 1:windowsize
    factor = (2*pi*(k-1))/windowsize;
    twiddle(k) = exp(1i*factor);
end
for k = 1:(length(signal1) - windowsize)
    
    %dft1 time domain values
    old_x = live_samples1(sample_index);
    live_samples1(sample_index) = signal1(k);
    
    % and for dft2
    old_y = live_samples2(sample_index);
    live_samples2(sample_index) = signal2(k);
    
    if k < 1025
    %update
        for i = 1:windowsize
            dft1(i) = twiddle(i) * (r * dft1(i) - rton*old_x + signal1(k));
            dft2(i) = twiddle(i) * (r * dft2(i) - rton*old_y + signal2(k));
        end

        %hanning window
        windowed1(1) = 0.5*dft1(1) - 0.25*(dft1(2) + dft1(windowsize));
        windowed2(1) = 0.5*dft2(1) - 0.25*(dft2(2) + dft2(windowsize));
        
        for i = 2:1023
            windowed1(i) = 0.5*dft1(i) - 0.25*(dft1(i-1) + dft1(i+1));
            windowed2(i) = 0.5*dft2(i) - 0.25*(dft2(i-1) + dft2(i+1));
        end
        
        windowed1(1024) = 0.5*dft1(1024) - 0.25*(dft1(1) + dft1(1023));
        windowed2(1024) = 0.5*dft2(1024) - 0.25*(dft2(1) + dft2(1023));
    end
    
    
    
    if k >=1025
        %update with narrower windowing
       for i = maxbin-2:maxbin+4  %it's offset because of matlab arrays starting at 1 (humph)
            dft1(i) = twiddle(i) * (r * dft1(i) - rton*old_x + signal1(k));
            dft2(i) = twiddle(i) * (r * dft2(i) - rton*old_y + signal2(k));
       end 
       %window
       
       for i = maxbin-2:maxbin+4
            windowed1(i) = 0.5*dft1(i) - 0.25*(dft1(i-1) + dft1(i+1));
            windowed2(i) = 0.5*dft2(i) - 0.25*(dft2(i-1) + dft2(i+1));
       end
        
       
    end
    
    

    
    %next element in circular buffer
    sample_index = sample_index + 1;
    if sample_index>windowsize
        sample_index=1;
    end    
    
    
    if k>windowsize
        
        if delta>0
            gapCentre = angle(windowed2(maxbin+1)) - angle(windowed1(maxbin+1));
            gapRight = angle(windowed2(maxbin+2)) - angle(windowed1(maxbin+2));
            
            if gapCentre < 0
                gapCentre = gapCentre + 2*pi;
            end
            if gapRight < 0
                gapRight = gapRight + 2*pi;
            end

            gaps(k) = gapCentre;% + delta*(gapRight - gapCentre); 
        end
        
        if delta<=0
            gapCentre = angle(windowed2(maxbin+1)) - angle(windowed1(maxbin+2));
            gapLeft = angle(windowed2(maxbin)) - angle(windowed1(maxbin));
            
            if gapCentre < 0
                gapCentre = gapCentre + 2*pi;
            end
            if gapLeft < 0
                gapLeft = gapLeft + 2*pi;
            end
            
            
            gaps(k) = gapCentre;% - delta*(gapLeft - gapCentre);
        end
    
        
        
        
        mag1 = abs(windowed1);
        mag2 = abs(windowed2);
        allph1 = angle(windowed1(1:windowsize-1, 1)); 
        allph2 = angle(windowed2(1:windowsize-1, 1));

        

        gaps(k) = gaps(k)*(180/pi);

        
        if gaps(k)<0
            gaps(k) = gaps(k) + 360;
        end

      
    end
   %%crunch
end
allph1 = wrapTo2Pi(allph1);
allph2 = wrapTo2Pi(allph2);
allgap = wrapTo2Pi(allph2-allph1);
for i=1:length(mag1)
    
    if mag1(i) < 20
        allgap(i) = 0;
        
    end
end


outCount=1;
avCount = 1;
inVals = gaps(1025:998975);


for i=1:length(inVals)
    
   box(avCount) = inVals(i);
   
   if avCount == 64
       output(outCount) = mean(box);
       outCount = outCount + 1;
       avCount = 1;
   end
   
   avCount = avCount + 1;
end



std(output)

%plot output
subplot(3,1,3)
plot(timeav(1:length(output)), output)

%plot spectra
subplot(3,1,2)
yyaxis left
plot(freqspace(1:30),mag1(1:30),freqspace(1:30),mag2(1:30))
yyaxis right
plot(freqspace(1:30), allgap(1:30))

%mean(gaps(2000:length(gaps)-2000))




