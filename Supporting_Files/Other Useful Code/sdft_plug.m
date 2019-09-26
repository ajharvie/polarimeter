function dft_out = sdft_plug(y0)

    a = zeros(1,1024);
    twiddle = complex(a,0);

    windowsize = 1024;


    %circular buffers for time domain values
    live_samples1 = zeros(windowsize,1);

    dft1 = zeros(windowsize,1);
    sample_index = 1;
    r = 0.99999999999999; %damping factor
    rton = power(r,windowsize);

    %compute twiddle factors
    for k = 1:windowsize
        factor = (2*pi*(k-1))/windowsize;
        twiddle(k) = exp(1i*factor);
    end
    
    for k = 1:length(y0)

        %dft1 time domain values
        old_x = live_samples1(sample_index);
        live_samples1(sample_index) = y0(k);



        %update
        for i = 1:windowsize
            dft1(i) = twiddle(i) * r * (dft1(i) - rton*old_x + y0(k));
        end

        %next element in circular buffer
        sample_index = sample_index + 1;
        if sample_index>windowsize
            sample_index=1;
        end    



    end
       %%crunch

    dft_out=transpose(dft1);
    
end





