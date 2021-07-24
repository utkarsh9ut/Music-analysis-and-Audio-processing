function helper_plot(y,Fs)
pausetym = .1;
num_sample = ceil(Fs*pausetym);       %used ceil function for obtaining the next greatest intger incase of decimals, sampling frequency into pausetym would give us the number of samples in that amount of time                              

figure('Name','Group2_EAD'), hold on, xlim([0 length(y)]), ylim([.9*min(y) 1.1*max(y)])
                                xlabel("SAMPLES(176000samples=1second)-->"); 
                                ylabel("NORMALIZED DATA (between -1.0 and 1.0)-->"); 
                                title('Window Sample')
sound(y,Fs)
    for i = 1:num_sample:length(y)-num_sample      %plotting the sound by plotting and correspondingly pausing for the amount of time which corresponds to the number of samples in one pause time
        plot(i:i+num_sample,y(i:i+num_sample),'b')
        pause(pausetym)                             %waiting for a constant each time after plotting it's samples in one go, so that the time and samples plotted are in a proper sync
    end