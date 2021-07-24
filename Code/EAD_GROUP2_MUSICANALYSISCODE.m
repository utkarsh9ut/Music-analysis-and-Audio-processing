%% BASIC SOUND EFFECTS GROUP2 EAD  PART-1 OF PROJECT
%PLEASE TRAVERSE THROUGH THE DIFFERENT SECTIONS AND RUN THE PROGRAM
%SECION-WISE
%DONOT RUN THE ENTIRE PROGRAM AT ONCE
%% FILTERING EFFECTS (USE THE SONG:- SHAPE OF YOU)
%program to implement a LPF with cut off frequency 8khz to denoise audio signal
%it is assumed that original signal is 'x'


Sf = 44100; %Sampling frequency
Pf = 8e3; %passband frequency
Stf = 8.4e3; %Stopband frequency
Ap = 1; %passband ripple in dB
Ast = 95; %stopband attenuation in dB

df = designfilt('lowpassfir','PassbandFrequency', Pf, 'StopbandFrequency',Stf,'PassbandRipple',Ap,'StopbandAttenuation',Ast,'SampleRate',Sf);

fvtool(df);

[locate,pathname] = uigetfile('*.*','Select the Input Audio');
[x,Fs] = audioread(num2str(locate));
xn = awgn(x,15,'measured'); %signal corrupted

y = filter(df,xn);

%plotting signals
figure('Name','Group2_EAD'),
subplot(3,1,1)
plot(x);
title('Original Signal');
subplot(3,1,2)
plot(xn);
title('Noisy Signal');
subplot(3,1,3)
plot(y);
title('Filtered Signal');
%% Graphic Equalizer (USE THE SONG:- SHAPE OF YOU)
%program to implement real time Graphic Equaliser
clear all;
close all;
clc;

deviceReader = dsp.AudioFileReader("Filename","shapeofyou.mpeg");

deviceWriter = audioDeviceWriter('SampleRate',deviceReader.SampleRate);

equalizer = graphicEQ('Bandwidth', '1 octave', 'Structure', 'Parallel', 'SampleRate', deviceReader.SampleRate);

%Simulating effects i.e. pop, rock etc.
equalizer.Gains = [10 10 10 10 10 0 0 0 0 0];

visualize(equalizer)
    
nUnderruns = 0;

tic;

while toc < 30 %30 seconds of simulation
    in = deviceReader();
    out = equalizer(in);
    nUnderruns = nUnderruns + deviceWriter(out);
end

%clean up
release(deviceReader);
release(deviceWriter);
%% echo effect (USE THE SONG:- SHAPE OF YOU)
[locate,pathname]=uigetfile('*.*','Select the Input Audio');   %to manually select the audio from computer system

[x,Fs]= audioread(num2str(locate));     %reading the audio file and storing the samples in vector 'x' and sampling frequency in 'Fs
sample_num=length(x);                   
delaypathgain=0.7;
delaysamples=1800;


xn=padarray(x,delaysamples,0,'pre');           %padding the memory array with 0's equal to the number of delaysamples in the beginning.
audioout=zeros(sample_num+delaysamples,1);         %creating an audioout array of appropriate size


for i=(delaysamples+1):1:sample_num                 %loop to iterate through samples and storing them in the audioout
    audioout(i-delaysamples,1)=x(i)+delaypathgain*xn(i-delaysamples);
end
  

%% need for synchronization of rythm (EXTENSION AND IMPROVEMENT OF ECHO EFFECT)(USE THE SONG:- SHAPE OF YOU)
clear all;
close all;
clc;
[locate,pathname]=uigetfile('*.*','Select the Input Audio');
[x,Fs]= audioread(num2str(locate));

bps=1.6;
secondperbeat=1/bps;
delaysampled = fix(secondperbeat*Fs);

out=zeros(size(x));

g=0.5;
for n=1:length(x)
    
    if(n-delaysampled)<1
        out(n,1) = x(n,1)+ 0;
    else
        out(n,1)=x(n,1)+g*x(n-delaysampled,1);
    end 
end 
%% flang (USE FLANG48Hz audio in the folder)
[locate,pathname]=uigetfile('*.*','select the input audio');
[x,fs]=audioread(num2str(locate));
len=length(x); %get the length of music file
ts=len/fs;
fc=0.5;       % frequency of sine wave
ti=linspace(0,ts,len);
delay=100;        %delay factor
a=0.75;                      %attentuation factor

sinewave=sin(2*pi*fc*ti);
sinedelay=round(delay.*sinewave')+delay; %variable delay
y=zeros(len+delay,1);              %initialize output music signal 
xn=padarray(x,[delay,0],0,'pre'); %padding the array of input samples
for i=(delay+1):1:len
    y(i-delay,1)=x(i)+a*xn(i-sinedelay(i-delay)); %producing difference equations
 end;
%% reverb (USE THE SONG:- SHAPE OF YOU)
clear all;
close all;
clc;

%Program to implement Reverb effect

%Reading audio input from an audio file
deviceReader = dsp.AudioFileReader("Filename","shapeofyou.mpeg")

%Creating a deviceWriter to write audio to the output speakers using
%predefined function

deviceWriter = audioDeviceWriter('SampleRate',deviceReader.SampleRate);

%reverb effect is created using predefined function in matlab which take
%inputs : predelay, the delay added with input signal is 0.5 seconds,
%wetdrymix is the ratio of delayed and orignal input being fed to mixer.
reverb_effect = reverberator('PreDelay',(deviceReader.SampleRate)*0.625 , 'WetDryMix',0.2);

underruns = 0;

tic;

%30 seconds simulation using while loop
while toc<30 
    %taking input from the audio file read 
    input = deviceReader();
    
    %Proccesing and applying reverb effect to input audio using
    %reverb_effect which uses predefined reverberator function.
    output = reverb_effect(input);
    
    %underrun buffer is updated so that program does not keep repeating the 
    %sound or not playing sound at all
    underruns = underruns + deviceWriter(output);
end

%Clean up of variables
release(reverb_effect);
release(deviceReader);
release(deviceWriter);
%% PROJECT PART-2 DECODING THE NOTES AND FREQUENCY DETECTION
%% 
clear;
clc; 
clf; 
close all

[song,sampling_freq] = audioread('songfile.mp3');
sampling_freq= sampling_freq*4;   % fastening the audio since the original file is relatively slow
figure('Name','Group2_EAD'), plot(song(:,1)), xlabel("SAMPLES(44000samples=1second)-->"); ylabel("NORMALIZED DATA (between -1.0 and 1.0)-->");title('SONG, entire song')

tstart = 2900000;        % setting the basic parameters for our song
tstop = 4900000;     


% beginning with analysis of a certain window
w = song(tstart:tstop);   %storing the samples in 'w' array
[~,n] = size(w);          
t = linspace(tstart,tstop,n);  %creating n equally spaced samples between tstart and tstop


helper_plot(w,sampling_freq);  %plotting the window of song using helper function

audiowrite('songfile_window.wav',w,sampling_freq);  %storing the window in matlab created file "songfile_window.wav"
%% reducing the number of samples
clc

down_slide = 20;                                   %downsample factor, we will reduce the number of samples by a factor of 20
sample_mod = floor(n/down_slide);
sampling_mod = round(sampling_freq/down_slide);    %modifying the sampling frequency accordingly

average_out = zeros(1,sample_mod);  %intializing a vector with a length of modified number of samples

for k = 1:sample_mod
    average_out(k) = mean(w(down_slide*(k-1)+1:down_slide*k));
end
figure('Name','Group2_EAD'), subplot(3,1,1)
                              plot(linspace(0,100,n),abs(w))
                              title('WAVEFORM BEFORE DOWN-SAMPLING')
                              subplot(3,1,2)
                              plot(linspace(0,100,sample_mod),abs(average_out))
                              title('WAVEFORM AFTER DOWN-SAMPLING')
                              subplot(3,1,3)
                             plot(linspace(0,100,n),abs(w)), hold on
                             plot(linspace(0,100,sample_mod),abs(average_out))
                             xlabel("SAMPLES (X-AXIS ADJUSTED ACCORDING TO THE NUMBER OF SAMPLES-->"); ylabel("NORMALIZED DATA -->");
                             title('COMPARATIVE ANALYSIS (DISCRETE NOTES OF SONG)')
legend('INTIAL WINDOW', 'AFTER DOWN SAMPLING AND AVERAGING')
 sound(average_out,sampling_mod);
 %% Finding the prominent frequency points of the song
 close all
k = 1;
thresh_out = zeros(1,sample_mod);

while (k <= sample_mod)
    thresh_comp = 5*median(abs(average_out(max(1,k-5000):k)));    %defining a threshold to find the chunks of song with prominent amplitude
   
    if (abs(average_out(k)) > thresh_comp)                        % if the average value is greater than threshold than we consider the following 500 points
        for j = 0:500                                             %we collect those average values which are greater than the defined threshold
            if (k + j <= sample_mod)
                thresh_out(k) = average_out(k);
                k = k + 1;
            end
        end
        k = k + 1400;              %once we are done with a certain chunk of prominent points we jump 1400 samples as they we would be useless for us and would correspond to the same note whose sampling has already been done
    end
    k=k+1;
end

figure('Name','Group2_EAD'), subplot(2,1,1), plot(abs(average_out)), title('INITIAL SONG WINDOW'), ylim([0 1.1*max(average_out)])
        subplot(2,1,2), plot(abs(thresh_out)), xlabel("SAMPLES-->");ylabel("NORMALIZED DATA -->");title('Detection of peaks(frequency points) using moving average threshold')
        
 sound(thresh_out,fix(sampling_mod)); 
 %% find frequencies corresponding to each prominent frequency point
 clc; 
close all

k = 1;note_initial = 0;


while k < sample_mod                 %iterating the loop through the number of samples 
    
    note_last = 0;                   %pivot element to check for end of our prominent frequency point
    
     j = 1;                          %initializing a variable to start counting the note content of detected chunk of prominent frequency point
     
    while (((thresh_out(k) ~= 0) || (note_last > 0)) && (k < sample_mod))   %iterating the loop throug the number of prominent frequency chunks  
        note(j) = thresh_out(k);                                            % storing the non-zero detected frequency points          
        k = k + 1;   j = j + 1;
        if (thresh_out(k) ~= 0)            % if threshold value not equal to zero then make the pivot equal to 1
            note_last = 20;
        else
            note_last = 0;                  % now, as soon as the threshold becomes 0, we need to analyze the frequency content of the collected points
        end
        if (note_last == 0)                 % calculating the frequency of the collected points using fast fourier methods
           if (j > 25)
               resized_note = padarray(note,[0 j],0,'post');    % padding the note with zeros to double size (N --> 2*N-1)
               fft_note = fft(resized_note);            % storing the result obtained after fast fourier transform in a variable
               
               len_note = length(note);
               f = linspace(0,(1+len_note/2),len_note);
               [~,index] = max(abs(fft_note(1:length(f))));
               if (f(index) > 20)
                   note_initial = note_initial + 1;          %iterating the note number
                   note_collec(note_initial) = f(index)*2;   
                   figure('Name','Group2_EAD'), 
                            subplot(2,1,1)
                            plot(f,abs(fft_note(1:length(f))))
                            title(['Fundamental frequency = ',num2str(note_collec(note_initial)),' Hz'])
                            subplot(2,1,2)
                            plot(f,angle(fft_note(1:length(f))))
               end
               k = k + 50;
           end
           clear note;
           break
        end
        
    end
    k = k + 1;
end
%% recreating the music to check the correctness of our decoded notes
fs = 20500;  % approx sampling frequency for recreating the song

duration = .6;  % the duration for which we play a particular note

song_reproduced = zeros(1,duration*fs*length(note_collec));    %intializing a array with size equal to the approximate sample number of our reproduced song

for k = 1:length(note_collec)                                  % loop to iterate through the number of notes collected after the entire decoding process
    [letter(k,1),freq(k)]= freq_helper(note_collec(k));
    values = 0:1/fs:duration;                                   %creating an array from 0 to duration with values at every 1/fs interval 
    a = sin(2*pi*freq(k)*values*2);                             % using sine wave as medium to express the frequency in audio form
    song_reproduced((k-1)*fs*duration+1:k*fs*duration+1) = a;    % adding the sine wave to our reproduced song vector
     sound(a,fs); pause(.6);
end

letter  
audiowrite('fur_elise_recreated.wav',song_reproduced,fs);
