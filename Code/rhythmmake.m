%%PROJECT PART-3 MAKING A RHYTHM FROM GIVEN FREQUENCY NOTATIONS
function rhythmmake()
%Fundamental Frequency
fs = 44100;

% rest
zero = rest(4,fs);     % 1.0 sec
zeroh = rest(8,fs);    % 0.5 sec
zerohh = rest(16,fs);  % 0.25 sec

% eighth note as 0.5 sec
C = key(52,8,fs);  %C 1
D = key(54,8,fs);  %D 2
E = key(56,8,fs);  %E 3
F = key(57,8,fs);  %F 4
G = key(59,8,fs);  %G 5
A = key(61,8,fs);  %A 6
B = key(63,8,fs);  %B 7

% quarter note as 1 sec
C_4 = key(52,4,fs);  %C 1
D_4 = key(54,4,fs);  %D 2
E_4 = key(56,4,fs);  %E 3
F_4 = key(57,4,fs);  %F 4
G_4 = key(59,4,fs);  %G 5
A_4 = key(61,4,fs);  %A 6
B_4 = key(63,4,fs);  %B 7

% make a song
line_1 = [ E E E_4 E E E_4];
line_2 = [ E G C D E_4];
line_3 = [ F F F_4 F F E E_4];
line_4 = [ E E D D E D_4 G_4];
line_5 = [ E E E_4 E E E_4 E G];
line_6 = [ C D_4 E F F F F F E_4];
line_7 = [ E E G G F D C_4];

song = [ line_1 line_2 line_3 line_4 line_5 line_6 line_7];
plot(song)
title('plot of rhythm')
sound(song,fs,24);
% audiowrite('star.wav',song,fs,'BitsPerSample',32);

function wave = key(p, n, fs)
    t = 0:1/fs:4/n;
    idx = 440*2^((p-49)/12);
    
%     method 1 - orginal
%     wave = (sin(2*pi*idx*t));

%     method 2 - exponential decreasing 
     tt = 4/n:-1/fs:0;
     wave = (sin(2*pi*idx*t)).*exp(tt);
     wave = wave./max(wave);
    
%     method 3 - triangle decreasing 
%     mid = (t(1)+t(end))/2;
%     tri = -(abs(t-mid)-mid);
%     tri = tri./max(tri);
%     wave = (sin(2*pi*idx*t)).*tri;
    
%     plot(wave)

function wave = rest(n,fs)
    t = 0:1/fs:4/n;
    tt = 4/n:-1/fs:0;
    wave = 0*sin(2*pi*t).*exp(tt);