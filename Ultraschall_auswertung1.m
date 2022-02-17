%--- dependencies --------------------------------------------------------------

% pkg install -forge control;
% pkg install -forge signal;
pkg load control; % mag2db
pkg load signal;  % hilbert

%--- loading data and global settings ------------------------------------------

load einzelpuls.mat;
set(0, 'DefaultLineLineWidth', 1); % bit thicker lines in the plots ... couldn't properly see them

%--- fig 1: taking a look at the raw data --------------------------------------

figure(1);
plot(S.data);
xlim([0 160]);
ylim([0 250]);

%--- fig 2: scaled to volts vs. microseconds -----------------------------------

V=(S.sample_offset-S.data)./S.sample_resolution*S.gain;
t=(S.starting_address:S.ending_address-1)/S.FS;

figure(2);
plot(t * 1e6, V); % times 1e6 so we get microseconds, like in the PDF
xlabel("TOF [microsecs]");
ylabel("Voltage [V]");
ylim([-0.5 0.5]);

%--- fig 3: normalized-amplitude-scaled-to-dB vs. frequency (from FFT) ---------

volt = V;
f = (0:2048 - 1) * S.FS / 2048;
yfft = fft(volt, 2048);
amplitude = 20*log10( abs(yfft) ); % it seems "20*log10" is all that mag2db does???

% "Normalization can be done in many different ways - depending on window, 
% number of samples, etc. Common trick: take FFT of known signal and normalize
% by the value of the peak. Say in the above example your peak is 123 - if you
% want it to be 1 , then divide it ( and all results obtained with this 
% algorithm) by 123" (https://stackoverflow.com/a/20170743/5354137)
normalized_amplitude = amplitude - max(amplitude); % is that what he means?!

figure(3);
plot(f / 1e6, normalized_amplitude);
title("Amplitude");
xlabel("Frequency [MHz]");
ylabel("Amplitude [dB]");
xlim([0 150]);
ylim([-70 0]);

%--- range for where signal goes above -6dB, and below again -------------------

normalized_amplitude = normalized_amplitude(1:2048/2); % get rid of that sort of mirror-image created from the 
                                                       % FFT (the whole figure 3 is from 0-500MHz, we just 
                                                       % don't see it because we give xlim of 0-150MHz)
x1 = find(normalized_amplitude > -6)(1);   % time point of first rise above 6dB
x2 = find(normalized_amplitude > -6)(end); % time point of first falling below 6dB again

frange6dB = [f(x1) f(x2)]; % Hz, not MHz! (it goes from 26.367MHz to 63.477MHz)

%--- fig 4: phase ("all the angles") vs. frequency (from same FFT) -------------

phase = angle(yfft); % so, it's not in decibels? what is it in, then? nevermind...

figure(4);
plot(f / 1e6, phase);
title("Phase");
xlabel("Frequency [MHz]");
xlim([0 150]);
ylim([-4 4]);

%--- fig 5: hilbert real & hilbert imag & hilbert absolute ---------------------

yhil160 = hilbert(V); % hilbert(xr) returns the analytic signal, x, from a real data sequence, xr. 
                      % If xr is a matrix, then hilbert finds the analytic signal corresponding to each column.

% https://en.wikipedia.org/wiki/Hilbert_transform
% - invented by David Hilbert (1862-1943), who liked to wear a very funny hat!
% - it's a "specific linear operator", he invented it because he was Friends with Bernhard Riemann.
% - takes a function, u(t) of a real variable and produces another function of a real variable H(u)(t).
% - has a particularly simple representation in the frequency domain:
%   It imparts a phase shift of ±90° (π⁄2 radians) to every frequency component of a function,
%   the sign of the shift depending on the sign of the frequency.
% - has the effect of shifting the phase of the negative frequency components of u(t) by +90° (π⁄2 radians)
%   and the phase of the positive frequency components by −90°,
%   and i·H(u)(t) has the effect of restoring the positive frequency components
%   while shifting the negative frequency ones an additional +90°,
%   resulting in their negation (i.e., a multiplication by −1).

microsecs = t * 1e6;

figure(5);
plot5 = plot(microsecs,real(yhil160), microsecs,imag(yhil160), microsecs,abs(yhil160));
set(plot5(1), 'Color', 'blue');
set(plot5(2), 'Color', 'red');
set(plot5(3), 'Color', 'black');
xlabel("TOF [microsecs]");
ylabel("Voltage [V]");
ylim([-0.4 0.6]);

%--- fig 6: phase-after-hilbert-shift vs. frequency  ---------------------------

yhil = hilbert(V, 2048); % hilbert(xr,n) uses an n-point fast Fourier transform (FFT) to compute 
                         % the Hilbert transform. The input data is zero-padded or truncated to length n, as appropriate
% yhil = abs(yhil); % (real: original time signal / imaginary: phase-shifted signal / abs: both together)
shifted_y = fft(yhil, 2048); % crossing back over into the frequency domain
shifted_phase = angle(shifted_y);

figure(6);
plot(f / 1e6, shifted_phase);
title("Phase");
xlabel("Frequency [MHz]");
xlim(frange6dB/1e6);

%--- fig 7: linearized-phase-after-hilbert-shift vs. frequency -----------------

xfromto = intersect(find(f>frange6dB(1)), find(f<frange6dB(2))); % all indexes from ~26MHz to ~63MHz
ffromto = f(xfromto);                                            % all frequencies from ~26MHz to ~63MHz
yfromto = shifted_phase(xfromto);                                % the part of the phase-shifted plot that goes from ~26MHz to ~63MHz
yfromto(yfromto < 0) = 0;                                        % on that one, get rid of all negative numbers, findpeaks() doesn't like them
[pval, pind] = findpeaks(yfromto);
xpoints = ffromto(pind) / 1e6; % the frequency (x axis) at each peak (= at each phase jump) 
                               % ... only doing this in MHz because the calculation takes too long otherwise
                               % ... also, doing it in MHz instead of Hz makes the resulting line look smoother
N = length(xpoints); % the number of peaks (= the number of phase jumps)
ypoints = (1:N)*-1;  % i.e., [-1, -2, -3, -4, -5], since N is 5
fit = polyfit(xpoints, ypoints, 1); % first order; from the return values, only polyf is needed
slope = fit(1)
intercept = fit(2)
x = [frange6dB(1)/1e6:frange6dB(2)/1e6];
y = intercept + (x .* slope);

figure(7);
plot(x, y);
title("Lin. Phase");
xlabel("Frequency [MHz]");
xlim(frange6dB/1e6);
ylim([-6 0]);

%--- time of flight (TOF = t0 + dN / df) ---------------------------------------

tinterp = interp1(1:161, t, linspace(1, 161, 2048)); % make it 2048 elements long instead of just 161 elements
t0 = tinterp(xfromto(1));
% trange6dB = [tinterp(x1) tinterp(x2)];
t_ph = N / (frange6dB(2) - frange6dB(1));
TOF = t0 + t_ph; % ABSOLUTELY NO CLUE IF THIS IS CORRECT!!!! 😰😰😰
disp(["TOF = " num2str(TOF*1e6) " microsecs"]);
