function [y, t] = rcpuls(a,tau,fs,span)
% This function is a modified version of the one given in the text, Program 2.2-1
% a: Roll off/ Excess bandwidth factor
% tau: Nyquist period or symbol time 
% fs: sampling frequency at which the continuous pulse is sampled
% span: Defines width for truncation, number of tau periods on either side of the peak of the pulse
% The pulse has a one sided bandwidth, BW = (1+alpha)/(2*tau);
% Choose fs > 2*BW
% Use all the provided functions with caution!

t_positive = eps:(1/fs):span*tau;  % Replace 0 with eps (smallest positive number MATLAB can produce) to prevent NANs
t = [-fliplr(t_positive(2:end)) t_positive];
tpi = pi/tau; 
atpi = tpi*a; 
at = 4*a^2/tau^2;
y = sin(tpi*t).*cos(atpi*t)./(tpi*t.*(1-at*t.^2));
%y = y/norm(y);
