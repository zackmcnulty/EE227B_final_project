function [Sc] = real2complex(Sr)
% Given a real signal vector Sr,
% output a complex signal vector Sc

s_real = Sr(1:end/2);
s_imag = Sr(end/2+1:end);

Sc = s_real + 1i*s_imag;
