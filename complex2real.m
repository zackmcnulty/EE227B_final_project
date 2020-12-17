function [Hr, Sr] = complex2real(Hc, Sc)
% Given a complex channel matrix Hc and signal vector Sc,
% output a real channel matrix Hr and signal vector Sr

H_real = real(Hc);
H_imag = imag(Hc);
Hr = [H_real -1*H_imag; H_imag H_real];

Sr = [real(Sc); imag(Sc)];
