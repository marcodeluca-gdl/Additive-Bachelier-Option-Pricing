function [I] = integralFFT(phi, M, queryPoints,a,s,t)
% This function computes the integral of the Lewis formula
% using the Fast Fourier Transform (FFT).
% INPUTS:
% phi: characteristic function of the forward price process
% M: number of points for the FFT
% dz: step size for the FFT
% queryPoints: points at which to evaluate the integral
% OUTPUTS:
% I: vector of integral values at the query points


% FFT parameters
N = 2^M;
xi_1 = - 5*sqrt(t-s);
d_xi = -2*xi_1/(N-1);
dz=2*pi/(N*d_xi);
z_1 = -dz*(N-1)/2;
xi = xi_1:d_xi:-xi_1;
z = z_1:dz:-z_1;

% integrand function
f = phi(xi - 1i*a) ./ (1i*xi + a);
f_tilde = f .* exp(-1i * z_1 * d_xi .* (0:N-1));

FFT = fft(f_tilde);

prefactor = d_xi * exp(-1i * xi_1 * z);
I = prefactor .* FFT;
I = real(I);

% We interpolate the integral in the relevant points
I = interp1(z, I, queryPoints);

end