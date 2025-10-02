function [cumulative] = lewisFormula(phi, a, Ra, x, s, t)
% Lewis formula for computing cumulative distribution function
%
% Inputs:
%   phi - Characteristic function handle
%   a, Ra - Parameters in Lewis formula
%   x   - Evaluation points
%   s, t - Time parameters
%
% Output:
%   cumulative - CDF values at x points

    M = 15;  % Number of points for FFT integration
    
    % Compute integral using FFT
    [intFFT] = integralFFT(phi, M, x, a, s, t);
    
    % Apply Lewis formula
    cumulative = Ra - exp(-a*x) / (2*pi) .* intFFT;
    
end