function lg = loggamma_complex(z)
% Lanczos approximation (g = 7, n = 9)
coeff = [ ...
    0.9999999999998099322768470047348,...
    676.5203681218851,...
   -1259.1392167224028,...
     771.32342877765313,...
    -176.61502916214059,...
      12.507343278686905,...
      -0.13857109526572012,...
       9.9843695780195716e-6,...
       1.5056327351493116e-7];

mask = real(z) < 0.5;
lg     = zeros(size(z),'like',z);

if any(mask)
    lg(mask) = log(pi) - log(sin(pi*z(mask))) - loggamma_complex(1 - z(mask));
end

if any(~mask)
    z2 = z(~mask) - 1;
    x  = coeff(1);
    for k = 2:numel(coeff)
        x = x + coeff(k)./(z2 + k - 1);
    end
    t  = z2 + 7.5;               % g + 0.5  con g = 7
    lg(~mask) = 0.5*log(2*pi) + (z2+0.5).*log(t) - t + log(x);
end
end