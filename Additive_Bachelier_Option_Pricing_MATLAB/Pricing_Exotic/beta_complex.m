function B = beta_complex(x,y)
% B(x,y) for complex inputs
B = gamma_complex(x).*gamma_complex(y) ./ gamma_complex(x+y);
end