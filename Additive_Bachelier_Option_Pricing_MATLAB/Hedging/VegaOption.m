function vega = VegaOption(strike, fwdPrice,sigmat, ttm, disc)
    y = (strike-fwdPrice)/(sigmat*sqrt(ttm));
    vega=disc*normpdf(-y/sigmat);
end