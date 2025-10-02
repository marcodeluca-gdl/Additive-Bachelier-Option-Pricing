function delta = deltaOption(strike, fwdPrice,sigmat, ttm, disc, flag)
    y = (strike-fwdPrice)/(sigmat*sqrt(ttm));
    delta = disc * normcdf(-y/sigmat);

    if flag == 'put'
        delta = delta-disc;
    end

end