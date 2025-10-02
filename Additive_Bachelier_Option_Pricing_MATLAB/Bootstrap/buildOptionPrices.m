function [optionPrices, strikes] = buildOptionPrices(callFolder, putFolder, target)

%  Inputs
%  -------
%  callFolder : char | string
%       Folder containing ONLY the *call.csv files
%  putFolder  : char | string
%       Folder containing ONLY the *put.csv files
%  target     : numeric (yyyymmdd)
%       Date that will be passed straight to readCSVData.
%
%  Output
%  -------
%  optionPrices : double(2*M-by-N)
%       Row-pairs (call, put) for each expiration
%  strikes      : double(1-by-N)
%       Strike vector extracted from the first call file

    %% Gather call files ---------------------------------------------------
    callFiles = dir(fullfile(callFolder, '*call.csv'));
    if isempty(callFiles)
        error('No "*call.csv" files found in %s', callFolder);
    end

    % Sort alphabetically – because YYYYMM in the name makes this
    % chronological as well.
    [~, idx] = sort({callFiles.name});
    callFiles = callFiles(idx);
    M = numel(callFiles);          % Number of expirations

    %% Strike vector (taken from first call file) --------------------------
    firstCallTbl = readtable(fullfile(callFolder, callFiles(1).name));
    firstCallMat = table2array(firstCallTbl);
    strikes      = firstCallMat(1, 2:end);   % Row-1, columns 2:end
    N            = numel(strikes);

    %% Pre-allocate final matrix ------------------------------------------
    optionPrices = zeros(2*M, N);

    %% Loop over each expiration ------------------------------------------
    for k = 1:M
        callPath = fullfile(callFolder, callFiles(k).name);

        % Build put-file name by replacing 'call' with 'put'
        putName = strrep(callFiles(k).name, 'call', 'put');
        putPath = fullfile(putFolder, putName);

        if ~isfile(putPath)
            error('Missing put file: %s', putPath);
        end

        % Fill rows (2k-1 : 2k) with the two-row block from readCSVData
        optionPrices(2*k-1 : 2*k, :) = readCSVData(callPath, putPath, target);
    end
end
