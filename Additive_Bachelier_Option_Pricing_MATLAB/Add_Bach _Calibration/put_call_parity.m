function OptionPrices = put_call_parity(OptionPrices, strikes, discounts, forward)
%PUT_CALL_PARITY  Fill ONLY missing call prices via put-call parity.
%
%
%   OptionPrices : (2*M)×N matrix
%                  odd rows  (1,3,5,…) = calls
%                  even rows (2,4,6,…) = puts
%                  set missing prices to NaN
%   strikes      : 1×N row vector of strikes (one per column)
%   discounts    : 1×M row vector or scalar  (discount factor per maturity)
%   forward      : 1×M row vector or scalar  (forward price per maturity)
%
%   The routine fills C when it is NaN and the corresponding P is available:
%       C = P + DF·(F0 − K)
%   It does NOT attempt to fill missing puts.

    % --- dimension checks --------------------------------------------------
    [nRows, nStrikes] = size(OptionPrices);
    if mod(nRows,2) ~= 0
        error('OptionPrices must have an even number of rows (call/put pairs).');
    end
    if numel(strikes) ~= nStrikes
        error('Length of strikes vector must equal number of columns in OptionPrices.');
    end

    M = nRows / 2;   % number of maturities / call-put pairs

    % --- main loop ---------------------------------------------------------
    for i = 1:M              % maturity index
        callRow = 2*i - 1;   % odd row = call
        putRow  = 2*i;       % even row = put

        DF = discounts(i+1);
        F0 = forward(i);

        for j = 1:nStrikes   % strike/column index
            K = strikes(j);

            C = OptionPrices(callRow, j);
            P = OptionPrices(putRow,  j);

            % Only fill CALL if it is missing and PUT is available
            if isnan(C) && ~isnan(P)
                call = P + DF * (F0 - K);
                if call>0
                    OptionPrices(callRow, j) = call;
                end
            end

            if ~isnan(C) && isnan(P)
                put = C - DF * (F0 - K);
                if put>0
                    OptionPrices(putRow, j) = put;
                end
            end
        end
    end
end
