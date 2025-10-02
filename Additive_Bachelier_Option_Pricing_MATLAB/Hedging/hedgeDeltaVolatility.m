function futureInfo = hedgeDeltaVolatility(callInfo, fwdPrices, sigmat, ...
    ttm, disc)
%--------------------------------------------------------------------------
% Hedge the delta of a portfolio of calls with futures.
%
% INPUT
% callInfo : 5 × nCall matrix
% row 1 → notional (number of option contracts) for each call
% row 2 → strike of each call
% row 3 → time-to-maturity (years) of each call
% row 4 → time-to-maturity indexes of each call
% row 5 → number of calls
%
% fwdPrices : vector of forward prices. Element k matches sigmat(k), ttm(k).
% sigmat : vector of implied volatilities (same ordering as ttm).
% ttm : vector of maturities (years).
% disc : vector of discount factors
%
% OUTPUT
% futureInfo : 4 × nCall matrix
% row 1 → futures notional amounts (number of futures to buy (+) or sell (−) to hedge)
% row 2 → maturity (years) of the corresponding future
% row 3 → time-to-maturity indexes
% row 4 → number of futures contracts

nCall = size(callInfo, 2);
futureNotional = zeros(1, nCall);
futureTTM = callInfo(3, :); % same maturity as each call
TtmIdx = callInfo(4,:);
nFutures = zeros(1,nCall);

for i = 1:nCall
    % ----------- extract data for the i-th call -----------
    N_opt = callInfo(1, i); % option notional (contracts)
    K = callInfo(2, i); % strike
    tau = callInfo(3, i); % maturity in years
    
    % Locate the market row that has the same maturity
    k = find(abs(ttm - tau) < 1e-10, 1, 'first');
    if isempty(k)
        error('Maturity %.6g not found in ttm vector.', tau);
    end
    
    F = fwdPrices(k); % forward price
    vol = sigmat(k); % implied volatility
    DF = disc(k); % discount factor
    
    % -------------- option delta (per contract) --------------
    delta_call = deltaOption(K, F, vol, tau, DF, 'Nan'); % ∂C/∂F
    
    % Total delta contributed by this call position
    portfolio_delta = N_opt * delta_call;
    
    % -------------- futures needed to neutralize delta ----------
    % One futures contract has delta +1 → need −portfolio_delta
    futureNotional(i) = -portfolio_delta;
    nFutures(i) = futureNotional(i) / fwdPrices(TtmIdx(i));
end

% Assemble output: 4 rows × nCall columns
futureInfo = [futureNotional;    % Row 1: Futures notional amounts
              futureTTM;         % Row 2: Futures maturities
              TtmIdx;            % Row 3: Time-to-maturity indexes
              ceil(nFutures)];   % Row 4: Number of futures contracts
end