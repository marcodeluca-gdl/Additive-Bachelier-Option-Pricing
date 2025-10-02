function [disc, fwd_prices] = bootstrap(optionPrices, strikes)
    % This function calculates discount factors and forward prices from option price data
    %
    % Inputs:
    %   optionPrices: Matrix of option prices (rows = maturities, columns = strikes)
    %   strikes:      Vector of strike prices corresponding to the columns of optionPrices
    %
    % Outputs:
    %   disc:         Vector of discount factors for each maturity
    %   fwd_prices:   Vector of forward prices for each maturity
    
    % Get dimensions of the input matrix
    N = size(optionPrices);
    
    % Calculate number of discount factors (maturities)
    % Each maturity has 2 rows
    disc_num = N(1)/2;
    
    % Initialize output vectors
    disc = ones(disc_num, 1);        % Discount factors (first one is 1.0 at t=0)
    fwd_prices = zeros(disc_num-1, 1); % Forward prices
    
    % Loop through each maturity (starting from the second one)
    for i = 2:disc_num
        % Initialize empty matrix to store valid option price differences and strikes
        matrix = [];
        
        % Loop through each strike price
        for j = 1:N(2)
            % Check if both option prices (call and put) for this maturity and strike are valid
            if ~isnan(optionPrices(2*i-1, j)) && ~isnan(optionPrices(2*i, j))
                % Store the difference between the two option prices and the strike
                matrix = [matrix; optionPrices(2*i-1, j) - optionPrices(2*i, j), strikes(j)];
            end
        end
        
        % Calculate means for both the sinthetic forwards and strikes
        G_mean = mean(matrix(:, 1));
        K_mean = mean(matrix(:, 2));
        
        % Calculate the slope of the regression line using covariance/variance formula
        % This estimates the discount factor as the negative of the slope
        numerator = (matrix(:, 1) - G_mean)' * (matrix(:, 2) - K_mean);
        denominator = sum((matrix(:, 2) - K_mean).^2);
        disc(i) = -numerator / denominator;
       
        par = [matrix(:, 2) ones(size(matrix(:, 2)))] \ matrix(:, 1);
        
        % Calculate forward price from regression intercept and discount factor
        fwd_prices(i-1) = par(2) / disc(i);
    end
end