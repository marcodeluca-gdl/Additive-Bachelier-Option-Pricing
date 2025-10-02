function selectedIndices = paretoSelection(values, threshold)
% PARETOSELECTION - Selects the most significant elements according to Pareto principle
% This function implements the Pareto principle (80/20 rule) to identify the subset
% of elements that contribute to a specified cumulative fraction of the total
%
% INPUTS:
%   values    - Vector of positive values to analyze
%   threshold - Desired cumulative fraction (e.g., 0.8 for 80% of total)
%               Default: 0.8 if not specified
%
% OUTPUTS:
%   selectedIndices - Original indices of the selected values that meet the threshold

% Set default threshold to 80% if not provided
if nargin < 2
    threshold = 0.8; % Default: 80% Pareto threshold
end

% Sort values in descending order to prioritize largest contributors
% sortedValues: values arranged from largest to smallest
% sortedIndices: original positions of the sorted values
[sortedValues, sortedIndices] = sort(values, 'descend');

% Calculate cumulative sum of sorted values
% This shows the running total as we add each element
cumulativeSum = cumsum(sortedValues);

% Calculate total sum of all values for normalization
total = sum(values);

% Convert cumulative sum to cumulative fraction (percentage of total)
% Each element shows what fraction of the total is covered up to that point
cumulativeFraction = cumulativeSum / total;

% Find the first position where cumulative fraction meets or exceeds threshold
% This identifies the minimal set of elements needed to reach the target percentage
idx = find(cumulativeFraction >= threshold, 1, 'first');

% Return the original indices of the selected elements
% These are the most significant contributors according to Pareto principle
selectedIndices = sortedIndices(1:idx);

end