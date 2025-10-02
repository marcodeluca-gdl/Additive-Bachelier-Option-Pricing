function dates = getDates(t0)
%   Input:
%     t0: Value Date
%   Output: 
%     Vector of dates in format datetime containing all futures expiries

dates=  [t0, datetime('17-Aug-2020'),datetime('17-Nov-2020'),...
        datetime('17-Feb-2021'), datetime('17-May-2021'), datetime('17-Aug-2021'),...
        datetime('16-Nov-2021'), datetime('17-May-2022'),datetime('16-Nov-2022')]';
end