
% -----------------------------------------------------------------
%  randvar_probval.m
%
%  This functions computes the probability of a random
%  variable be less than or equal to a given value.
%
%  input:
%  data_bins - (Nbins x Ndt) bins matrix
%  data_freq - (Nbins x Ndt) frequency matrix
%  value     - (    1 x Ndt) reference values vector
%
%  output:
%  prob - (    1 x Ndt) probability values vector
% ----------------------------------------------------------------- 
%  programmer: Americo Barbosa da Cunha Junior
%              americo.cunhajr@gmail.com
%
%  last update: Oct 6, 2016
% -----------------------------------------------------------------

% -----------------------------------------------------------------
function prob = randvar_probval(data_bins,data_freq,value)

    % check number of arguments
    if nargin < 3
        error('Too few inputs.')
    elseif nargin > 3
        error('Too many inputs.')
    end
    
    % check arguments
    if size(data_bins) ~= size(data_freq)
        error('data_bins and data_freq must be matrices with same dimensions')
    end
    
    % compute matrices dimensions
	[Nbins,Ndt] = size(data_freq);
    
    % preallocate memory for values matrix
	prob = zeros(1,Ndt);
    
    % loop over time instants
    for n=1:Ndt
        
        % indicator of value
        ind_value = data_bins(:,n) - value(1,n) <= 0.0;

        % compute probability
        prob(1,n) = trapz(data_bins(:,n).*ind_value,data_freq(:,n));
        
    end

return
% -----------------------------------------------------------------
