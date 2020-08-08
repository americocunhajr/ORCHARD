
% -----------------------------------------------------------------
%  shannon_entropy.m
%
%  This function computes the Shannon entropy of pX, the PDF of a
%  random variable X: Omega -> R (with finite support).
%  
%  Shannon entropy of pdf_X is defined as
%
%                   --
%                  |
%  S(pdf_X) :=   - | pdf_X(x) ln(pdf_X(x)) dx
%                  |
%                -- R
%  
%  where pdf_X is the PDF of X.
%
%  Input:
%  bins_X - (numbins x Ndt) random variable X bins
%  freq_X - (numbins x Ndt) random variable X frequencies
%
%  Output:
%  SpX    - (1 x Ndt) Shannon entropy of X's PDF
% ----------------------------------------------------------------- 
%  programmer: Americo Barbosa da Cunha Junior
%              americo.cunhajr@gmail.com
%
%  last update: July 16, 2016
% -----------------------------------------------------------------

% -----------------------------------------------------------------
function H_pdf_X = shannon_entropy(bins_X,freq_X)
    
    % check number of arguments
    if nargin < 2
        error(' Too few inputs.')
    elseif nargin > 2
        error(' Too many inputs.')
    end
    
    % check input arguments
    if size(bins_X) ~= size(freq_X)
        error('bins_X and freq_X must have the same size.')
    end
    
    % compute data matrix dimensions
	[numbins,Ndt] = size(bins_X);
        
    % preallocate memory for entropy vector
	H_pdf_X = zeros(1,Ndt);
    
    % loop over time instants
    for n=1:Ndt
        
        % compute the width of the bins
        binwidth = bins_X(2:end,n) - bins_X(1:end-1,n);
    
        % aproximate the probability over a bin
        binprob = binwidth.*freq_X(1:end-1,n);
        
        % compute Shannon entropy of freq_X
        shannon = - cumsum(binprob.*log(binprob./binwidth));
        
        % save Shannon entropy of freq_X
        H_pdf_X(1,n) = shannon(end);
    end
    
return
% -----------------------------------------------------------------
