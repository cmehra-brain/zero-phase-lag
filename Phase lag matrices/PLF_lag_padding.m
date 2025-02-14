function [plf,lag] = PLF_lag_padding(S1,S2)
% This function computes the phase locking value (PLV) (or phase locking
% factor, PLF) between two signals and the phase difference between them.
% 
% -------------------------------------------------------------------------
%   [plf,lag] = PLF_lag(S1,S2)
%
% INPUTS:
%   S1, S2 - two 1D signals [1 x time points]
%
% OUTPUT:
%   plf - phase locking factor/value UNdirected index
%   lag - lag between the two signals
%
% Note: output elements ( i , j ) represent a connection from i to j. 
%
% -------------------------------------------------------------------------
% References:
%  Lachaux et al., (1999); 
%  Petkov, G., Goodfellow, M., Richardson, M. P., & Terry, J. R. (2014). 
%    A critical role for network structure in seizure onset: 
%      a computational modeling approach. Front Neurol, 5, 261. 
%      http://doi.org/10.3389/fneur.2014.00261
% -------------------------------------------------------------------------
% Version v0: GPetkov, 02.06.2013, K&P lab
% Version v1: M.Lopes, 2017
%
% Summary of all changes implemented in this function (an item per version)
%   
%  - Marinho Lopes, m.lopes@exeter.ac.uk, 2017
%    I just added the 'lag' as an output of this function.
%  - Chirag Mehra, chirag.1.mehra@kcl.ac.uk - Dec 2022
%       edited to include 2 seconds of data padding before the hilbert transform
%       (included in the input signal), but to crop this out after
% -------------------------------------------------------------------------
% This function was tested by:
% - Leandro Junges L.L.L.Junges@exeter.ac.uk, Nov/2018
%   - no modification
% -------------------------------------------------------------------------


plf = [];
    if nargin<2 %returns the number of function input arguments given in the call to the currently executing function. Use this syntax in the body of a function only.
        disp( 'usage: [plf,lag] = PLF_lag(S1,S2);' );             
        return;     
    end
    
    %doing this has the same effect of the more efficient cropping method
    %below
%     hilb_1 = hilbert(S1');
%     hilb_1_cropped = hilb_1(:,500:5500-1); %discard first 2 and last 2 seconds
%     
%     hilb_2 = hilbert(S2');
%     hilb_2_cropped = hilb_2(:,500:5500-1); %discard first 2 and last 2 seconds
%     Cm = angle(hilb_1_cropped)-angle(hilb_2_cropped);

    
    Cm = angle(hilbert(S1'))-angle(hilbert(S2')); 
    Cm = Cm(:,500:5500-1); %crop out first 2 and last 2 seconds
    arg = mean(exp(1i*Cm));    
    plf = abs(arg);
    lag = angle(arg);     
return
end
