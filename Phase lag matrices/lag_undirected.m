function [lag_matrix]=lag_undirected(data)
% This function calculates the mean phase difference between pairs of ROIs
% uses code for PLV, but focusing on the lag value.
%
%
% Input: 'data': matrix containing scalp EEG data
%
% Output: 'lag_matrix'
%
% Note: output elements ( i , j ) represent a connection from i to j.
%
% This function uses the following function developed by others:
%    1. generate_iAAFT_it : to generate univariate surrogates
%    2. PLF_lag           : to compute the PLV
% -----------------------------------------
%
% References:
% [1] Schmidt, Helmut, et al. "Dynamics on networks: the role of local
% dynamics and global networks on the emergence of hypersynchronous neural
% activity." PLoS Comput Biol 10.11 (2014): e1003947.
%
% ------------------------------------------
% % Created by Marinho Lopes and Petroula Laiou
% m.lopes@exeter.ac.uk, petroula.laiou@kcl.ac.uk
% ------------------------------------------

sf = 250; % WE HAVE TO SPECIFY THE SAMPLING FREQUENCY

n_ch = min(size(data)); % number of channels
n_samp=max(size(data)); % number of sample points

res = 2*pi/sf;

if size(data,1)>size(data,2)
    data=data'; % Rearrange data as channels x time
end

lag_matrix=zeros(n_ch,n_ch); % connectivity matrix

data=data';
DataMean = mean(data);
DataStd = std(data);
data = (data-repmat(DataMean,[n_samp,1]))./repmat(DataStd,[n_samp,1]);
%data=(data-mean(data))./std(data); % normalization
for ch1=1:n_ch-1
    for ch2=ch1+1:n_ch
        [plf,lag]=PLF_lag_padding(data(:,ch1),data(:,ch2)); %contains data padding
        
        %for undirected with zero lag
        lag_matrix(ch1,ch2)=lag;
        lag_matrix(ch2,ch1)=lag;
        
    end
end


end