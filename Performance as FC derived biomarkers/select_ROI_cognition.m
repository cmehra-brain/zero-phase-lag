%% script to select ROIs for spatial working memory analysis
% this script takes the whole FC adjacency matrix and only keeps ROIs at
% are based on SMW literature. ROI numbers and literature references
% included

% -----------------------------------------------------------------------
% This script was produced and tested by Chirag Mehra, for the work found in the manuscript: 
% Mehra et al., (2025): "Zero-phase-delay synchrony between interacting neural populations: implications for functional connectivity derived biomarkers"
% Please cite the most up to date version of the manuscript when using this script
% -----------------------------------------------------------------------

%% All cortical ROIs from DK atlas

% ROI numbers involved in SWM: 3	5	7	11	13	14	15	17	18	19	26	27	28	31	37	39	41	45	47	48	49	51	52	53	60	61	62	65

% ROI label	    ROI number	    refs to support relatioship with spatial WM
% 'ctx-lh-bankssts'	1	
% 'ctx-lh-caudalanteriorcingulate'	2	
% 'ctx-lh-caudalmiddlefrontal'	3	https://link.springer.com/article/10.3758/CABN.3.4.255; https://www.mdpi.com/2076-3425/7/4/43
% 'ctx-lh-cuneus'	4	
% 'ctx-lh-entorhinal'	5	https://jps.biomedcentral.com/articles/10.1186/s12576-021-00805-1; https://pubmed.ncbi.nlm.nih.gov/26180122/
% 'ctx-lh-fusiform'	6	
% 'ctx-lh-inferiorparietal'	7	https://www.nature.com/articles/s41598-017-06293-x#Sec4 / https://www.biorxiv.org/content/10.1101/650077v1.full.pdf ;https://pubmed.ncbi.nlm.nih.gov/26180122/
% 'ctx-lh-inferiortemporal'	8	
% 'ctx-lh-isthmuscingulate'	9	
% ctx-lh-lateraloccipital'	10	
% 'ctx-lh-lateralorbitofrontal'	11	https://link.springer.com/article/10.3758/CABN.3.4.255; https://www.mdpi.com/2076-3425/7/4/43
% 'ctx-lh-lingual'	12	
% ctx-lh-medialorbitofrontal'	13	https://link.springer.com/article/10.3758/CABN.3.4.255; https://www.mdpi.com/2076-3425/7/4/43
% ctx-lh-middletemporal'	14	https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3246590/
% 'ctx-lh-parahippocampal'	15	https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4522562/ (spatial configuration)
% 'ctx-lh-paracentral'	16	
% 'ctx-lh-parsopercularis'	17	https://link.springer.com/article/10.3758/CABN.3.4.255; https://www.mdpi.com/2076-3425/7/4/43; https://link.springer.com/article/10.1007/s10339-016-0772-7
% 'ctx-lh-parsorbitalis'	18	https://link.springer.com/article/10.3758/CABN.3.4.255; https://www.mdpi.com/2076-3425/7/4/43; https://link.springer.com/article/10.1007/s10339-016-0772-7
% 'ctx-lh-parstriangularis'	19	https://link.springer.com/article/10.3758/CABN.3.4.255; https://www.mdpi.com/2076-3425/7/4/43; https://link.springer.com/article/10.1007/s10339-016-0772-7
% 'ctx-lh-pericalcarine'	20	
% ctx-lh-postcentral'	21	
% 'ctx-lh-posteriorcingulate'	22	
% ctx-lh-precentral'	23	
% 'ctx-lh-precuneus'	24	
% 'ctx-lh-rostralanteriorcingulate'	25	
% ctx-lh-rostralmiddlefrontal'	26	https://link.springer.com/article/10.3758/CABN.3.4.255; https://www.mdpi.com/2076-3425/7/4/43; https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6659551/
% 'ctx-lh-superiorfrontal'	27	https://link.springer.com/article/10.3758/CABN.3.4.255; https://www.mdpi.com/2076-3425/7/4/43; https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6659551/
% 'ctx-lh-superiorparietal'	28	https://www.nature.com/articles/s41598-017-06293-x#Sec4 / https://www.biorxiv.org/content/10.1101/650077v1.full.pdf ; https://link.springer.com/article/10.3758/CABN.3.4.255; https://pubmed.ncbi.nlm.nih.gov/26180122/; https://link.springer.com/article/10.1007/s10339-016-0772-7
% 'ctx-lh-superiortemporal'	29	
% 'ctx-lh-supramarginal'	30	
% 'ctx-lh-frontalpole'	31	https://link.springer.com/article/10.3758/CABN.3.4.255; https://www.mdpi.com/2076-3425/7/4/43
% 'ctx-lh-temporalpole'	32	
% 'ctx-lh-transversetemporal'	33	
% 'ctx-lh-insula'	34	
% 'ctx-rh-bankssts'	35	
% 'ctx-rh-caudalanteriorcingulate'	36	
% 'ctx-rh-caudalmiddlefrontal'	37	https://link.springer.com/article/10.3758/CABN.3.4.255; https://www.mdpi.com/2076-3425/7/4/43
% 'ctx-rh-cuneus'	38	
% 'ctx-rh-entorhinal'	39	https://jps.biomedcentral.com/articles/10.1186/s12576-021-00805-1; https://pubmed.ncbi.nlm.nih.gov/26180122/
% 'ctx-rh-fusiform'	40	
% 'ctx-rh-inferiorparietal'	41	https://www.nature.com/articles/s41598-017-06293-x#Sec4 / https://www.biorxiv.org/content/10.1101/650077v1.full.pdf ;https://pubmed.ncbi.nlm.nih.gov/26180122/
% 'ctx-rh-inferiortemporal'	42	
% 'ctx-rh-isthmuscingulate'	43	
% 'ctx-rh-lateraloccipital'	44	
% 'ctx-rh-lateralorbitofrontal'	45	https://link.springer.com/article/10.3758/CABN.3.4.255; https://www.mdpi.com/2076-3425/7/4/43
% 'ctx-rh-lingual'	46	
% 'ctx-rh-medialorbitofrontal'	47	https://link.springer.com/article/10.3758/CABN.3.4.255; https://www.mdpi.com/2076-3425/7/4/43
% 'ctx-rh-middletemporal'	48	https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3246590/
% 'ctx-rh-parahippocampal'	49	https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4522562/ (spatial configuration)
% 'ctx-rh-paracentral'	50	
% 'ctx-rh-parsopercularis'	51	https://link.springer.com/article/10.3758/CABN.3.4.255; https://www.mdpi.com/2076-3425/7/4/43; https://link.springer.com/article/10.1007/s10339-016-0772-7
% 'ctx-rh-parsorbitalis'	52	https://link.springer.com/article/10.3758/CABN.3.4.255; https://www.mdpi.com/2076-3425/7/4/43; https://link.springer.com/article/10.1007/s10339-016-0772-7
% 'ctx-rh-parstriangularis'	53	https://link.springer.com/article/10.3758/CABN.3.4.255; https://www.mdpi.com/2076-3425/7/4/43; https://link.springer.com/article/10.1007/s10339-016-0772-7
% ctx-rh-pericalcarine'	54	
% 'ctx-rh-postcentral'	55	
% 'ctx-rh-posteriorcingulate'	56	
% 'ctx-rh-precentral'	57	
% 'ctx-rh-precuneus'	58	
% 'ctx-rh-rostralanteriorcingulate'	59	
% 'ctx-rh-rostralmiddlefrontal'	60	https://link.springer.com/article/10.3758/CABN.3.4.255; https://www.mdpi.com/2076-3425/7/4/43; https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6659551/
% 'ctx-rh-superiorfrontal'	61	https://link.springer.com/article/10.3758/CABN.3.4.255; https://www.mdpi.com/2076-3425/7/4/43; https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6659551/ 
% 'ctx-rh-superiorparietal'	62	https://www.nature.com/articles/s41598-017-06293-x#Sec4 / https://www.biorxiv.org/content/10.1101/650077v1.full.pdf ; https://link.springer.com/article/10.3758/CABN.3.4.255; https://pubmed.ncbi.nlm.nih.gov/26180122/; https://link.springer.com/article/10.1007/s10339-016-0772-7
% 'ctx-rh-superiortemporal'	63	
% 'ctx-rh-supramarginal'	64	
% 'ctx-rh-frontalpole'	65	https://link.springer.com/article/10.3758/CABN.3.4.255; https://www.mdpi.com/2076-3425/7/4/43
% 'ctx-rh-temporalpole'	66	
% 'ctx-rh-transversetemporal'	67	
% 'ctx-rh-insula'	68	



%% clear variables and add paths as required
clear
clc

addpath('...') %add path to scripts

%EEGs paths 
MainEEG = '...EEG_networks_DK_1200_methodspaper/full_nets'; % define a main directory
EEG_list = dir(fullfile(MainEEG,'**/nets_ROIdata*'));

%sort rows of EEG
EEG_list = struct2table(EEG_list); % convert the struct array to a table
EEG_list = sortrows(EEG_list, {'name', 'folder'}); % sort by name followed by folder type
EEG_list = table2struct(EEG_list); % change it back to struct array if necessary

%% output each of the SWI specific networks
path2output_COH_1_4 = '...SWM_nets\COH\1_4Hz';
path2output_COH_4_8 = '...SWM_nets\COH\4_8Hz';
path2output_COH_8_13 = '...SWM_nets\COH\8_13Hz';
path2output_COH_13_20 = '...SWM_nets\COH\13_20Hz';
path2output_COH_20_32 = '...SWM_nets\COH\20_32Hz';

path2output_img_COH_1_4 = '...SWM_nets\img_COH\1_4Hz';
path2output_img_COH_4_8 = '...SWM_nets\img_COH\4_8Hz';
path2output_img_COH_8_13 = '...SWM_nets\img_COH\8_13Hz';
path2output_img_COH_13_20 = '...SWM_nets\img_COH\13_20Hz';
path2output_img_COH_20_32 = '...SWM_nets\img_COH\20_32Hz';

path2output_PLV_1_4 = '...SWM_nets\PLV\1_4Hz';
path2output_PLV_4_8 = '...SWM_nets\PLV\4_8Hz';
path2output_PLV_8_13 = '...SWM_nets\PLV\8_13Hz';
path2output_PLV_13_20 = '...SWM_nets\PLV\13_20Hz';
path2output_PLV_20_32 = '...SWM_nets\PLV\20_32Hz';

path2output_wPLI_1_4 = '...SWM_nets\wPLI\1_4Hz';
path2output_wPLI_4_8 = '...SWM_nets\wPLI\4_8Hz';
path2output_wPLI_8_13 = '...SWM_nets\wPLI\8_13Hz';
path2output_wPLI_13_20 = '...SWM_nets\wPLI\13_20Hz';
path2output_wPLI_20_32 = '...SWM_nets\wPLI\20_32Hz';

path2output_AEC_1_4 = '...SWM_nets\AEC\1_4Hz';
path2output_AEC_4_8 = '...SWM_nets\AEC\4_8Hz';
path2output_AEC_8_13 = '...SWM_nets\AEC\8_13Hz';
path2output_AEC_13_20 = '...SWM_nets\AEC\13_20Hz';
path2output_AEC_20_32 = '...SWM_nets\AEC\20_32Hz';

path2output_orth_AEC_1_4 = '...SWM_nets\orth_AEC\1_4Hz';
path2output_orth_AEC_4_8 = '...SWM_nets\orth_AEC\4_8Hz';
path2output_orth_AEC_8_13 = '...SWM_nets\orth_AEC\8_13Hz';
path2output_orth_AEC_13_20 = '...SWM_nets\orth_AEC\13_20Hz';
path2output_orth_AEC_20_32 = '...SWM_nets\orth_AEC\20_32Hz';

%%

for i = 1:length(EEG_list) %bring up each set of networks one at a time
    disp(i)
    cd(EEG_list(i).folder);
    data = load(EEG_list(i).name);
    fnames = fieldnames(data);
    EEG_network = data.(fnames{1});
    
    %get rid of the temporarily needed variables
    clear data
    clear fnames
    
    %Ensure only the first network is chosen.
    EEG_network = EEG_network(:,:,1); %only keep the first, get rid of the rest
    
    %specify the ROIs you want - these are based on the DK atlas 
    SWM_ROI = [3	5	7	11	13	14	15	17	18	19	26	27	28	31	37	39	41	45	47	48	49	51	52	53	60	61	62	65];

    EEG_network = EEG_network(SWM_ROI, SWM_ROI);
    
    FC_freq = erase(string(EEG_list(i).folder), "...\EEG_networks_DK_1200_methodspaper\full_nets\");
    
    
    %%
    
    if FC_freq == "COH\1_4Hz"
        cd(path2output_COH_1_4)
        save([EEG_list(i).name],'EEG_network');
        
    elseif FC_freq == "COH\4_8Hz"
        cd(path2output_COH_4_8)
        save([EEG_list(i).name],'EEG_network');
        
    elseif FC_freq == "COH\8_13Hz"
        cd(path2output_COH_8_13)
        save([EEG_list(i).name],'EEG_network');
        
    elseif FC_freq == "COH\13_20Hz"
        cd(path2output_COH_13_20)
        save([EEG_list(i).name],'EEG_network');
        
    elseif FC_freq == "COH\20_32Hz"
        cd(path2output_COH_20_32)
        save([EEG_list(i).name],'EEG_network');
        
    elseif FC_freq == "img_COH\1_4Hz"
        cd(path2output_img_COH_1_4)
        save([EEG_list(i).name],'EEG_network');
        
    elseif FC_freq == "img_COH\4_8Hz"
        cd(path2output_img_COH_4_8)
        save([EEG_list(i).name],'EEG_network');
        
    elseif FC_freq == "img_COH\8_13Hz"
        cd(path2output_img_COH_8_13)
        save([EEG_list(i).name],'EEG_network');
        
    elseif FC_freq == "img_COH\13_20Hz"
        cd(path2output_img_COH_13_20)
        save([EEG_list(i).name],'EEG_network');
        
    elseif FC_freq == "img_COH\20_32Hz"
        cd(path2output_img_COH_20_32)
        save([EEG_list(i).name],'EEG_network');
        
    elseif FC_freq == "PLV\1_4Hz"
        cd(path2output_PLV_1_4)
        save([EEG_list(i).name],'EEG_network');
        
    elseif FC_freq == "PLV\4_8Hz"
        cd(path2output_PLV_4_8)
        save([EEG_list(i).name],'EEG_network');
        
    elseif FC_freq == "PLV\8_13Hz"
        cd(path2output_PLV_8_13)
        save([EEG_list(i).name],'EEG_network');
        
    elseif FC_freq == "PLV\13_20Hz"
        cd(path2output_PLV_13_20)
        save([EEG_list(i).name],'EEG_network');
        
    elseif FC_freq == "PLV\20_32Hz"
        cd(path2output_PLV_20_32)
        save([EEG_list(i).name],'EEG_network');
        
    elseif FC_freq == "wPLI\1_4Hz"
        cd(path2output_wPLI_1_4)
        save([EEG_list(i).name],'EEG_network');
        
    elseif FC_freq == "wPLI\4_8Hz"
        cd(path2output_wPLI_4_8)
        save([EEG_list(i).name],'EEG_network');
        
    elseif FC_freq == "wPLI\8_13Hz"
        cd(path2output_wPLI_8_13)
        save([EEG_list(i).name],'EEG_network');
        
    elseif FC_freq == "wPLI\13_20Hz"
        cd(path2output_wPLI_13_20)
        save([EEG_list(i).name],'EEG_network');
        
    elseif FC_freq == "wPLI\20_32Hz"
        cd(path2output_wPLI_20_32)
        save([EEG_list(i).name],'EEG_network');
        
    elseif FC_freq == "AEC\1_4Hz"
        cd(path2output_AEC_1_4)
        save([EEG_list(i).name],'EEG_network');
        
    elseif FC_freq == "AEC\4_8Hz"
        cd(path2output_AEC_4_8)
        save([EEG_list(i).name],'EEG_network');
        
    elseif FC_freq == "AEC\8_13Hz"
        cd(path2output_AEC_8_13)
        save([EEG_list(i).name],'EEG_network');
        
    elseif FC_freq == "AEC\13_20Hz"
        cd(path2output_AEC_13_20)
        save([EEG_list(i).name],'EEG_network');
        
    elseif FC_freq == "AEC\20_32Hz"
        cd(path2output_AEC_20_32)
        save([EEG_list(i).name],'EEG_network');
        
    elseif FC_freq == "orth_AEC\1_4Hz"
        cd(path2output_orth_AEC_1_4)
        save([EEG_list(i).name],'EEG_network');
        
    elseif FC_freq == "orth_AEC\4_8Hz"
        cd(path2output_orth_AEC_4_8)
        save([EEG_list(i).name],'EEG_network');
        
    elseif FC_freq == "orth_AEC\8_13Hz"
        cd(path2output_orth_AEC_8_13)
        save([EEG_list(i).name],'EEG_network');
        
    elseif FC_freq == "orth_AEC\13_20Hz"
        cd(path2output_orth_AEC_13_20)
        save([EEG_list(i).name],'EEG_network');
        
    elseif FC_freq == "orth_AEC\20_32Hz"
        cd(path2output_orth_AEC_20_32)
        save([EEG_list(i).name],'EEG_network');
        
    end
    
end

