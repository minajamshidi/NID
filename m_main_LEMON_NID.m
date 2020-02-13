%% main script for extraction of coupled paris of oscillation from LEMON data
% This scripts runs on resting-state EEG of young, right-handed men of
% LEMON dataset. You should download the data beforehand.
%
% LEMON paper:
% Babayan Anahit, Erbey Miray, Kumral Deniz, et al., Scientific Data 6, 180308 (2019).
%
% 
%%
% This script is a part of the supplementary material for the following paper:
% 
% -  Jamshidi Idaji et al,  "Nonlinear Interaction Decomposition (NID): A Method 
% for Separation of Cross-frequency Coupled Sources in Human Brain". doi: <https://doi.org/10.1016/j.neuroimage.2020.116599>
% 
% 
% (C) Mina Jamshidi Idaji, Oct 2019, @ MPI CBS, Leipzig, Germany
% https://github.com/minajamshidi/NID
% minajamshidi91@gmail.com
%%
% *Please cite the above paper in case of 
% any significant usage of this script or the simulation pipeline.
%% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS
% OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY
% AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER
% OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
% DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
% DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER
% IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
% OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%%

%%
clc; clear; close all;
%% set path
addpath('./NID'); % path of NID function

% add the path of EEGLAB

%the path of preprocessed data
path_to_data1 = '/data/pt_nro109/Share/EEG_MPILMBB_LEMON/EEG_Preprocessed_BIDS_ID/EEG_Preprocessed/';

% add the path+name of the metafile including subjects information
metafile_path = '/data/pt_nro109/Share/Behavioural_Data_MPILMBB_LEMON/META_File_IDs_Age_Gender_Education_Drug_Smoke_SKID_LEMON.csv';

path_save = './save_data';

%startup_bbci_toolbox;
%% parameters related to NID

% frequencies ~~~~~~~~~~~~~~~~~~~~~~~~~~
F_base = 10;
FrP = 1;
FrQ = 2;
fp = FrP*F_base;
fq = FrQ*F_base;
Fs = 250; % sampling frequency
SynchType = 1;

% number of pairs ~~~~~~~~~~~~~~~~~~~~~~~~~~
PairsNum1 = 5; % # of synchronized sources
Ncomp = PairsNum1; % # of components for SSD
Ncomp_ICA = 2*Ncomp; % # of components for ICA

% SSD options ~~~~~~~~~~~~~~~~~~~~~~~~~~
opt_SSD = struct;
opt_SSD.DoSSD = [1 1];
opt_SSD.Ncomp = Ncomp;
opt_SSD.fc = [fp fq];

%  ICA options ~~~~~~~~~~~~~~~~~~~~~~~~~~
opt_ICA = struct;
opt_ICA.ICA_Method = 'fastICA_E5';%'JADE'; %or 'fastICA'
opt_ICA.Ncomp_ICA = Ncomp_ICA;
%% Read data Info

% Here I first find the right handed young subjects
% -------------------------------------------------------------------------
T =  readtable(metafile_path); % the IDs are in this csv file
ID = table2array(T(:,1)); % IDs
hand = table2array(T(:,4)); % hand info
age = table2array(T(:,3)); % age info
gender = table2array(T(:,2)); %gender info

ind_young_male_right = find(strcmp(age,'20-25') | strcmp(age,'25-30') & ...
     ismember(gender,1) &  strcmp(hand,'right')); % Ind of young_male_righ-handed subjects
SubjN_young_male_right = numel(ind_young_male_right);%number of young_male_righ-handed subjects

subj_ID = ID(ind_young_male_right); % IDs of right handed young subjects

% -------------------------------------------------------------------------
%Then look at the available preprocessed subjects' data
% -------------------------------------------------------------------------
names = dir([path_to_data1,'*_EC.set']);
names = {names.name}; % names of all files in preprocessed file
SubjN_preproc = length(names);%number of available preprocessed subjects
ID_read = cell(SubjN_preproc,1);%ID of available preprocessed subjects
for subj = 1:SubjN_preproc
    ID_read{subj} = names{subj}(1:end-7);
end
%% load data and process it
ID_computed = intersect(ID_read,subj_ID);
SubjN = length(ID_computed);%# of young_male_righ-handed subjects, for whom we have the preprocessed data
fprintf('~~~~~~~~~~~~~~~~~~~~~~~~~\n~~~~~~~~~~~~~~~~~~~~~~~~~\n');
fprintf('There are totally %d young male right-handed subjects\n',SubjN);
fprintf('---------------------------\n---------------------------\n---------------------------\n');
PermNum = 1000; %# of permutation iterations for sattistical testing of synchronization

% -------------------------------------------------------------------------
% Build information containers:D:P
% -------------------------------------------------------------------------
load('LEMON_chanlocs61.mat'); %this chanlocs file is saved from a subject with 61 channels
CC_All = NaN(PairsNum1,SubjN);
PLV_sig = cell(SubjN,1);

% filter coefficients -------------------------------------------------
[butter_base_b,butter_base_a] = butter(2,[F_base-1 F_base+1]/(Fs/2));
[butter_p_b,butter_p_a] = butter(2,[FrP*F_base-1 FrP*F_base+1]/(Fs/2));
[butter_q_b,butter_q_a] = butter(2,[FrQ*F_base-1 FrQ*F_base+1]/(Fs/2));
    
for subj = 1:SubjN
    ID_subj = ID_computed{subj};
    fprintf('~~~~~~~~~~~~~~~~~~~~~~~~~\n~~~~~~~~~~~~~~~~~~~~~~~~~\n');
    fprintf('Subject number %d/%d with ID %s\n',subj,SubjN,ID_subj);
    fprintf('~~~~~~~~~~~~~~~~~~~~~~~~~\n');
    % load the data of current subject ------------------------
    SubjName = [ID_computed{subj},'_EC.set'];
    file_name = sprintf([path_to_data1,'/%s'], SubjName);
    eeg = pop_loadset('filename', file_name);
    
    % if any channel is missing: interpolate the channel --------------
    
    [~,bad_elec] = setdiff({chanlocs61.labels},{eeg.chanlocs.labels}); % check the existance of a bad channel
    if numel(bad_elec) % if there is a bad channel
        EEG = pop_interp(eeg, chanlocs61, 'spherical'); % interpolate it
    else
        EEG = eeg;
    end
    
    % set variables for permutation ---------------------------------
    
    % in the permutation part, we segment the signal and permute the
    % segments. Therefore there should be integer number of segments. Here
    % we remove the part of X that cannot be in a complete segment.
    
    X = double(EEG.data)';
    L = size(X,1); % the length of the signal
    Ls = 1*Fs; % 1s segment length
    Ns = floor(L/Ls); % # of segments
    X = X(1:Ns*Ls,:); % remove the parts that cannnot be segmented
    
    
    %% METHOD   
   
    % X -> Time x Channel
    [out_ssd] = SSD_BP(X,Fs,opt_SSD,1,[]);
    
    %perm------------------
    %the first iteration is non-permuted
    CC_pairsNum = NaN(PairsNum1,PermNum+1);
    
    ind_perm = 0;
    while ind_perm <= PermNum
        fprintf('***Subject %d/%d, perm %d try*****\n',subj,SubjN,ind_perm+1);
        X_p_ssd = out_ssd.X_NB{1};
        X_q_ssd = out_ssd.X_NB{2};
        if ind_perm>0 % the first iteration is the unpermuted data results
            p2 = randperm(Ns);
            for k = 1:PairsNum1
                x2 = reshape(X_p_ssd(:,k),Ls,Ns);
                y2 = reshape(X_q_ssd(:,k),Ls,Ns);
                y2 = y2(:,p2);
                X_q_ssd(:,k) = y2(:);
            end
        end
        out_ssd.X_Aug = [X_p_ssd'; X_q_ssd'];
        out_ssd.X_NB{2} = X_q_ssd;
        
        try
                        
           [Success,A_final_p,A_final_q,out_ICA,ICA_Method,D_synch,SynchFac_mean,out_ssd,IDX_src,Xp,Xq] = ...
                NID(X,FrP,FrQ,F_base,'source_set_num',PairsNum1,'ssd',out_ssd,'fs',Fs,'synchtype',1);
            if Success
                ind_perm = ind_perm + 1;
                fprintf('***Subject %d/%d, perm %d succeeded*****\n',subj,SubjN,ind_perm);
                CC_pairsNum(:,ind_perm) = D_synch';
                A_final1 = cell(1,2);
                A_final1{1} = A_final_p;
                A_final1{2} = A_final_q;
                %~~~~~~~~~~~~~~~~
                if ind_perm==1
                    out_ssd1 = out_ssd;
                    out_ICA1 = out_ICA;
                    SynchFac_pq1 = D_synch;
                    Xp1 = Xp;
                    Xq1 = Xq;
                    IDX_src1 = IDX_src;
                    ICA_Method1 = ICA_Method;
                    CC_All(:,subj) = D_synch';
                end
                
            end
            
            
        catch
            fprintf('ICA failed in the first iter\n');
        end
    end % end ind_perm
   
    % finding the significant pairs ---------------------------------------
    PLV1 = CC_pairsNum(:,1);
    PLV_perm = max(CC_pairsNum(:,2:end));
    Perc = NaN(5,1);
    for p = 1:5
        Perc(p) = sum(PLV_perm>=PLV1(p))/(PermNum/100);
    end
    ind1 = find(Perc<5); % significance level 0.05
    PLV_sig{subj} = PLV1(ind1);
    
    % Saving --------------------------------------------------------------
    %Name_save = sprintf('%s/%s_NID_fp_%d_fq_%d_perm.mat',path_save,ID_subj,FrP*F_base,FrQ*F_base);
    %save(Name_save, 'out_ssd1','out_ICA1','A_final1',...
    %    'SynchFac_pq1','Xp1','Xq1','CC_pairsNum','ICA_Method1');
   
    
    
    %%
end
%%

