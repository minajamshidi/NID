%% main script for running NID on synthetic data
% this script generates simulated EEG with METH toolbox, and then applies
% NID on it, in order to extract the coupling components. Then calculate 
% the errors in Monte Carlo runs.

% Note that in this script the results are computed but there is no graph
% in the output.
% --------------------------
% there are four loops:
% 1.loop on SynchType
% 2. loop on frequencies of fast and slow oscillations,
% 3. loop on different SNRs,
% 4.loop for doing the analysis in multiple iteration, in order to have box
% plots.

%To control the number of iterations, look at variables FrP1 and FrQ1 for
%the first loop, and SNRdb for the second loop. iterN variable determines
% the number of iterations of third loop.
% --------------------------
% by changing the variable SynchType, the type of coupling can be
% controlled.

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%% other toolboxes:
% * fastICA function of fastICA toolbox is modified to include other contrast fuctions.

% * METH toolbox by Guido Nolte is used in this script
% https://www.uke.de/english/departments-institutes/institutes/neurophysiology-and-pathophysiology/research/research-groups/index.html 
% (as of 06/11/2019)

% JADE

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
% +++++++++++++++
% The simulation pipeline code is based on the code used in:

% Nikulin, V.V., Nolte, G., Curio, G., 2012. Cross-frequency decomposi-
% tion: A novel technique for studying interactions between neuronal oscilla-
% tions with different frequencies. Clinical Neurophysiology 123, 1353–1360.
% doi:10.1016/j.clinph.2011.12.004.
%%
% *Please cite the above papers in case of 
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
clc; clear; close all;
%% set path

% add the path of METH toolbox here:
METH_path = '/data/pt_02076/MATLAB_Toolboxes/meth/';

addpath('./simulate_EEG')
addpath(genpath('./NID'));
addpath(genpath(METH_path))
%% parameters
% frequencies ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
F_base = 10;
FrP1 = 1;
FrQ1 = 2;
pq_N = length(FrP1);
Fs = 500; % sampling frequency
T_sec = 5*60; % length of each signal segment in seconds.

% Coupling Properties ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SynchType = 1; % parameter setting the type of coupling
kappa1 = 10;
kappa = kappa1*0.1; % variable controling the strength of coupling

%       - SynchType: the type of synchronization 
%               0: same amplitude, 100% locked phase
%               1: independent amplitude, 100% locked phase
%               2: independent amplitude, locked phase (strength of synch controled with kappa von mesis)
                % for this option you will need CircStat toolbox by Philipp Berens & Marc J. Velasco
                % https://de.mathworks.com/matlabcentral/fileexchange/10676-circular-statistics-toolbox-directional-statistics
                % (as of 06/11/2019)
%                 kappa is the kappa parameter in von mises distribution.
%                 larger kappa -> more peaky von mises -> less jitter -> stronger coupling
%               3: Correlated amplitudes, independent phases (strength of synch controled with kappa)
%                  A2 = kappa*A1 + (1-kappa)*(white Noise)
%                  Thus, kappa larger (near 1) means that the coupling strength is higher
%               4: phase-amplitude coupling
%               5: 2kinds, 1st phase locked, 2nd amp locked
%               6: independent amplitude, independent phases
%               7: nothing. only noise
%       - kappa: the parameter related to the strength of synchronization. 
%       (only used for SynchType=2,3, otherwise put it as 0). 

DiffTopo = 1; % indicates if the topographies of the sources should be generated randomly (1) or should be the same (0).
PhaseLag_q = 0; %phase lag of source at fast frequency

% number of pairs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
PairsNum = 2; % number of of synchronized sources
NoiseNum = 125; % number of noise sources
% SNR ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SNRdb = -5; % in db
SNR1 = 10.^(SNRdb/10);
SNR_N = length(SNRdb);

% Save and Plot ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
DoSave = 0; % save?
save_path = './'; %where?

% Iteration settings ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
iterN = 100; % number of iterations

tic
%%
for ind_pq = 1: pq_N % loop on fp and fq
    
    % Settings for each (fp,fq) couple ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    FrP = FrP1(ind_pq);
    FrQ = FrQ1(ind_pq);
    fp = FrP*F_base; % slow frequency
    fq = FrQ*F_base; % fast frequency
    rng('shuffle')
    seed1 = randperm(9999999,4*iterN);
    seed1 = reshape(seed1,4,iterN);
    
    
    for ind_snr = 1:SNR_N % loop on SNR
        
        % Settings for each SNR value ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if SNRdb(ind_snr)<0 % this is for the purpose of save name
            MinPos = 'minus';
        else
            MinPos = 'pos';
        end
        ep_NID = NaN(iterN,PairsNum); % errors of slow oscillations
        eq_NID = NaN(iterN,PairsNum); % errors of fast oscillations
        ep_SSD = NaN(iterN,PairsNum);
        eq_SSD = NaN(iterN,PairsNum);
        D_synch = NaN(iterN,PairsNum);
        SynchFac_mean = NaN(iterN,1);
        ICA_Method = cell(iterN,1);
        dipoles = cell(iterN,1);
        n_catch = 0;
        for ind_iter = 1:iterN
            try
                fprintf('***** SNR=%d (db) **** iter num.=%d\n',SNRdb(ind_snr),ind_iter)
                %% make synth. signal
                disp('Making the signals....');
                rng(seed1(ind_iter),'twister')% we make the same signal for different SNRs, just modify the signal powers with regards to SNR
                
                [X, src_p, src_q, pattern_p, pattern_q, electrode_locs, dipoles{ind_iter},~, ~] = ...
                    Produce_EEG_CFC(PairsNum, NoiseNum, F_base, FrP, FrQ, Fs,SNR1(ind_snr),DiffTopo,...
                    SynchType,kappa,PhaseLag_q,T_sec, [], []);
                % ------------------------------
                % X is Time*channel
                % X is bandpass-filtered in [0.5,50]
                % ------------------------------
                
                % ----------------------------------
                % to plot the PSD of the signal, you can use pwelch or
                % plot_sepec from https://github.com/minajamshidi/EEG-Functions-MATLAB
                % figure
                % plot_spec(X, Fs, 'frequency_mark', [10,20], 'freq_res', 0.5, 'noverlap', Fs, 'f_max',50)
                % set(gca, 'fontweight','bold', 'fontsize', 14)
                
                % ----------------------------------
                % to plot the mixing patterns, use this function of METH:
                % figure, showfield(pattern_p(:,1), electrode_locs);
                
                % ----------------------------------
                % Test the synchronization coefficients: 
                % PLValue = Phase_Locking(src_p,src_q,FrP,FrQ)
                % CorrAbs = corr(abs(hilbert(src_p)),abs(hilbert(src_q)))
                % Corr_PAC = corr(abs(angle(hilbert(src_p))),abs(hilbert(src_q)))

                
                %% NID
                [Success,A_final_p,A_final_q,out_ICA,ICA_Method{ind_iter},D_synch(ind_iter,:),SynchFac_mean(ind_iter),out_ssd] = ...
                    NID(X,FrP,FrQ,F_base,'source_set_num',PairsNum,'ssd_or_bp',[1 1],'fs',Fs,'synchtype',SynchType);
                
                
                Error_p = Pattern_Angle( pattern_p, A_final_p);
                Error_q = Pattern_Angle(pattern_q, A_final_q);
                % --------
                
                Perm1 = perms(1:PairsNum); % all permutations
                
                P1 = Error_p; % errors of selected components at fp
                Q1 = Error_q; % errors of selected components at fq
                PQ1 = P1.^2 + Q1.^2;
                Er = NaN(size(Perm1,1),1);
                for k = 1:size(Perm1,1) %for all permutations
                    E = 0;
                    for j = 1:PairsNum
                        E = E + PQ1(Perm1(k,j),j); %sum of square errors
                    end
                    Er(k) = E/PairsNum; %MSE
                end
                [~,IDX] = min(Er); % permutation with min MSE
                IDX_patt2 = Perm1(IDX,:);
                
                % save the errors of selected components ~~~~~~~~~~~~~~~~~~~~~~~~~~
                % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                idx_error = sub2ind(size(Error_p), IDX_patt2(1:PairsNum), 1:PairsNum);
                ep_NID(ind_iter,:) = Error_p(idx_error);
                eq_NID(ind_iter,:) = Error_q(idx_error);
                
                             
                
                %% plot the topo
                disp('~~~~~~~~');
                for j = 1:PairsNum
                    disp(['Source P, Pair ',num2str(j),', error=',num2str(Error_p(IDX_patt2(j),j))])
                    disp(['Source Q, Pair ',num2str(j),', error=',num2str(Error_q(IDX_patt2(j),j))])
                    
                end
                %%
                disp('~~~~~~~~');
            catch
                n_catch = n_catch + 1;
            end
        end %for ind_iter
        
        %% save
        if DoSave
            Name = sprintf('Data_iterN_%d_fp_%d_fq_%d_SynchType_%d_%dpairs_Fs_%d_SNR_%s%d_Tsec_%d',...
                iterN,fp,fq,SynchType,PairsNum,Fs,MinPos,abs(SNRdb(ind_snr)),T_sec);
            Name_save = [save_path,Name,'_Result','.mat'];
            save(Name_save,'ep_NID','eq_NID','seed1',...
                'SynchFac_mean','D_synch','ICA_Method', 'dipoles');
        end
    end % end ind_snr
end % ind_pq
toc
