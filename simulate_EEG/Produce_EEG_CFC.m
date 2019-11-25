function [ eeg, src_sig_p, src_sig_q, pattern_sig_p, pattern_sig_q, ...
    electrode_locs, dipoles, FilterCoeff, source_analysis ] = ...
    Produce_EEG_CFC(PairsNum, NoiseNum, FrBase, FrP, FrQ, Fs,SNR,DiffTopo,...
    SynchType,kappa,PhaseLag_q,T_sec, src_sig_p, src_sig_q)
%% help
% -------------------------------------------------------------------------
% INPUT: 
%       *PairsNum: number of coupled pairs
%       * NoiseNum: number of noise sources
%       *FrBase: base frequency f_b
%       *FrP and FrQ: frequency ratios 
%           --> the oscillations are at f_p = FrP*f_b and f_q = FrQ*f_b
%       *Fs : sampling frequnecy
%       *SNR: signal-to-noise ratio
%       *DiffTopo: the binary parameter specifying wether the two oscillations are
%       allowed to have different topographies: 0 (similar), 1 (random)
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
%           
%       *PhaseLag_q: phase lag of the singla at f_q
%       *T_sec: the length of the simulated EEG in seconds
%       *Seed: a 1x3 matrix containing the seeds used for the random number
%           generator. This is used if you wanna reproduce your results. In case you
%           don't need Seed, just pass random numbers to it, e.g. randperm(999999,3).
% OUTPUTS:
%       *eeg: multi-channel simulated EEG (Time x channel)
%       *src_sig_p and src_sig_q: simulated oscillation at f_p and fq (T x PairsNum)
%       *pattern_sig_p and pattern_sig_q: mixing patterns of oscillations at 
%       f_p and f_q (channel x PairsNum)
%       *electrode_locs: electrode locations for plotting purposes
%       * dipoles: structure with the information about the dipoles
%       * FilterCoeff: structure containing the filter coefficients
%       * source_analysis: the structure containing the information of the
%       head model
%%
% METH toolbox by Guido Nolte is used in this script
% https://www.uke.de/english/departments-institutes/institutes/neurophysiology-and-pathophysiology/research/research-groups/index.html 
% (as of 06/11/2019)
%%
% This script is a part of the supplementary material for the following paper:
% 
% -  Jamshidi Idaji et al,  "Nonlinear Interaction Decomposition (NID): A Method 
% for Separation of Cross-frequency Coupled Sources in Human Brain". doi: <https://doi.org/10.1101/680397 
% https://doi.org/10.1101/680397>
% 
% 
% (C) Mina Jamshidi Idaji, Oct 2019, @ MPI CBS, Leipzig, Germany
% https://github.com/minajamshidi/NID
% minajamshidi91@gmail.com
% +++++++++++++++
% This code is based on the code used in:
% 
% Nikulin, V.V., Nolte, G., Curio, G., 2012. Cross-frequency decomposi-
% tion: A novel technique for studying interactions between neuronal oscilla-
% tions with different frequencies. Clinical Neurophysiology 123, 1353–1360.
% doi:10.1016/j.clinph.2011.12.004.
%%
% *Please cite the above papers (or a future peer-reviewed version) in case of 
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
load('sa_eeg.mat');
% load('/data/pt_02076/MATLAB_Toolboxes/Guido_toolbox/sa_montreal.mat')
source_analysis = sa; clear sa;
electrode_locs = source_analysis.locs_2D;
%% Parameters
% putting all the needed parameters into one structure
param = struct;
param.PairsNum = PairsNum;
param.FrBase = FrBase;
param. FrP = FrP;
param. FrQ = FrQ;
param. Fs = Fs;
param. SNR = SNR;
param.T = Fs * T_sec; %number of points in the signal
param.NoiseNum = NoiseNum;% number of (pink) noise sources
param.SourceNum = PairsNum * 2;

%% create band-pass filters
FilterCoeff = struct;
[FilterCoeff.butter_base_b,FilterCoeff.butter_base_a] = butter(2,[FrBase-1 FrBase+1]/(Fs/2));          
[FilterCoeff.butter_p_b,FilterCoeff.butter_p_a] = butter(2,[FrP*FrBase-1 FrP*FrBase+1]/(Fs/2));    
[FilterCoeff.butter_q_b,FilterCoeff.butter_q_a] = butter(2,[FrQ*FrBase-1 FrQ*FrBase+1]/(Fs/2));    

%% Make signals
[eeg, src_sig_p, src_sig_q, pattern_sig_p, pattern_sig_q, dipoles] = ...
        mk_elec_CFC (source_analysis, param, FilterCoeff, SynchType, kappa, ...
        DiffTopo, PhaseLag_q, src_sig_p, src_sig_q);
       
%%
[butter_preproc_b, butter_preproc_a] = butter(2,[0.5 50]/(Fs/2));
eeg = filtfilt (butter_preproc_b, butter_preproc_a, eeg);
end

