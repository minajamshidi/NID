function [eeg, src_sig_p, src_sig_q, pattern_sig_p, pattern_sig_q, dipoles] = ...
    mk_elec_CFC( source_analysis, param, FilterCoeff, SynchType, kappa, ...
    DiffTopo, PhaseLag_q, src_sig_p, src_sig_q)
% The function makes the simulated EEG with cross-frequency coupled sources(CFC)
% INPUT:
%       - source_analysis: the head model from METH toolbox
%       - param: includes the simulation parameters (go to Produce_EEG)
%       - FilterCoeff: filter coefficients for the two frequency bands and
%       the base frequency
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
%       - diffTopo: 0 or 1. determines whether the two sources of each pair 
%       are spatially distinct.
%       - PhaseLag: the amount of phase lag between the two sources of one
%       pair.
%       - src_sig_p and src_sig_q: the sources to be mixed with noise. if
%       you want random generation of these sources, put them as [].
% OUTPUTS:
%         - eeg: the simulated multichannel EEG: time x channel
%         - src_sig_p, src_sig_q: the coupled oscillations: time x PairsNum
%         - pattern_sig_p, pattern_sig_q: the mixing patterns of the
%         simulated coupled oscillations: channel x PairsNum
%         -dipoles: structure with the information of the dipoles
%%
% METH toolbox by Guido Nolte is used in this script
% https://www.uke.de/english/departments-institutes/institutes/neurophysiology-and-pathophysiology/research/research-groups/index.html 
% (as of 06/11/2019)
%%
% This script is a part of the supplementary material for the following paper:
% 
% -  Jamshidi Idaji et al,  "Nonlinear Interaction Decomposition (NID): A Method 
% for Separation of Cross-frequency Coupled Sources in Human Brain". doi: https://doi.org/10.1016/j.neuroimage.2020.116599>
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
% %
% *Please cite the above papers in case of 
% any significant usage of this script or the simulation pipeline.
% % 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS
% OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY
% AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER
% OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
% DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
% DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER
% IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
% OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

%%
unpack_param(param, 'param');
unpack_param(FilterCoeff, 'FilterCoeff');
%% building the mixing patterns
[pattern_sig, pattern_noise, dipoles] = mk_model_mni (source_analysis, SourceNum, NoiseNum);
if ~DiffTopo % if the sources of the two frequency come from the same location
    pattern_sig = [pattern_sig(:,1:PairsNum),pattern_sig(:,1:PairsNum)];
end

pattern_sig_p = pattern_sig(:,1:PairsNum);
pattern_sig_q = pattern_sig(:,PairsNum+1:SourceNum);

NoiseNum = size(pattern_noise, 2);
%% generate PairsNum synchronized pairs
% rng(Seed(3),'twister');
if ~numel(src_sig_p) && ~numel(src_sig_q)
    [src_sig_p, src_sig_q] = mk_sync_sig (T, PairsNum, butter_base_b, butter_base_a,butter_p_b,butter_p_a,butter_q_b,butter_q_a, FrP, FrQ, 0, PhaseLag_q,SynchType,kappa);
end
%% prepare and mix the noise
% rng(Seed(4),'twister');
src_noise = pinknoise(T,NoiseNum);
elec_noise = src_noise * pattern_noise';

elec_noise_p = filtfilt(butter_p_b,butter_p_a,elec_noise);
elec_noise_q = filtfilt(butter_q_b,butter_q_a,elec_noise);

var_elec_noise_p = sum(var(elec_noise_p));
var_elec_noise_q = sum(var(elec_noise_q));

%% adjust SNR if necessary
[src_sig_p, src_sig_q] = snr_adjust (SNR, src_sig_p, src_sig_q, pattern_sig, var_elec_noise_p, var_elec_noise_q);
%% mix the sources (Src) into channels
src_sig = [src_sig_p src_sig_q];
elec_sig = src_sig * pattern_sig';

%% Finally, we have signals in the sensors space:
eeg = elec_sig + elec_noise;

end
%% ------------------------------------------------------------------------------------------------------------------------------------------
% -------------------------------------------------------------------------------------------------------------------------------------------
%%
function unpack_param (var, str)
% Given a structure with fields, spawns variables with the same names in the 'caller' workspace
%
% INPUT:
%     var - structure variable
%     str - string name for the variable in the caller workspace
%
fns = fieldnames(var);
for i = 1:length(fns)
    evalin('caller', sprintf('%s = %s.%s;', fns{i}, str, fns{i}));
end
end
%% ------------------------------------------------------------------------------------------------------------------------------------------
% -------------------------------------------------------------------------------------------------------------------------------------------
%%
function [pattern_sig, pattern_noise, dipoles] = mk_model_mni (source_analysis, SourceNum, NoiseNum)
% INPUT:
%       source_analysis: the head model from METH toolbox
%       SourceNum: number of desired sources
%       NoiseNum: number of noise sources
% OUTPUT:
%       pattern_sig (chNxSourceNum): the mixing patterns of desired sources 
%       pattern_sig (chNxNoiseNum): the mixing patterns of noise sources 
%%
% this function is using METH toolbox
% https://www.uke.de/english/departments-institutes/institutes/neurophysiology-and-pathophysiology/research/working-groups/index.html
%%
% rng(seed(1),'twister');
[dipole_pos_3d_noise, NoiseNum] = noise_dipoles_position(source_analysis, NoiseNum);
%
%% generate random dipole orientation vectors for synced sources and noise
% rng(seed(2),'twister');
theta_phi = rand(NoiseNum+SourceNum,2);
theta = theta_phi(:,1)*2*pi - pi;
phi = theta_phi(:,2)*pi - pi/2;
[X, Y, Z] = sph2cart(theta, phi, ones(NoiseNum+SourceNum, 1));
dipole_orient_3d = [X, Y, Z];

dipole_orient_3d_sig = dipole_orient_3d(1:SourceNum,:);
dipole_orient_3d_noise = dipole_orient_3d(SourceNum+1:SourceNum+NoiseNum,:);

dipoles = struct;
dipoles.sig_ori = dipole_orient_3d_sig;
dipoles.noise_ori = dipole_orient_3d_noise;
%%

%% generate random dipole positions for synced sources and noise
% rng(seed(2),'twister');
dipole_idx = randperm(size(source_analysis.cortex.vc,1));

dipole_pos_3d_sig = source_analysis.cortex.vc(dipole_idx(1:SourceNum),:);
% dipole_pos_3d_noise = source_analysis.cortex.vc(dipole_idx(SourceNum+1:SourceNum+NoiseNum),:);
dipoles.sig_pos_idx = dipole_idx(1:SourceNum);
% dipoles.noise_pos_idx = dipole_idx(SourceNum+1:SourceNum+NoiseNum);
dipoles.sig_pos = dipole_pos_3d_sig;
dipoles.noise_pos = dipole_pos_3d_noise;
%% use METH to generate mixing patterns
pattern_sig = forward_general([dipole_pos_3d_sig dipole_orient_3d_sig],source_analysis.fp);
pattern_noise = forward_general([dipole_pos_3d_noise dipole_orient_3d_noise],source_analysis.fp);

end
%% ------------------------------------------------------------------------------------------------------------------------------------------
function [dipole_pos_3d_noise, N_dipoles_f] = noise_dipoles_position(source_analysis, N_dipoles)
N = ceil(nthroot(N_dipoles,3));
vc = source_analysis.cortex.vc;
vc_x = vc(:,1);
vc_y = vc(:,2);
vc_z = vc(:,3);
% x_range = max(vc_x) - min(vc_x);
x_bins = linspace(min(vc_x)-0.1,  max(vc_x)+0.1, N+1);
y_bins = linspace(min(vc_y)-0.1,  max(vc_y)+0.1, N+1);
z_bins = linspace(min(vc_z)-0.1,  max(vc_z)+0.1, N+1);
n = 0;
dipole_pos_3d_noise = NaN(N^3,3);
for ix = 1:N
    for iy = 1:N
        for iz = 1:N
            n = n + 1;
            IND = vc_x>=x_bins(ix) & vc_x<x_bins(ix+1) & ...
                vc_y>=y_bins(iy) & vc_y<y_bins(iy+1) & ...
                vc_z>=z_bins(iz) & vc_z<z_bins(iz+1);
            if ~sum(IND)
                dipole_pos_3d_noise(n,:) = [];
                n = n - 1;
            else
                vc_ind = vc(IND,:);
                idx_vox = randi(sum(IND), 1);
                dipole_pos_3d_noise(n, :) = vc_ind(idx_vox,:);
            end
            
        end
    end
end
N_dipoles_f = size(dipole_pos_3d_noise, 1);
end

%% ------------------------------------------------------------------------------------------------------------------------------------------
% -------------------------------------------------------------------------------------------------------------------------------------------
%%
function [src_sig_p, src_sig_q] = snr_adjust (snr, src_sig_p, src_sig_q, pattern_sig, var_elec_noise_p, var_elec_noise_q)
% INPUT:
%       snr:the target snr
%       src_sig_p,src_sig_q : source signals in the two frequency band
%       pattern_sig (chNxSourceNum): the mixing patterns of desired sources 
%       var_elec_noise_p, var_elec_noise_q: the variance of the noise in the two frequency
%       bands
%       
% OUTPUT:
%       the source signals with adjusted SNR       
%%

%
%%
PairsNum2 = size(src_sig_p,2);
SourceNum = PairsNum2*2;
T = size(src_sig_p,1);

snr_cur = NaN(1,SourceNum);
%% FrP signals
for k = 1:PairsNum2
    source = zeros(T,SourceNum);
    source(:,k) = src_sig_p(:,k);
    elec_sig_p = source * pattern_sig';
    var_elec_sig_p = var(elec_sig_p);
    snr_cur(k) = sum(var_elec_sig_p) / var_elec_noise_p;
end
%% FrQ signals
for k = 1:PairsNum2
    source = zeros(T,SourceNum);
    source(:,PairsNum2+k) = src_sig_q(:,k);
    elec_sig_q = source * pattern_sig';
    var_elec_sig_q = var(elec_sig_q);
    
    snr_cur(PairsNum2+k) = sum(var_elec_sig_q) / var_elec_noise_q;
end
%% Normalize SNRs to a given value
src_sig = [src_sig_p src_sig_q];
snr_factor = sqrt(snr_cur / snr);
src_sig = bsxfun(@rdivide,src_sig,snr_factor);
src_sig_p = src_sig(:,1:PairsNum2);
src_sig_q = src_sig(:,PairsNum2+1:SourceNum);
end
