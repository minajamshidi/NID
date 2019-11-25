function [eeg, pattern_sig] = mk_elec_focal ( source_analysis, param, src_sig)
% The function makes the simulated EEG wuth the given focal sources.
% for each of the focal sources a dipole will be generated randomely and
% then mixed with pinknoise.
% INPUT:
%       - source_analysis: the head model from METH toolbox
%       - param: includes the simulation parameters (go to Produce_EEG)
%       - FilterCoeff: filter coefficients for the two frequency bands and
%       the base frequency
%       - src_sig: the the focal sources to be mixed with noise    
% OUTPUTS:
%       - eeg: the simulated multichannel EEG: time x channel
%       - pattern_sig: the mixing pattern of the focal sources. 
%       channel x n_sources
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
unpack_param(param, 'param');
%% building the mixing patterns
[pattern_sig, pattern_noise] = mk_model_mni (source_analysis, SourceNum, NoiseNum);

%% prepare and mix the noise
src_noise = pinknoise(T,NoiseNum);
elec_noise = src_noise * pattern_noise';


var_elec_noise = sum(var(elec_noise));

%% adjust SNR if necessary
src_sig = snr_adjust_broadband (SNR, src_sig, pattern_sig, var_elec_noise);
%% mix the sources (Src) into channels
elec_sig = src_sig * pattern_sig';

%% Finally, we have signals in the sensors space:
eeg = elec_sig + elec_noise;

%% END of PREPARING the DATA (real or simulated)


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
function [pattern_sig, pattern_noise] = mk_model_mni (source_analysis, SourceNum, NoiseNum)
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

%
%% generate random dipole orientation vectors for synced sources and noise
dipole_orient_3d = randn(NoiseNum+SourceNum,3);  
dipole_orient_3d = bsxfun(@rdivide,dipole_orient_3d,sqrt(sum(dipole_orient_3d.^2,2))); 

dipole_orient_3d_sig = dipole_orient_3d(1:SourceNum,:);
dipole_orient_3d_noise = dipole_orient_3d(SourceNum+1:SourceNum+NoiseNum,:);

%% generate random dipole positions for synced sources and noise
dipole_idx = randperm(size(source_analysis.cortex.vc,1));

dipole_pos_3d_sig = source_analysis.cortex.vc(dipole_idx(1:SourceNum),:);
dipole_pos_3d_noise = source_analysis.cortex.vc(dipole_idx(SourceNum+1:SourceNum+NoiseNum),:);

%% use METH to generate mixing patterns
pattern_sig = forward_general([dipole_pos_3d_sig dipole_orient_3d_sig],source_analysis.fp);
pattern_noise = forward_general([dipole_pos_3d_noise dipole_orient_3d_noise],source_analysis.fp);

end
%% ------------------------------------------------------------------------------------------------------------------------------------------
% -------------------------------------------------------------------------------------------------------------------------------------------
%%
function src_sig = snr_adjust_broadband (snr, src_sig, pattern_sig, var_elec_noise)
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
SourceNum = size(src_sig,2);
snr_cur = NaN(1,SourceNum);
%% FrP signals
for k = 1:SourceNum
    source = zeros(size(src_sig));
    source(:,k) = src_sig(:,k);
    elec_sig_p = source * pattern_sig';
    var_elec_sig = var(elec_sig_p);
    snr_cur(k) = sum(var_elec_sig) / var_elec_noise;
end

%% Normalize SNRs to a given value
snr_factor = sqrt(snr_cur / snr);
src_sig = bsxfun(@rdivide,src_sig,snr_factor);
end
