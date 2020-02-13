function [ eeg, pattern_sig, electrode_locs ] = Produce_EEG(SourceNum, Fs,SNR, src_sig)
%% 
% The function makes the simulated EEG wuth the given focal sources.
% for each of the focal sources a dipole will be generated randomely and
% then mixed with pinknoise.
% -------------------------------------------------------------------------
% INPUT: 
%       *PairsNum: number of coupled pairs
%       *FrBase: base frequency f_b
%       *FrP and FrQ: frequency ratios 
%           --> the oscillations are at f_p = FrP*f_b and f_q = FrQ*f_b
%       *Fs : sampling frequnecy
%       *SNR: signal-to-noise ratio
%       *DiffTopo: the binary parameter specifying wether the two oscillations are
%       allowed to have different topographies: 0 (similar), 1 (random)
%      * src_sig: the given source signals
% -------------------------------------------------------------------------
% OUTPUTS:
%       *eeg: multi-channel simulated EEG (Time x channel)
%       *pattern_sig: mixing patterns of  the sources (channel x n_sources)
%       *electrode_locs: electrode locations for plotting purposes
%%
% METH toolbox by Guido Nolte is used in this script
% https://www.uke.de/english/departments-institutes/institutes/neurophysiology-and-pathophysiology/research/research-groups/index.html 
% (as of 06/11/2019)
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
% This code is based on the code used in:

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
load('sa_eeg.mat');
source_analysis = sa; clear sa;
electrode_locs = source_analysis.locs_2D;
%% Parameters
% putting all the needed parameters into one structure
param = struct;
param. Fs = Fs;
param. SNR = SNR;
param.T = length(src_sig); %number of points in the signal
param.NoiseNum = 100;% number of (pink) noise sources
param.SourceNum = SourceNum;


%% Make signals
[eeg, pattern_sig] = mk_elec_focal ( source_analysis, param, src_sig);
%%
[butter_preproc_b, butter_preproc_a] = butter(2,[0.5 50]/(Fs/2));
eeg = filtfilt (butter_preproc_b, butter_preproc_a, eeg);
end

