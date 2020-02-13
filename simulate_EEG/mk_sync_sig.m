function [src_sig_p, src_sig_q] = ...
    mk_sync_sig (T, PairsNum, butter_base_b, butter_base_a,...
    butter_p_b,butter_p_a, butter_q_b,butter_q_a,FrP, FrQ,...
    Delay_p, Delay_q,SynchType,kappa)
%% 
% Function for producing synchronized oscillations.
% INPUT:
%       -T: number of time samples
%       - PairsNum: number of pairs of coupled oscillations
%       - butter_base_... : filter coefficients for base frequency
%       - butter_p_... : filter coefficients for slower oscillation
%       - butter_p_... : filter coefficients for faster oscillation
%       - FrP and FrQ: frequency ratios
%       - Delay_p, Delay_q (rad): phase lags of slow and fast oscillations
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
s1 = filtfilt(butter_base_b,butter_base_a,randn(T,PairsNum));
s1 = bsxfun(@rdivide,s1,std(s1,0,1));
sp = filtfilt(butter_p_b,butter_p_a,randn(T,PairsNum));  % for using its envelope
sq = filtfilt(butter_q_b,butter_q_a,randn(T,PairsNum)); % for using its envelope
sp_abs = abs(hilbert(sp)); % amplitude envelope of a signal at fp
sq_abs = abs(hilbert(sq)); % amplitude envelope of a signal at fq

s1_cplx = hilbert(s1);
sp_abs = abs(s1_cplx); % ****
src_sig_p = real(sp_abs.*exp(1i.*(FrP*angle(s1_cplx) + Delay_p)));
switch SynchType
    case 0 % 100% phase-phase coupling and amp-amp coupling
        sq_abs = sp_abs;
        src_sig_q = real(sq_abs.*exp(1i.*(FrQ*angle(s1_cplx) + Delay_q)));
    case 1 % 100% phase-phase + independent amps
        %         abs_q = sq_abs;% ****
        abs_q = flip(sp_abs);% ****
        src_sig_q = real(abs_q.*exp(1i.*(FrQ*angle(s1_cplx) + Delay_q)));
    case 2 % von mises controled phase-phase + indep amps
        dPhi = NaN(T,PairsNum);
        for k = 1:PairsNum
            dPhi(:,k) = circ_vmrnd(0, kappa, T);
        end
        abs_q = sq_abs;
        src_sig_q = real(abs_q.*exp(1i.*(FrQ*angle(s1_cplx) + dPhi)));
        src_sig_q = filtfilt(butter_q_b,butter_q_a,src_sig_q);
        src_sig_q = real(abs_q.*exp(1i.*angle(hilbert(src_sig_q))));
    case 3 % amp-amp coupling + phase independent
        src_sig_q = filtfilt(butter_q_b,butter_q_a,randn(T,PairsNum));
        Phi_q = angle(hilbert(src_sig_q));
        abs_q = kappa*sp_abs + (1-kappa)*randn(T,PairsNum);
        src_sig_q = real(abs_q.*exp(1i.*(Phi_q)));
        src_sig_q = filtfilt(butter_q_b,butter_q_a,src_sig_q);
    case 4 % phase-amp
        abs_q = abs(angle(hilbert(src_sig_p)));
        Phi_q = angle(hilbert(sq));
        src_sig_q = real(abs_q.*exp(1i.*(Phi_q)));
    case 5 %indep amp + indep phase
        src_sig_q = sq;
    case 6 % nothing
        src_sig_p = zeros(size(src_sig_p));
        src_sig_q = zeros(size(src_sig_p));
end

end
