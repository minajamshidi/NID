function Synch_ind = Phase_Locking(x1,x2,p,q)
% computes the phase locking value of two signals
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
%%
% *Please cite the above paper (or a future peer-reviewed version) in case of 
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
phi_p = unwrap(angle(hilbert(x1)));
phi_q = unwrap(angle(hilbert(x2)));

N1 = size(x1,2);
N2 = size(x2,2);
if N1>1 && N2>1
    Psi_pq = cell(N1,N2);
%     Psi_peak = NaN(N,N);
    Synch_ind = zeros(N1, N2);
    for k = 1: N1
        for j = 1:N2
            Psi_pq{k,j} = mod(q*phi_p(:,k) - p*phi_q(:,j),2*pi);
            Synch_ind(k,j) = abs(mean(exp(1i*Psi_pq{k,j})));

        end
    end

else
    Psi_pq = mod(q*phi_p - p*phi_q,2*pi);
    Synch_ind = abs(mean(exp(1i*Psi_pq)));

end
end

