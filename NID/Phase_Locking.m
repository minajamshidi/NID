function Synch_ind = Phase_Locking(x1,x2,p,q)
% computes the phase locking value of two signals
% INPUTs: 
%       *x1, x2 [Time x N]: the signals for which the plv should be calculated
%       * p, q: the frequency ratios (x1 and x2 are p:q coupled) 
% OUTPUT:
%        * Synch_ind [N x N]: a matrix, for which the element (i,j) is the plv of
%       x1(:,i) and x2(:,j)
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% this code is part of the PhD thesis of Mina Jamshidi Idaji
% Neurology Dept @ MPI CBS, Leipzig, Germany + Machine learning group @ TU
% Berlin, Germany
% reference paper:
%  .....
%
% In case you use this script as a significant part of your coding, please
% cite the paper appropriately.
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS 
% OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY
% AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER 
% OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
% DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
% DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER 
% IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT 
% OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%
phi_p = unwrap(angle(hilbert(x1)));
phi_q = unwrap(angle(hilbert(x2)));

N = size(x1,2);
if N>1
    Psi_pq = cell(N,N);
%     Psi_peak = NaN(N,N);
    Synch_ind = zeros(N);
    for k = 1: N
        for j = 1:N
            Psi_pq{k,j} = mod(q*phi_p(:,k) - p*phi_q(:,j),2*pi);
            Synch_ind(k,j) = abs(mean(exp(1i*Psi_pq{k,j})));

        end
    end

else
    Psi_pq = mod(q*phi_p - p*phi_q,2*pi);
    Synch_ind = abs(mean(exp(1i*Psi_pq)));

end
end

