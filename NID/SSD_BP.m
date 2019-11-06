function [out] = SSD_BP(X,Fs,opt,DoAugment,varargin)
% extracts the SSD components of input multichannel signal
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
% !!!!!! bandpass mode does not work yet
%%
%X must be (time x channels)
% default values
Prune = [];

% read varargin
if ~isempty(varargin)
    for j = 1:2:(length(varargin)-1)
        switch lower (varargin{j})
            case 'prune'
                Prune = varargin{j+1};
            case 'filter_b'
                butter_b = varargin{j+1};
            case 'filter_a'
                butter_a = varargin{j+1};
        end
    end
end
%%
out = struct;

DoSSD = opt.DoSSD ;
Ncomp = opt.Ncomp;
fc = opt.fc;
out.DoSSD = DoSSD;
out.Ncomp = Ncomp;

Nf = length(fc);
A_ssd = cell(Nf,1);
W_ssd = cell(Nf,1);
X_NB = cell(Nf,1);
T = size(X,1);
if DoAugment
    X_Aug = NaN(Nf*Ncomp,T);
end
for n = 1:Nf
    if DoSSD(n)
        fpass = 1*(fc(n)/10);
        fNoise1 = fpass +2;
        fNoise2 = fpass+1;
        
        FiltFreq = [fc(n)- fpass fc(n)+fpass;
            fc(n)-fNoise1 fc(n)+fNoise1;
            fc(n)-fNoise2 fc(n)+fNoise2];

        [W_ssd{n}, A_ssd{n}, ~, ~, X_NB{n}] = ssd2(X, FiltFreq, Fs, 2, Prune,Ncomp);
        if DoAugment
            X_Aug((n-1)*Ncomp+1:n*Ncomp,:) = X_NB{n}';
        end
    else
         X_NB{n} = filtfilt(butter_b{n},butter_a{n},X); % Time*channel
    end
end

if DoAugment
    out.X_Aug = X_Aug;
end
out.X_NB = X_NB;
out.A_ssd = A_ssd;
out.W_ssd = W_ssd;

end

