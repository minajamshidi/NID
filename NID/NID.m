function [Success,A_final_p,A_final_q,out_NGMD,ICA_Method,D_synch,SynchFac_mean,out_ssd,IDX_src,Xp,Xq] = ...
    NID(X,FrP,FrQ,F_base,varargin)
% The function implimenting nonlinear interaction decomposition

% INPUT
%       * X [time x channel]: multi-channel signal to be decomposed
%       *  FrP, FrQ (integer): frequency ratios
%       * F_base (Hz): f_b ~~> f1 = FrP*f_b and f2 = FrQ*f_b
% Inputs in varargin (those marked with ** are mandatory):
%           ** 'source_set_num' (integer): number of sets of coupled sources --> PairsNum
%           * 'ssd_or_bp'[1 x 2](binary): each element = 1 if SSD should be applied
%           in that frequency band
%           ** 'fs' (Hz): sampling frequnecy 
%           ** 'synchtype': type of synchronization
%               = 1: phase coupling
%               = 3: amplitude-amplitude coupling
%               = 4: phase-amplitude coupling
%           * 'ngmd_method' ('E5' or 'jade'): one can decide beforehand which method 
%           of non-gaussianity maximization to use
%           * 'ssd': computed ssd struction. if this option is used SSD is not
%           applied and the input is used. SSD structures's fields are:
%                   DoSSD: [1 1]
%                   Ncomp: # of ssd components in each frequency
%                   X_Aug: [2Ncomp×T double]
%                   X_NB: {2×1 cell}: including the ssd component in both
%                   frequencies
%                   A_ssd: {2×1 cell}: including the mixing matrices of the
%                   two frequencies.
%           * 'prune': prune variabe of EEG in X
%
% OUTPUTs:
%       * Success: 0/1: if one of NGMD methods converged
%       * A_final_p and A_final_q [channel x PairsNum]: mixing patterns of extracted 
%       sources. column i of each is coupled to the same column of the other.
%       * out_NGMD: output of NGMD: a structure including: 
%           * X_Aug_ICA: [2PairsNum × Time]: the output of NGMD
%           * A_ICA: [2PairsNum×2PairsNum]: the mixing patterns of ICA
%           * W_ICA: [2PairsNum×2PairsNum]: the filters of ICA
%       * ICA_Method: the method of NGMD that is finally selected
%       * D_synch: the synchronization indeces
%       * SynchFac_mean: the mean of synchronization indeces
%       * out_ssd: SSD structure
%       * IDX_src: index of the selected ICA components
%       * Xp,Xq: the extracted oscillations
%
%Example:
%   [Success,A_final_p,A_final_q,out_ICA,ICA_Method,D_synch,SynchFac_mean,out_ssd,~,Xp,Xq] = ...
%                NID(X,FrP,FrQ,F_base,'source_set_num',PairsNum1,'ssd',out_ssd,'fs',Fs,'synchtype',1);
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Note that two methods are applied for NGMD, namely JADE (implimeted by Cardoso)
% and the contrast function introduced in the paper (we call it E5). Both
% the methods are appliied on the augmented matrix and then we decide for
% the method which has maximum non-gaussianity based on negentropy. 
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Note that this function is general in the sense that it can be used for
% NID in its generalized case (coupling in >2 frequency bands). However it
% is not tested extensively for possible errors.
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
%%%% default values
Prune = [];
SSD_out = 0;
E5_flag = 1;
JADE_flag = 1;
DoSSD = [1 1];
butter_a = [];
butter_b = [];
%% input check ------------------------------------
% ------------------------------------------------

% X 
if size(X,1)<size(X,2)
    warning('NID accepts raw data as time*channel arrangement. It seems the raw data is not this format. Therefore it is transposed')
    X = X';
end

if (rem(length(varargin),2) == 1)
    error('Optional parameters should always go by pairs');
else
    for j = 1:2:(length(varargin)-1)
        if ~ischar (varargin{j})
            error (['Unknown type of optional parameter name (parameter' ...
                ' names must be strings).']);
        end
        switch lower (varargin{j})
            case 'source_set_num'
                SourceSetNum = varargin{j+1};
                if ischar (SourceSetNum)
                    error('The number of ICA components should be a number, not a character.')
                else
                    if SourceSetNum<=0 || floor(SourceSetNum)~=SourceSetNum
                        error('The number of source sets must be positive integer.')
                    end
                end
                Ncomp_ICA = 2*SourceSetNum;
            case 'ssd_or_bp'
                DoSSD = varargin{j+1};
                if ~sum(ismember(DoSSD,1)) && ~sum(ismember(DoSSD,0))
                    error('The SSD_or_BP variable in accepts only binary values (0 or 1).');
                end
            case 'fs'
                Fs = varargin{j+1};     
            case 'prune'
                Prune = varargin{j+1};   
            case 'synchtype'
                SynchType = varargin{j+1};   
            case 'ssd'
                out_ssd = varargin{j+1};
                SSD_out = 1;
                Ncomp_ICA = size(out_ssd.X_Aug,1);
            case 'ngmd_method'
                if strcmp(varargin{j+1},'E5')
                    E5_flag = 1;
                   JADE_flag = 0;
                elseif strcmp(varargin{j+1},'JADE')
                    E5_flag = 0;
                   JADE_flag = 1;
                end
            case 'filter_b'
                butter_b = varargin{j+1};
            case 'filter_a'
                butter_a = varargin{j+1};
            
        end
    end
end
if ~exist('Fs')
    error('No input for sampling frequency.')   
end
if ~exist('SourceSetNum')
    error('No input for number of source sets.')   
end


%% SSD
if ~SSD_out
    opt_SSD = struct;
    opt_SSD.DoSSD = DoSSD;
    opt_SSD.Ncomp = SourceSetNum;
    fp = FrP*F_base;
    fq = FrQ*F_base;
    opt_SSD.fc = [fp fq];
    [out_ssd] = SSD_BP(X,Fs,opt_SSD,1,'prune',Prune,'filter_b',butter_b,'filter_a',butter_a);
end
%%  Non-gaussianity maximization
[out_E5,out_jade] =  DoICA_2(out_ssd.X_Aug,Ncomp_ICA,'flags',[E5_flag,JADE_flag]);
%% back-projecting to sensor-space
if ~isempty(out_E5)
    [A_final_E5] = Compute_Afinal(out_ssd,out_E5);
end
if ~isempty(out_jade)
    [A_final_jade] = Compute_Afinal(out_ssd,out_jade);
end
%~~~~~~
X_p_ssd = out_ssd.X_NB{1};
X_q_ssd = out_ssd.X_NB{2};
if ~isempty(out_E5)
    W_ICA_E5 = out_E5.W_ICA;
    X_Aug_ICA_E5 = out_E5.X_Aug_ICA;
    A_final_p_E5 = A_final_E5{1};
    A_final_q_E5 = A_final_E5{2};
end
if ~isempty(out_jade)
    W_ICA_jade = out_jade.W_ICA;
    X_Aug_ICA_jade = out_jade.X_Aug_ICA;
    A_final_p_jade = A_final_jade{1};
    A_final_q_jade = A_final_jade{2};
end
%% cmputing the synchronization index
if ~isempty(out_E5)
    Xp_E5 = X_p_ssd*W_ICA_E5(:,1:SourceSetNum)'; % extracted sources of slow oscillation
    Xq_E5 = X_q_ssd*W_ICA_E5(:,SourceSetNum+1:end)'; % extracted sources of fast oscillation
    
   
    if ismember(SynchType,[0,1,2]) % if phase coupling
        SynchFac_pq_E5 = Phase_Locking(Xp_E5,Xq_E5,FrP,FrQ);
    elseif ismember(SynchType,3) % if power coupling
        SynchFac_pq_E5 =  corr(abs(hilbert(Xp_E5)),abs(hilbert(Xq_E5)));
    elseif ismember(SynchType,4) % if power coupling
        SynchFac_pq_E5 =  corr(abs(angle(hilbert(Xp_E5))),abs(hilbert(Xq_E5)));
    end
end
% ~~~~~~~~~~~~~~+
if ~isempty(out_jade)
    Xp_jade = X_p_ssd*W_ICA_jade(:,1:SourceSetNum)'; % extracted sources of slow oscillation
    Xq_jade = X_q_ssd*W_ICA_jade(:,SourceSetNum+1:end)'; % extracted sources of fast oscillation
    

    if ismember(SynchType,[0,1,2]) % if phase coupling
        SynchFac_pq_jade = Phase_Locking(Xp_jade,Xq_jade,FrP,FrQ);
    elseif ismember(SynchType,3) % if power coupling
        SynchFac_pq_jade =  corr(abs(hilbert(Xp_jade)),abs(hilbert(Xq_jade)));
    elseif ismember(SynchType,4) % if power coupling
        SynchFac_pq_jade =  corr(abs(angle(hilbert(Xp_jade))),abs(hilbert(Xq_jade)));
    end
end

%% selecting PairsNum components which have the maximum non-gaussianity
f = @(x) (kurtosis(x) - 3).^2/48 + mean(x.^3).^2/12; %negentropy
if ~isempty(out_E5)
    K1 = [];
    K2 = [];
    K3 = 1:Ncomp_ICA;
    for k = 1:Ncomp_ICA-1
        for j = k+1:Ncomp_ICA
            err1 = Pattern_Angle( A_final_p_E5(:,k), A_final_p_E5(:,j));
            err2 = Pattern_Angle( A_final_q_E5(:,k), A_final_q_E5(:,j));
            if err1<0.1 && err2<0.1
                K1 = [K1,k];
                K2 = [K2,j];
            end
        end
    end
    %K = [K1;K2];
    f1 = f(X_Aug_ICA_E5(K1,:)'); f2 = f(X_Aug_ICA_E5(K2,:)');
    [~,ii]=min([f1;f2]);
    ii = ii-1;
    Ks1 = unique( K1.*(~ii) + K2.*ii);
    K3(ismember(K3,Ks1)) = []; %delete the min
    CC = diag(SynchFac_pq_E5);
    
    [~,ii] = sort(f(X_Aug_ICA_E5(K3,:)')'.*CC(K3),'descend');
    IDX_src_E5 = K3(ii(1:SourceSetNum)); % index of selected components
end
% ~~~~~~~~~~~~~~+
if ~isempty(out_jade)
    K1 = [];
    K2 = [];
    K3 = 1:Ncomp_ICA;
    for k = 1:Ncomp_ICA-1
        for j = k+1:Ncomp_ICA
            err1 = Pattern_Angle( A_final_p_jade(:,k), A_final_p_jade(:,j));
            err2 = Pattern_Angle( A_final_q_jade(:,k), A_final_q_jade(:,j));
            if err1<0.1 && err2<0.1
                K1 = [K1,k];
                K2 = [K2,j];
            end
        end
    end
    f1 = f(X_Aug_ICA_jade(K1,:)'); f2 = f(X_Aug_ICA_jade(K2,:)');
    [~,ii]=min([f1;f2]);
    ii = ii-1;
    Ks1 = unique( K1.*(~ii) + K2.*ii);
    K3(ismember(K3,Ks1)) = []; %delete the min
    CC = diag(SynchFac_pq_jade);
    
    [~,ii] = sort(f(X_Aug_ICA_jade(K3,:)')'.*CC(K3),'descend');
    IDX_src_jade = K3(ii(1:SourceSetNum)); % index of selected components
end
%% Jade or E5?
St0 = 2;
Success = 1;
if ~isempty(out_jade) && ~isempty(out_E5)
    negent_E5 = f(X_Aug_ICA_E5(IDX_src_E5,:)');
    negent_jade = f(X_Aug_ICA_jade(IDX_src_jade,:)');
    St0 = mean(negent_jade)>=mean(negent_E5);
elseif isempty(out_jade) && isempty(out_E5)
    Success = 0;
end

if Success
    St1 = ~isempty(out_jade) && isempty(out_E5);
    St2 = isempty(out_jade) && ~isempty(out_E5);
    
    
    if St0==1 || St1
        ICA_Method = 'JADE';
        out_NGMD = out_jade;
        IDX_src = IDX_src_jade;
        SynchFac_pq = SynchFac_pq_jade;
        A_final_p = A_final_p_jade(:,IDX_src);
        A_final_q = A_final_q_jade(:,IDX_src);
        Xp = Xp_jade;
        Xq = Xq_jade;
    elseif St0==0 || St2
        ICA_Method = 'fastICA_E5';
        out_NGMD = out_E5;
        IDX_src = IDX_src_E5;
        SynchFac_pq = SynchFac_pq_E5;
        A_final_p = A_final_p_E5(:,IDX_src);
        A_final_q = A_final_q_E5(:,IDX_src);
        Xp = Xp_E5;
        Xq = Xq_E5;
    end
    
    D_synch = diag(SynchFac_pq);
    D_synch = D_synch(IDX_src)';
    SynchFac_mean = mean(D_synch);
else
    A_final_p = [];
    A_final_q = [];
    out_NGMD = [];
    ICA_Method = [];
    D_synch = [];
    SynchFac_mean = [];
    out_ssd = [];
    IDX_src = [];
    Xp = [];
    Xq = [];
    
end
    

end %end function
%% Sub-functions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function [out_E5,out_jade] =  DoICA_2(X,Ncomp_ICA,varargin)
% this function does the non-gaussianity maximization decomposition
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

E5_flag = varargin{2}(1);
JADE_flag = varargin{2}(2);
disp('Do ICA~~~~~~~~~~~~~~~~~~~~~~');
disp('~~~~~~~~~~~~~~~~~~~~~~');
notConverge = zeros(2,1);
%% E5
if E5_flag
    try
        [out_E5.X_Aug_ICA, out_E5.A_ICA, out_E5.W_ICA,out_E5.W_2,out_E5.deMat,out_E5.WMat,out_E5.WSig] = ...
            fastica_NID(X,'g','E5','LastEig',Ncomp_ICA);
        if size(out_E5.X_Aug_ICA,1) ~= Ncomp_ICA
            warning('fastICA_E5 could not find enough components!');
            notConverge(1) = 1;
            out_E5 = [];
        end
    catch
        warning('fastICA_E5 did not converge!');
        notConverge(1) = 1;
        out_E5 = [];
    end
else
    notConverge(1) = 1;
    out_E5 = [];
end
%% jade
if JADE_flag
    try
        out_jade.W_ICA = jadeR(X,Ncomp_ICA);
        out_jade.X_Aug_ICA = out_jade.W_ICA*X;
        out_jade.A_ICA = cov(X')*out_jade.W_ICA'*inv(cov(out_jade.X_Aug_ICA'));%pinv(out_jade.W_ICA);
        if size(out_jade.X_Aug_ICA,1) ~= Ncomp_ICA
            warning('JADE could not find enough components!');
            notConverge(2) = 1;
            out_jade = [];
        end
    catch
        warning('JADE did not converge!');
        notConverge(2) = 1;
        out_jade = [];
    end
else
    notConverge(2) = 1;
    out_jade = [];
end
%%
if sum(notConverge)==2
    error('Non of the ICA algorithms converged!')
end

disp('ICA Finished~~~~~~~~~~~~~~~~~~~~~~');
disp('~~~~~~~~~~~~~~~~~~~~~~');
end
%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function [A_final] = Compute_Afinal(out_SSD,out_ICA)
%this function back-projects the ICA weights to SSD space
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
Nf = length(out_SSD.A_ssd);
A_ICA = out_ICA.A_ICA;
Ncomp = out_SSD.Ncomp;
A_final = cell(Nf,1);
for n = 1:Nf
    if out_SSD.DoSSD(n)
        A_ssd_n = out_SSD.A_ssd{n};
        A_final{n} = A_ssd_n*A_ICA((n-1)*Ncomp+1:n*Ncomp,:);
    else
        A_final{n}= A_ICA((n-1)*Ncomp+1:n*Ncomp,:);
    end
end
end
%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
