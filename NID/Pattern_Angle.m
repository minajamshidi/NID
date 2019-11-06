function [mismatch_Mat] = Pattern_Angle( p1, p2)
% computes the dissimilarity of two spatial patterns
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
%%
vec_mismatch = @(V1,V2) 1-abs(dot(V1,V2,1))./(sqrt(sum(V1.^2,1)).*sqrt(sum(V2.^2,1)));
N1 = size(p1,2);
N2 = size(p2,2);
mismatch_Mat = NaN (N1, N2);
for i = 1:N1
    for j = 1:N2
        mismatch_Mat(i,j) = vec_mismatch (p1(:,i), p2(:,j));
    end
end

end