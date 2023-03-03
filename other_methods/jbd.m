function [B, RR, RC, phi, sweep]=jbd(RR, W, threshold, max_sweep, c_dim)
% function [B, RR, RC, phi, sweep]=jbd(RR, W, threshold, max_sweep, c_dim)
%
% Approximate BLOCK joint diagonalization of positive symmetric matrices
% by a quasi-Newton technique
%
%    Copyright (C) 2011 Dana Lahat and Jean-François Cardoso
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
% This function calculates and returns the matrix B which best jointly
% block-diagonalizes the matrices R1...RQ.
% The function also returns the block-diagonalized matrices in RR.
% The stopping criterion is based on the Kullback-Leibler divergence.
% The criterion is weighted with the positive weights W with size 1xQ
% The transformation B is computed with an accuracy equal to threshold.
%
%
% Input:
% RR = [ R1 R2 ... RQ] has size [m,m*Q] and is the concatenation
%           of Q symmetric real positive-definite m*m matrices
% W         is a length-Q vector of weights for the Q matrices in RR
% threshold is the stopping criterion. Typically 1e-6 -- 1e-9
% max_sweep is the maximal number of sweeps, that is, iterations over all blocks.
% c_dim     is the vector holding the block-pattern: sum(c_dim)=m
%
% Output:
% B     is the estimated demixing matrix, invertible. Size [m,m]
% RR    is the set of block-diagonalized matrices. Size [m,m*Q]
% RC =  return code. 1 if success, 0 if did not converge (that is, max_sweep
%       reached and threshold not achieved).
% phi   is a vector containing all scalar values of the Kullback-Leibler
%       based criterion which we try to minimize.
% sweep is the number of sweeps used.
%
%
% This algorithm is an extension to Pham's work, "Joint Approximate
% Diagonalization of Positive Definite Hermitian Matrices", 2001.
% This is an exact coding of the algorithm which was published in ICA/LVA
% 2012, see citation below.
%
% Authors: Dana Lahat, Jean-François Cardoso, 
%          Dana@Lahat.org.il, cardoso@tsi.enst.fr
% First version: February 2009
%
% Publication to cite:
% @INPROCEEDINGS{Lahat12_jbd,
%   author = {D. Lahat and J.-F. Cardoso and H. Messer},
%   title = {Joint Block Diagonalization Algorithms for Optimal Separation of Multidimensional Components},
%   booktitle = {Latent Variable Analysis and Signal Separation},
%   year = {2012},
%   editor = {Fabian Theis and Andrzej Cichocki and Arie Yeredor and Michael Zibulevsky},
%   volume = {7191},
%   series = LNCS,
%   pages = {155--162},
%   address = {Heidelberg},
%   publisher = {Springer}
% }

% Parameter Initializations
m = sum(c_dim); % Total dimension of all components
n = length(c_dim); % Number of components
Q = length(W); % Number of matrices 
W = W / sum(W); % normalize. W is used only as a relative weight, not as a number of samples.
RR = reshape(RR, m, m, Q) ;
%------------------------------------------------
% A binary matrix of indices on the main block-diagonal
BlockPattern = blockpattern(c_dim) ;
%------------------------------------------------
% Calculate indices for block start
BlockStart=zeros(n,1); % Block Start index
for i=1:n % loop over component index
    BlockStart(i)=sum(c_dim(1:i-1))+1; % Block "i" starts at...
end
%--------------------------------------------------
% Calculate block indices
BlockIndex=blockindex(c_dim, BlockStart);
%------------------------------------------------
% Calculate commutation matrix (cell array, for all required dimensions)
T_commutation=cell(n,n);
for ii=1:n
    for jj=1:ii-1
        T_commutation{ii,jj}=T_operator(c_dim(ii),c_dim(jj));
    end
end
%------------------------------------------------
sweep = 1; % init number of sweeps
inv_B = eye(m); % init for A
RC=1; % if failed, RC=0;
lambda=1; % Normally, this is one

E_ =cell(n,n); % init temporary cell-array E_
%------------------------------------------------


% Begin iterations

while (sweep <max_sweep)
    
    phi(sweep)=calc_mismatch(RR,W,BlockPattern); % mismatch for previous sweep
    
    nabla= calc_g(RR, W, BlockPattern);
    ng=norm(nabla,'fro') ;
    %fprintf('jbd: sweep= %i   phi= %e norm(g)=%e\n', sweep, phi(sweep), ng);
    if (ng < threshold), break; end;
    ng_vec(sweep)=ng;
    %------------------------------------------------
    
    H=calc_Hess(RR, W, c_dim, T_commutation, BlockIndex);
    for ii=1:n
        for jj=1:ii-1
            gij = nabla(BlockIndex{ii},BlockIndex{jj});
            gji = nabla(BlockIndex{jj},BlockIndex{ii});
            gij = gij(:);
            gji = gji(:);
            Eij_ji = H{ii,jj}\[gij;gji];% inversion of small matrices instead of the entire Hessian
            E_{ii,jj} = reshape(Eij_ji(1:c_dim(ii)*c_dim(jj)) , c_dim(ii), c_dim(jj));
            E_{jj,ii} = reshape(Eij_ji(c_dim(ii)*c_dim(jj)+1:end), c_dim(jj), c_dim(ii) );
        end
        E_{ii,ii}= zeros(c_dim(ii), c_dim(ii));
    end
    E = cell2mat(E_);
    T = eye(m)-lambda*E ;
    inv_B = inv_B*T ;
    RR=update_P(RR,inv(T),Q);
    
    sweep=sweep +1;
end

B=inv(inv_B);
if 0;%verbose
    fprintf('jbd.m: Total number of sweeps =%3i. Criterion =%10.5e\n', sweep, phi(end));
end %if

end % function jbd

%----------------------------------

function g=calc_g(P, W, BlockPattern)
% function g=calc_g(P, W, BlockPattern)
%
% Calculate relative gradient

g=0;
[m,~,Q]=size(P);
% -----------------------------
% Relative Gradient g
% -----------------------------
for q=1:Q
    g=g + W(q)*((P(:,:,q).*BlockPattern)\P(:,:,q));
end
g=eye(m)-g;% zero on main diagonal

end % function calc_g

%---------------------------------------------

function H=calc_Hess(P, W, c_dim, T_commutation, BlockIndex)
% function H=calc_Hess(P, W, c_dim, T_commutation, BlockIndex)
%
% Calculate Hessian

n=length(c_dim);
H=cell(n,n);
for ii=1:n
    for jj=1:ii-1
        T=T_commutation{ii,jj};
        Riiq = extract_comp(P,c_dim,ii,ii, BlockIndex);
        Rjjq = extract_comp(P,c_dim,jj,jj, BlockIndex);
        iRiiq = inv_q(Riiq);
        iRjjq = inv_q(Rjjq);
        Hij = kron_nq(Rjjq, iRiiq,W);
        Hji = kron_nq(Riiq, iRjjq,W);
        Hessij = [Hij, T'; T, Hji];
        H{ii,jj} = Hessij;
        H{jj,ii} = Hessij';
    end
end

end % function calc_Hess

%-----------------------------------------------

function kld=calc_mismatch(P,W,BlockPattern)
% function kld=calc_mismatch(P,W,BlockPattern)
%
% kld is the weighted Kullback-Leibler induced divergence (that is, the mismatch) between the set
% of matrices P and their block-diagonal counterparts, with block-pattern
% BlockPattern

Q=length(W);
temp=0;
for q=1:Q
    temp = temp + W(q)*log(det(P(:,:,q).*BlockPattern) / det(P(:,:,q)));
end
kld=temp/2;

end % function calc_mismatch

%-----------------------------------------------

function Pout=update_P(Pin, T, Q)
% function Pout=update_P(Pin, T, Q)
%
% Apply the linear transformation T to P:
% P(:,:,q)=T*Pin(:,:,q)*T', q=1,...,Q
%
% Input:    Pin [n,n,Q], T [n,n]
% Output:   Pout [m,m,Q]
[nrows, ~]=size(T);
Pout=zeros(nrows, nrows, Q);
for q=1:Q
    Pout(:,:,q)=T*Pin(:,:,q)*T';
end

end % function update_P

%-----------------------------------------------

function iPq=inv_q(Pq)
% function iPq=inv_q(Pq)
%
%   Purpose: Pq is [ni,nj,Q]
%       iPq(:,:,q)= inv(Pq(:,:,q))

[ni,~,Q]=size(Pq);
iPq=zeros(ni,ni,Q);% init
for q=1:Q
    iPq(:,:,q)= inv(Pq(:,:,q));
end

end % function inv_q
%-----------------------------------------------
function g=kron_nq(Pi,Pj,W)
%
% g= <kron(Pi, Pj)>
% Weighted average of kron(Piq,Pjq)

Q=length(W);
g=0;
for q=1:Q
      g=g+W(q)*kron(Pi(:,:,q),Pj(:,:,q));
end

end % function kron_nq
%-----------------------------------------------

function Pq_c=extract_comp(Pq, c_dim, ii, jj, BlockIndex)
% function Pq_c=extract_comp(Pq, c_dim, ii, jj, BlockIndex)
%
% Purpose:
%       Extract the (ii,jj) block- component from a set of Q matrices
% Output:
%   Pq_c   [c_dim(ii),c_dim(jj),Q]

[~,~,Q]=size(Pq);
Pq_c=zeros(c_dim(ii), c_dim(jj), Q);
for q=1:Q
    Pq_c(:,:,q)=Pq(BlockIndex{ii},BlockIndex{jj},q);
end

end % function extract_comp


%----------------------------------------------------
function bi=blockindex(c_dim, BlockStart)
% function bi=blockindex(c_dim, BlockStart)
%
% bi is an [n,1] cell-array. The i'th cell is a [c_dim(i),1] vector which holds the
% indices of the i'th block in the block-pattern c_dim

n=length(c_dim);
bi=cell(n,1);
bs=BlockStart;
for ni=1:n
    bi{ni}=bs(ni):bs(ni)+c_dim(ni)-1;
end;

end % function blockindex

%----------------------------------------------------

function B = blockpattern(c_dim)
% function B = blockpattern(c_dim)
%
% Build a 0/1 matrix B of booleans for components
% of dimension c_dim(1), ... c_dim(end)
% B(i,j) is true if indices i and j belong to the same block


ncomp = length(c_dim) ;
N = sum(c_dim) ;

if (ncomp == 1) %% a single block !
    B = true(N) ;
    return
end%if

B = false(N) ;
istart = 1 ;
for icomp = 1:ncomp
    istop  = istart+ c_dim(icomp)-1 ;
    B(istart:istop, istart:istop) = true ;
    istart = istop+1 ;
end%for

end% function blockpattern

%-----------------------------------------------

function T=T_operator(ii,jj)
% function T=T_operator(ii,jj)
%
% T is the commutation matrix such that vec(X')=T*vec(X), where X is an
% [ii,jj] matrix

Z=zeros(ii,jj);
T=zeros(ii*jj);
for r=1:ii
    for s=1:jj
        Ers=Z;
        Ers(r,s)=1;
        Ers_Ers=kron(Ers,Ers');
        T=T+Ers_Ers;
    end
end

end % function T_operator
%-----------------------------------------------
