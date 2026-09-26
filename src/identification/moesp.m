function [sigma, ssfun] = moesp(y, u, ell)
% MOESP  Multivariable Output-Error State-Space Subspace Identification
% Author: Guanru Pan, TUHH ICS
%
%   [SIGMA, SSFUN] = MOESP(Y, U, ELL) identifies the observable subspace 
%   from measured data using the MOESP algorithm.
%
%   Inputs:
%     Y   - An N-by-ny output data matrix, where N is the number of 
%           sampling points and ny is the number of output variables.
%     U   - An N-by-nu input data matrix corresponding to Y.
%     ELL - The embedding dimension (block size), which should be greater 
%           than the expected system order.
%
%   Outputs:
%     SIGMA - A vector of singular values (subspace scores), which can be 
%             used to determine the system order. The system order n can be 
%             chosen such that sum(SIGMA(1:n)) ≈ sum(SIGMA).
%
%     SSFUN - A function handle for identifying the state-space matrices. 
%             Once the system order n is determined, call:
%
%               [A, B, C, D] = SSFUN(n)
%
%             to obtain the corresponding state-space realization of the 
%             identified system.

% Input and output check
error(nargchk(1,3,nargin));
error(nargoutchk(0,4,nargout));
[ndat,ny]=size(y);
[mdat,nu]=size(u);
if ndat~=mdat
    error('Y and U have different length.')
end

% block Hankel matrix
N = ndat-ell+1;
Hy = zeros(ell*ny,N);
Hu = zeros(ell*nu,N);
sN=sqrt(N);
sy=y'/sN;
su=u'/sN;
for s=1:ell
    Hy((s-1)*ny+1:s*ny,:)=sy(:,s:s+N-1);
    Hu((s-1)*nu+1:s*nu,:)=su(:,s:s+N-1);
end

% SVD
[U,S,V] = svd(Hy*(eye(size(Hu,2))-pinv(Hu)*Hu));

sigma = diag(S);
% n=find(cumsum(ss)>0.85*sum(ss),1);
ssfun = @ssmat;

    function [A,B,C,D]=ssmat(n)
        % C and A
        Ol = U(:,1:n)*diag(sqrt(sigma(1:n)));
        C=Ol(1:ny,:);
        A=pinv(Ol(1:ny*(ell-1),:))*Ol(ny+1:ell*ny,:);
        % B and D
        U2_t = U(:,n+1:end)';
        M1 = U2_t*Hy*pinv(Hu);
        m = ny*ell-n;
        M = zeros(m*ell,nu);
        L = zeros(m*ell,ny+n);
        for k=1:ell
            M((k-1)*m+1:k*m,:)=M1(:,(k-1)*nu+1:k*nu);
            L((k-1)*m+1:k*m,:)=[U2_t(:,(k-1)*ny+1:k*ny) U2_t(:,k*ny+1:end)*Ol(1:end-k*ny,:)];
        end
        DB=pinv(L)*M;
        D=DB(1:ny,:);
        B=DB(ny+1:end,:);
    end
end