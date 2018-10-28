function [X,I] = asgn(alpha,R,n,varargin)

% Creates random variates for aSGN(m) from the acceptance-rejection method.
% For the latter, a Student's t pdf is employed as a majorizing function
% and the associated alpha-sub-Gaussian (aSG) pdfs are evaluated from
% pre-computed tabulated .mat files of v(r;alpha,d). The algorithm is based
% on the outline provided in [1].
%
% The program computes the required aSG pdfs from the tabulated
% v(r;alpha,d) .mat files located at "/tab_files/vr_repo". The
% expression for v(r;alpha,d) is given by eq. (7) of [1]. The tabulated
% files allow generating distributions for 'alpha' in [1.1:0.01:1.98] and
% 'd' in [1:1:10], where d=length(R). The program employs persistent
% variables and linear interpolation between sampled points of
% v(r;alpha,d).
%
% ******** Inputs *********
% alpha:        The characteristic exponent associated with the aSGN(m)
%               process. This is a scalar input and should lie within
%               (1.1:0.01:1.98)
%
% R:            The covariance matrix of any adjacent 'm+1' samples in an
%               aSGN(m) process. The dimension of 'R' is equal to
%               'm+1'.It should be a symmetric toeplitz matrix. The maximum
%               acceptable size of 'R' is 10x10.
%
% n:            The number of samples required.
%
% init_samples: An optional argument, which forces the first 'm' samples of
%               the noise process to 'init_samples'. The length
%               of the vector should be equal to 'm'.
%
% ******** Output *********
% X:    The 1 x n vector of aSGN(m) samples
%
% I:    The innovation (driving) samples of the corresponding white noise
%       process. For the alpha=1, these are IID samples stemming from a
%       t-distribution with degrees-of-freedom m+1 and scale=1.
%
%---------------------------------------
%
% There are a few things to be kept in mind. 'R' is a symmetric toeplitz
% matrix and must be of full-rank with (m+1) rows/columns. If the diagonal
% elements are delta^2, then each sample is S(alpha,delta) distributed. See
% [1] for more information. Moreover, based on the tabulation of
% v(r;alpha,d), 'R' can be of a maximum possible dimension 10 x 10.
%
% 'init_samples' is useful when a continuous realization of aSGN(m) is
% required, but the realization is too large to be evaluated in one go and
% needs to be constructed over several iterations. In such a case the
% following minimalistic example highlights its application
%
% for j=1:N
%   if j==1
%       noise =asgn(alpha,R,n);
%   else
%       noise =asgn_read(alpha,R,noise(end-length(R)+2:end));
%   end
% end
%
%---------------------------------------
%
% References:
%
% [1] A. Mahmood and M. Chitre, "Generating random variates for stable
%     sub-Gaussian processes with memory", Signal Processing, Volume 131,
%     Pages 271-279, 2017. (https://doi.org/10.1016/j.sigpro.2016.08.016.)
%
%------------------------
% Author: Ahmed Mahmood
% Year: 2016


persistent c vJoint InvR SigRootX1 InvRX1 kappa k1 k2 ...
    kmarg step Nmin Nmax rind ModeFactor alphaPrev RPrev

if nargin>4
    error('Exceeding the number of input paramters')
end

m=size(R,1)-1;
    
if m~=0
    funk1=@(d) (2^alpha)*sin(pi*alpha/2)*gamma((alpha+2)/2)*gamma((alpha+d)/2)./(gamma(d/2)*(pi*alpha)/2); 
    funk2=@(d) 4*gamma(d./alpha)./(alpha*(2.^d).*((gamma(d/2)).^2));
    funkmarg=@(d) gamma(d/2)/(gamma((d-1)/2)*sqrt(pi));
    c=1.2;  % is the acceptance-rejection factor. See Table-I in [1] for
                % appropriate values c=1.2 is sufficient for 'alpha' in [1.1:0.01:1.98]
                % and 'd=m+1' in [2:1:10]. One can reduce 'c' for faster computation
                % This depends on 'alpha' and 'd' and can be ascertained from the
                % aforementioned table
                
    %******* 'R' Operations ********
    
    if ~isequal(RPrev,R)
        if ~isequal(R,R.') || sum(diag(R)~=R(1,1))>0
            error('R should be a symmteric Toeplitz matrix of size m+1')
        end
        k1=funk1(m:m+1);
        k2=funk2(m:m+1);
        kmarg=funkmarg(m+1);
        kappa=det(R)/det(R(1:end-1,1:end-1));
        InvR=inv(R);
        InvRX1=inv(R(1:end-1,1:end-1));
        SigRootX1=chol(R(1:end-1,1:end-1),'lower');
        RPrev=R;
        ModeFactor=(R(end,1:end-1)/(R(1:end-1,1:end-1)));
        disp('*Change in R -- Evaluating Cholesky decompositions, inverses and the conditional location parameter*')
    end
    
    %******* Generating Constants  & Loading Tabulated Data ********
    
    if ~isequal(alphaPrev,alpha)
        k1=funk1(m:m+1);
        k2=funk2(m:m+1);
        %fpath=mfilename('fullpath');
        %[fpath,~,~] = fileparts(fpath); 
        %load([fpath,'\tab_files\vr_repo\vr_alpha=',num2str(alpha),'.mat']) % Gets the variables: 'vJoint','Nmax','Nmin','res'
        load(['vr_alpha=',num2str(alpha),'.mat'])
        alphaPrev=alpha;
        step=(log10(Nmax)-log10(Nmin))/res;
        disp('*Loading v(r;alpha,d)*')
    end
    
    if m>size(vJoint,1)-1
        error(char("The dimensions of 'R' exceed the maximum possible 10x10"))
    end 
    
    %****************************************************************
    
    if n<=m
        if length(varargin)==1
            if length(varargin{1})==m
                X=varargin{1};
            else
                error(char("'init_samples' should be of length 'm'"))
            end
        else
            A=stabrnd(alpha/2,1,2*(cos(pi*alpha/4))^(2/alpha),0,1,1);
            T= chi2rnd(m);
            S=randn(m,1);
            S=S/sqrt(sum(S.^2));
            X=(SigRootX1*sqrt(A*T)*S).'; % from pg. 4, of "Multivariate elliptically contoured stable
                                         % distributions: theory and estimation" by J. P. Nolan
        end
        X=X(1:n);
    else
        X = zeros(1,n); % Preallocate m
        I=zeros(1,n);
        if length(varargin)==1
            if length(varargin{1})==m
                X(1:m)=varargin{1};
            else
                error(char("'varargin{1}' should be of length 'm'"))
            end
        else
            A=stabrnd(alpha/2,1,2*(cos(pi*alpha/4))^(2/alpha),0,1,1);
            T= chi2rnd(m);
            S=randn(m,1);
            S=S/sqrt(sum(S.^2));
            X(1:m)=(SigRootX1*sqrt(A*T)*S).';
        end
        vstud=alpha+m;
        norms=tpdf(0,vstud);
        %k=m;
        for i = m+1:n
            
            X1=X(1,i-m:i-1).';
            mode=ModeFactor*X1;
            
            norm1=SubGaussCondProbTabulate(alpha,X1,mode,InvRX1,InvR,vJoint,Nmin,Nmax,step,rind,kappa,k1,k2,kmarg);
            
            accept = false;
            while accept == false
                
                u = rand();
                t_innov=trnd(vstud,1,1);
                v=(norms/norm1)*t_innov+mode;
                g_v=(norm1/norms)*tpdf((v-mode)*(norm1/norms),vstud);
                f_v=SubGaussCondProbTabulate(alpha,X1,v,InvRX1,InvR,vJoint,Nmin,Nmax,step,rind,kappa,k1,k2,kmarg);
                
                %f_v=SubGaussCondProb(alpha,R,X1,v);
                if c*u <= f_v/g_v
                    X(i) = v;
                    I(i) = t_innov;
                    accept = true;
                end
                %      k=k+1;
            end
            %disp(['i=',num2str(i),' k=',num2str(k)])
        end
    end
else
    if isscalar(R)
        X=stabrnd(alpha,0,sqrt(R),0,1,n);
    else
        error(char("'R' should be a symmteric Toeplitz matrix of size 'm+1'"))
    end
end
end

function f_X2X1=SubGaussCondProbTabulate(alpha,X1,X2_ind,InvRX1,InvR,vJoint,Nmin,Nmax,step,rind,kappa,k1,k2,kmarg)

% Generates the conditional probability f(X2|X1) if [X1, X2] is a
% sub-Gaussian stable random vector such that X1(i)~X2~S(alpha,delta) and
% rho is the correlation coefficient of the underlying Gaussian vector. We
% assume the joint-probabiluty is given by f(X1,X2).

m=length(X1);
% if size(InvSigRoot,1)*size(InvSigRoot,2)~=dim^2;
%     error('R should be a length(X1)+1 square matrix')
% end

%r1=sqrt((InvSigRootX1*X1).'*(InvSigRootX1*X1));
r1=sqrt(X1.'*InvRX1*X1);
X=[X1; X2_ind];
%r=sqrt((InvSigRoot*X).'*(InvSigRoot*X));
r=sqrt(X.'*InvR*X);

if r1<Nmin
    %vJointR1=k2(1);
    grad=(vJoint(m,1)-k2(1))/Nmin;
    cons=k2(1);
    vJointR1=grad*r1+cons;
elseif r1>Nmax
    vJointR1=alpha*k1(1)*(r1^(-alpha-m));
else
    tempInd=(log10(r1)-log10(Nmin))/step+1;
    tempInd=[floor(tempInd), ceil(tempInd)];
    grad=(vJoint(m,tempInd(1))-vJoint(m,tempInd(2)))/(rind(tempInd(1))-rind(tempInd(2)));
    cons=vJoint(m,tempInd(1))-grad*rind(tempInd(1));
    vJointR1=grad*r1+cons;
    %vJointR1=interp1([rind(tempInd(1)),rind(tempInd(2))],[vJoint(m,tempInd(1)),vJoint(m,tempInd(2))],r1);
end

if r<Nmin
    %vJointR=k2(2);
    grad=(vJoint(m+1,1)-k2(2))/Nmin;
    cons=k2(2);
    vJointR=grad*r+cons;
elseif r>Nmax
    vJointR=alpha*k1(2)*(r^(-alpha-m-1));
else
    tempInd=(log10(r)-log10(Nmin))/step+1;
    tempInd=[floor(tempInd), ceil(tempInd)];
    grad=(vJoint(m+1,tempInd(1))-vJoint(m+1,tempInd(2)))/(rind(tempInd(1))-rind(tempInd(2)));
    cons=vJoint(m+1,tempInd(1))-grad*rind(tempInd(1));
    vJointR=grad*r+cons;
    %vJointR=interp1([rind(tempInd(1)),rind(tempInd(2))],[vJoint(m+1,tempInd(1)),vJoint(m+1,tempInd(2))],r);
end

f_X2X1=(1/sqrt(kappa))*kmarg*vJointR/vJointR1;

end