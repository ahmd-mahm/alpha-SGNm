function logf_XN = asgnpdf(X,alpha,R)

% Accepts the N x L real matrix 'X'. Each column of 'X' represents an
% N-dimensional outcome of a random vector consitituting N immediately
% adjacent samples of an aSGN(m) process with characteristic exponent
% 'alpha' and covariance matrix 'R'. The program computes the joint-PDF of
% X(:,i) and outputs its natural log as 'f_XN(1,i)'. Note that 'f_XN' is a
% L-dimensional row vector. Moreover if 'N<=length(R)', then 'f_XN(1,i)' is
% the log of the pdf of an N-dimensional alpha-sub-Gaussian random vector
% at outcome 'X(:,i). We employ persistent variables and spline
% interpolation between sampled points in vJoint(r;alpha,d).
%
% The program computes the required pdf from the tabulated v(r;alpha,d)
% .mat files located at "..\Tabulated Files\File Repo". The expression for
% v(r;alpha,d) is given by eq. (7) of [1]. The tabulated files allow
% generating distributions for 'alpha' in [1.1:0.01:1.98] and 'd' in
% [1:1:10], where d=length(R). The program employs persistent variables and
% spline interpolation between sampled points of v(r;alpha,d).
%
% ******** Inputs *********
% X     : a matrix of real numbers
%
% alpha :The characteristic exponent of the required pdf, should lie within
%        the range (1.1:0.01:1.98)
%
% R     :The covariance matrix between any consecutive samples of a aSGN(m)
%        process. The dimension of 'R' is equal to 'm+1', where m is the
%        order of the aSGN(m) process. It should be a symmetric toeplitz
%        matrix. The maximum acceptable size of 'R' is 10x10.
%
% ******** Output *********

% logf_XN  :The 1 x L vector, such that 'f_XN(1,i)' is the (natural log) of
%           the N-dimensional pdf at 'X(:,i)', where the latter is assumed
%           to be an outcome of 'N' immediately adjacent samples of an
%           aSGN(m) process parameterized by 'alpha' and 'R'
%
%---------------------------------------
%
% The actual pdf may be acquired by typing 'exp(log_f_XN)' in the main
% program. However, this may be a really small value for large 'N' and thus
% prone to error.
%
% References:
%
% [1] A. Mahmood and M. Chitre, "Generating random variates for stable
%     sub-Gaussian processes with memory", Signal Processing, Volume 131,
%     Pages 271-279, 2017. (https://doi.org/10.1016/j.sigpro.2016.08.016.)
%
%------------------------
% Author: Ahmed Mahmood
% Year: 2015

persistent vJoint k1 k2 Nmin Nmax alphaPrev pp pp1 RPrev kdet...
    InvR InvRX1 NPrev rind


m=size(R,1)-1;
N= size(X,1);

funk1=@(d) (2^alpha)*sin(pi*alpha/2)*gamma((alpha+2)/2)*gamma((alpha+d)/2)./(gamma(d/2)*(pi*alpha)/2);
funk2=@(d) 4*gamma(d./alpha)./(alpha*(2.^d).*((gamma(d/2)).^2));

if ~isequal(RPrev,R) || ~isequal(alphaPrev,alpha) || ~isequal(NPrev,N)
    if ~isequal(RPrev,R)
        if ~isequal(R,R.') || sum(diag(R)~=R(1,1))>0
            error('R should be a symmteric Toeplitz matrix of size m+1')
        end
        kdet=[det(R(1:end-1,1:end-1)), det(R)];
        InvR=inv(R);
        InvRX1=inv(R(1:end-1,1:end-1));
        disp('*Change in R -- Evaluating inverses and determinants*')
    end
    if ~isequal(alphaPrev,alpha)
        fpath=mfilename('fullpath');
        [fpath,~,~] = fileparts(fpath);
        load([fpath,'\..\tab_files\vr_repo\vr_alpha=',num2str(alpha),'.mat']);
        disp('*Loading v(r;alpha,d)*')
    end
    if m>size(vJoint,1)-1
        error(char("The dimensions of 'R' exceed the maximum possible 10x10"))
    end
    
    log_rind=log10(rind);
    if N > m
        k1=funk1(m:m+1);
        k2=funk2(m:m+1);
        log_vJoint1=log10(vJoint(m,:));
        log_vJoint=log10(vJoint(m+1,:));
        pp1=spline(log_rind,log_vJoint1);
        pp=spline(log_rind,log_vJoint);
    else
        k1=funk1(N);
        k2=funk2(N);
        log_vJoint=log10(vJoint(N,:));
        pp=spline(log_rind,log_vJoint);
    end
    NPrev=N;
    alphaPrev=alpha;
    RPrev=R;
end


%****************************************************************

if N<=m
    logf_XN=log(SubGaussJointPDFTabulate(alpha,X,inv(R(1:N,1:N)),vJoint,Nmin,Nmax,pp,k1,k2,det(R(1:N,1:N))));
else
    logf_XN=log(SubGaussJointPDFTabulate(alpha,X(1:m,:),InvRX1,vJoint,Nmin,Nmax,pp1,k1(1),k2(1),kdet(1)));
    for i=1:N-m
        logf_XN=logf_XN+log(SubGaussJointPDFTabulate(alpha,X(i:m+i,:),InvR,vJoint,Nmin,Nmax,pp,k1(2),k2(2),kdet(2)))...
            -log(SubGaussJointPDFTabulate(alpha,X(i:m+i-1,:),InvRX1,vJoint,Nmin,Nmax,pp1,k1(1),k2(1),kdet(1)));
    end
end
end

function f_X=SubGaussJointPDFTabulate(alpha,X,InvR,vJoint,Nmin,Nmax,pp,k1,k2,kdet)

dim=size(X,1);

r=sqrt(sum((InvR.'*X).*X,1));
vJointR=zeros(1,length(r));

temp=r(r<Nmin);
grad=(vJoint(dim,1)-k2)/Nmin;
cons=k2;
vJointR(r<Nmin)=grad*temp+cons;


temp=(r(and(r>=Nmin,r<Nmax)));
if ~isempty(temp)
    vJointR(and(r>=Nmin,r<Nmax))=10.^ppval(pp,log10(temp));
end

temp=r(r>=Nmax);
vJointR(r>=Nmax)=alpha*k1*(temp.^(-alpha-dim));

f_X=(1/sqrt(kdet))*gamma(dim/2)/(2*(pi^(dim/2)))*vJointR;

end