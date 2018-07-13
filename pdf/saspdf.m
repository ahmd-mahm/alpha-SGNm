function f_X=saspdf(X,alpha,delta)

% Numerically computes the SaS univariate pdf corresponding to
% S(alpha,delta)=S(alpha,0,delta,0) for each element (representing random
% outcomes) in the matrix 'X' and returns 'f_X', where 'f_X(i,j)' is the
% probability density of the S(alpha,delta) distributed random variable at
% outcome 'X(i,j)'
%
% The program computes the required SaS pdf from the tabulated v(r;alpha,d)
% .mat files located at "../tab_files/vr_repo". The expression for
% v(r;alpha,d) is given by eq. (7) of [1]. The tabulated files allow
% generating distributions for 'alpha' in [1.1:0.01:1.98]. The program
% employs persistent variables and spline interpolation between sampled
% points of v(r;alpha,d).
%
% ******** Inputs *********
% X    : a matrix of real numbers
%
% alpha:The characteristic exponent of the required pdf, should lie within
%       the range (1.1:0.01:1.98)
%
% delta:the scale of the required pdf, should be a positive real number
%
% ******** Output *********
% f_X: The S(alpha,delta) pdf evaluate at 'X'
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
% Year: 2015

persistent vJoint k1 k2 Nmin Nmax rind alphaPrev pp

if  ~isequal(alphaPrev,alpha)
    %fpath=mfilename('fullpath');
    %[fpath,~,~] = fileparts(fpath);
    %load([fpath,'/../tab_files/vr_repo/vr_alpha=',num2str(alpha),'.mat']);
    load(['vr_alpha=',num2str(alpha),'.mat']);
    alphaPrev=alpha;
    k1=(2^alpha)*sin(pi*alpha/2)*gamma((alpha+2)/2)*gamma((alpha+1)/2)./(gamma(1/2)*(pi*alpha)/2);
    k2=4*gamma(1/alpha)/(alpha*2*((gamma(1/2)).^2));
    log_vJoint=log10(vJoint(1,:));
    log_rind=log10(rind);
    pp=spline(log_rind,log_vJoint);
    disp('*Loading v(r;alpha,d)*')
end

if delta<=0
    error("'delta' must be a positive real number")
end

X=abs(X)/delta;

vJointR=zeros(size(X,1),size(X,2));

temp=X(X<Nmin);
grad=(vJoint(1,1)-k2)/Nmin;
cons=k2;
vJointR(X<Nmin)=grad*temp+cons;

temp=X(and(X>=Nmin,X<Nmax));
vJointR(and(X>=Nmin,X<Nmax))=10.^ppval(pp,log10(temp));

temp=X(X>=Nmax);
vJointR(X>=Nmax)=alpha*k1*(temp.^(-alpha-1));

f_X=gamma(1/2)/(2*delta*sqrt(pi))*vJointR;

end