function [alpha, delta, Cov] = asgnfit(X,m,varargin)

% Call as
% [alpha, delta, Cov] = asgnfit(X,m) or
% [alpha, delta, Cov] = asgnfit(X,m,p,'covariation')
%
% This function trains the aSGN(m) model to the data vector 'X'. It
% evaluates the characteristic exponent 'alpha' and the scale 'delta' via
% McCulloch's method using 'sstabfit()'. The underlying normalized
% covariance matrix 'Cov' can be evaluated via either least-squares (LS) or
% the covariation method. In both cases the dimension of 'Cov' is  m x m.
% Revert to [1] for details on these methods.
%
% The default method to evaluate 'Cov' is that of LS. This offers an
% effective estimate of 'Cov' which is guaranteed to be positive
% semi-definite. Call the function as 'asgnfit(X,m)' to tune the model via
% the LS method.
% 
% The covariation method requires an additional parameter 'p'. Ideally, 1 <
% p < alpha. In most practical impulsive scenarios p=1 is sufficient. Call
% the function as 'asgnfit(X,m,p,'covariation')' to tune the model via the
% covariation method. This does not always guarantee a positive definite
% 'Cov' but offers a robust estimate. For larger 'm', we have observed this
% method to increasingly produce non-positive definite estimates.
%
%
% If the actual covariance matrix of 'X' is required, it can be obtained as
% (delta^2)*Cov in the main program.
%
%---------------------------------------
%
% References:
%
% [1] C. L. Nikias and M. Shao, Signal processing with Alpha-Stable Distributions
%     and Applications. New York: Chapman-Hall, 1996.
%
%------------------------
% Author: Ahmed Mahmood
% Year: 2015

if ~isvector(X)
    error('X should be a vector')
end

[alpha,delta,~]=sstabfit(X); % McCulloch's estimator of 'alpha' and 'delta'.
Cov=zeros(m+1,m+1);
Xlen=length(X);

switch nargin
    case 2
        del=sqrt(sum(X.^2)/(2*Xlen));
        for i=1:m
            r=sum(X(1:end-i).*X(1+i:end))/(Xlen-i);
            Cov=Cov+diag(r*ones(m+1-i,1),i);
        end
        Cov=(Cov+Cov.')+2*(del^2)*eye(m+1);
        Cov=Cov/(2*del^2);
        
    case 4
        p=varargin{1};
        est=varargin{2};
        if ~isscalar(p) || ~strcmp(est,'covariation')
            error(char("Inputs are incorrect: Should be of the form asgnfit(X,m,p,'covariation')"))
        else
            C=((sum(abs(X).^p)/Xlen)^(1/p))/delta; % From...
            for i=1:m
                tempXlen=Xlen-mod(Xlen,i);
                Xtemp=reshape(X(1:end-mod(Xlen,i)),i,tempXlen/i);
                if mod(tempXlen/i,2)~=0
                    Xtemp=Xtemp(:,1:end-1);
                    tempXlen=size(Xtemp,1)*size(Xtemp,2);
                end
                Xtemp=reshape(Xtemp.',2,tempXlen/2);
                r= (2/(C^p))*(delta^(2-p))*(Xtemp(1,:)*((sign(Xtemp(2,:)).*(abs(Xtemp(2,:)).^(p-1))).'))/(tempXlen/2);
                Cov=Cov+diag(r*ones(m+1-i,1),i);
            end
            Cov=(Cov+Cov.')+2*(delta^2)*eye(m+1);
            Cov=Cov/(2*delta^2);
        end
    otherwise
        error(char("Either use asgnfit(X,m) for a least-squares computation of Cov or asgnfit(X,m,p,'covariation') for a covariation based computation of Cov"))
end

end