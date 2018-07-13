clc; clear; close all

% Checks the validity of all tabulated v(r;alpha,d) .mat-files, by plotting
% the v(r;alpha,d) curve against 'r' and comparing it to its closed-form
% convergent expressions for 'r=0' and 'r->Inf'. The files are located in
% '.\File Repo'.
%
% The expression for v(r;alpha,d) is given by eq. (7) of [1]. One can
% derive any d-dimensional elliptic symmetric alpha-stable distribution
% from v(r;alpha,d). The tabulated files allow generating distributions for
% 'alpha' in [1.1:0.01:1.98] and 'd' in [1:1:10].
%
%---------------------------------------
%
% For each tabulated .mat file, the corresponding 'alpha' is explicitly
% highlighted in its title. The following variables are stored in each
% file:
%
% Nmax = 10^3 
% Nmin = 10^-3
% res  = 2000
% rind = 10.^(log10(Nmin):(log10(Nmax)-log10(Nmin))/res:log10(Nmax));
% vJoint = v(rind;alpha,d)
%
% The numerically computed values of v(r;alpha,d) are stored as the matrix
% vJoint, which has 10 rows and length(rind) columns. The rows of vJoint
% corresponds to 'd' in(1:1:10). The range of 'rind' is meant to
% encapsulate the structure of v(rind;alpha,d) within its convergent
% regimes between r=0 and r->Inf. The sample points 'rind' that run between
% 'Nmin' and 'Nmax' is adequate in this regard
%
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

alphaInd=1.1:0.01:1.98;

fpath=mfilename('fullpath');
[fpath,~,~] = fileparts(fpath);

for j=1:length(alphaInd)
    load([fpath,'\vr_repo\vr_alpha=',num2str(alphaInd(j)),'.mat'])
    dimVec=1:size(vJoint,1);
    grid on
    rindLarge=rind(rind>10);
    
    K0=(4*gamma(dimVec/alphaInd(j)))./(alphaInd(j).*(2.^dimVec).*...
                 (gamma(dimVec/2)).^2);                 % Theoretical value of v(r) as r->0 
    K1=(2^alphaInd(j))*sin(pi*alphaInd(j)/2)*gamma((alphaInd(j)+2)/2)*...
        (gamma((alphaInd(j)+dimVec)/2)./...
        ((gamma(dimVec/2))*pi/2));                      % Theoretical value of v(r) as r->Inf
    
    vJointNearZero=K0;
    vJointLargeCorrected=zeros(length(dimVec),length(rindLarge));
    
    for i=1:length(dimVec)
        vJointLargeCorrected(i,:)=K1(i)*rindLarge.^(-alphaInd(j)-dimVec(i));
    end
    
    % overall comparison
    sfigure(1)
    subplot(2,1,1)
    loglog(rind,vJoint.','-k')
    hold on
    loglog(rindLarge,vJointLargeCorrected.','--r','LineWidth',2)
    loglog(ones(1,length(dimVec))*Nmin,vJointNearZero,'ob');
    title(['{alpha}= ',num2str(alphaInd(j)),', tails'])
    xlabel('$r$','Interpreter','Latex')
    ylabel('$v(r,\alpha,d)$','Interpreter','Latex')
    ylim([10^-35 10^0])
    hold off
    grid on
    
    % zoomed-in visual
    subplot(2,1,2)
    loglog(rind,vJoint.','-k')
    hold on
    loglog(ones(1,length(dimVec))*Nmin,vJointNearZero,'ob');
    title(['{alpha}= ',num2str(alphaInd(j)),', near r=0'])
    xlabel('$r$','Interpreter','Latex')
    ylabel('$v(r,\alpha,d)$','Interpreter','Latex')
    ylim([10^-5 10^0])
    hold off
    grid on
    drawnow  
end