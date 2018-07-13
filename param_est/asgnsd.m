function [pxx,pasg,f]=asgnsd(X,fs)

% This function takes in a vector (of noise/data samples) 'X' and the
% sampling frequency 'fs' and computes its spectral density 'pxx' via the
% pwelch method. Moreover, the aSGN(m) model is tuned to 'X' by employing
% 'asgnfit()'. The resulting paramters are then used to evaluate the
% closed-form spectral density (see eq. (7) of [1]) within the aSGN(m)
% framework. This results in 'pasg', which is computed for 'm' in (2:2:8).
% Both 'pxx' and 'pasg' are plotted against 'f' and are such that the areas
% under each curve are normalized to 1.
%
% 'asgnsd(X,fs)' shows the impact of varying 'm' on the spectrum
% (statistics) tracking ability of the aSGN(m) model on 'X'.
%
%
% ******** Inputs *********
% X : a data vector
%
% fs: the sampling frequency in (kHz)
%
% ******** Outputs *********
% pxx : the normalized spectral density of 'X' computed by the Welch
%       method.
%
% pasg: the normalized spectral density of aSGN(m) tuned to 'X'. 'pasg' has
%       four coulmns corresponding to m in (2:2:8).
%
% f   : the frequency vector (in kHz) associated with the samples in 'pxx'
%       and 'pasg'
%
% ******** References *********
% [1] A. Mahmood, M. Chitre and H. Vishnu, "Locally Optimal Inspired
% Detection in Snapping Shrimp Noise," in IEEE Journal of Oceanic
% Engineering, vol. 42, no. 4, pp. 1049-1062, Oct. 2017.
%
%------------------------
% Author: Ahmed Mahmood
% Year: 2016

if ~isvector(X)
    error("'X' should be a vector")
end


N=512;
f=(0:1/N:1)/2;

m=2:2:8;
dim=m+1;
max_dim=length(dim);
rho=zeros(1,max_dim);
rho(1)=1;

pasg=zeros(max_dim,length(f));
[~, delta, R_max] = asgnfit(X,max(m));

for i=1:max_dim
    R=R_max(1:m(i)+1,1:m(i)+1);
    rho(i)=R(1,dim(i));
    a=(R(end,1:end-1)/(R(1:end-1,1:end-1))).';
    sigma2=2*(delta^2)*(R(end,end)-R(end,1:end-1)*a);
    temp=fliplr(a.')*(exp(-1i*2*pi*((1:dim(i)-1).'*f)));
    pasg(i,:)=sigma2./(abs(1-temp).^2);
end
pasg=pasg.';
pasg=pasg./(sum(pasg,1)*(f(2)-f(1))*fs);
plot(f*fs,10*log10(pasg),'LineWidth',1)
xlabel('frequency (kHz)')
ylabel('spectral density')
grid on
hold on

[pxx,~] = pwelch(X,N,[],N*2,1);
pxx=pxx./(sum(pxx)*fs/(2*N));
plot(f*fs,10*log10(pxx),'LineWidth',1)
hold off
legend(['m=',num2str(m(1))],['m=',num2str(m(2))],['m=',num2str(m(3))],['m=',num2str(m(4))],'pwelch(X)')