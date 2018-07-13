function [alpha,c,delta] = sstabfit(x)

% SSTABFIT Fit a symmetric alpha stable distribution to data
%
% [alpha,c,delta] = sstabfit(x)
%
% alpha, c and delta are the characteristic exponent, scale parameter
% (dispersion^1/alpha) and location parameter respectively.
%
% alpha is computed based on McCulloch (1986) fractile
% c is computed based on Fama & Roll (1971) fractile
% delta is the 50% trimmed mean of the sample
%
% Author: Mandar Chitre
% Last Modified: July 30, 2004
% $Revision: 1.1 $

%% lookup table from McCulloch (1986)
ena = [    2.4388     2.4388     2.4388     2.4388     2.4388;
           2.5120     2.5117     2.5125     2.5129     2.5148;
           2.6080     2.6093     2.6101     2.6131     2.6174;
           2.7369     2.7376     2.7387     2.7420     2.7464;
           2.9115     2.9090     2.9037     2.8998     2.9016;
           3.1480     3.1363     3.1119     3.0919     3.0888;
           3.4635     3.4361     3.3778     3.3306     3.3161;
           3.8824     3.8337     3.7199     3.6257     3.5997;
           4.4468     4.3651     4.1713     4.0052     3.9635;
           5.2172     5.0840     4.7778     4.5122     4.4506;
           6.3140     6.0978     5.6241     5.2195     5.1256;
           7.9098     7.5900     6.8606     6.2598     6.1239;
          10.4480     9.9336     8.7790     7.9005     7.6874;
          14.8378    13.9540    12.0419    10.7219    10.3704;
          23.4831    21.7682    18.3320    16.2163    15.5841;
          44.2813    40.1367    33.0018    29.1399    27.7822 ];

%% compute
delta = trimmean(x,50);
p = prctile(x,[5 25 28 72 75 95]);
c = (p(4)-p(3))/1.654;
an = (p(6)-p(1))/(p(5)-p(2));
if an < 2.4388
  alpha = 2;
else
  j = min(find(an <= ena(:,1)));
  if isempty(j)
    error('Unable to determine alpha');
  end
  t = (an-ena(j-1,1))/(ena(j,1)-ena(j-1,1));
  alpha = (22-j-t)/10;
  if alpha < 0.5
    alpha = 0.5
  end
end
