%We generated 3 aSGN(m) files in '\a1_5__m_8\' as an example. This was done
%by

alpha=1.57;
Cov= toeplitz([1.0000, 0.6873, 0.3420, 0.2021, 0.0663, 0.0451, 0.0240,...
    -0.0338, -0.0960]); % note that Cov is 9 x 9, thus m=8. 
FileSize=0.5; % 0.5 MB
NumFiles=3;
fs=180; % This has no impact on generating the realization
asgn_write(alpha,Cov,FileSize,NumFiles,fs);


%To read this data, equate 'fpath' to the absolute path string of the data
%files. If 10000 samples need to be retrieved, use the following snippet:

fpath=mfilename('fullpath');
[fpath,~,~] = fileparts(fpath);
fpath=[fpath,'\a1_57__m_8\'];
samps=10000;
delta=5;
x = delta*asgn_read(fpath,samps);
figure
plot(x)

% running x = asgn_read(fpath,samps) again would call in the next 10000
% samples and so on untill all the data in all the files is exhausted. Note
% that 'delta' is required scale, i.e., each sample in x is S(alpha,delta)
% distributed.