function asgn_write(alpha,Cov,FileSize,NumFiles,fs)

% Generating aSGN(m) realizations is time taking. It is perhaps pertinent
% to generate realizations offline and save them as data files. These can
% be read later to significantly reduce simulation time for Monte Carlo
% trials. This program generates an aSGN(m) realization and saves it as
% 'NumFiles' data files, each of which is sized 'FileSize' MB. The program
% also displays periodic status updates about its completion.
%
% The data files are saved in a folder within the 'asgn_write()' root
% directory, and are titled on the basis of 'alpha' and 'm=length(Cov)-1'.
% As an example, if alpha=1.57 and m=4, the files are saved in
% '/a1_57__m_4/' and are labelled as 'asgn_1.bin', 'asgn_2.bin' and so on.
% 'asgn_write()' will produce a warning if the directory already exists
%
% 'asgn_write()' is intended to be used once for a set of parameters, after
% which 'asgn_read()' can be used to retrieve the data.
%
%
% ******** Inputs *********
% alpha   : the characteristic exponent associated with the aSGN process
%
% Cov     : the (m+1)-dimensional covariance matrix. This needs to be a
%           symmetric Toeplitz matrix wth unit diagonal entries
%
% FileSize: the size of each data file in megabytes (MBs). We recommend a
%           'FileSize' up until 2048 MB (2GB)
%
% NumFiles: the number of files required
%
% fs      : the sampling frequency in (kHz)
%
%
% example inputs:
%
% alpha=1.5;
% Cov= toeplitz([1.0000, 0.6873, 0.3420, 0.2021, 0.0663, 0.0451, 0.0240, -0.0338, -0.0960]);
% FileSize=2048;
% NumFiles= 10;
% fs=180;
% asgn_write(alpha,Cov,FileSize,NumFiles,fs)
%
%---------------------------------------
%
% An important note: fs is the supposed sampling frequency. Though fs has
% no impact on generating the aSGN(m) samples, it is necessary for
% associating the aSGN(m) noise samples to practical scenarios. For
% example, let's say we were to estimate the parameters of aSGN(m) from a
% practical dataset sampled at 'fs' and we generate an aSGN(m) realization
% (based on these estimates) for later use. If an algorithm's performance
% is to be tested on these samples, the results would be meaningless from a
% practical perspective if 'fs' is not mentioned. This stems from the fact
% that impulsive noise attributes depend heavily on the frequency regimes
% one plans working in.
%
% To generate aSGN(m) with arbitrary scale 'delta', retrieve the samples
% via 'asgn_read()' and multiply the samples by 'delta'.
%
% The 'FileSize' and 'NumFiles' values depend on the simulations you
% eventually intend to accomplish. A back-of-the-hand calculation can be
% performed as follows: Say you want to compute the bit error rate
% performance of a communication system up to 10^-5, with atleast 100 bits
% in error. Let the number of passband samples per baseband symbol be k.
% Then the number of noise samples required for such a simulation would be
% approximately k*100/10^-5 = k*10^7. As each sample is stored as an 8-byte
% float, the total space required would be (k*8*10^7)/(1024^2) MBs.

%------------------------
% Author: Ahmed Mahmood
% Year: 2015


if ~isequal(Cov,Cov.') || sum(diag(Cov)~=1)>0
    error("Ensure that 'Cov' is a symmetric Teoplitz matrix with unit diagonal elements")
end

d=size(Cov,1);
samps=(FileSize*(1024^2))/8; % Dividing by 8 as each sample is stored as 8-byte floating point.
samps_blk=10000;             % The code generates 'samps_blk' samples of aSGN(m) at a time.

Iter= ceil(samps/samps_blk); % Number of calls to the 'aSGN' function. The number of samples is rounded of to the nearest multiple of 'samps_blk' 
%samps=Iter*samps_blk;

a_str=num2str(alpha);
m_str=num2str(d-1);
fpath=mfilename('fullpath');
[fpath,~,~] = fileparts(fpath); 



disp(['Final file size: ',num2str(Iter*samps_blk*8/(1024^2)),' MB'])    % Displays the (corrected) size of a file.
disp(['Samples per file: ',num2str(Iter*samps_blk)])                    % Displays the number of samples in a file.
disp(['fs = ',num2str(fs),' kHz'])                                      % Displays 'fs'
disp(['File duration: ',num2str(Iter*samps_blk/(fs*1000)),'secs ('...
    ,num2str(Iter*samps_blk/(fs*1000*3600)),'hrs)'])                    % The time (w.r.t 'fs') that each data file spans.
disp(['Number of files: ',num2str(NumFiles)])                           % Displays the total number of files to be saved
disp(['Total duration: ',num2str(NumFiles*Iter*samps_blk/(fs*1000)),'secs or ('...
    ,num2str(NumFiles*Iter*samps_blk/(fs*1000*3600)),'hrs)'])           % The total duration (w.r.t 'fs') ofthe data files span.

disp('===========================')
disp(['alpha = ',a_str])    % Display 'alpha'.
disp(['m = ',m_str])        % Displays the order of the aSGN process.
disp('===========================')
fpath=[fpath,'/a',a_str(1),'_',a_str(3:end),'__m_',m_str];

[~,ex,~]=mkdir(fpath);
if ~isempty(ex)
    warning('on')
    warning(['The directory already exists: ', fpath])
    disp('---------------------------')
end

for i=1:NumFiles
    
    disp(['File Number: ',num2str(i)])
    for j=1:Iter
        
        if i==1 && j==1
            noise=asgn(alpha,Cov,samps_blk);
        else
            noise=asgn(alpha,Cov,samps_blk+d-1,noise(end-d+2:end));
            noise=noise(d:end);
        end
        
        if j==1
            fid = fopen([fpath,'/asgn_',num2str(i),'.bin'],'w');
        else
            fid = fopen([fpath,'/asgn_',num2str(i),'.bin'],'a');
        end
        fwrite(fid, noise,'double');
        fclose(fid);
        
        if mod(j,round(Iter/5))==0
            disp(['Percentage Completed: ',num2str(j/Iter*100)])
        end
    end
    disp('---------------------------')
end
