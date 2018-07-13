# The $\alpha$SGN($m$) MATLAB toolbox

**Ahmed Mahmood  (July 2018)**



## Brief Intro.

The aim of this toolbox is to incorporate emerging techniques related to $\alpha$-stable processes within Matlab's functionality. Our topic of interest is the $\alpha$SGN($m$) model. As of now, literature on the topic can be found primarily at <https://arl.nus.edu.sg/twiki6/bin/view/ARL/Publications > and consists of (1)--(6) amongst others.

This document is offers a quick overview of the toolbox. For more details, revert to the actual scripts. Before listing down the functions within this toolbox, we introduce a few concepts to aid the discussion: 

The heavy-tailed $\alpha$SGN($m$) process is derived from an $(m+1)$-dimensional $\alpha$-sub-Gaussian ($\alpha$SG) distribution, which in turn is parameterized by the characteristic exponent $\alpha\in(0,2)$ and the covariance matrix $\mathbf{R}$. The $\alpha$SG distribution is elliptic and can be denoted by $\alpha$SG($\alpha,\mathbf{R}$) or equivalently by $\alpha$SG($\alpha,\delta,\mathbf{Cov}$), where $\mathbf{Cov}$ is the normalized covariance matrix and $\delta\in(0,\infty)$ is the scale, i.e., $\mathbf{R}=\delta^2\mathbf{Cov}$. The $\alpha$SGN($m$) process essentially constrains any of its $(m+1)$ adjacent samples to follow the aforementioned $\alpha$SG distribution. This ensures $\mathbf{R}$ (and thus $\mathbf{Cov}$) is a symmetric Toeplitz matrix. We note that each sample of an $\alpha$SGN($m$) process is a symmetric $\alpha$-stable (S$\alpha$S) random variable. The latter's distribution is completely characterized by the tuple $(\alpha,\delta)$ and can be expressed as $\mathcal{S}(\alpha,\delta)$.

For more information about the employed parameterization, revert to (3). For a great introduction to $\alpha$SG distributions, do go through Nolan's discussion on the topic at <http://fs2.american.edu/jpnolan/www/stable/EllipticalStable.pdf>.

##Functions

Add all directories to Matlab's path before use.

###Random variates
- ${\tt X=asgn(alpha,R,N,\_')}$: Generates ${\tt N}$ samples of $\alpha$SGN($m$) with underlying distribution $\alpha$SG(${\tt alpha,R}$). The function accepts optional inputs as well (highlighted by ${\tt\_}$).

  

- ${\tt x=stabrnd(alpha,beta,delta,mu,N,M)}$: McCulloch's original script that returns independent outcomes of the stable distribution $\mathcal{S}({\tt alpha,beta,delta,mu})$, where ${\tt beta}$ is the skew parameter and ${\tt mu}$ is the location. Note that $\mathcal{S}({\tt alpha,0,delta,0})$ is statistically equivalent to $\mathcal{S}({\tt alpha,delta})$. The outcomes are returned as the ${\tt N\times M}$ matrix ${\tt x}$.

	​	

###Fitting

- ${\tt [alpha,mu,delta]=sstabfit(x)}$:  Mandar Chitre's script that fits an S$\alpha$S distribution to the data vector ${\tt x}$.

  

- ${\tt [alpha,delta,Cov]=asgnfit(X,m,\_)}$: Estimates $\alpha$SGN($m$) parameters from the data vector ${\tt X}$ for a given ${\tt m}$. Also, ${\tt \_}$ signifies optional inputs.

  

- ${\tt [pxx,pasg,f]=asgnsd(X,fs)}$: Determines ${\tt pxx}$, the normalized non-parameterized spectral density, for the data vector ${\tt X}$ sampled at ${\tt fs}$ kHz. The associated frequency vector ${\tt f}$ (in kHz) is also returned. Moreover, ${\tt pasg}$ is the closed-form *parametric* spectral density of ${\tt X}$ evaluated under the assumption that ${\tt X}$ holds consecutive samples of $\alpha$SGN($m$). Comparing ${\tt pxx}$ to ${\tt pasg}$ allows discerning a suitable $m$ to characterize ${\tt X}$ within the $\alpha$SGN($m$) framework.

  

### Computing PDFs

- ${\tt logf\_XN=asgnpdf(X,alpha,R)}$: Samples the joint-pdf and returns its logarithm for each column (outcome) of ${\tt X}$, under the assumption that the columns hold consecutive samples of $\alpha$SGN($m$) with underlying distribution $\alpha$SG(${\tt alpha, R}$).

- ${\tt f\_X= saspdf(X,alpha,delta)}$: Samples and returns the pdf associated with $\mathcal{S}({\tt alpha},{\tt delta})$ at each element (outcome) in the matrix ${\tt X}$. Note that ${\tt f\_X}$ is the same size of ${\tt X}$.


### Storing and retrieving $\alpha$SGN($m$) data

As $\alpha$SGN($m$) is a GARCH process, samples have to be sequentially computed. Consequently, ${\tt asgn()}$ can become a potential bottleneck for computation time when running intensive performance tests/simulations. To circumvent this, realizations may be pre-computed and stored offline, only to be retrieved for later use. The following functions are helpful in this regard: 

- 
  ${\tt asgn\_write(alpha,Cov,FileSize,NumFiles,fs)}$: Computes and stores a realization of $\alpha$SGN($m$) with underlying distribution $\alpha$SG(${\tt alpha},{1},{\tt Cov}$)= $\alpha$SG(${\tt alpha},{\tt Cov}$). The samples are stored as ${\tt NumFiles}$ files, each of which is sized ${\tt FileSize}$ MBs. Also, ${\tt fs}$ is the associated sampling frequency. The files are stored in a sub-directory of the function's root folder. As an example, if ${\tt alpha}=1.57$ and $m=$4 (implying ${\tt Cov}$ is a $5\times5$ matrix), the files are saved in "\a1_57__m_4\\" and are labeled as "asgn\_1.bin", "asgn\_2.bin" and so on. 
  	

- $ {\tt [x,\_] = asgn\_read(fpath,samps,\_)}$: Reads ${\tt samps}$ samples from the files written by ${\tt asgn\_write()}$. The absolute path directory of these files is passed as the input string ${\tt fpath}$.


### Miscellaneous

- ${\tt vr\_validate}$: refer to the script for details.




## References

1. A. Mahmood and M. Chitre,  *Temporal Analysis of Stationary Markov $\alpha$-sub-Gaussian Noise,* in OCEANS 2016 MTS/IEEE, (Monterey, CA, USA), September 2016. 
2. A. Mahmood, M. Chitre, and V. Hari, *Locally optimal inspired detection in snapping shrimp noise,* IEEE Journal of Oceanic Engineering, vol. 42, pp. 1049--1062, October 2017.
3. A. Mahmood and M. Chitre, *Generating Random Variates for Stable Sub-Gaussian Processes with Memory,* Signal Processing, vol. 131, pp. 271--279, February 2017.
4. A. Mahmood and M. Chitre, *Optimal and Near-Optimal Detection in Bursty Impulsive Noise,* IEEE Journal of Oceanic Engineering, vol. 42, pp. 639--653, October 2016.
5. A. Mahmood and M. Chitre, *Uncoded Acoustic Communication in Shallow Waters with Bursty Impulsive Noise,* in Underwater Communications Networking (Ucomms 2016), (Lerici, Italy), September 2016. (Invited).
6. A. Mahmood, V. Hari, and M. Chitre, *Model-Based Signal Detection in Snapping Shrimp Noise,* in Underwater Communications Networking (Ucomms 2016), (Lerici, Italy), September 2016. (Invited).