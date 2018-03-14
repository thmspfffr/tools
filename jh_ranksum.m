%  * jh_ranksum.c
%  * July 2012, (c) Joerg Hipp
%  *
%  * This function computes statistical z-values for the 
%  * Mann-Whitney U test (also called the Mann-Whitney-Wilcoxon (MWW) 
%  * or Wilcoxon rank-sum test)
%  *
%  * http://en.wikipedia.org/wiki/Mann-Whitney_U
%  *
%  * To rebuild type "mex jh_ranksum.c"
%  * 
%  * Usage:
%  * outMatrix = jh_ranksum(inMatrix);
%  * 
%  * inMatrix .. matrix with single precision input data [n_vox,n_vox,n_subjects]
%  *             equal amount of subjects are assumed, while first all data from the
%  *             one group and then all data from the other group.
%  * outMatrix .. single precision output (statistical z values [n_vox,n_vox]
%  *
%  * To derive p-values in matlab: p = 2*normcdf(-abs(outMatrix));
