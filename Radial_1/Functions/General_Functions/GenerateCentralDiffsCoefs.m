function coefs = GenerateCentralDiffsCoefs(m,n)
% Build the coefficients for a central finite-difference stencil that
% corresponds to the m-th derivative with n-th order of accuracy. 
%
% Written by Manuel A. Diaz, ENSMA 2020.
% 
% There are 2p+1 or 2*floor((m+1)/2)-1+n central coefficients 
%                  a_{-p}, ..., 0, ..., a_{+p}.
% These are given by the solution of the linear equation system A * w = b :
%
%  [   1        1      .. 1 ..   1      1    ][a_-p]   [ 0 ]
%  [  -p      -p+1     .. 0 ..  p-1     p    ][  : ]   [ : ]
%  [ (-p)^2  (-p+1)^2  .. 0 .. (p-1)^2  p^2  ][  : ] = [ m!]
%  [   :        :         :      :       :   ][  : ]   [ : ]
%  [ (-p)^2p (-p+1)^2p .. 0 .. (p-1)^2p p^2p ][a_+p]   [ 0 ]
%
% where the only non-zero value on the R.H.S is in the (m+1)-th row.
% 
%   Inputs:  
%       m : m-th derivative
%       n : order of accuracy (must be an even value!)
%
%   Outputs:
%       coefs: column vector with central finite difference coefficients.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Notes: This formulation is stable (at most) for m<10 and n<20 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n_coefs=2*floor((m+1)/2)-1+n; p=(n_coefs-1)/2;

% Solve system A*w = b
A=power(-p:p,(0:2*p)'); b=zeros(2*p+1,1); b(m+1)=factorial(m); coefs=A\b;

% Round elements near values close to machine-epsilon to zero
coefs = coefs.*not(abs(coefs)<2000*eps);



