function C = ip(A,B)
% function C = ip(A,B)
% 
% Assuming A and B are 3D matrices this function multiplies the "slices" of
% A and B together. For example if A is a (4x5x10) matrix while B is a
% (5x7x10) matrix then multiplying the "slices" of A and B together
% involves multiplying 10 4x5 matrices by 10 5x7 matrices. The resulting 
% matrix C will then be a 3D matrix with the dimensions (4x7x10). Matlab 
% does not automatically handle this kind of operation with 3D matrices, it
% can only perform elementwise multiplication or division of 3D matrices. 
% Therefore, this function was written to perform the desired operation. In
% the context of the collocation problem this function enables the
% polynomials governing a trajectory of many segments to be rapidly
% computed, this is much faster than looping through segments and
% performing each subsequent matrix multiplication.
%
% INPUTS:
%    A      full time matrix for nodes that are variables
%           (n+1 x n+1 x m)
%    B      time matrix for nodes that are constraints (assume does not
%           contain endpoints of segment)  (n+1 x (n-1)/2 x m)
%
% OUTPUTS:
%    C      matrix of polynomial coefficients
%
% Originally Written by: D. Grebow, 07/15/2015
% Last Update: R. Pritchett, 02/17/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%-------------------------------------------------------------------------%
% Original code implemented by Daniel Grebow
%-------------------------------------------------------------------------%

% sa=size(A);sb=size(B);ii=ones(sa(1),1)*[1:sb(2)];
% C=reshape(sum(repmat(permute(A,[2 1 3:numel(sa)]),[1 sb(2)]).*B(:,ii(:)',:)),[sa(1) sb(2) sa(3:end)]);

%-------------------------------------------------------------------------%
% Original code broken into distinct steps for clarity by Robert Pritchett
%-------------------------------------------------------------------------%

% Calculate size of matrices A and B
sa = size(A);
sb = size(B);

% Create matrix of column indices for sub-matrices of inner product result
ii=ones(sa(1),1)*[1:sb(2)];

% Step 1: Use permute to transpose each sub-matrix of A, i.e. each 2D that exists along the third dimension of A
step1 = permute(A,[2 1 3:numel(sa)]);

% Step 2: Replicate the result of step 1 along the second dimension of the matrix,
% the number of replications should equal the value of the second dimension of B
step2 = repmat(step1,[1 sb(2)]);

% Step 3: This step essentially calculates the dot product of each row of A with
% each column of B, however it does so in a manner that makes really
% efficient use of Matlab's matrix math capabilities. The step's prior to
% this one were simply rearranging A and B in a manner that facilitated the
% quick computation of these dot products. For insight onto why this is
% done look up the explanation provided in the inner product section of the
% Wikipedia page on matrix multiplication.
step3 = sum(step2.*B(:,ii(:)',:));

% Step 4: Reshape the result of the dot product operations such that each 2D
% sub-matrix of step 3 is sa(1) x sb(2) matrix. In other words the
% sub-matrices will have the same number of rows as A and the same number
% of columns of B. Note: A and B must have the same value along the third
% dimension.
C = reshape(step3,[sa(1) sb(2) sa(3:end)]);

return
