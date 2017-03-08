function [iRowC,jColC] = JacIndC_LT(colt)
% function [iRowC,jColC] = JacIndC(N,n_seg,n_state)
% 
% This function calculates the row and column indices of the nonzero 
% values of the collocation constraint Jacobian. These indices are later 
% used in the construction of the sparse DF matrix. Note that the output 
% column vectors do not include indices for the boundary constraints.
% Indices for the boundary constraints are calculated in JacIndB.
%
% INPUTS:
%    colt     structure containing collocation and optimization parameters
%
% OUTPUTS:
%    iRowC    row indices of the collocation constraint jacobian 
%    jColC    column indices of the collocation constraint jacobian 
%
% Originally Written by: R. Pritchett, 06/07/2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Extract necessary parameters from colt stucture
N = colt.N; % degree of interpolating polynomial
n_seg = colt.n_seg; % number of segments
n_state = colt.n_state; % number of states, same as number of equations of motion
n_cntrl = colt.n_cntrl; % number of control variables
n_slack = colt.n_slack; % number of slack variables
n_coast = colt.n_coast; % number of coasting parameters

%% Preallocate Matrices %%

% Total design variables per variable node
l = (N+1)*n_state/2 + n_cntrl + n_slack;

Dx0_iRow=zeros(n_state,l,n_seg);
Dx0_jCol=zeros(n_state,l,n_seg);

Dxf_iRow=zeros(n_state,l,n_seg);
Dxf_jCol=zeros(n_state,l,n_seg);

Ddelta_iRow=zeros(n_state*(N-1)/2,l,n_seg);
Ddelta_jCol=zeros(n_state*(N-1)/2,l,n_seg);

%% Fill Matrices %%

for ii=1:n_seg
    
    for jj=1:n_state
        Dx0_iRow(jj,:,ii)=jj+n_state*(ii-1);
    end
    
    for jj=1:l
        Dx0_jCol(:,jj,ii)=n_coast+l+jj+l*(ii-1);
    end
 
    for jj=1:n_state
        Dxf_iRow(jj,:,ii)=jj+n_state*(ii-1);
    end
    
    for jj=1:l
        Dxf_jCol(:,jj,ii)=n_coast+jj+l*(ii-1);
    end    
    
    for jj=1:n_state*(N-1)/2
        Ddelta_iRow(jj,:,ii)=n_state*(n_seg-1)+jj+n_state*(ii-1)*(N-1)/2;
    end    
    
    for jj=1:l
        Ddelta_jCol(:,jj,ii)=n_coast+jj+l*(ii-1);
    end
    
end

%% Assemble Column Vectors %%

iRowC=[ reshape(Dx0_iRow(:,:,1:end-1),[n_state*l*(n_seg-1) 1]);
        reshape(Dxf_iRow(:,:,1:end-1),[n_state*l*(n_seg-1) 1]);
        reshape(Ddelta_iRow,[n_state*l*(N-1)*n_seg/2 1])  ];
        
jColC=[ reshape(Dx0_jCol(:,:,1:end-1),[n_state*l*(n_seg-1) 1]);
        reshape(Dxf_jCol(:,:,1:end-1),[n_state*l*(n_seg-1) 1]);
        reshape(Ddelta_jCol,[n_state*l*(N-1)*n_seg/2 1])  ];