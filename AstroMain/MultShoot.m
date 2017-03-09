function [Z,DF_fin,mult] = MultShoot(Z,mult)
% function [Z,DF_fin,mult] = MultShoot(Z,mult)
% 
% Given an initial design variable vector and a structure of multiple 
% shooting parameters this function solves a multiple shooting problem.
%
% INPUTS:
%    Z          design variable vector (n_state*n+1 x 1)
%    mult       structure containing multiple shooting parameters
%
% OUTPUTS:
%    Z          final design variable vector (n_state*n+1 x 1)
%    DF_fin     DF matrix used for final update needed for use in Pseudo-Arclength Continuation method
%    mult       structure containing multiple shooting parameters
%
% Written by R. Pritchett, 02/09/17
% Last Update: R. Pritchett, 02/09/17
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Setup %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Extract necessary parameters from mult stucture
BndCase = mult.BndCase;
newt_tol = mult.newt_tol;
atten_tol = mult.atten_tol;
atten = mult.atten;
n = mult.n;

% Initialize Newton's method iteration counter
newti = 0;

% Print header for output table
fprintf('\n                         Sparsity                                                          \n')
fprintf(' Segs  len(X)  len(F)  Non-0s    (%%%%)        max(|dX|)         max(|F|)       Iter  \n')
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Sparse Matrix Setup %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate indices associated with nonzeros of the continuity constraints in the jacobian
[iRowC,jColC] = JacIndC(mult);

% Calculate indices associated with boundary constraints in the jacobian
switch BndCase
    case 'FixEndPtOnly'
        [iRowB,jColB] = JacIndB_FixEnd(mult);
    case 'Periodicity'
        [iRowB,jColB] = JacIndB_Period(mult);
end

% Assemble row and column indices into individual column vectors
iRow = [iRowC; iRowB]; 
jCol = [jColC; jColB];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Form Initial Constraint Vector %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate initial constraint vector
[F,STMi] = MakeF(Z,mult);

% Calculate error
error = max(abs(F));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate Partial Derivatives of the Jacobian %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate lengths of design and constraint vectors
Zl = length(Z);
Fl = length(F);

% Compute full Jacobian with either Forward or Central Difference Methods
% [DF_fwrd] = MakeDF_Num(Z,Zl,Fl,mult);
DF_fwrd = zeros(Fl,Zl);

% Compute only nonzeros of the Jacobian
switch BndCase
    case 'FixEndPtOnly'
        [DF] = MakeDF_FixEnd(Z,STMi,mult);
    case 'Periodicity'
        [DF] = MakeDF_Period(Z,STMi,mult);
end

% Convert DF into a sparse matrix
DF_sparse_fwrd = sparse(DF_fwrd); % for forward step case
DF_sparse_new = sparse(iRow,jCol,DF,Fl,Zl); % for numerically integrated case

% Calculate length of DF_sparse and DF sparsity percentage
DFl_fwrd = nnz(DF_sparse_fwrd); % for forward step case
DFl = length(DF_sparse_new); % for numerically integrated case
sparsity = 100*(1-DFl/(Fl*Zl));

%-------------------------------------------------------------------------%
% Code for Checking the Jacobian %
%-------------------------------------------------------------------------%

% % Compare forward step and sparse numericially integrated results
% DF_comp = full(DF_sparse_new) - DF_fwrd; % difference
% DF_comp_pcnt = (DF_comp./DF_fwrd)*100; % percent difference
% max_diff = max(DF_comp(:));
% max_pcnt_diff = max(DF_comp_pcnt(:));
% 
% % Plot finite difference and numerically integrated Jacobians (for debugging only)
% figure (1)
% hold on
% grid on
% spy(DF_sparse_new,'b')
% spy(DF_fwrd,'or')
% hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Update Design Variables %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Apply update equation to calculate a new design variable vector
[dZ] = Update(F,DF_sparse_new,Zl,Fl);
maxdZ = max(abs(dZ)); % calculate max change in X for print output

newti = newti + 1; % iterate counter by one

%Print output for one iteration of Newton's method
fprintf(' %4.0f  %6.0f  %6.0f  %6.0f  %8.2f  %16.8e  %16.8e  %5.0f \n',n-1,Zl,Fl,DFl,sparsity,maxdZ,error,newti)    

% Reset DF_sparse_new to DF_sparse (placed here in case Newton's Method while loop is not entered)
DF_sparse = DF_sparse_new;
      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Begin Newton's Method %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

while error > newt_tol
    
    % Reset DF_sparse_new to DF_sparse
    DF_sparse = DF_sparse_new;
    
    if max(dZ) > atten_tol
        Z = Z + (1/atten).*dZ; % take quarter step
    else
        Z = Z + dZ; % take a full step
    end

    % Calculate initial constraint vector
    [F,STMi] = MakeF(Z,mult);

    % Calculate error
    error = max(abs(F));

    % Calculate lengths of design and constraint vectors
    Zl = length(Z);
    Fl = length(F);

    % Compute full Jacobian with either Forward or Central Difference Methods
%     [DF_fwrd] = MakeDF_Num(Z,Zl,Fl,mult);
    DF_fwrd = zeros(Fl,Zl);

    % Compute only nonzeros of the Jacobian
    switch BndCase
        case 'FixEndPtOnly'
            [DF] = MakeDF_FixEnd(Z,STMi,mult);
        case 'Periodicity'
            [DF] = MakeDF_Period(Z,STMi,mult);
    end
    
    % Convert DF into a sparse matrix
    DF_sparse_fwrd = sparse(DF_fwrd); % for forward step case
    DF_sparse_new = sparse(iRow,jCol,DF,Fl,Zl); % for numerically integrated case

    % Calculate length of DF_sparse and DF sparsity percentage
    DFl_fwrd = nnz(DF_sparse_new); % for forward step case
    DFl = length(DF); % for numerically integrated case
    sparsity = 100*(1-DFl/(Fl*Zl));

    % Apply update equation to calculate a new design variable vector
    [dZ] = Update(F,DF_sparse_new,Zl,Fl);
    maxdZ = max(abs(dZ)); % calculate max change in X for print output

    newti = newti + 1; % iterate counter by one

    %Print output for one iteration of Newton's method
    fprintf(' %4.0f  %6.0f  %6.0f  %6.0f  %8.2f  %16.8e  %16.8e  %5.0f \n',n-1,Zl,Fl,DFl,sparsity,maxdZ,error,newti)

end

%% Output %%

% Save DF matrix used for final update for use in Pseudo-Arclength Continuation method
DF_fin = DF_sparse;

% Save number of Newton's method iterations
mult.newti = newti;

