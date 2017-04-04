function [Z,x_bnd,tau_nodes,C,colt] = CollSolve_LT(Z,t_bnd,header,colt)
% function [Z,x_bnd,tau_nodes,C,colt] = CollSolve_LT(Z,t_bnd,header,colt)
% 
% Given an initial design variable vector and a set of boundary node times 
% this function solves the collocation problem when low-thrust is included.
%
% INPUTS:
%    Z          design variable vector (n_coast+(n_state+n_cntrl+n_slack)*n_seg*(N+1)/2 x 1)
%    t_bnd      non-normalized times at  boundary nodes (n_seg x 1)
%    header     string that indicates whether to print output of each Newton's method iteration 
%    colt       structure containing collocation and optimization parameters
%
% OUTPUTS:
%    Z          final design variable vector (n_coast+(n_state+n_cntrl+n_slack)*n_seg*(N+1)/2 x 1)
%    slack      column vector containing final slack variable values (n_slack x 1)
%    x_bnd      matrix of boundary node states (l x 1 x n_seg+1)
%    tau_nodes  vector of normalized variable and defect node times (N x 1) (same for all segments)
%    C          matrix of polynomial coefficients (l x (N+1) x n_seg)
%    colt       structure containing collocation and optimization parameters
%
% Written by R. Pritchett, 6/07/16
% Last Update: R. Pritchett, 01/31/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Setup %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Extract necessary parameters from colt stucture
NodeSpace = colt.NodeSpace;
Mesh = colt.Mesh;
N = colt.N;
n_seg = colt.n_seg;
n_state = colt.n_state;
n_coast = colt.n_coast;
newt_tol = colt.newt_tol;
newti = colt.newti;
max_iter = colt.max_iter;
max_iter_chk = colt.max_iter_chk;
atten_tol = colt.atten_tol;
atten = colt.atten;

% Parameters unique to deBoor mesh refinement
if strcmp(Mesh,'deBoor') == 1
    mini = colt.mini;
    maji = colt.maji;
    maxe = colt.maxe;
    maxdiffe = colt.maxdiffe;
end

% Parameters unique to CEP mesh refinement
if strcmp(Mesh,'CEP') == 1
    remi = colt.remi;
    addi = colt.addi;
    maxe = colt.maxe;
end

% Print header for output table if switched on
if strcmp('On',header) == 1
    fprintf('\n                                    Sparsity                                                          \n')
    fprintf(' Deg  Segs  len(X)  len(F)  Non-0s    (%%%%)        max(|dX|)         max(|F|)       Iter  \n')
end
    
% Calculate constants
[tau_nodes,A,A_inv,B,D,W] = CollSetup(N,NodeSpace); % LG is interpolation method

% Calculate segment time intervals, divide by 2, and vectorize
t_seg = repmat(reshape(diff(t_bnd),[1 1 n_seg]),[n_state,(N+1)/2])/2; % for variable nodes
t_seg_d = repmat(reshape(diff(t_bnd),[1 1 n_seg]),[n_state,(N-1)/2])/2; % for defect nodes

% Store collocation matrices in structure
collmat = struct;

% Copy constant matrices into three dimensional matrices
collmat.A = repmat(A,[1 1 n_seg]);
collmat.Ainv = repmat(A_inv,[1 1 n_seg]);
collmat.Bnew = repmat(B(:,2:end-1),[1 1 n_seg]);
collmat.B0 = repmat(B(:,1),[1 1 n_seg]);
collmat.Bf = repmat(B(:,end),[1 1 n_seg]);
collmat.Dnew = repmat(D,[1 1 n_seg]);
collmat.Wnew = repmat(diag(W).',[n_state 1 n_seg]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Sparse Matrix Setup %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate indices associated with nonzeros of the defect constraints in the jacobian
[iRowC,jColC] = JacIndC_LT(colt);

% Calculate indices associated with boundary constraints in the jacobian
if n_coast > 0 % if coast parameters are included
    [iRowB,jColB] = JacIndB_LT_Coast(colt); % NOTE: This function has not been modified to match the new JacIndB_LT structure
else
    [iRowB,jColB] = JacIndB_LT(colt);
end

% Assemble row and column indices into individual column vectors
iRow = [iRowC; iRowB]; 
jCol = [jColC; jColB];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Form Initial Constraint Vector %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate initial constraint vector
[F,x0,xf,C] = MakeF_LT(Z,t_seg,t_seg_d,collmat,colt);

% Calculate error
error = max(abs(F));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate Partial Derivatives of the Jacobian %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate lengths of design and constraint vectors
Zl = length(Z);
Fl = length(F);

% Compute full Jacobian with either Forward or Central Difference Methods
% [DF_fwrd] = MakeDF_Num(Z,Fl,t_seg,t_seg_d,collmat,colt);
% DF_fwrd = zeros(Fl,Zl);

% Compute only nonzeros of the Jacobian with Complex Step Differentiation
if n_coast > 0 % if coast parameters are included
    [DF] = MakeDF_LT_Coast(Z,t_seg,t_seg_d,collmat,colt); % NOTE: This function has not been modified to match the new MakeDF_LT structure
else % if no coast parameters are included
    [DF] = MakeDF_LT(Z,t_seg,t_seg_d,collmat,colt);
end

% Convert DF into a sparse matrix
% DF_sparse_fwrd = sparse(DF_fwrd); % for forward step case
DF_sparse=sparse(iRow,jCol,DF,Fl,Zl); % for complex step case

% Calculate length of DF_sparse and DF sparsity percentage
% DFl_fwrd = nnz(DF_sparse_fwrd); % for forward step case
DFl = length(DF); % for complex step case
sparsity = 100*(1-DFl/(Fl*Zl));

% % Compare forward step and sparse complex step results
% DF_comp = full(DF_sparse) - DF_fwrd; % difference
% DF_comp_pcnt = (DF_comp./DF_fwrd)*100; % percent difference
% max_diff = max(DF_comp(:));
% max_pcnt_diff = max(DF_comp_pcnt(:));
% 
% % Plot forward step and complex step Jacobians (for debugging only)
% figure (1)
% hold on
% grid on
% spy(DF_sparse,'b')
% spy(DF_sparse_fwrd,'or')
% hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Update Design Variables %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Apply update equation to calculate a new design variable vector
[dZ] = Update(F,DF_sparse,Zl,Fl);
maxdZ = max(abs(dZ)); % calculate max change in X for print output

newti = newti + 1; % iterate counter by one

%Print output for one iteration of Newton's method
if strcmp('NoMesh',Mesh) == 1 % if no mesh
    fprintf(' %2.0f  %4.0f  %6.0f  %6.0f  %6.0f  %8.2f  %16.8e  %16.8e  %5.0f \n',N,n_seg,Zl,Fl,DFl,sparsity,maxdZ,error,newti)
elseif strcmp('deBoor',Mesh) == 1 % if de Boor mesh
    fprintf(' %2.0f  %4.0f  %6.0f  %6.0f  %6.0f  %8.2f  %16.8e  %16.8e  %5.0f  %5.0f  %5.0f  %16.8e  %5.0f \n',N,n_seg,Zl,Fl,DFl,sparsity,maxdZ,error,newti,mini,maji,maxe,maxdiffe)
elseif strcmp('CEP',Mesh) == 1 % if CEP mesh
    fprintf(' %2.0f  %4.0f  %6.0f  %6.0f  %6.0f  %8.2f  %16.8e  %16.8e  %5.0f  %5.0f  %5.0f  %16.8e \n',N,n_seg,Zl,Fl,DFl,sparsity,maxdZ,error,newti,remi,addi,maxe)
end
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Begin Newton's Method %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

while error > newt_tol && ~max_iter_chk
    
    if max(dZ) > atten_tol
        Z = Z + (1/atten).*dZ; % take quarter step
    else
        Z = Z + dZ; % take a full step
    end

    % Calculate initial constraint vector
    [F,x0,xf,C] = MakeF_LT(Z,t_seg,t_seg_d,collmat,colt);

    % Calculate error
    error = max(abs(F));

    % Calculate lengths of design and constraint vectors
    Zl = length(Z);
    Fl = length(F);

    % Compute full Jacobian with either Forward or Central Difference Methods
%     [DF_fwrd] = MakeDF_Num(Z,Fl,t_seg,t_seg_d,collmat,colt);
%     DF_fwrd = zeros(Fl,Zl);

    % Compute only nonzeros of the Jacobian with Complex Step Differentiation
    if n_coast > 0 % if coast parameters are included
        [DF] = MakeDF_LT_Coast(Z,t_seg,t_seg_d,collmat,colt); % NOTE: This function has not been modified to match the new MakeDF_LT structure
    else % if no coast parameters are included
        [DF] = MakeDF_LT(Z,t_seg,t_seg_d,collmat,colt);
    end

    % Convert DF into a sparse matrix
%     DF_sparse_fwrd = sparse(DF_fwrd); % for forward step case
    DF_sparse=sparse(iRow,jCol,DF,Fl,Zl); % for complex step case

    % Calculate length of DF_sparse and DF sparsity percentage
%     DFl_fwrd = nnz(DF_sparse); % for forward step case
    DFl = length(DF); % for complex step case
    sparsity = 100*(1-DFl/(Fl*Zl));

    % Apply update equation to calculate a new design variable vector
    [dZ] = Update(F,DF_sparse,Zl,Fl);
    maxdZ = max(abs(dZ)); % calculate max change in X for print output

    newti = newti + 1; % iterate counter by one

    %Print output for one iteration of Newton's method
    if strcmp('NoMesh',Mesh) == 1 % if no mesh
        fprintf(' %2.0f  %4.0f  %6.0f  %6.0f  %6.0f  %8.2f  %16.8e  %16.8e  %5.0f \n',N,n_seg,Zl,Fl,DFl,sparsity,maxdZ,error,newti)
    elseif strcmp('deBoor',Mesh) == 1 % if de Boor mesh
        fprintf(' %2.0f  %4.0f  %6.0f  %6.0f  %6.0f  %8.2f  %16.8e  %16.8e  %5.0f  %5.0f  %5.0f  %16.8e  %5.0f \n',N,n_seg,Zl,Fl,DFl,sparsity,maxdZ,error,newti,mini,maji,maxe,maxdiffe)
    elseif strcmp('CEP',Mesh) == 1 % if CEP mesh
        fprintf(' %2.0f  %4.0f  %6.0f  %6.0f  %6.0f  %8.2f  %16.8e  %16.8e  %5.0f  %5.0f  %5.0f  %16.8e \n',N,n_seg,Zl,Fl,DFl,sparsity,maxdZ,error,newti,remi,addi,maxe)
    end    
    
    % Check if max iteration limit has been reached
    if newti == max_iter
        max_iter_chk = true;
        warning('The solver stopped because the Newton''s method maximum iteration limit was reached')
    end
    
end

%% Output %%

%Create vector of boundary nodes 
x_bnd = cat(3,x0,xf(:,:,end));

% Save newti and max_ter_chk in colt structure
colt.newti = newti;
colt.max_iter_chk = max_iter_chk;
