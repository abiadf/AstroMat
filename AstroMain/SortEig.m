function [eig_sorted,evec_sorted] = SortEig(eigs,evec)
% SORTEIG Sort a set of eigenvalues to acheive a consistent order
%
%   [OUT,ILS] = sorteig(EIGS,X) ; 
%       EIGS is the m-by-n array of n eigenvalues of the monodromy matrix
%       X (optional) is the x coordinate of the crossing.  If x is not an 
%       input, the code will use the indices in interpolation 
%
%       OUT is the sorted set of eigenvalues
%       ILS is the indicies of the IN matrix corresponding to OUT matrix
%
%   Authors:
%       Geoff Wawrzyniak, Purdue University, 2006
%       Andrew Cox, Purdue University, 2015
%       Robert Pritchett, Purdue University, 2017

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Loop Through Each Set of Eigenvalues %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Permutations of eigenvalue and eigenvector columns indices
eig_ind_perm = perms([1:6]); % permutations of column indices

% Define total number of permutations for a set of six eigenvalues
n_perm = size(eig_ind_perm,1);

% Calculate number of eigenvalue sets
n_eig = size(eigs,1);

% Preallocate storage vectors for eigenvalues and eigenvectors
eig_sorted = zeros(n_eig,6);
evec_sorted = zeros(6,6,n_eig);

for ii = 1:n_eig  % for each row of eigs  
    
    % Reset total eigenvalue cost with each iteration
    tot_cost = zeros(n_perm,1);

    %---------------------------------------------------------------------%
    % Generate all 720 Permutations of Eigenvalues (and Vectors) in Current Row
    %---------------------------------------------------------------------%

    % Define eigenvalues and eigenvectors for current iteration
    eig_i = eigs(ii,:);
    evec_i = evec(:,:,ii); % matrix of eigenvectors for current iteration
    
    % Permutations of eigenvalue and eigenvector columns
    eig_perm = zeros(n_perm,6); % preallocate matrix of eigenvalue permutations
    evec_perm = zeros(6,6,n_perm); % preallocate matrix of eigenvector permutations

    % Permute column indices
    for jj = 1:n_perm
        eig_perm(jj,:) = eig_i(eig_ind_perm(jj,:)); 
        evec_perm(:,:,jj) = evec_i(:,eig_ind_perm(jj,:)); 
    end   
    
    %---------------------------------------------------------------------%
    % Check Sequential Pairs
    %---------------------------------------------------------------------%
    
    % Calculate difference between 1 and products of pairs 
    pair_chk = abs(ones(n_perm,3) - eig_perm(:,[1 3 5]).*eig_perm(:,[2 4 6]));
    
    % Assign cost based on sum of differences calculated in previous step
    pair_cost = sum(pair_chk,2);
    
    % Sort cost
%     [pair_cost_sort,pair_cost_ind] = sort(pair_cost);
    
    % Add pair cost to total eig_cost sum
    tot_cost = tot_cost + pair_cost;
    
    %---------------------------------------------------------------------%
    % Check if Subsequent Eigenvectors are Parallel
    %---------------------------------------------------------------------%
    
    if ii == 1
        evec_cost = zeros(n_perm,1); % No previous eigenvectors to compare with
    else
        % Create 3D matrix using eigenvectors of previous sort iteration
        evec_prev = repmat(evec_sorted(:,:,ii-1),[1 1 n_perm]);
        
        % Calculate dot products 
        evec_chk = ones(1,6,n_perm) - abs(dot(evec_prev,evec_perm,1));
        
        % Sum and reshape dot product results into column vector
        evec_cost = reshape(sum(evec_chk,2),[n_perm 1]);
        
        % Sort cost
%         [evec_cost_sort,evec_cost_ind] = sort(evec_cost);
     
    end

    % Add eigenvector cost to total eig_cost sum
    tot_cost = tot_cost + evec_cost;
    
    %---------------------------------------------------------------------%
    % Check if Subsequent Eigenvalues are Very Different
    %---------------------------------------------------------------------%
    
    if ii == 1
        eig_cost = zeros(n_perm,1); % No previous eigenvectors to compare with
    else
        % Create matrix using eigenvalues of previous sort iteration
        eig_prev = repmat(eig_sorted(ii-1,:),[n_perm 1]);
        
        % Calculate difference 
        eig_chk = abs(eig_prev - eig_perm);
        
        % Sum and reshape dot product results into column vector
        eig_cost = sum(eig_chk,2);
        
        % Sort cost
%         [eig_cost_sort,eig_cost_ind] = sort(eig_cost);
     
    end
    
    % Add eigenvector cost to total eig_cost sum
    tot_cost = tot_cost + eig_cost;
    
    %---------------------------------------------------------------------%
    % Determine Permutation with Minimum Cost and Store Result
    %---------------------------------------------------------------------%
    
    [~,min_ind] = min(tot_cost);
    eig_sorted(ii,:) = eig_perm(min_ind,:);
    evec_sorted(:,:,ii) = evec_perm(:,:,min_ind);
    
end
    
