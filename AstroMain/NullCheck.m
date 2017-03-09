function [delt_Z,null_chk] = NullCheck(delt_Z,delt_Z_old,null_chk)
% function [delt_Z] = NullCheck(delt_Z,delt_Z_old,null_chk)
% 
% This function checks the sign and dimension of the null vector computed 
% in a pseudo-arclength continuation scheme. After performing these checks
% it outputs the null vector needed to take the next step forward in the 
% family being computed via the continuation scheme.
%
% INPUTS:
%    delt_Z        null vector of current pseudo-arclength step 
%    delt_Z_old    null vector of previous pseudo-arclength step
%    null_chk      structure containing parameters relevant to the null_chk
%                  function
%
% OUTPUTS:
%    delt_Z        null vector selected for next pseudo-arclength step 
%                  following the sign and dimension checks of this function
%    null_chk      structure containing parameters relevant to the null_chk
%                  function
% 
% Written by R. Pritchett, 02/14/17
% Last Update: R. Pritchett, 02/14/17
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Check for Bifurcation %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Explanation: A bifurcation occurs when the dimension of the null space
% is greater than 1, i.e. the null space is no longer a single column 
% vector, but rather a matrix of two or more column vectors. In this case 
% it is desirable to remain in the same family as pseudo-arclength 
% continuation is performed. To ensure this, whenever a bifurcation occurs 
% the dot product of each of the null vectors is calculated with respect to
% the previous null vector. The smallest resulting dot product value 
% corresponds to the null space vectors with the smallest angle between 
% them. It is assumed that the new null space vector that has the smallest
% angle with respect to the old null space vector is the null space vector 
% that remains in the originally computed family. Therefore, the new null 
% space vector corresponding to the smallest angle is used to calculate 
% the next guess in the continuation method.

% Calculate the number of columns in the null space
null_dim = size(delt_Z,2);

% Check if the number of columns is greater than 1 
if null_dim > 1

    % Iterate counter for each bifurcation recorded
    null_chk.bif_cnt = null_chk.bif_cnt + 1;

    % Calculate angle between null vectors and previous null vector
    null_ang = zeros(null_dim,1);
    for jj = 1:null_dim
        null_ang(jj,1) = acos(dot(delt_Z(:,jj),delt_Z_old));
    end

    % Determine smallest angle between null vectors
    [~,min_ind] = min(abs(null_ang));

    % Save relevant information for bifurcation
    null_chk.bif_Z{null_chk.bif_cnt} = null_chk.Zpls;
    null_chk.bif_null{null_chk.bif_cnt} = delt_Z;
    null_chk.bif_ind{null_chk.bif_cnt} = null_chk.ps_ind;

    % Select null vector with smallest angle wrt previous null vector
    delt_Z = delt_Z(:,min_ind);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Ensure Null Vector Maintains the Same Sign %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Explanation: The sign of the null vector is arbitrary, therefore when
% Matlab calculates the null vector with the "null" function it can
% sometimes flip the sign arbitrarily. In pseudo-arclength continuation the
% null vector is used to step along the solution space of the orbital
% family, and it is desirable that each subsequent step be in the same
% direction. Therefore, the signs of the elements of the null vector with
% respect to the signs of the same elements of the previous null vector are
% compared. If more than 5 elements (# arbitrarily chosen) of the null
% vector change sign between iterations then it is determined that the sign
% has been flipped by Matlab. This changing of the sign is then countered
% by multiplying the entire null space vector by -1 to ensure the current
% null space vector steps the algorithm in the same direction as the
% previous null space vector. Note, that elements of the null vector very
% near machine precision can arbitrarily change sign, therefore these 
% elements are ignored when determining whether the sign of the null vector
% has flipped.

% Calculate angle between null vectors and previous null vector
null_ang = zeros(2,1);
delt_Z_dir = [-delt_Z delt_Z];
for jj = 1:2
%     null_ang(jj,1) = acos(dot(delt_Z_dir(:,jj),delt_Z_old));
    null_ang(jj,1) = 1-dot(delt_Z_dir(:,jj),delt_Z_old);
end

% Determine smallest angle between null vectors
% min_ind = find(abs(null_ang) < pi/2);
[~,min_ind] = min(abs(null_ang));

if min_ind == 1
    null_chk.flip_cnt = null_chk.flip_cnt + 1;
    null_chk.flip_iter = [null_chk.flip_iter; null_chk.ps_ind];
end

% Select null vector with smallest angle wrt previous null vector
delt_Z = delt_Z_dir(:,min_ind);
