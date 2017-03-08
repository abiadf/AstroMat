function [dX] = Update(F,DF,Xl,Fl)
% function [dX] = Update(F,DF,Xl,Fl)
% 
% This function contains update equations for Newton's method. Note that if 
% the Jacobian is not a square matrix then the minimum norm rule is used.
%
% INPUTS:
%    F      constraint vector, column vector
%    DF     jacobian, sparse or full matrix data type
%    Xl     length of design vector
%    Fl     length of constraint vector
%
% OUTPUTS:
%    dX     values to update design vector, column vector
%
% Originally Written by: R. Pritchett, 07/15/2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Apply Newton's Method or Minimum-Norm depending on dimensions
if Fl == Xl
    %Newton's Method
    dX = -DF\F;
else
    %Minimum-Norm Equation
    dX = -DF'*((DF*DF')\F);
end