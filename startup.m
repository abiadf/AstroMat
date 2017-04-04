% startup.m
%
% Description: This script runs when a new Matlab session is first started.
% It can be used to automatically define variables, run scripts, or add
% folders to the search path.

% Add library of core astrodynamics functions to the Matlab search path at startup
addpath C:\Users\Robert\Documents\GitRepo_AstroMat\AstroMain
addpath C:\Users\Robert\Documents\GitRepo_AstroMat\AstroMain\include
addpath C:\Users\Robert\Documents\GitRepo_AstroMat\AstroMain\Textures

% Add COLT library to the Matlab search path at startup
addpath C:\Users\Robert\Documents\GitRepo_AstroMat\COLT

% Add SNOPT library to the Matlab search path at startup
addpath C:\Users\Robert\Documents\snopt-matlab-2.3.0\matlab
addpath C:\Users\Robert\Documents\snopt-matlab-2.3.0\matlab\examples
addpath C:\Users\Robert\Documents\snopt-matlab-2.3.0\matlab\examples\fmincon
addpath C:\Users\Robert\Documents\snopt-matlab-2.3.0\matlab\examples\hs76
addpath C:\Users\Robert\Documents\snopt-matlab-2.3.0\matlab\examples\hs116
addpath C:\Users\Robert\Documents\snopt-matlab-2.3.0\matlab\examples\hsmain
addpath C:\Users\Robert\Documents\snopt-matlab-2.3.0\matlab\examples\snmain
addpath C:\Users\Robert\Documents\snopt-matlab-2.3.0\matlab\examples\sntoy
addpath C:\Users\Robert\Documents\snopt-matlab-2.3.0\matlab\examples\spring
addpath C:\Users\Robert\Documents\snopt-matlab-2.3.0\matlab\examples\t1diet
