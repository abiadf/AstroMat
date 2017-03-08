function bodyplot(imdir,name,r,x0,y0,z0,opaq,rot_ax,rot_deg)
%BODYPLOT  Plot a spherical surface map from given image.
%   BODYPLOT(IMDIR,NAME,R,X0,Y0,Z0,OPAQ) reads an image from the directory 
%   IMDIR beginning with the string NAME.  The body's radius is input as R,
%   and the body's origin is input with cartesian coordinates X0, Y0, Z0.
%   Lighting is added with a surface opaqueness of OPAQ percent.  If
%   BODYPLOT is called without the OPAQ argument, the surface is plotted
%   with zero transparency.
%
%   Example:
%       bodyplot('Textures\','Earth',6378,0,0,0,0.9);
%
%   See also SURFACE, IMREAD, IMFINFO.

%   Mark C Jesick
%   University of Texas
%   May 3, 2008

%   Updates:
%   July 8, 2009    -   Removed any need for image processing software in
%                   	image reorientation section.

%-------------------------------------------------------------------------%
% Input Arguments Check
%-------------------------------------------------------------------------%
if nargin == 6 % if the surface transparency isn't specified
    opaq = 1; % make surface opaque
elseif nargin < 6 % if incorrect number of input arguments
    error('Not enough input arguments.'); % terminate programme
end % end input arguments block


%-------------------------------------------------------------------------%
% Get Texture Image File
%-------------------------------------------------------------------------%
ext = {'bmp' 'jpg' 'png' 'tif'}; % possible image extensions
next = length(ext); % number of image extensions
for i = 1:next % for each image extension
    imname = [imdir name 'Texture.' ext{i}]; % concatenate image name
    if exist(imname,'file') % if the image exists
        break % exit image extensions loop
    end % end extension block
end % end image extensions loop
body = imread(imname); % load texture map
inform = imfinfo(imname); % get texture image info


%-------------------------------------------------------------------------%
% Properly Orient Image
%-------------------------------------------------------------------------%
dims = ndims(body); % number of dimensions of body
if dims == 2 % if image has a bit depth of two
    body = flipud(body); % mirror image vertically
elseif dims == 3 % if image has a bit depth of three  
    for i = 1:dims % for each dimension
        body(:,:,i) = flipud(body(:,:,i)); % mirror dimension vertically
    end % end dimension loop
end % end surface orientation block

%-------------------------------------------------------------------------%
% Plot Image on Sphere
%-------------------------------------------------------------------------%
n = 50; % number of longitudinal and latitudinal sphere faces
[x y z] = sphere(n); % calculate coordinates for sphere
s = surface(r*x + x0,r*y + y0,r*z + z0,... % create spherical surface
    'facecolor','texturemap',... % specify an image will be used
    'cdata',body); % add image data to sphere surface

% Prevent image display irregularities due to the color map
if strcmp(ext{i},'bmp') && ... % if this is a bitmap image type
        ~isempty(inform.Colormap) % and if image specifies a colormap
    colormap(inform.Colormap); % use the image's colormap
end % end colormap block

% Specify surface properties and add lighting
set(s,'edgecolor','none',... % delete edges from sphere faces
    'FaceLighting','gouraud',... % body lighting scheme
    'DiffuseStrength',1,... % diffusivity of light on body
    'AmbientStrength',1,... % strength of overall scene lighting
    'FaceAlpha',opaq); % surface transparency (0 = invisible, 1 = opaque)
light; % add lighting source

%-------------------------------------------------------------------------%
% Rotate Sphere
%-------------------------------------------------------------------------%

rotate(s,rot_ax,rot_deg);

