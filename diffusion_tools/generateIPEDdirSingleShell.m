% 
% IPED - Improved B0-distortion correction in diffusion MRI
% Copyright (C) 2013-2023 C Bhushan, D Varadarajan, AA Joshi, RM Leahy, and JP Haldar.
% 
% This work is released under either of Apache-2.0 OR GPL-2.0 licenses. 
% Please see https://github.com/cbhushan/IPED for details.
% 
% SPDX-License-Identifier: Apache-2.0 OR GPL-2.0-only
% 


function [pos, theta, phi] = generateIPEDdirSingleShell(nDirections, nPED, opts)
% Generate uniformly distributed diffusion encoding directions for interlaced PED. nDirections and
% nPED are scalar numbers such that nDirections is a multiple of nPED. POS is a 3D matrix of size
% (nDirections/nPED x nPED x 3). See default options below. 
%

defaultoptions = struct(...
   'dt', 0.01, ...
   'pos_init', 'spiral', ... archimedes/spiral/rand/init_matrix
   'maxIter', 10000, ...
   'thetaTol', 0.00003, ... in radians
   'alpha', 0.5 ...
   );

if rem(nDirections, nPED)~=0
   error('nDirections must be multiple of nPED')
end

if(~exist('opts','var')),
   opts = defaultoptions;
else
   tags = fieldnames(defaultoptions);
   for i = 1:length(tags)
      if ~isfield(opts, tags{i})
         opts.(tags{i}) = defaultoptions.(tags{i});
      end
   end
end

% initialize positions
N = nDirections/nPED;
if ischar(opts.pos_init)   
   % alternating PED index
   indp = reshape(1:nDirections, nPED, [])';
   indp = indp(:);
   
   switch opts.pos_init
      case 'archimedes'
         pos = archimedes_init(nDirections);
         pos = pos(randperm(nDirections), :);
      case 'spiral'
         pos = spiral_init(nDirections);
         pos = pos(indp, :);
      case 'rand'
         pos = rand_uniform_init(nDirections);         
         pos = pos(indp, :);
      otherwise
         error('opts.pos_init must be either archimedes/spiral/rand')
   end
   
elseif isequal(size(opts.pos_init), [N, nPED, 3])   
   pos = reshape(opts.pos_init, nDirections, 3);
   
else
   error('size of opts.pos_init is not correct!')
   
end

% visualize initialization
h = displayDiffusionDirections(pos(1:N,:));
for p = 2:nPED
   st = (p-1)*N + 1;
   ed = (p-1)*N + N;
   displayDiffusionDirections(pos(st:ed,:), rand(1, 3), h);
end

% Minimize electrostatic repulsion
exitflag = false;
nIter = 0;
a1 = 1-opts.alpha;
a2 = opts.alpha;

while ~exitflag
   % due all direction
   [F_all, F_all_tangent] = netSymmetricESforce(pos);
   F_all = F_all./size(pos,1);
   F_all_tangent = F_all_tangent./size(pos,1);
   
   % due to direction from same PED
   if nPED>1
      F_PED = zeros(size(pos));
      F_PED_tangent = zeros(size(pos));
      for p = 1:nPED
         st = (p-1)*N + 1;
         ed = (p-1)*N + N;
         ptemp = squeeze(pos(st:ed,:));
         [F, Ft] = netSymmetricESforce(ptemp);
         F_PED(st:ed,:) = F./size(ptemp,1);
         F_PED_tangent(st:ed,:) = Ft./size(ptemp,1);
      end
   else
      F_PED = F_all;
      F_PED_tangent = F_all_tangent;
   end
   
   F_total = a1*F_PED + a2*F_all;
   F_total_tangent = a1*F_PED_tangent + a2*F_all_tangent;
   
   axisOfRot = cross(pos, F_total, 2);
   F_total_tangent_mag = sqrt(sum(F_total_tangent.^2, 2));
   angleOfRot = F_total_tangent_mag*(opts.dt^2);
   
   nextPos = rotateVec3DAxis(pos, axisOfRot, angleOfRot);
   pos = nextPos;
   
   nIter = nIter+1;   
   theta_change = max(abs(angleOfRot));
   if theta_change<opts.thetaTol
      exitflag = 1;
      fprintf('\nMax theta is smaller than thetaTol.\n');
   elseif nIter>opts.maxIter
      exitflag = 2;
      fprintf('\nReached maximum number of allowed iterations.\n');
   end
end

% Rotate the found solution so that first one point to z-axis [0,0,1]
axisOfRot = cross(pos(1,:),[0 0 1]);
angleOfRot =  acos(pos(1,3));
pos = rotateVec3DAxis(pos, axisOfRot(ones(1,nDirections),:), angleOfRot*ones(nDirections,1));
pos = reshape(pos, N, nPED, 3); % reshape to make 3D matrix

if nargout>1
   [theta, phi, ~] = cart2sph(pos(:,:,1), pos(:,:,2), pos(:,:,3));
end

end

%
% Some related references (not implemented here)
%   http://www.math.niu.edu/~rusin/known-math/95/sphere.faq
%   http://eqsp.sourceforge.net/
%   http://www.mathworks.com/matlabcentral/fileexchange/37004-uniform-sampling-of-a-sphere


function pos = spiral_init(N, perturb)
% Returns N points which are approx. evenly distributed one hemishphere. 2*N points are sampled
% from spiral using algorithm described by Saff and Kuijlaars and only first half is returned.
% A small pertubation is optionally added to avoid huge energy when using netSymmetricESforce().
%
% This arranges the nodes along a spiral in such a way that the distance between nodes along the
% spiral is approximately equal to the distance between coils of the spiral. 
% Another way to understand this method: Divide the sphere along lines of "latitude" into parallel
% bands of equal area and place a node at some "longitude" in the middle of each band.
%
% Reference: 
%    Saff and Kuijlaars, "Distributing many points on a sphere", The Mathematical Intelligencer,
%    1997, DOI: 10.1007/BF03024331   

if ~exist('perturb', 'var')
   perturb = true;
end

N = 2*N; % Choose 2N points on sphere

k = 1:N;
hk = (2*(k-1)/(N-1)) - 1;
theta = acos(hk); % [0 pi]

phi = zeros(1, N); % [0 2*pi]
temp = sqrt(N*(1-(hk.^2)));
for k = 2:(N-1)
   phi(k) = mod(phi(k-1)+(3.6/temp(k)), 2*pi);
end

[X, Y, Z] = sph2cart(phi, theta-(pi/2), 1);
pos = [X(:) Y(:) Z(:)];
pos = pos(1:N/2, :); % first half

if perturb
   temp = pos + 0.15*(rand(size(pos))-0.5); % small random perturbation 
   temp_m = sqrt(sum(temp.^2, 2));
   pos = temp./temp_m(:, [1 1 1]);
end

end


function pos = spiral2_init(N)
% OLDERE implementation does not match paper! Was taken from:  
%  http://web.archive.org/web/20120421191837/http://www.cgafaq.info/wiki/Evenly_distributed_points_on_sphere
%
% Returns 2N points which are approx. evenly distributed on surface of sphere. 2N points are sampled
% from spiral using algorithm described by Saff and Kuijlaars.
%
% This arranges the nodes along a spiral in such a way that the distance between nodes along the
% spiral is approximately equal to the distance between coils of the spiral. 
% Another way to understand this method: Divide the sphere along lines of "latitude" into parallel
% bands of equal area and place a node at some "longitude" in the middle of each band.
% This implementation chooses successive longitudes according to the "most irrational number" (known
% as the golden section) so that no two nodes in nearby bands come too near each other in longitude.
%
% Reference: 
%    Saff and Kuijlaars, "Distributing many points on a sphere", The Mathematical Intelligencer,
%    1997, DOI: 10.1007/BF03024331   

dlong = pi*(3-sqrt(5)); %  /* ~2.39996323 */
dz = 2.0/N;
long = 0;
z = 1-dz/2;

X = zeros(1, N*2);
Y = zeros(1, N*2);
Z = zeros(1, N*2);

for k = 1:N
  r  = sqrt(1-z*z);
  X(2*k-1) = cos(long)*r;
  X(2*k) = -X(2*k-1);
  Y(2*k-1) = sin(long)*r;
  Y(2*k) = -Y(2*k-1);
  Z(2*k-1) = z;
  Z(2*k) = -z;
  z = z - dz;
  long = long + dlong;
end
pos = [X(:) Y(:) Z(:)];

end


function pos = archimedes_init(N)
% Returns 2*N points which are approx. evenly distributed on surface of sphere. N points are sampled
% using stratified sampling using Archimedes' Theorem only on top hemisphere (algorithm described
% below) and then for each of the sampled point the diamterically opposite point is also included.
%
% Archimedes' Theorem (Local): The axial projection of any measurable region on a sphere on the
%     right circular cylinder circumscribed about the sphere preserves area.
%
% By above theorem, for any two measurable regions S1 and S2 on the unit sphere with equal areas,
% their axial projections S1i and S2i on the circumscribed cylinder will have equal areas. The reverse
% is also true because the axial projection is a bijection (except for the two poles and their
% corresponding base circles which all have zero area). This naturally leads to a spherical sampling
% algorithm: Generate a random point on the cylinder [-1, 1]x[0, 2*pi] and then find its inverse
% axial projection on the unit sphere. If a random point is uniformly distributed on the cylinder,
% by the above argument, its inverse axial projection will be uniformly distributed on the sphere.
%
% This algorithm can be enhanced by applying stratified sampling over cylinder i.e. first sub-divide 
% cylinder into homogeneous areas before sampling. In this implementation, stratified sampling is
% modified so that points are sampled only on the top half of the sphere, 
% i.e. cylinder [0, 1]x[0, 2*pi] & then diameterically opposite points are added later.
%
% Reference:
%    Shao & Badler, 1996, Spherical Sampling by Archimedes' Theorem

% Sub divide the unfolded-cylinder into pxq sub-rectangles such that p:q = 2:(2*pi)
p = ceil(sqrt(2*N/pi)); % 2*N because we really want to select 2N points
q = ceil(sqrt(2*N*pi));
dp=2/p;
dq=2/q;
[Xc, Yc] = ndgrid(linspace(dp/4, (1-dp/2), ceil(p/2)), (-1+dq/2):dq:(1-dq/2));

% Random sample in sub-rectangle
x = Xc(:) + dp*(rand(ceil(p/2)*q,1)-0.5);
y = Yc(:) + dq*(rand(ceil(p/2)*q,1)-0.5);
clear Xc Yc

% Remove excess samples
R = ceil(p/2)*q-N;
if R>0
   idx=randperm(ceil(p/2)*q,R);
   x(idx)=[];
   y(idx)=[];
end

% Convert z to latitude
lon = (y+1)*pi;
z = x;
lat = acos(z);

% Convert spherical to rectangular co-ords
X = cos(lon).*sin(lat);
Y = sin(lon).*sin(lat);
pos = [X Y z];

end


function pos = rand_uniform_init(N)
% Returns N points which are approx. evenly distributed on surface of sphere. N points are sampled
% from a distribution which is uniformly distribute over surface of sphere. This method gives
% 'good' results with large N. 
%
% Multivariate normal distribution with identity covariance matrix is rotationally symmetric around
% the origin. Hence, if samples are picked up from this distribution and normalized to be unit
% length then they will be samples from a distribution which is uniformly distributed on surface of
% unit n-sphere. 
%
% Reference: (Muller 1959, Marsaglia 1972).
%     http://mathworld.wolfram.com/SpherePointPicking.html

p = randn([N, 3]); % standard normal distribution
mag_p = sqrt(sum(p.^2, 2));
pos = p ./ mag_p(:,[1 1 1]);

end


