% 
% IPED - Improved B0-distortion correction in diffusion MRI
% Copyright (C) 2013-2023 C Bhushan, D Varadarajan, AA Joshi, RM Leahy, and JP Haldar.
% 
% This work is released under either of Apache-2.0 OR GPL-2.0 licenses. 
% Please see https://github.com/cbhushan/IPED for details.
% 
% SPDX-License-Identifier: Apache-2.0 OR GPL-2.0-only
% 


function [pos, theta, phi] = generateIPEDdirMultiShell(nDirections, nPED, nS, bval, opts)
% Generate uniformly distributed diffusion encoding directions for interlaced PED. nDirections and
% nPED are scalar numbers such that nDirections is a multiple of nPED. POS is a 3D matrix of size
% (nDirections/nPED x nPED x 3). See default options below. 
%
% nDirections = sum(nS)
% all nS should be multiple of nPED
%

defaultoptions = struct(...
   'dt', 0.01, ...
   'pos_init', 'spiral', ... archimedes/spiral/rand/init_matrix
   'maxIter', 20000, ...
   'thetaTol', 5e-6, ... in radians
   'alpha', [1 1 1 1]./4, ...
   'indp', []...
   );

if nDirections ~= sum(nS)
   error('nDirections must follow nDirections = sum(nS)')
end

if max(abs(rem(nS, nPED))) > 0
   error('All of nS must be multiple of nPED')
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
% Shells can have different number of directions.
for ind = 1:size(nS,2),
    N(ind) = nS(ind)/nPED;
end;
Ntot = nDirections/nPED;  % Total number for 1 PE

if ischar(opts.pos_init)
    indp = opts.indp;
%    indp = [];
%    sz=0;
%    sy=0;
%    for ind = 1:size(N,2),
%        for p = 1:N(ind)
%           indp = [indp; [sy+sz + p:N(ind): sy+sz + nS(ind)]'];
%           if(size(indp,1) == N(ind))
%               sy = Ntot;
%           end;
%        end
%        sz = sz + nS(ind);
%    end;
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
elseif isequal(size(opts.pos_init), [sum(N), nPED, 3])
   pos = reshape(opts.pos_init, nDirections, 3);   
else
   error('size of opts.pos_init is not correct!')
   
end

% Visualize initialization
h = displayDiffusionDirections(pos);
title('All directions');
sz =0; %size of each PE for calculation simplification
for ind = 1:size(N,2),
    h = displayDiffusionDirections(pos(sz+1:sz+N(ind),:));
    for p = 2:nPED
       st = (p-1)*Ntot + sz+1;
       ed = (p-1)*Ntot + sz+ N(ind);
       displayDiffusionDirections(pos(st:ed,:), [1 0 0], h);
    end
    sz=sz + N(ind);
end;


% Minimize electrostatic repulsion
opts.alpha = opts.alpha./norm(opts.alpha);
a1 = opts.alpha(1);
a2 = opts.alpha(2);
a3 = opts.alpha(3);
a4 = opts.alpha(4);
nShells = length(nS);
% nS2 = [0; nS(:); nDirections]./nPED; % for easy indexing 
nS2 = [0; cumsum(nS)']./nPED;
exitflag = false;
nIter = 0;
while ~exitflag
   
   % due to directions from same PED on each shell seperately
   F_PED_shell = zeros(size(pos));
   F_PED_shell_tangent = zeros(size(pos));
   for p = 1:nPED
      for s = 1:nShells
         st = (p-1)*Ntot + nS2(s) + 1;
         ed = (p-1)*Ntot + nS2(s+1);
         ptemp = squeeze(pos(st:ed,:));
         [F, Ft] = netSymmetricESforce(ptemp);
         F_PED_shell(st:ed,:) = F./(size(ptemp,1).^2);
         F_PED_shell_tangent(st:ed,:) = Ft./(size(ptemp,1).^2);
      end
   end
   
   % due to directions from same PED across all shells
   F_PED = zeros(size(pos));
   F_PED_tangent = zeros(size(pos));
   for p = 1:nPED
      st = (p-1)*Ntot + 1;
      ed = (p-1)*Ntot + Ntot;
      ptemp = squeeze(pos(st:ed,:));
      [F, Ft] = netSymmetricESforce(ptemp);
      F_PED(st:ed,:) = F./(size(ptemp,1).^2);
      F_PED_tangent(st:ed,:) = Ft./(size(ptemp,1).^2);
   end
   
   % due to all directions on each shell seperately 
   F_shell = zeros(size(pos));
   F_shell_tangent = zeros(size(pos));
   for s = 1:nShells
      temp = (nS2(s)+1):nS2(s+1);
      shind = [];
      for p = 0:(nPED-1)
         shind = [shind temp+p*Ntot];
      end
      ptemp = squeeze(pos(shind,:));
      [F, Ft] = netSymmetricESforce(ptemp);
      F_shell(shind,:) = F./(size(ptemp,1).^2);
      F_shell_tangent(shind,:) = Ft./(size(ptemp,1).^2);
   end
   
   % due all direction across all shells
   [F_all, F_all_tangent] = netSymmetricESforce(pos);
   F_all = F_all./(size(pos,1).^2);
   F_all_tangent = F_all_tangent./(size(pos,1).^2);
   
   F_total = a1*F_PED_shell + a2*F_PED  + a3*F_shell + a4*F_all;
   F_total_tangent =  a1*F_PED_shell_tangent  + a2*F_PED_tangent + a3*F_shell_tangent + a4*F_all_tangent;
   
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

nIter
theta_change
sz=0;
for ind = 1:size(N,2),
h = displayDiffusionDirections(pos(sz+1:sz+N(ind),:));
for p = 2:nPED
   st = (p-1)*Ntot + 1 + sz;
   ed = (p-1)*Ntot + sz+ N(ind);
   displayDiffusionDirections(pos(st:ed,:), [1 0 0], h);
end
sz=sz+N(ind);
end;
% Rotate the found solution so that first one point to z-axis [0,0,1]
axisOfRot = cross(pos(1,:),[0 0 1]);
angleOfRot =  acos(pos(1,3));
pos = rotateVec3DAxis(pos, axisOfRot(ones(1,nDirections),:), angleOfRot*ones(nDirections,1));
pos = reshape(pos, Ntot, nPED, 3); % reshape to make 3D matrix

% Normalize the directions with b-values
maxb = max(bval);
bv=[];
for ind = 1:size(N,2)
    bv = [bv; bval(ind)*ones(N(ind),1)];
end;
bfac = sqrt(bv./maxb);
pos = repmat(bfac, [1  nPED 3]).*pos;
if nargout>1
   [theta, phi, ~] = cart2sph(pos(:,:,1), pos(:,:,2), pos(:,:,3));
end

end


% function pos = spiral_init(N)
% % Returns 2*N points which are approx. evenly distributed on surface of sphere. N points are sampled
% % from spiral using algorithm described below and for each of the point the diamterically opposite
% % point is also included.
% %
% % This arranges the nodes along a spiral in such a way that the distance between nodes along the
% % spiral is approximately equal to the distance between coils of the spiral. 
% % Another way to understand this method: Divide the sphere along lines of "latitude" into parallel
% % bands of equal area and place a node at some "longitude" in the middle of each band.
% % This implementation chooses successive longitudes according to the "most irrational number" (known
% % as the golden section) so that no two nodes in nearby bands come too near each other in longitude.
% %
% % Reference: 
% %    Saff and Kuijlaars, "Distributing many points on a sphere", The Mathematical Intelligencer,
% %    1997, DOI: 10.1007/BF03024331   
% 
% dlong = pi*(3-sqrt(5)); %  /* ~2.39996323 */
% dz = 2.0/N;
% long = 0;
% z = 1-dz/2;
% 
% X = zeros(1, N*2);
% Y = zeros(1, N*2);
% Z = zeros(1, N*2);
% 
% for k = 1:N
%   r  = sqrt(1-z*z);
%   X(2*k-1) = cos(long)*r;
%   X(2*k) = -X(2*k-1);
%   Y(2*k-1) = sin(long)*r;
%   Y(2*k) = -Y(2*k-1);
%   Z(2*k-1) = z;
%   Z(2*k) = -z;
%   z = z - dz;
%   long = long + dlong;
% end
% pos = [X(:) Y(:) Z(:)];
% 
% end


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
