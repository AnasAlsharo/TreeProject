% This file is part of TREEQSM.
%
% TREEQSM is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% TREEQSM is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with TREEQSM.  If not, see <http://www.gnu.org/licenses/>.

function inputs = define_input(Clouds,nPD1,nPD2Min,nPD2Max)

% ---------------------------------------------------------------------
% DEFINE_INPUT.M       Defines the required inputs (PatchDiam and BallRad 
%                        parameters) for TreeQSM based in estimated tree
%                        radius.
%
% Version 1.0.0
% Latest update     4 May 2022
%
% Copyright (C) 2013-2022 Pasi Raumonen
% ---------------------------------------------------------------------

% Takes in a single tree point clouds, that preferably contains only points 
% from the tree and not e.g. from groung. User defines the number of
% PatchDiam1, PatchDiam2Min, PatchDiam2Max parameter values needed. Then
% the code estimates automatically these parameter values based on the 
% tree stem radius and tree height. Thus this code can be used to generate
% the inputs needed for QSM reconstruction with TreeQSM.
%
% Inputs:
% P         Point cloud of a tree OR string specifying the name of the .mat
%             file where multiple point clouds are saved              
% nPD1      Number of parameter values estimated for PatchDiam1
% nPD2Min   Number of parameter values estimated for PatchDiam2Min
% nPD2Max   Number of parameter values estimated for PatchDiam2Max
%
% Output:
% inputs    Input structure with the estimated parameter values
% ---------------------------------------------------------------------


% Create inputs-structure
create_input
Inputs = inputs;

% If given multiple clouds, extract the names
if ischar(Clouds) || isstring(Clouds)
    matobj = matfile(Clouds);
    if isprop(matobj, 'trees') % Check for 'trees' instead of 'P'
        trees = matobj.trees; % Assign the variable 'trees' instead of 'P'
        nt = numel(trees); % Determine the number of trees
    else
        error('The .mat file does not contain trees variable.');
    end
else
  P = Clouds; % Assuming Clouds is the variable with your point clouds
  nt = numel(P); % Determine the number of point clouds
end
% Now, inputs should be an array of structures, one for each tree
for i = 1:nt
    inputs(i).PatchDiam1 = 0; % Initialize PatchDiam1 for each tree

end



%% Estimate the PatchDiam and BallRad parameters
for i = 1:nt
  if nt > 1
    % Select point cloud
    P = matobj.(names{i});
    inputs(i) = Inputs;
    inputs(i).name = names{i};
    inputs(i).tree = i;
    inputs(i).plot = 0;
    inputs(i).savetxt = 0;
    inputs(i).savemat = 0;
    inputs(i).disp = 0;
  end

  %% Estimate the stem diameter close to bottom
  % Define height
  Hb = min(P(:,3));
  Ht = max(P(:,3));
  TreeHeight = double(Ht-Hb);
  Hei = P(:,3)-Hb;

  % Select a section (0.02-0.1*tree_height) from the bottom of the tree
  hSecTop = min(4,0.1*TreeHeight);
  hSecBot = 0.02*TreeHeight;
  hSec = hSecTop-hSecBot;
  Sec = Hei > hSecBot & Hei < hSecTop;
  StemBot = P(Sec,1:3);

  % Estimate stem axis (point and direction)
  AxisPoint = mean(StemBot);
  V = StemBot-AxisPoint;
  V = normalize(V);
  AxisDir = optimal_parallel_vector(V);

  % Estimate stem diameter
  d = distances_to_line(StemBot,AxisDir,AxisPoint);
  Rstem = double(median(d));

  % Point resolution (distance between points)
  Res = sqrt((2*pi*Rstem*hSec)/size(StemBot,1));

 %% Define the PatchDiam parameters
% User-defined parameters
userPatchDiam1 = [0.01 0.03]; 
userPatchDiam2Min = [0.02 0.03 0.04];
userPatchDiam2Max = [0.03 0.05];

% PatchDiam1
if nPD1 == 1
  inputs(i).PatchDiam1 = userPatchDiam1;
else
  n = nPD1;
  inputs(i).PatchDiam1 = linspace(min(userPatchDiam1), max(userPatchDiam1), n);
end

% PatchDiam2Min
if nPD2Min == 1
  inputs(i).PatchDiam2Min = userPatchDiam2Min;
else
  n = nPD2Min;
  inputs(i).PatchDiam2Min = linspace(min(userPatchDiam2Min), max(userPatchDiam2Min), n);
end

% PatchDiam2Max
if nPD2Max == 1
  inputs(i).PatchDiam2Max = userPatchDiam2Max;
else
  n = nPD2Max;
  inputs(i).PatchDiam2Max = linspace(min(userPatchDiam2Max), max(userPatchDiam2Max), n);
end

  % Define the BallRad parameters:
  inputs(i).BallRad1 = max([inputs(i).PatchDiam1+1.5*Res;
    min(1.25*inputs(i).PatchDiam1,inputs(i).PatchDiam1+0.025)]);
  inputs(i).BallRad2 = max([inputs(i).PatchDiam2Max+1.25*Res;
    min(1.2*inputs(i).PatchDiam2Max,inputs(i).PatchDiam2Max+0.025)]);

  %plot_point_cloud(P,1,1)
end
