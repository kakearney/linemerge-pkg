function [latout, lonout, segindices] = linemerge(lat, lon)
%LINEMERGE Merges line segments that share endpoints
%
% [latout, lonout] = linemerge(lat, lon)
% [latout, lonout, segindices] = linemerge(lat, lon)
%
% Merges line segments that have matching endpoints.  This function
% performs the same function claimed by the Mapping Toolbox function
% polymerge; however, polymerge returns incorrect results for many
% situations.  
%
% Input variables:
%
%   lat:        cell array or NaN-delimited vector of line segment
%               latitudes 
%
%   lon:        cell array or NaN-delimited vector of line segment
%               longitudes 
%
% Output variables:
%
%   latout:     cell array or NaN-delimited vector (matches input) of
%               merged line segment latitudes.  These may overlap each
%               other if branches lead to mutliple path possibilites 
%
%   lonout:     cell array or NaN-delimited vector (matches input) of
%               merged line segment longitudes.  These may overlap each
%               other if branches lead to mutliple path possibilites.
%
%   segindices: cell array holding indices of the original line segments
%               that make up each output segment

% Copyright 2006 Kelly Kearney

%---------------------------
% Check input
%---------------------------

inputisvector = false;
if ~iscell(lat)
    [lat, lon] = polysplit(lat, lon);
    inputisvector = true;
end

nseg = length(lat);

%---------------------------
% Find all unique vertices
%---------------------------

alllat = cat(1, lat{:});
alllat = alllat(:);
alllon = cat(1, lon{:});
alllon = alllon(:);
coords = [alllat alllon];
verts = unique(coords, 'rows');

for iseg = 1:nseg
    iverts(iseg,1) = find((lat{iseg}(1) == verts(:,1)) & (lon{iseg}(1) == verts(:,2)));
    iverts(iseg,2) = find((lat{iseg}(2) == verts(:,1)) & (lon{iseg}(2) == verts(:,2)));
end

%---------------------------
% Find connected segments
%---------------------------

% Two-point runs

runs = num2cell(iverts,2);

% Look for longer runs

count = 3;
while count <= size(verts,1)
    
    fprintf('Pass %d of %d\n', count, size(verts,1));
    
    % See if you can add on one point to each existing run
    tic;
    newruns = cell(0);
    for irun = 1:length(runs)
        
        %fprintf('  Looking for additions to run %d of %d\n', irun, length(runs));
        
        ivertend = runs{irun}(end);
        
        inextseg = find(ivertend == iverts(:,1));
        for icnect = 1:length(inextseg)
            newruns = [newruns; {[runs{irun} iverts(inextseg(icnect),2)]}];
        end
    end
    
    t1 = toc;
    fprintf('  Search for continuations: %f s\n', t1); 
    tic;
    
    % Remove short runs that are subsets of new runs
    
    if isempty(newruns)
        fprintf('No new runs, exiting\n');
        break
    else
        
        
        runsText = cellfun(@(a) [',' sprintf('%d,',a)], runs, 'uni', 0);
        newrunsText = cellfun(@(a) [',' sprintf('%d,',a)], newruns, 'uni', 0);
        printtextarray(runsText, 'runs.temp');
        printtextarray(newrunsText, 'newruns.temp');

        %fprintf('  Checking for subsets\n');  
        perl('linemerge.pl');

        issubset = logical(load('issub.temp'));
        delete('issub.temp', 'runs.temp', 'newruns.temp');
    %     issubset = false(size(runs));
    %     for irun = 1:length(runs)
    %         
    %         %fprintf('  Checking new runs as subset of run %d of %d\n', irun, length(runs));
    %         
    %         for inewrun = 1:length(newruns)
    %             [tf, loc] = ismember(runs{irun}, newruns{inewrun});
    %             if all(tf)
    %                 if all(diff(loc) == 1)
    %                     issubset(irun) = true;
    %                 end
    %             end
    %         end
    %     end
        runs(issubset) = [];

        runs = [runs; newruns];
        count = count + 1;
        
        t2 = toc;
        fprintf('  Check for subsets: %f s\n', t2); 
    end
    
    
end

%---------------------------
% Format output
%---------------------------

% Get segment indices

if nargout == 3
    
    segindices = cell(size(runs));
    
    for irun = 1:length(runs)
        npoints = length(runs{irun});
        segindices{irun} = zeros(npoints-1,1);
        for ipt = 1:npoints-1
            [tf, segindices{irun}(ipt)] = ismember(runs{irun}(ipt:ipt+1), iverts, 'rows');
        end
    end
    
end
        
% Longitude and latitude

latout = cell(size(runs));
lonout = cell(size(runs));
for irun = 1:length(runs)
    latout{irun} = verts(runs{irun}, 1);
    lonout{irun} = verts(runs{irun}, 2);
end
   
if inputisvector
    [latout, lonout] = polyjoin(latout, lonout);
end
