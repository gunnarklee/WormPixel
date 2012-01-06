function [Totalcelldata, cellnums] = cloadfinal2(ims, varargin)
%CLOAD Load a single image file and return cell data.
%   DATA = CLOAD loads a single image file, identifies objects via the CONTOUR
%   function, and returns a data structure which includes the original file
%   name and loaded image file, the parameters used in the analysis, a
%   count of the number of cells found, and a substructure containing the
%   cells borders and internal coordinate system data.
%
%   If called without any parameters the user is prompted graphically to
%   choose a file. If called with parameters, if the first parameter is a
%   valid file, that file is analyzed. Otherwise the user is prompted to
%   choose a file.
%
%   DATA = CLOAD(...,'PARAMETERNAME1', PARAMETERVALUE1, ...) sets
%   PARAMETERNAME1 to PARAMETERVALUE1. If a parameter is not specified its
%   default value is used. Valid parameters and their default values are:
%
%   contourlevel (100): the pixel value for the contour of the adjusted
%       image, from 0 (dark) to 255 (light).
%   mincellperim (30): the minimum number of points in the perimeter for a
%       cell to be counted.
%   maxcellperim (500): the maximum number of points in the perimeter for a
%       cell to be counted.
%   minarea (50): the minimum area for a cell to be counted.
%   maxarea (Inf): the maximum area for a cell to be counted.
%   npoints (200): the number of points in the middle line of the cell.
%   optmid (true): find the optimal midline.
%   filterwidth (true): filter cells by width.
%   wdfilter (true): filter cells by the standard deviation of the second
%       derivative with respect to length of the cell width (trust me, it
%       works really well).

% First, set default arguments; then, load new values if they exist

param.mincellperim = 10;
param.maxcellperim = Inf;
param.minarea = 10;
param.maxarea = Inf;
param.npoints = 30;
param.optmid = false;
param.filterwidth = false;
param.contourlevel =1; %10
param.channel = 1;
param.wdfilter = false;
filename = [];
% 
%     function paramset(startpoint)
%         if length(varargin) < startpoint + 1
%             %error(['Parameter ' varargin{startpoint} ' specified without value!'])
%             disp(['WARNING: Parameter "' varargin{startpoint} '" was specified without a value! It is being ignored.'])
%         elseif ischar(varargin{startpoint})==0
%             %error(['Parameter ' num2str(startpoint) ' is not a string: invalid syntax'])
%             paramset(startpoint + 1)
%         else
%             if isfield(param, lower(varargin{startpoint}))
%                 if isa(varargin{startpoint+1},class(param.(lower(varargin{startpoint}))))
%                     param.(lower(varargin{startpoint})) = varargin{startpoint+1};
%                     startpoint = startpoint +1;
%                 else
%                     disp(['WARNING: data specified is incorrect type for parameter "' varargin{startpoint} '" and is being ignored.'])
%                 end
%             else
%                 disp(['WARNING: "' varargin{startpoint} '" is not a valid parameter. It is being ignored.'])
%             end
%             if length(varargin) >= startpoint+1
%                 paramset(startpoint + 1)
%             end
%         end
%     end
% 
%     function linesearch(endshift)
%         shiftedpole = indmod(polepos + endshift,length(celldata.cells(cellcount).border));
%         if length(min(shiftedpole):max(shiftedpole)) < length(celldata.cells(cellcount).border)/4
%             return
%         end;
%         top = evenint(celldata.cells(cellcount).border(:, min(shiftedpole):max(shiftedpole)) ,param.npoints);
%         bot = fliplr(evenint([celldata.cells(cellcount).border(:, max(shiftedpole):end) celldata.cells(cellcount).border(:, 1:min(shiftedpole))],param.npoints));
%         mid = (top+bot)./2;
%         if linelen(mid) > linelen(celldata.cells(cellcount).mid)
%             celldata.cells(cellcount).top = top;
%             celldata.cells(cellcount).bot = bot;
%             celldata.cells(cellcount).mid = mid;
%             polepos = shiftedpole;
%             celldata.cells(cellcount).poles = [celldata.cells(cellcount).border(:,polepos(1)) celldata.cells(cellcount).border(:,polepos(2))];
%             linesearch(endshift)
%         end
%     end
% 
% % Check arguments
% if nargin == 0 % If no arguments are given, prompt for a filename
%     [filename, pathname] = uigetfile('*.tif','Pick an image file...');
%     cd(pathname)
%     if isequal(filename,0)
%         error('No file specified');
%     else
%         filename = [pathname filename];
%     end
% else
%     % Check if first argument is a string
%     if isnumeric(varargin{1}) && ~isscalar(varargin{1}) && ~isvector(varargin{1})
%         image = varargin{1};
%         if nargin >= 2
%             paramset(2);
%         end
%     elseif ischar(varargin{1}) == 0
%         error('First parameter is not an image or filename: invalid syntax')
%     else
%         switch exist(varargin{1},'file')
%             case 2 % If the first argument is a valid filename, use it
%                 filename = varargin{1};
%                 if nargin >= 2
%                     paramset(2);
%                 end
%             case 7 % If the first argument is a directory name, open the file selector in that directory
%                 [filename, pathname] = uigetfile(fullfile(varargin{1}, '*.tif'),'Pick an image file...');
%                 if isequal(filename,0)
%                     error('No file specified');
%                 else
%                     filename = [pathname filename];
%                     if nargin >= 2
%                         paramset(2);
%                     end
%                 end
%             otherwise % Otherwise, open up the file selector UI
%                 if isfield(param, varargin{1})
%                     [filename, pathname] = uigetfile('*.tif','Pick an image file...');
%                     if isequal(filename,0)
%                         error('No file specified');
%                     else
%                         filename = [pathname filename];
%                     end
%                     if nargin >= 1
%                         paramset(1);
%                     end
%                 else
%                     error([varargin{1} ' is neither a valid filename nor a valid parameter name; please check your spelling and current working directory and try again'])
%                 end
%         end
%         image = imread(filename); % Read in the specified file
%     end
% end
%  imsubs = makeSubStack(filename, pathname);
imsubs = ims;
% for i = 1:size(imsubs, 3)
    
% image = imread(filename);
for i = 1:size(imsubs,3)
image = imsubs(:,:,i);
celldata.SourceName = filename;
celldata.SourceImage = image;
if size(image,3) > 1, image = image(:,:,param.channel); end; % If the image is multi-channel, choose the first channel
celldata.parameters = param;
%NOTE: removed "imadjust" in front of "image" in next line

contourdata = contourc(double((image)),[param.contourlevel param.contourlevel]); % Get contour matrix with 1 level from image
%contourdata = contourc(double((image)),1); % Get contour matrix with 1 level from image

% Split contour matrix (from CONTOURC, CONTOUR, etc.) into individual cell data
ind = 1;
celldata.CellCount = 0;
celldata.cells = [];
while (ind+1 <= length(contourdata)) && (contourdata(2,ind)+ind <= length(contourdata))
    % Filter for perimeter size and cells touching the edge of the image
    if (contourdata(2,ind) > param.mincellperim && contourdata(2,ind) < param.maxcellperim) ...
            && ~any(any(contourdata(:,(ind+1):(contourdata(2,ind)+ind)) == 1)) ...
            && ~any(contourdata(1,(ind+1):(contourdata(2,ind)+ind)) == size(celldata.SourceImage,2)) ...
            && ~any(contourdata(2,(ind+1):(contourdata(2,ind)+ind)) == size(celldata.SourceImage,1))
        celldata.CellCount = celldata.CellCount+1;
        celldata.cells(celldata.CellCount).border = contourdata(:,(ind+1):(contourdata(2,ind)+ind));
        % Calculate cell area and perimeter
        celldata.cells(celldata.CellCount).area = polyarea(celldata.cells(celldata.CellCount).border(1,:),celldata.cells(celldata.CellCount).border(2,:));
        celldata.cells(celldata.CellCount).perimeter = sum(sqrt(sum(diff(celldata.cells(celldata.CellCount).border,1,2).^2)));
        % Calculate cell "center"
        celldata.cells(celldata.CellCount).center = mean(celldata.cells(celldata.CellCount).border');
        % Filter by area
        if celldata.cells(celldata.CellCount).area < param.minarea || celldata.cells(celldata.CellCount).area > param.maxarea
            celldata.cells(celldata.CellCount) = [];
            celldata.CellCount = celldata.CellCount - 1;
        end
    end
    ind = ind+1+contourdata(2,ind);
end

% Create internal coordinate system

for cellcount = 1:length(celldata.cells)
    % Identify poles
    distsx = repmat(celldata.cells(cellcount).border(1,:),length(celldata.cells(cellcount).border(1,:))-1,1) - fliplr(reshape(repmat(celldata.cells(cellcount).border(1,:)', 1, length(celldata.cells(cellcount).border(1,:))-1),length(celldata.cells(cellcount).border(1,:))-1,length(celldata.cells(cellcount).border(1,:))));
    distsy = repmat(celldata.cells(cellcount).border(2,:),length(celldata.cells(cellcount).border(2,:))-1,1) - fliplr(reshape(repmat(celldata.cells(cellcount).border(2,:)', 1, length(celldata.cells(cellcount).border(2,:))-1),length(celldata.cells(cellcount).border(2,:))-1,length(celldata.cells(cellcount).border(2,:))));
    dists = (distsx.^2 + distsy.^2)';

    [polepos(1), polepos(2)] = find(dists == max(max(dists)),1);
    polepos(2) = indmod(polepos(1)+polepos(2),length(celldata.cells(cellcount).border));
    celldata.cells(cellcount).poles = [celldata.cells(cellcount).border(:,polepos(1)) celldata.cells(cellcount).border(:,polepos(2))];
    % Split into two halves in order to find midlines
    celldata.cells(cellcount).top = evenint(celldata.cells(cellcount).border(:, min(polepos):max(polepos)) ,param.npoints);
    celldata.cells(cellcount).bot = fliplr(evenint([celldata.cells(cellcount).border(:, max(polepos):end) celldata.cells(cellcount).border(:, 1:min(polepos))],param.npoints));
    celldata.cells(cellcount).mid = (celldata.cells(cellcount).top+celldata.cells(cellcount).bot)./2;
    celldata.cells(cellcount).centermid = mean(celldata.cells(cellcount).mid,2);
    % Find "best" midline
    if param.optmid
        linesearch([1 0]);
        linesearch([-1 0]);
        linesearch([0 1]);
        linesearch([0 -1]);
    end

    clear dists poles top bot mid newpole1

    % Find partial lengths of midline to use as basis for internal axis
    celldata.cells(cellcount).internalx(2:param.npoints) = cumsum(sqrt(sum(diff(celldata.cells(cellcount).mid,1,2).^2)));

    % Get length of midline
    celldata.cells(cellcount).midlength = celldata.cells(cellcount).internalx(param.npoints);

    % Find thicknesses
    celldata.cells(cellcount).thickness = sqrt(sum((celldata.cells(cellcount).top-celldata.cells(cellcount).bot).^2));
    
    % Find second derivatives for filtering
    celldata.cells(cellcount).wdfilter = std(diff(celldata.cells(cellcount).thickness,2));
end

% Filter by thickness
if param.filterwidth && isfield(celldata.cells, 'thickness')
    maxwidth = max(reshape([celldata.cells.thickness],param.npoints,celldata.CellCount));
    celldata.cells((maxwidth > (mean(maxwidth)+std(maxwidth))))=[];%celldata.cells(find(maxwidth > (mean(maxwidth)+std(maxwidth))))=[];
end
% Filter by second derivative of thickness
if param.wdfilter && isfield(celldata.cells, 'thickness')
    celldata.cells([celldata.cells.wdfilter] > 2)=[];
    celldata.cells = rmfield(celldata.cells,'wdfilter');
end
celldata.CellCount = length(celldata.cells);
Totalcelldata(i) = celldata;
cellnums(i) = celldata.CellCount;
end

end 
%---------------
% Subfunctions

%---------------
function z = indmod(x,y)
% Same as MOD, except indmod(x,x) = x, not 0. Useful for circular shifted
% indices.
z = mod(x,y);
z(z==0)=y;
end

%---------------
function points = evenint(curve, n)
% Returns n equally spaced points distributed along curve

del = linelen(curve)/(n-1);
if del == 0
    error('Can not subdivide a curve with length zero')
end

partials(2:length(curve)) = cumsum(sqrt(sum(diff(curve,1,2).^2)));

partmat = repmat(partials,n-2,1)';
cdel = del*(1:n-2);
cdelrep = repmat(cdel,length(curve),1);
indexmat = length(curve)+1-sum(partmat >= cdelrep);

points(:,2:n-1) = curve(:,indexmat-1) + repmat((cdel-partmat(indexmat-1))./(partmat(indexmat)-partmat(indexmat-1)),2,1).*(curve(:,indexmat)-curve(:,indexmat-1));
points(:,[1 n]) = curve(:,[1 end]);
end

%---------------
function leng = linelen(points)
% Returns the length of a line. The input format is 2 by n matrix, where the top row is x and bottom is y

if size(points,2) < 2
    leng = 0;
else
    leng = sum(sqrt(sum(diff(points,1,2).^2)));
end
end
