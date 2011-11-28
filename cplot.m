function cplot(celldata, str, feat,cellnums,  width)
% Plot cell borders on adjusted original image
%
% CPLOT(DATA,STR,FEAT,CELLNUMS,WIDTH)
% DATA is the only required argument, and refers to the variable output from CLOAD, CELLPARAM, or SICIC
% STR is a string describing the color, linestyle, and markerstyle to plot. Syntax is the same as for the PLOT command
% FEAT is a string describing what to plot. If empty, it defaults to 'border'. 'mid' plots the midline,
% CELLNUMS is a vector containing which cells to plot, useful for highlighting cells with certain
%  characteristics from the output of a FIND command. An empty vector plots all cells
%  'cross' plots a subset of the lateral lines, and 'top' and 'bot' plot opposite halves of the cell
% WIDTH is the width of the line to draw

if ~exist('cellnums') || isempty(cellnums), cellnums = [1:length(celldata.cells)]; end
if ~exist('feat') || isempty(feat), feat = 'border'; end
if ~exist('str') || isempty(str), str = 'b';end
if ~exist('width') || isempty(width), width = 1;end

%figure
%imshow(imadjust(celldata.SourceImage(:,:,1)))
hold on
for count = cellnums
    if strcmp(feat,'cross')
		for c=1:5:celldata.parameters.npoints
		     plot([celldata.cells(count).top(1,c) celldata.cells(count).bot(1,c)], [celldata.cells(count).top(2,c) celldata.cells(count).bot(2,c)],str,'LineWidth', width)
		end
	else
		plot(celldata.cells(count).(feat)(1,:),celldata.cells(count).(feat)(2,:),str,'LineWidth',width)
	end
end