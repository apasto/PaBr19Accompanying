function TOUT_sectionCall()
%TOUT_sectionCall
% Wrapper for the 'TOUT_section' function.
% Plot a map of the study area, wait for 2 clicks on section extremes,
% then plot a section (slice) between them.
% TOUT_section will then  wait for 2 clicks on the section, to define columns.
%
% 2019, Alberto Pastorutti and Carla Braitenberg

narginchk(0,0)
nargoutchk(0,0)

% ETOPO1 in area, 0.05 deg step
load('../topo/ETOPO1_005d_crop.mat','ETOPO1_005d');

% get UTM zone and extent polygon from grid reference
load('../thermal/Tgrid.mat','Tgrid');

% plot topography map
[MapFigure, MapAxes, MapImage, MapColorbar] = ...
    mapimagesc(ETOPO1_005d.lon,ETOPO1_005d.lat,ETOPO1_005d.z); %#ok<ASGLU>
colormap(MapAxes,parula);
MapAxes.YLim = [44, 56];
MapAxes.XLim = [14, 36];
MapAxes.CLim = [-5,1250]; % static, since area stays the same
MapAxes.Colormap(1,:) = [0.75 0.75 0.75]; % gray water bodies
MapAxes.XGrid = 'on';
MapAxes.XMinorGrid = 'on';
MapAxes.YGrid = 'on';
MapAxes.YMinorGrid = 'on';
MapAxes.Layer = 'top';
MapAxes.XLabel.String = 'Lon. [deg]';
MapAxes.YLabel.String = 'Lat. [deg]';
MapColorbar.Label.String = 'Topography [m]';

% overlay extents polygon, after un-projection to WGS84
% 'explode' sides, reproject one by one
% doing so to make truly un-projected lines
NSegm = 50; % number of segement for each side
[WSide_Lat,WSide_Lon] = ...
    minvtran(...
        Tgrid.UTMstruct,...
        linspace(Tgrid.Extents(1,1),Tgrid.Extents(2,1),NSegm),... % Easting
        linspace(Tgrid.Extents(1,2),Tgrid.Extents(2,2),NSegm)); % Northing
[NSide_Lat,NSide_Lon] = ...
    minvtran(...
        Tgrid.UTMstruct,...
        linspace(Tgrid.Extents(2,1),Tgrid.Extents(3,1),NSegm),... % Easting
        linspace(Tgrid.Extents(2,2),Tgrid.Extents(3,2),NSegm)); % Northing
[ESide_Lat,ESide_Lon] = ...
    minvtran(...
        Tgrid.UTMstruct,...
        linspace(Tgrid.Extents(3,1),Tgrid.Extents(4,1),NSegm),... % Easting
        linspace(Tgrid.Extents(3,2),Tgrid.Extents(4,2),NSegm)); % Northing
[SSide_Lat,SSide_Lon] = ...
    minvtran(...
        Tgrid.UTMstruct,...
        linspace(Tgrid.Extents(4,1),Tgrid.Extents(1,1),NSegm),... % Easting
        linspace(Tgrid.Extents(4,2),Tgrid.Extents(1,2),NSegm)); % Northing
Pl.WSide = plot(WSide_Lon,WSide_Lat,...
    'Color','red','LineWidth',2,'Parent',MapAxes);
Pl.NSide = plot(NSide_Lon,NSide_Lat,...
    'Color','red','LineWidth',2,'Parent',MapAxes);
Pl.ESide = plot(ESide_Lon,ESide_Lat,...
    'Color','red','LineWidth',2,'Parent',MapAxes);
Pl.SSide = plot(SSide_Lon,SSide_Lat,...
    'Color','red','LineWidth',2,'Parent',MapAxes);

% pick points
[Input_lon,Input_lat] = ginput(2);
% project extremes to UTM, unproject segments to WGS84
[Input_x,Input_y] = mfwdtran(Tgrid.UTMstruct,Input_lat,Input_lon);
[Section_Lat,Section_Lon] = minvtran(...
    Tgrid.UTMstruct,...
    linspace(Input_x(1),Input_x(2),NSegm),... % Easting
    linspace(Input_y(1),Input_y(2),NSegm)); % Northing

%plot on map
Pl.Section = plot(Section_Lon,Section_Lat,...
    'Color','red','LineWidth',2,'Parent',MapAxes); %#ok<STRNU>

% call section, points in the form [x1,y1; x2,y2]
Section_points = [Input_x,Input_y];
[Sect.FigSect, Sect.Ax, Sect.Pl, Sect.ColorbarsH, Sect.PatchH] = ...
    TOUT_section(Section_points,'all',1); %#ok<STRNU>

end

function varargout = mapimagesc(x, y, c)
%mapimagesc wrap 'imagesc' to plot a 1:1 aspect, x,y-axis image
%
% Syntax: [figure, axes, image, colorbar](optional) = mapimagesc(x, y, c)
%
%  Plot a scaled-colour image, directly with 1:1 aspect ratio
%  and origin placed in the bottom left corner.
%  Commonly used as a simple and fast way to check map data
%  on rectangular grids.
%  This is equivalent to calling imagesc(x,y,c), opening plot tools
%  then clicking 'axis image' and 'axis xy' in the image property editor.
%
%  Input arguments:
%    x,y : axis vectors
%    c   : array of image data, with x-axis along rows (dimension 2)
%
%  Output arguments: handles to [figure, axes, image, colorbar]
%   either none or exactly 4 output arguments are allowed

% Create figure
figureMAP = figure;

% Create axes
axesMAP = axes('Parent',figureMAP);
hold(axesMAP,'on');

% to get correct origin placement, use 'low level' syntax for imagesc
imageMAP = imagesc('XData',x,'YData',y,'CData',c,'Parent',axesMAP);
box(axesMAP,'on');
axis(axesMAP,'tight');
set(axesMAP,'DataAspectRatio',[1 1 1],'Layer','top');
colorbarMAP = colorbar('peer',axesMAP);

% write optional argouts
if nargout==4
    varargout{1} = figureMAP;
    varargout{2} = axesMAP;
    varargout{3} = imageMAP;
    varargout{4} = colorbarMAP;
end

end
