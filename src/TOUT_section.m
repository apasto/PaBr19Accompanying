function varargout = TOUT_section(varargin)
%TOUT_section plot section (vertical slice, 1 segment) of thermal model (and input data)
%
% Syntax: [FigSect, Ax, Pl, ColorbarsH, PatchH] = HF_section(Points,ContourType,DoColumns)
%
%  Input arguments: all optional
%    Points : Section profile, a line defined by a 2x2 matrix in the form [x1,y1; x2,y2].
%             If no arguments are provided, or if this is empty "[]",
%             the default section is loaded - same section plotted in the
%             manuscript.
%    ContourType : either 'none','notext','all'
%                  no contours, non labelled contours, labelled contours
%    DoColumns : if this argument is provided, plot a second figure with columns
%                (curves of T, k, A, Q vertically trough the model)
%                - if set to 1, it waits for the user to click on the section twice,
%                  to define the location of two columns
%                - if set to 2, 4 columns are plotted at the same distance
%                  along section used in the manuscript figure
%                - it can be set to a vector of column positions, along section
%
%  Output arguments: graphic handles to [figure, axes, plots, colorbars, section fill]
%   either none or exactly 5 output arguments are allowed
%
% 2019, Alberto Pastorutti and Carla Braitenberg

%% manage output arguments
% before performing anything
if ~or(nargout==0,nargout==5)
    error('Number of output arguments must be either 0 (none) or 5 (graphic handles)')
end

%% manage input arguments
narginchk(0,3)
if nargin==0
    LoadDefaultSection = true;
    % no arguments, 'ContourType' default is 'all'
    ContourType = 'notext';
else
    if isempty(varargin{1})
        LoadDefaultSection = true;
    else
        LoadDefaultSection = false;
        % check if input is a 2x2 matrix
        assert(and(size(varargin{1},1)==2,size(varargin{1},1)==2),...
            'First input argument must be a 2x2 matrix, [x1,y1; x2,y2]')
        SectV.X = varargin{1}(:,1);
        SectV.Y = varargin{1}(:,2);
    end
    if nargin==1
        % no optional argument 'ContourType', default is 'all'
        ContourType = 'all';
    else % implies nargin==2
        ContourType = varargin{2};
        % check if is legit case
        assert(any(strcmpi(ContourType,{'none','notext','all'})),...
            'Second input argument must be one of ''none'',''notext'',''all''')
    end
end

if nargin==3
    DoColumns = varargin{3};
else
    DoColumns = false;
end

if LoadDefaultSection
    % use default vertices
    SectV.X = [-60000, 820000];
    SectV.Y = [5180000, 5880000];
end

%% load data: output volumes
load('../../mohoHFpaper/data/2018-09_paper/heatflow/Iter.mat','Iter');
% get UTM zone
load('../thermal/Tgrid.mat','Tgrid');
% ETOPO1 in area, 0.01 deg step
load('../topo/ETOPO1_005d_crop.mat','ETOPO1_005d');

%% define section interpolation parameters
% in metres
SectI_horStep = 10e3;
SectI_verStep = 250;

%% build grids

% build section grid
% distance along profile
Sect_d = ...
    0:...
    SectI_horStep:...
    sqrt((SectV.X(2)-SectV.X(1))^2+(SectV.Y(2)-SectV.Y(1))^2);
SectX = interp1([Sect_d(1), Sect_d(end)],SectV.X,Sect_d,'linear');
SectY = interp1([Sect_d(1), Sect_d(end)],SectV.Y,Sect_d,'linear');
SectZ = ...
    Iter.Layers.DefGrid.zshift:...
    SectI_verStep:...
    Iter.Layers.DefGrid.z(end)+Iter.Layers.DefGrid.zshift;

% for interp3, build meshgrid of slice-coordinates
[SectXmeshX,SectXmeshZ] = ...
    meshgrid(SectX,SectZ);
[SectYmeshY,~] = ...
    meshgrid(SectY,SectZ);

% build meshgrids from 'Iter' struct - this includes edge padding
[PX,PY] = ...
    meshgrid(...
        Iter.Layers.DefGrid.y,...
        Iter.Layers.DefGrid.x);

[VX,VY,VZ] = ...
    meshgrid(...
    Iter.Layers.DefGrid.y,...
    Iter.Layers.DefGrid.x,...
    Iter.Layers.DefGrid.z+Iter.Layers.DefGrid.zshift);

%% un-project profile points to WGS84 and extract topography profile
[SectLat,SectLon] = minvtran(Tgrid.UTMstruct,SectX,SectY);

InT = interp2(...
    ETOPO1_005d.lat',ETOPO1_005d.lon,...
    ETOPO1_005d.z',...
    SectLat,SectLon,'makima');

%% extract section-horizons from depth maps
% using interp2
% sequence of query as Y,X keeps coords convention

Interp2Method = 'linear';

InP.AIR = interp2(...
    PX,PY,...
    Iter.Layers.AIR.RoundedDepthMap+Iter.Layers.DefGrid.zshift,...
    SectY,SectX,Interp2Method)*1e-3;
InP.SED = interp2(...
    PX,PY,...
    Iter.Layers.SEDS.RoundedDepthMap+Iter.Layers.DefGrid.zshift,...
    SectY,SectX,Interp2Method)*1e-3;
InP.UCR = interp2(...
    PX,PY,...
    Iter.Layers.UCRUST.RoundedDepthMap+Iter.Layers.DefGrid.zshift,...
    SectY,SectX,Interp2Method)*1e-3;
InP.LCR = interp2(...
    PX,PY,...
    Iter.Layers.LCRUST.RoundedDepthMap+Iter.Layers.DefGrid.zshift,...
    SectY,SectX,Interp2Method)*1e-3;
InP.LAB = interp2(...
    PX,PY,...
    Iter.Layers.LID.RoundedDepthMap+Iter.Layers.DefGrid.zshift,...
    SectY,SectX,Interp2Method)*1e-3;

%% extract heat flow profiles from maps

Interp2QMethod = 'linear';

InQ.Q0 = interp2(...
    PX,PY,...
    Iter.IterData.Q0{end},...
    SectY,SectX,Interp2QMethod);
InQ.Qm = interp2(...
    PX,PY,...
    Iter.IterData.Qm{end},...
    SectY,SectX,Interp2QMethod);

%% extract section-slices from volumes
% using interp3, and permuting volumes from [z,x,y] to [x,y,z]

VolPermuteOrder = [2 3 1]; % [z,x,y] to [x,y,z]

% graphical workaround for k and A interpolation:
%  fill areas outside boundaries with values at boundary
%  this prevent the linear interpolation from interpolating
%  the unused values in 'air' and 'asthenosphere'
Iter.Vol.Ai = Iter.Vol.A;
Iter.Vol.ki = Iter.Vol.k;
Iter.Vol.Ai(Iter.Vol.A==99) = 0.02e-6; % under LAB
Iter.Vol.Ai(Iter.Vol.A==88) = 0; % over TOPO
Iter.Vol.ki(Iter.Vol.k==0) = 2.45; % approx value of k(T) at LAB
% 2.45 will show slightly at TOPO, when zooming a lot
% these volumes ('i' suffix) are NOT used later
% when (and if) we are extracting columns

InV.T = interp3(VX,VY,VZ,...
    permute(Iter.Vol.T,VolPermuteOrder),...
    SectYmeshY',SectXmeshX',SectXmeshZ','linear')' - 273.15;
% InV.L = interp3(VX,VY,VZ,...
%     permute(Iter.Vol.L,VolPermuteOrder),...
%     SectYmeshY',SectXmeshX',SectXmeshZ','linear')';
InV.k = interp3(VX,VY,VZ,...
    permute(Iter.Vol.ki,VolPermuteOrder),...
    SectYmeshY',SectXmeshX',SectXmeshZ','linear')';
InV.A = interp3(VX,VY,VZ,...
    permute(Iter.Vol.Ai,VolPermuteOrder),...
    SectYmeshY',SectXmeshX',SectXmeshZ','linear')' * 1e6;

%% declare embedded functions for repeated plot calls

function BoundPlotsStruct = PlotBoundariesOnSection(SectD,InP,ParentAxes,BoundaryColor,BoundaryLineWidth)
InPList = fieldnames(InP);
for pn=2:(length(InPList)-1) % skip topography and bottom
    BoundPlotsStruct.(InPList{pn}) = ...
        stairs(...
            SectD,InP.(InPList{pn}),'Parent',ParentAxes,...
            'Color',BoundaryColor,'LineWidth',BoundaryLineWidth);
end
end

function PatchPlotsStruct = PlotOutOfBoundariesOnSection(SectD,InP,ParentAxes,BoundaryColor,BoundaryLineWidth,FillColor)
% each point is sampled twice (thus the reshaped inner concatenation)
% but y is lagging by one sample
% this is the same strategy used in the 'stairs' function
InPList = fieldnames(InP);
for pn=[1, length(InPList)] % only topography and bottom
    if pn==1
        YLimIndex = 1; % topography: get first YLim, bottom: get second YLim
        OverEdge = -5; % go further than Ylim
    else
        YLimIndex = 2;
        OverEdge = 5;
    end
    PatchPlotsStruct.(InPList{pn}) = ...
        patch(...
            [...
                reshape([SectD; SectD],1,length(SectD)*2),...
                SectD(end), SectD(end), SectD(1)],...
            [...
                InP.(InPList{pn})(1),reshape([InP.(InPList{pn}); InP.(InPList{pn})],1,length(SectD)*2),...
                ParentAxes.YLim(YLimIndex)+OverEdge, ParentAxes.YLim(YLimIndex)+OverEdge],...
            FillColor,...
            'EdgeColor',BoundaryColor,...
            'LineWidth',BoundaryLineWidth,...
            'Parent',ParentAxes);
end
end

%% figure

BoundaryColor = 'white';
BoundaryLineWidth = 0.5;

Sect_d_km = Sect_d*1e-3;
Sect_d_km_steps = (Sect_d-SectI_horStep/2)*1e-3; % shift to make steps fit
SectZ_km = SectZ*1e-3;

FigSect = figure();
FigSect.Colormap = parula;
FigSectsWidth = 5; % width of sections (in subplot cols)
FigPlotsWidth = 1; % width of col on the right (leave space for colorbars)
FigCols = FigSectsWidth + FigPlotsWidth;
FigRows = 9; % topo + Q0 + Qm + 2*(T,k,A)

FigSect.Units = 'normalized';
FigSect.Position = [0.05,0.05,0.4,0.8];

% create subplots axes
Ax.S_topo = subplot(FigRows,FigCols,1:FigSectsWidth);
Ax.S_Q0 = subplot(FigRows,FigCols,FigCols+(1:FigSectsWidth));
Ax.S_Qm = subplot(FigRows,FigCols,2*FigCols+(1:FigSectsWidth));
Ax.S_T = subplot(FigRows,FigCols,[3*FigCols+(1:FigSectsWidth), 4*FigCols+(1:FigSectsWidth)]);
Ax.S_k = subplot(FigRows,FigCols,[5*FigCols+(1:FigSectsWidth), 6*FigCols+(1:FigSectsWidth)]);
Ax.S_A = subplot(FigRows,FigCols,[7*FigCols+(1:FigSectsWidth), 8*FigCols+(1:FigSectsWidth)]);

% set common properties
Ax_names = fieldnames(Ax);
for n=1:length(Ax_names)
    Ax.(Ax_names{n}).Box = 'on';
    Ax.(Ax_names{n}).XGrid = 'on';
    Ax.(Ax_names{n}).YGrid = 'on';
    hold(Ax.(Ax_names{n}),'on');
    Ax.(Ax_names{n}).Layer = 'top';
end

% specific properties
Ax.S_topo.XAxisLocation = 'top';
Ax.S_Q0.XTickLabel = [];
Ax.S_Qm.XTickLabel = [];
Ax.S_T.XTickLabel = [];
Ax.S_k.XTickLabel = [];

% y direction
Ax.S_T.YDir = 'reverse';
Ax.S_k.YDir = Ax.S_T.YDir;
Ax.S_A.YDir = Ax.S_T.YDir;

% link axes and set limits
linkaxes([Ax.S_topo,Ax.S_Q0,Ax.S_Qm,Ax.S_T,Ax.S_k,Ax.S_A],'x')
linkaxes([Ax.S_T,Ax.S_k,Ax.S_A],'y')
Ax.S_topo.XLim = ...
    [min(Sect_d_km), max(Sect_d_km)] + ...
    [SectI_horStep*1e-3, -SectI_horStep*1e-3]; % trim away one step at the edges
Ax.S_T.YLim = [min(InP.AIR)-10, max(InP.LAB)+10];

% axes labels
Ax.S_topo.XAxis.Label.String = 'Distance along section [km]';
Ax.S_A.XAxis.Label.String = Ax.S_topo.XAxis.Label.String;
Ax.S_topo.YAxis.Label.String = 'Topo. [m]';
Ax.S_Q0.YAxis.Label.String = 'Q_0 [mW/m^3]';
Ax.S_Qm.YAxis.Label.String = 'Q_M [mW/m^3]';
Ax.S_T.YAxis.Label.String = 'Depth [km]';
Ax.S_k.YAxis.Label.String = Ax.S_T.YAxis.Label.String;
Ax.S_A.YAxis.Label.String = Ax.S_T.YAxis.Label.String;

% plots
% topography
Pl.S_topo = plot(Sect_d_km,InT,...
    'Parent',Ax.S_topo,'Color','blue','LineWidth',1);
Ax.S_topo.YLim(2) = Ax.S_topo.YLim(2)+250;

% heat flow, Q0
Pl.S_Q.Q0 = plot(Sect_d_km,InQ.Q0,...
    'Parent',Ax.S_Q0,'Color','blue','LineWidth',1,...
    'DisplayName','Q0');
Ax.S_Q0.YLim(2) = Ax.S_Q0.YLim(2)+5;
Ax.S_Q0.YLim(1) = Ax.S_Q0.YLim(1)-5;
Ax.S_Q0.YMinorGrid = 'on';
Ax.S_Q0.YTick = ...
    RoundToStep(5,Ax.S_Q0.YLim(1),'ceil'):...
    25:...
    RoundToStep(5,Ax.S_Q0.YLim(2),'floor');

% heat flow, Qm
Pl.S_Q.Qm = plot(Sect_d_km,InQ.Qm,...
    'Parent',Ax.S_Qm,'Color','blue','LineWidth',1,...
    'DisplayName','QM');
Ax.S_Qm.YLim(2) = Ax.S_Qm.YLim(2)+5;
Ax.S_Qm.YLim(1) = Ax.S_Qm.YLim(1)-5;
Ax.S_Qm.YMinorGrid = 'on';

% T section
Pl.S_T.img = imagesc(Sect_d_km,SectZ_km,InV.T,'Parent',Ax.S_T);
Pl.S_T.bound = PlotBoundariesOnSection(Sect_d_km_steps,InP,Ax.S_T,BoundaryColor,BoundaryLineWidth);

% k section
Pl.S_k.img = imagesc(Sect_d_km,SectZ_km,InV.k,'Parent',Ax.S_k);
Pl.S_k.bound = PlotBoundariesOnSection(Sect_d_km_steps,InP,Ax.S_k,BoundaryColor,BoundaryLineWidth);

% A section
Pl.S_A.img = imagesc(Sect_d_km,SectZ_km,InV.A,'Parent',Ax.S_A);
Pl.S_A.bound = PlotBoundariesOnSection(Sect_d_km_steps,InP,Ax.S_A,BoundaryColor,BoundaryLineWidth);

% correct and round extremes of colormaps
Ax.S_k.CLim(1) = RoundToStep(0.1,min(InV.k(InV.k~=0)),'ceil');
Ax.S_k.CLim(2) = RoundToStep(0.1,Ax.S_k.CLim(2),'floor');
Ax.S_A.CLim(2) = RoundToStep(0.1,max(InV.A(InV.A<88)),'floor'); % 88 and 99 are boundary tags

% colormaps 
ColorbarsH.S_T = colorbar(Ax.S_T,'Location','east');
ColorbarsH.S_k = colorbar(Ax.S_k,'Location','east');
ColorbarsH.S_A = colorbar(Ax.S_A,'Location','east');
ColorBarsList = fieldnames(ColorbarsH);
for cn=1:length(ColorBarsList)
    % move to empty space on the right
    ColorbarsH.(ColorBarsList{cn}).AxisLocation = 'out';
    ColorbarsH.(ColorBarsList{cn}).Position(1) = ...
        ColorbarsH.(ColorBarsList{cn}).Position(1)+0.1;
    % larger ticks, both directions
    ColorbarsH.(ColorBarsList{cn}).TickDirection = 'both';
    ColorbarsH.(ColorBarsList{cn}).TickLength = 0.025;
end
ColorbarsH.S_T.Label.String = 'T [\circC]';
ColorbarsH.S_k.Label.String = 'k [W/(m\cdotK)]';
ColorbarsH.S_A.Label.String = 'A [\muW/m^3]';

% k and A, closer ticks
ColorbarsH.S_k.Ticks = Ax.S_k.CLim(1):0.25:Ax.S_k.CLim(2);
ColorbarsH.S_A.Ticks = ...
    RoundToStep(0.5,Ax.S_A.CLim(1),'ceil'):0.5:Ax.S_A.CLim(2);

% T contours
if any(strcmpi(ContourType,{'notext','all'}))
        [Pl.S_t.cont.C,Pl.S_t.cont.h] = ... % contour provides contours and handle
            contour(Sect_d_km,SectZ_km,InV.T,200:200:1000,...
            'Parent',Ax.S_T,...
            'LineColor','black');
        if strcmpi(ContourType,'all')
            clabel(Pl.S_t.cont.C,Pl.S_t.cont.h,'FontSize',9)
        end
end

% white patch in areas out of boundaries
PatchH.S_T = PlotOutOfBoundariesOnSection(...
    Sect_d_km_steps,InP,Ax.S_T,'black',BoundaryLineWidth,'white');
PatchH.S_k = PlotOutOfBoundariesOnSection(...
    Sect_d_km_steps,InP,Ax.S_k,'black',BoundaryLineWidth,'white');
PatchH.S_A = PlotOutOfBoundariesOnSection(...
    Sect_d_km_steps,InP,Ax.S_A,'black',BoundaryLineWidth,'white');

%% plot a figure with T, k, A on selected columns
if DoColumns==1 % select 2 columns on screen
    [Col_d,~] = ginput(2);
elseif DoColumns==2
    Col_d = [110, 415, 735, 1000]; % default sections
elseif DoColumns~=0
    Col_d = DoColumns;
end

if DoColumns~=0
    % plot columns on plots
    for n=1:length(Ax_names)
        Ax.(Ax_names{n}).ColorOrderIndex = 1; % reset plot color order
                                              % to get consistent colors
        for ColN=1:length(Col_d)
            Pl.Cols.(Ax_names{n}){ColN} = ...
                plot(...
                    [Col_d(ColN), Col_d(ColN)],...
                    [Ax.(Ax_names{n}).YLim(1), Ax.(Ax_names{n}).YLim(2)],...
                    'LineWidth',2,...
                    'Parent',Ax.(Ax_names{n}),...
                    'DisplayName',['Col. ',num2str(ColN,'%.0f')]);
        end
    end
    % extract columns from volumes
    % adopted interpolation methods change
    % therefore we do not extract them from slices
    Col_x = zeros(1,length(Col_d));
    Col_y = Col_x;
    ExtractedCols.T = zeros(length(SectZ),ColN);
    ExtractedCols.k = ExtractedCols.T;
    ExtractedCols.linear_k = ExtractedCols.T; % used to calculate Q
    ExtractedCols.A = ExtractedCols.T;
    ExtractedCols.Q = NaN(size(ExtractedCols.T));
    Iter.Vol.k(Iter.Vol.k==0) = NaN;
    Iter.Vol.A(Iter.Vol.A>=88) = NaN;
    minT = zeros(1,length(Col_d));
    maxT = minT;
    minT_i = minT;
    maxT_i = maxT;
    for ColN=1:length(Col_d)
        % get x,y of 'distance along section' of columns
        % round to nearest model node
        Col_x(ColN) = RoundToStep(...
            Iter.Layers.DefGrid.xstep,...
            interp1(Sect_d,SectX,Col_d(ColN)*1e3,'linear'),...
            'round');
        Col_y(ColN) = RoundToStep(...
            Iter.Layers.DefGrid.ystep,...
            interp1(Sect_d,SectY,Col_d(ColN)*1e3,'linear'),...
            'round');
        % if model origin was not a multiple of step, shift accordingly
        Col_x(ColN) = Col_x(ColN) + mod(Iter.Layers.DefGrid.xmin,40e3);
        Col_y(ColN) = Col_y(ColN) + mod(Iter.Layers.DefGrid.ymin,40e3);
        % extract columns
        ExtractedCols.T(:,ColN) = interp3(VX,VY,VZ,...
            permute(Iter.Vol.T,VolPermuteOrder),...
            Col_y(ColN),Col_x(ColN),SectZ','linear') - 273.15;
        ExtractedCols.k(:,ColN) = interp3(VX,VY,VZ,...
            permute(Iter.Vol.k,VolPermuteOrder),...
            Col_y(ColN),Col_x(ColN),SectZ','nearest');
        ExtractedCols.linear_k(:,ColN) = interp3(VX,VY,VZ,...
            permute(Iter.Vol.k,VolPermuteOrder),...
            Col_y(ColN),Col_x(ColN),SectZ','linear');
        ExtractedCols.A(:,ColN) = interp3(VX,VY,VZ,...
            permute(Iter.Vol.A,VolPermuteOrder),...
            Col_y(ColN),Col_x(ColN),SectZ','nearest') *1e6;
        % set to NaN the values above surface and below bottom
        minT(ColN) = min(ExtractedCols.T(:,ColN));
        maxT(ColN) = max(ExtractedCols.T(:,ColN));
        minT_i(ColN) = find(ExtractedCols.T(:,ColN)==minT(ColN),1,'last');
        maxT_i(ColN) = find(ExtractedCols.T(:,ColN)==maxT(ColN),1,'first');
        if minT_i(ColN)~=1
            ExtractedCols.T(1:minT_i(ColN)-1,ColN) = NaN;
            ExtractedCols.k(1:minT_i(ColN)-1,ColN) = NaN;
            ExtractedCols.linear_k(1:minT_i(ColN)-1,ColN) = NaN;
            ExtractedCols.A(1:minT_i(ColN)-1,ColN) = NaN;
        else
            minT_i(ColN) = 2; % for subsequent calls
        end
        % remove last 2 elements i any case, noisy
        ExtractedCols.T(maxT_i(ColN)-1:end,ColN) = NaN;
        ExtractedCols.k(maxT_i(ColN)-1:end,ColN) = NaN;
        ExtractedCols.linear_k(maxT_i(ColN)-1:end,ColN) = NaN;
        ExtractedCols.A(maxT_i(ColN)-1:end,ColN) = NaN;
        maxT_i(ColN) = maxT_i(ColN)-1; % for subsequent calls
    end
    % calculate Q_z heat flow in columns
    for ColN=1:length(Col_d)
        ExtractedCols.Q(minT_i(ColN)+2:maxT_i(ColN),ColN) = ...
            (...
                (diff(ExtractedCols.T(minT_i(ColN)+1:maxT_i(ColN),ColN))/SectI_verStep) ...
                .* (ExtractedCols.linear_k(minT_i(ColN)+2:maxT_i(ColN),ColN))...
            )*1e3;
    end
    % create figure and plots
    ColsFigure.Fig = figure();
    ColsFigure.Ax.T = subplot(1,4,1,'Parent',ColsFigure.Fig);
    ColsFigure.Ax.k = subplot(1,4,2,'Parent',ColsFigure.Fig);
    ColsFigure.Ax.A = subplot(1,4,3,'Parent',ColsFigure.Fig);
    ColsFigure.Ax.Q = subplot(1,4,4,'Parent',ColsFigure.Fig);
    ColsFigure_Ax_names = fieldnames(ColsFigure.Ax);
    for n=1:length(ColsFigure_Ax_names)
        ColsFigure.Ax.(ColsFigure_Ax_names{n}).Box = 'on';
        ColsFigure.Ax.(ColsFigure_Ax_names{n}).XGrid = 'on';
        ColsFigure.Ax.(ColsFigure_Ax_names{n}).YGrid = 'on';
        ColsFigure.Ax.(ColsFigure_Ax_names{n}).YDir = 'reverse';
        ColsFigure.Ax.(ColsFigure_Ax_names{n}).XAxisLocation = 'top';
        hold(ColsFigure.Ax.(ColsFigure_Ax_names{n}),'on');
    end
    ColsFigure.Ax.T.YAxis.Label.String = Ax.S_T.YAxis.Label.String;
    ColsFigure.Ax.T.XAxis.Label.String = ColorbarsH.S_T.Label.String;
    ColsFigure.Ax.k.XAxis.Label.String = ColorbarsH.S_k.Label.String;
    ColsFigure.Ax.A.XAxis.Label.String = ColorbarsH.S_A.Label.String;
    ColsFigure.Ax.Q.XAxis.Label.String = 'Q_{z} [mW/m^3]';
    linkaxes([ColsFigure.Ax.T,ColsFigure.Ax.k,ColsFigure.Ax.A,ColsFigure.Ax.Q],'y');
    ColsFigure.Ax.T.YLim = [SectZ_km(min(minT_i)), SectZ_km(max(maxT_i))];
    % plots
    for ColN=1:length(Col_d)
        % T
        ColsFigure.Pl.T{ColN} = ...
            plot(ExtractedCols.T(:,ColN),SectZ_km,...
            'Parent',ColsFigure.Ax.T,...
            'LineWidth',1,...
            'DisplayName',['Col. ',num2str(ColN,'%.0f')]);
        % k
        ColsFigure.Pl.T{ColN} = ...
            plot(ExtractedCols.k(:,ColN),SectZ_km,...
            'Parent',ColsFigure.Ax.k,...
            'LineWidth',1,...
            'DisplayName',['Col. ',num2str(ColN,'%.0f')]);
        % A
        ColsFigure.Pl.T{ColN} = ...
            plot(ExtractedCols.A(:,ColN),SectZ_km,...
            'Parent',ColsFigure.Ax.A,...
            'LineWidth',1,...
            'DisplayName',['Col. ',num2str(ColN,'%.0f')]);
        % Qz
        ColsFigure.Pl.Q{ColN} = ...
            plot(ExtractedCols.Q(:,ColN),SectZ_km,...
            'Parent',ColsFigure.Ax.Q,...
            'LineWidth',1,...
            'DisplayName',['Col. ',num2str(ColN,'%.0f')]);
    end
end

%% manage output arguments
if nargout==5 % output handles to graphic objects
    varargout{1} = FigSect;
    varargout{2} = Ax;
    varargout{3} = Pl;
    varargout{4} = ColorbarsH;
    varargout{5} = PatchH;
end

end

function out = RoundToStep(step,in,direction)
%RoundToStep round(/floor/ceil/fix) input vector to nearest arbitrary 'step' unit
% digit-counting code comes from this answer by user Jaymin on matlabcentral:
% https://mathworks.com/matlabcentral/answers/10795-counting-the-number-of-digits#answer_68482
% 
% Syntax: out = RoundToStep(step,in,[direction])
%
% Inputs:
%    step        : scalar, unit value
%    in          : vector of values to be rounded
%    [direction] : char vector, optional, direction of rounding
%                  'round', 'floor', 'ceil', 'fix'
%                  'round' is default behaviour
%                  case insensitive
%
% Outputs:
%    out         : rounded 'in' vector
%

% check and manage input
narginchk(2,3)
if nargin==2
    direction='round';
end
assert(isscalar(step),'Argument ''step'' must be a scalar.')
assert(isvector(in),'Argument ''in'' must be a vector.')
% check if 'direction' is an allowed round-like function
allowed_directions = {'round','floor','ceil','fix'};
assert(any(strcmpi(allowed_directions,direction)),...
    ['''',direction,''' is not an allowed rounding function.'])

% count the number of digits
stepDIG = numel(num2str(step))-numel(strfind(num2str(step),'-'))-numel(strfind(num2str(step),'.'));

% number of zsteps in next power of 10
stepdiv = (10^(stepDIG))/step; 

% round/floor/ceil/fix to step
roundfunc = str2func(lower(direction)); % create handle to selected rounding function
out = ((roundfunc((in/(10^(stepDIG)))*stepdiv))/stepdiv)*(10^(stepDIG));

end