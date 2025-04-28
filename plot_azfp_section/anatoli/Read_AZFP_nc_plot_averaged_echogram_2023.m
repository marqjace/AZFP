clearvars

%----------------------------------------------------------------------------------------------------------
gmt_local_offset = 0.29166667; % add 7 hours to local time to get GMT (in Julianday)

month_days = [31;28;31;30;31;30;31;31;30;31;30;31];

% Load Sunset/Sunrise times for year 2022
%=========================================
sunrise_sunset = load('D:\gliders\processing\AZFP_Chris\Sunrise-Sunset_2023n.dat');
YDays = sum(month_days);
JD_sunrise = NaN(YDays,1);
JD_sunset = NaN(YDays,1);
yd=1;
for mm=1:12
    for dd=1:month_days(mm)
        sunrise=sunrise_sunset(dd,mm*2);
        sunset=sunrise_sunset(dd,mm*2+1);
        % Sunrise times
        JD_sunrise(yd) = tmtonumb_d( 23, mm, dd, fix(sunrise/100), ((sunrise/100)-fix(sunrise/100))*100, 0 ); 
        % Sunset times
        JD_sunset(yd) = tmtonumb_d( 23, mm, dd, fix(sunset/100), ((sunset/100)-fix(sunset/100))*100, 0 ); 
        yd=yd+1;
    end
end

% Load glider data
Gl=load('D:\gliders\processing\AZFP_Chris\WA_202305241820-deployment_osu592_pass3.mat');

% Interpolate Lat, Lon to eliminate NaNs
Gl.Lon_int=interp1(Gl.Julday(~isnan(Gl.Lon)),Gl.Lon(~isnan(Gl.Lon)), Gl.Julday,'linear','extrap');
Gl.Lat_int=interp1(Gl.Julday(~isnan(Gl.Lat)),Gl.Lat(~isnan(Gl.Lat)), Gl.Julday,'linear','extrap');

% Get indexes of zig-zag lines
Ln=get_min_max(Gl.Lon_int,0.07);


% zig-zag line to process
Line=5;
LName='Line #5';

% Find indexes of beginning and end of selected transect
istart = Ln(Line);
iend = Ln(Line+1);

% Time frame for selected transect (convert to local time)
T1 = Gl.Julday(istart);
T2 = Gl.Julday(iend);
Lon1 = Gl.Lon_int(istart);
Lon2 = Gl.Lon_int(iend);
Line_text=sprintf('Time (UTC):   %02d:%02d, %02d/%02d/%4d - %02d:%02d, %02d/%02d/%4d', ...
                  Gl.Cur_Time(istart,4), Gl.Cur_Time(istart,5), Gl.Cur_Time(istart,2),Gl.Cur_Time(istart,3),Gl.Cur_Time(istart,1), ...
                  Gl.Cur_Time(iend,4), Gl.Cur_Time(iend,5), Gl.Cur_Time(iend,2),Gl.Cur_Time(iend,3),Gl.Cur_Time(iend,1) );

% If "heading" <0 glider is heading West, if heading is > 0 glider is hading East
heading = Lon2-Lon1;
%==========================================================================
% Find Sunset start time for selected transect (Local Time)
istart_bar = find( (JD_sunset >= T1 - gmt_local_offset) & (JD_sunset <= T2 - gmt_local_offset) );
if(~isempty(istart_bar))
    TM_SunSet = JD_sunset(istart_bar);
else
    TM_SunSet=T1;
end
% Find Sunrise start time for selected transect (Local Time)
istop_bar = find( (JD_sunrise <= T2 - gmt_local_offset) & (JD_sunrise >= T1 - gmt_local_offset));
TM_SunRise = JD_sunrise(istop_bar);

if( (length(TM_SunSet)==1) && (TM_SunSet == T1) )
    Lon_SunSet_bar=Lon1;
else
    for k=1:length(TM_SunSet)
        tmp = find((Gl.Julday - gmt_local_offset) >= TM_SunSet(k));
        Lon_SunSet_bar(k) = Gl.Lon_int(tmp(1));
    end
end

for k=1:length(TM_SunRise)
    tmp = find((Gl.Julday - gmt_local_offset) <= TM_SunRise(k));
    Lon_SunRise_bar(k) = Gl.Lon_int(tmp(length(tmp)));
end
if(length(Lon_SunRise_bar) < length(Lon_SunSet_bar))
    Lon_SunRise_bar=[Lon_SunRise_bar,Lon2];
elseif(length(Lon_SunRise_bar) > length(Lon_SunSet_bar))
    Lon_SunSet_bar=[Lon1,Lon_SunSet_bar];
end

%===============2023===========
dirin='D:\gliders\processing\AZFP_Chris\processed_2023\proc_nc';
dr = dir(fullfile('D:\gliders\processing\AZFP_Chris\processed_2023\proc_nc', '*.proc.nc'));
%===============2023===========

nf={dr.name}';

date_base = datenum(1900, 1, 1, 0, 0, 0); % 1970-01-01 00:00:00

% Find AZFP start and end file number for selected transect
FN1 = 0;
FN2 = 0;
for f=1:length(nf)
    Ping_time=ncread([ dirin '\' char(nf(f))], 'ping_time');
    TM = datevec(datestr((Ping_time)/3600/24+date_base));
    JD = tmtonumb_d( TM(:,1)-2000, TM(:,2), TM(:,3), TM(:,4), TM(:,5), TM(:,6) );

    if( (( T1 <= JD(1) ) || ( T1 <= JD(length(JD)) )) && FN1 == 0 )
        FN1=f;
        JD1=(JD(1));
    end
    if( ((( T2 >= JD(1) ) && ( T2 <= JD(length(JD)) )) || (( T2 <= JD(1) ) && ( T2 <= JD(length(JD)) ))) && FN2 == 0 )
        FN2=f;
        JD2=JD(length(JD));
    end
    if(FN1>0 && FN2>0)
        break
    end
end

if(~exist('JD2','var'))
    JD2=JD(length(JD));
    FN2=length(nf);
end
%======================================
disp(['T1=' num2str(T1) '; T2=' num2str(T2)])
disp(['JD1=' num2str(JD1) '; JD2=' num2str(JD2)])
disp(['FN1=' num2str(FN1) '; FN2=' num2str(FN2)])

%%
echo_range=ncread([ dirin '\' char(nf(1))], 'echo_range');

Echo_range_67=echo_range(1,:,1)';
Echo_range_120=echo_range(1,:,2)';
Echo_range_200=echo_range(1,:,3)';

Sv_67=[];
Sv_120=[];
Sv_200=[];
Ping_time=[];
Ln_PT=[]; % number of vertical samples in a dive
for f=FN1:FN2
    PT=ncread([ dirin '\' char(nf(f))], 'ping_time');
    PT=PT;
    Ln_PT=[Ln_PT;length(PT)];
    Ping_time=[Ping_time;PT];
    Sv=ncread([ dirin '\' char(nf(f))], 'Sv');
    Sv_67=cat(2,Sv_67, Sv(:,:,1));
    Sv_120=cat(2,Sv_120,Sv(:,:,2));
    Sv_200=cat(2,Sv_200,Sv(:,:,3));
end
%%
% Find indexes of last sample in a dive
PT_ind(1)=1;
for i=2:length(Ln_PT)+1
    PT_ind(i,1)=sum(Ln_PT(1:i-1))+1;
end
%%
% Time lag between Glider and AZFP for each zig-zag line
TM_lags =[505;505;505;505;505;505]; %Lags for 2023 zig-zag lines

TM_lag = TM_lags(Line);

% AZFP Julian day time stamps
Cur_Time = datevec(datestr((Ping_time+TM_lag)/3600/24+date_base));
Julday = tmtonumb_d( Cur_Time(:,1)-2000, Cur_Time(:,2), Cur_Time(:,3), Cur_Time(:,4), Cur_Time(:,5), Cur_Time(:,6) );

% Interpolate Glider data to AZFP time
Gl_Depth = interp1(Gl.Julday(~isnan(Gl.Depth)), Gl.Depth(~isnan(Gl.Depth)), Julday,'linear','extrap');
Gl_Depth(Gl_Depth<0) = 0;
Gl_Bottom_Depth = interp1(Gl.Julday(~isnan(Gl.Bottom_depth)), Gl.Bottom_depth(~isnan(Gl.Bottom_depth)), Julday,'linear','extrap');
Gl_Lat = interp1(Gl.Julday(~isnan(Gl.Lat)), Gl.Lat(~isnan(Gl.Lat)), Julday,'linear','extrap');
Gl_Lon = interp1(Gl.Julday(~isnan(Gl.Lon)), Gl.Lon(~isnan(Gl.Lon)), Julday,'linear','extrap');

%=========================================================================
% Modify matrix of Sv to remove AZFP data above glider depth for each ping
% time
[ix,iy]=size(Sv_67);
Sv_67_new=nan(ix,iy);
Sv_120_new=nan(ix,iy);
Sv_200_new=nan(ix,iy);

updn=get_min_max(Gl_Depth,10);
%=========================================================================
% 67 kHz
ishift_67(1,1:iy) = 1;
icut_67(1,1:iy) = 1;
ishift_67 = ishift_67';
icut_67 = icut_67';

for icol=1:iy
    ind = find(Gl_Depth(icol) >= Echo_range_67);
    ind_b = find(Gl_Bottom_Depth(icol)+2 <= Echo_range_67);
    if(~isempty(ind_b))
        icut_67(icol)=ind_b(1);
    else
        icut_67(icol)=ix;
    end

    if(~isempty(ind) && length(ind)>1)
        ishift_67(icol)=ind(length(ind));
        Sv_67_new(ishift_67(icol):ix,icol) = Sv_67(1:ix-ishift_67(icol)+1,icol);
        Sv_67_new(icut_67(icol):ix,icol) = NaN;
    else
        Sv_67_new(1:icut_67(icol),icol) = Sv_67(1:icut_67(icol),icol);
    end
end
clear ind ind_b

Sv_67_aver=nan(ix,length(PT_ind)-1);

for k=1:length(PT_ind)-1
    tmp=mean(Sv_67_new(:,PT_ind(k):PT_ind(k+1)-1),2,'omitnan');
    Gl_Julday(k)=mean(Julday(PT_ind(k):PT_ind(k+1)-1));
    Gl_Lat_aver(k)=mean(Gl_Lat(PT_ind(k):PT_ind(k+1)-1));
    Gl_Lon_aver(k)=mean(Gl_Lon(PT_ind(k):PT_ind(k+1)-1));
    Gl_Depth_aver(k)=mean(Gl_Depth(PT_ind(k):PT_ind(k+1)-1));
    Gl_Bottom_Depth_aver(k)=mean(Gl_Bottom_Depth(PT_ind(k):PT_ind(k+1)-1));
    Sv_67_aver(:,k) = tmp;
end
%%
%=========================================================================
% 120 kHz
ishift_120(1:length(iy)) = 0;
icut_120(1:length(iy)) = 0;
for icol=1:iy
    ind = find(Gl_Depth(icol) >= Echo_range_120);
    ind_b = find(Gl_Bottom_Depth(icol)+2 <= Echo_range_120);
    if(~isempty(ind_b))
        icut_120(icol)=ind_b(1);
    else
        icut_120(icol)=ix;
    end
    if(~isempty(ind) && length(ind)>1)
        ishift_120(icol)=ind(length(ind));
        Sv_120_new(ishift_120(icol):ix,icol) = Sv_120(1:ix-ishift_120(icol)+1,icol);
        Sv_120_new(icut_120(icol):ix,icol) = NaN;
    else
        Sv_120_new(1:icut_120(icol),icol) = Sv_120(1:icut_120(icol),icol);
    end
end
clear ind ind_b

Sv_120_aver=nan(ix,length(PT_ind)-1);

for k=1:length(PT_ind)-1
    tmp=mean(Sv_120_new(:,PT_ind(k):PT_ind(k+1)-1),2,'omitnan');
    Sv_120_aver(:,k) = tmp;
end

%=========================================================================
% 200 kHz
ishift_200(1:length(iy)) = 0;
icut_200(1:length(iy)) = 0;
for icol=1:iy
    ind = find(Gl_Depth(icol) >= Echo_range_200);
    ind_b = find(Gl_Bottom_Depth(icol)+2 <= Echo_range_200);
    if(~isempty(ind_b))
        icut_200(icol)=ind_b(1);
    else
        icut_200(icol)=ix;
    end
    if(~isempty(ind) && length(ind)>1)
        ishift_200(icol)=ind(length(ind));
        Sv_200_new(ishift_200(icol):ix,icol) = Sv_200(1:ix-ishift_200(icol)+1,icol);
        Sv_200_new(icut_200(icol):ix,icol) = NaN;
    else
        Sv_120_new(1:icut_200(icol),icol) = Sv_200(1:icut_200(icol),icol);
    end
end
clear ind ind_b

Sv_200_aver=nan(ix,length(PT_ind)-1);

for k=1:length(PT_ind)-1
    tmp=mean(Sv_200_new(:,PT_ind(k):PT_ind(k+1)-1),2,'omitnan');
    Sv_200_aver(:,k) = tmp;
end

%%
XX=0.13;
YY=[0.91,0.67,0.4,0.13];
WW=0.71;
HH=0.23;

xlim = [-125.0 -124.17]; % 2023

ymax=130;
Eco_lim=[-90 -60];

figure
set(gcf,'Position',[100,50,900,1024]);

subplot(4,1,1, 'Color',[240/255 240/255 240/255])
patch([Lon1;Lon1;Lon2;Lon2;Lon1],[0;1;1;0;0],'k-', 'linewidth', 1,'FaceColor','white','EdgeColor','black');
hold on

for k=1:length(Lon_SunSet_bar)
   patch([Lon_SunSet_bar(k);Lon_SunSet_bar(k);Lon_SunRise_bar(k);Lon_SunRise_bar(k);Lon_SunSet_bar(k)],[0;1;1;0;0],'k-', 'linewidth',1,'FaceColor','black');
end
dl = abs(Lon1-Lon2)/3;
y = [1.5 1.5];
if(heading<0)
    x = [Lon2+dl Lon1-dl];
    drawArrow = @(x,y,varargin) quiver( x(1),y(1),x(2)-x(1),y(2)-y(1),0, varargin{:} );    
    drawArrow(x,y,'linewidth',3,'color','k', 'ShowArrowHead', 0, 'Marker','<');
else
    x = [Lon2-dl Lon1+dl];
%    x = [Lon1-dl Lon2+dl];
    drawArrow = @(x,y,varargin) quiver( x(1),y(1),x(2)-x(1),y(2)-y(1),0, varargin{:} );    
    drawArrow(x,y,'linewidth',3,'color','k', 'ShowArrowHead', 0, 'Marker','>');
end

axbar = gca;
set(axbar,'XColor',[1.0 1.0 1.0],'YColor',[1.0 1.0 1.0], 'xlim', xlim, 'ylim', [0 2]);
set(axbar,'XTickLabel',' ','YTickLabel',' ', 'XTick', [], 'YTick',[], 'XColor', [240/255 240/255 240/255], 'YColor', [240/255 240/255 240/255]);


title([LName ', ' Line_text]);

PosBar = get(axbar,'Position');
set(axbar, 'Position', [XX YY(1) WW 0.04]);

%===================================================================================================================
subplot(4,1,2)

h1=pcolor(Gl_Lon_aver, Echo_range_67(1:1:ix),Sv_67_aver(1:1:ix,:)); shading flat;

colormap jet;
ylabel('Depth range, (m)');
cb67=colorbar('FontSize',10);
cb67.Label.String = ['Volume backscattering strength'; '67 kHz, (Sv re 1 m-1) [dB]    '];
clim(Eco_lim);
hold on;

% add 2m to glider altimeter on the plot to better visualize AZFP bottom detection
plot(Gl.Lon(istart:iend), Gl.Bottom_depth(istart:iend)+2,'k-','linewidth',2, 'markersize',6);

ax1=gca;
set(ax1,'YColor',[0.0 0.0 0.0],'YAxisLocation','left', 'ydir','reverse', 'xlim', xlim, 'ylim', [0 ymax]);
set(ax1,'XTickLabel',' ', 'Ylim', [0 ymax], 'YTick',0:20:ymax,'Fontsize',12, 'Fontweight','bold');
ylabel('Depth, (m)', 'FontSize',12, 'FontWeight','bold')    % y axis label

ax1b=axes('Position',get(ax1,'Position'),'YAxisLocation','right', ...
    'Color','none','YColor','k', 'Fontsize',12, 'Fontweight', ...
    'bold','xlim', xlim,  'ylim', [0 ymax], 'YTick',0:20:ymax, 'ydir','reverse','YTickLabel',' ');

X=Gl.Lon(istart:iend);
Y=Gl.Press(istart:iend);
Z=Gl.SigT(istart:iend);

% Plot SigT contour lines (used in plot_realtime_sec.m)
Sigt_grid_X=min(Gl_Lon):0.001:max(Gl_Lon);
Sigt_grid_Y=0:0.1:110;
[XI,YI] = meshgrid(Sigt_grid_X,Sigt_grid_Y);
ind=find(~isnan(X) & ~isnan(Y) & ~isnan(Z));
hold on
if(length(ind)>1)
    ZI = griddata(X(ind),Y(ind),Z(ind),XI,YI);
    ZIs=smoothdata(ZI,1, 'movmean');
    for c=19:0.5:28
        [C,h] = contour(XI,YI,ZIs, [c c], 'LineColor', [0.99 0.99 0.99], 'LineWidth', 1,'ShowText','on','Parent',ax1b);
        clabel (C,h, 'FontSize',10, 'Color', [0.99 0.99 0.99]);
    end
end
Pos1=get(ax1,'Position');

%===================================================================================================================
subplot(4,1,3)
pcolor(Gl_Lon_aver, Echo_range_120(1:1:ix),Sv_120_aver(1:1:ix,:)); shading flat;

colormap jet;
ylabel('Depth range, (m)');
cb120=colorbar('FontSize',10);
cb120.Label.String = ['Volume backscattering strength'; ...
    '120 kHz, (Sv re 1 m-1) [dB]   '];
clim(Eco_lim);

hold on;
% add 2m to glider altimeter on the plot to better visualize AZFP bottom detection
plot(Gl.Lon(istart:iend), Gl.Bottom_depth(istart:iend)+2,'k-','linewidth',2, 'markersize',6);

ax2=gca;
set(ax2,'YColor',[0.0 0.0 0.0],'YAxisLocation','left', 'ydir', 'reverse', 'xlim', xlim, 'ylim', [0 ymax]);
set(ax2,'XTickLabel',' ', 'Ylim', [0 ymax], 'YTick',0:20:ymax, 'Fontsize',12, 'Fontweight','bold');
ylabel('Depth, (m)', 'FontSize',12, 'FontWeight','bold')    % y axis label

ax2b=axes('Position',get(ax2,'Position'),'YAxisLocation','right', ...
    'Color','none','YColor','k', 'Fontsize',12, 'Fontweight', ...
    'bold','xlim', xlim,  'ylim', [0 ymax], 'YTick',0:20:ymax, 'ydir','reverse','YTickLabel',' ');

% Plot SigT contour lines (used in plot_realtime_sec.m)
Sigt_grid_X=min(Gl_Lon):0.001:max(Gl_Lon);
Sigt_grid_Y=0:0.1:110;
[XI,YI] = meshgrid(Sigt_grid_X,Sigt_grid_Y);
ind=find(~isnan(X) & ~isnan(Y) & ~isnan(Z));
hold on
if(length(ind)>1)
    ZI = griddata(X(ind),Y(ind),Z(ind),XI,YI);
    ZIs=smoothdata(ZI,1, 'movmean');
    for c=19:0.5:28
        [C,h] = contour(XI,YI,ZIs, [c c], 'LineColor', [0.99 0.99 0.99], 'LineWidth', 1,'ShowText','on','Parent',ax2b);
        clabel (C,h, 'FontSize',10, 'Color', [0.99 0.99 0.99]);
    end
end
Pos2=get(ax1,'Position');

%===================================================================================================================
subplot(4,1,4)
pcolor(Gl_Lon_aver, Echo_range_200(1:1:ix),Sv_200_aver(1:1:ix,:)); shading flat;

colormap jet;
ylabel('Depth range, (m)');
xlabel('Longitude');
cb200=colorbar('FontSize',10);
cb200.Label.String = ['Volume backscattering strength'; ...
    '200 kHz, (Sv re 1 m-1) [dB]   '];
clim(Eco_lim);

hold on;
% add 2m to glider altimeter on the plot to better visualize AZFP bottom detection
plot(Gl.Lon(istart:iend), Gl.Bottom_depth(istart:iend)+2,'k-','linewidth',2, 'markersize',6);

ax3=gca;
set(ax3,'YColor',[0.0 0.0 0.0],'YAxisLocation','left', 'ydir', ...
    'reverse', 'xlim', xlim, 'ylim', [0 ymax]);
set(ax3,'XTickLabel',' ', 'Ylim', [0 ymax], 'YTick',0:20:ymax, ...
    'Fontsize',12, 'Fontweight','bold');
ylabel('Depth, (m)', 'FontSize',12, 'FontWeight','bold')    % y axis label

ax3b=axes('Position',get(ax3,'Position'),'YAxisLocation','right', ...
    'Color','none','YColor','k', 'Fontsize',12, 'Fontweight', ...
    'bold','xlim', xlim,  'ylim', [0 ymax], 'YTick',0:20:ymax, 'ydir','reverse','YTickLabel',' ');

% Plot SigT contour lines (used in plot_realtime_sec.m)
Sigt_grid_X=min(Gl_Lon):0.001:max(Gl_Lon);
Sigt_grid_Y=0:0.1:110;
[XI,YI] = meshgrid(Sigt_grid_X,Sigt_grid_Y);
ind=find(~isnan(X) & ~isnan(Y) & ~isnan(Z));
hold on
if(length(ind)>1)
    ZI = griddata(X(ind),Y(ind),Z(ind),XI,YI);
    ZIs=smoothdata(ZI,1, 'movmean');
    for c=19:0.5:28
        [C,h] = contour(XI,YI,ZIs, [c c], 'LineColor', [0.99 0.99 0.99], 'LineWidth', 1,'ShowText','on','Parent',ax3b);
        clabel (C,h, 'FontSize',10, 'Color', [0.99 0.99 0.99]);
    end
end
Pos3=get(ax3,'Position');

% Adjust subplots positions
set(axbar, 'Position', [XX YY(1) WW 0.04]);
set(ax1, 'Position', [XX YY(2) WW    HH]);
set(ax1b, 'Position', [XX YY(2) WW    HH]);
set(ax2, 'Position', [XX YY(3) WW    HH]);
set(ax2b, 'Position', [XX YY(3) WW    HH]);
set(ax3, 'Position', [XX YY(4) WW    HH]);
set(ax3b, 'Position', [XX YY(4) WW    HH]);
