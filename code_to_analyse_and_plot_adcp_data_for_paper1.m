



addpath(genpath('C:\Users\a5278\Documents\MATLAB\matlab_functions'))


latlim = [77.9 82.3];
lonlim = [2 25];
ibcaofile='C:\Users\a5278\Documents\MATLAB\matlab_functions\ibcao\IBCAO_V3_30arcsec_SM.grd';

in=ncinfo(ibcaofile);

x=ncread(ibcaofile,'x');
y=ncread(ibcaofile,'y');
[ibcao.lon,ibcao.lat]=meshgrid(x,y);
ilat=ibcao.lat';
ilon=ibcao.lon';
idepth=ncread(ibcaofile,'z');

[ibcao.lon,ibcao.lat]=meshgrid(x(x>lonlim(1)&x<lonlim(2)),y(y>latlim(1)&y<latlim(2)));

ibcao.depth=idepth( ilon>lonlim(1)&ilon<lonlim(2) & ilat>latlim(1)&ilat<latlim(2) );
ibcao.depth=reshape(ibcao.depth,size(ibcao.lon,2),size(ibcao.lon,1));
ibcao.lat=ibcao.lat';
ibcao.lon=ibcao.lon';

 [Z, refvec] = geoloc2grid( ibcao.lat,  ibcao.lon, ibcao.depth, mean(diff(x)));

 ma=ones(20); % smooth out variation smaller then 18.46 km!
 h = 1/numel(ma)*ma;
   z_smooth = filter2(h,Z);
   [aspect, slope, gradN, gradE] = gradientm(z_smooth, refvec);
   
   
clear x y ilat ilon idepth

load('adcp_all_2014_2017.mat')
load('ctd.mat')
dv=datevec(adcp.time);
load('p_boxes.mat')

%%
figure(2)
clf
hold on
set(gcf,'color',[1 1 1])

m_proj('lambert','long',lonlim,'lat',latlim); 
 m_contour(ibcao.lon,ibcao.lat,ibcao.depth,[-6000:300:0],'color',[.5 .5 .5]);

% [C,H]=m_contour(ibcao.lon,ibcao.lat,ibcao.depth,[-4000,-3000,-2000,-1000,-400,-200],'color',[.5 .5 .5]);
% clabel(C,H,'color',[.5 .5 .5])


m_gshhs_h('patch',[.9 .9 .9]);
m_grid('xlabeldir','end','fontsize',10);

 ix=adcp.dist>0.5 &  dv(:,2)>6 &  dv(:,2)<10 & dv(:,1)==2014;
 a1=m_plot(adcp.lon(ix),adcp.lat(ix),'.','color',[255, 102, 102]./255)
 ix=adcp.dist>0.5 &  dv(:,2)>6 &  dv(:,2)<10 & dv(:,1)==2015;
 a2=m_plot(adcp.lon(ix),adcp.lat(ix),'.','color',[153, 255, 153]./255 )
 ix=adcp.dist>0.5 &  dv(:,2)>6 &  dv(:,2)<10 & dv(:,1)==2016;
 a3=m_plot(adcp.lon(ix),adcp.lat(ix),'.','color',[102, 153, 255]./255 )
 ix=adcp.dist>0.5 &  dv(:,2)>6 &  dv(:,2)<10 & dv(:,1)==2017;
 a4=m_plot(adcp.lon(ix),adcp.lat(ix),'.','color',[255, 153, 255]./255 )
 
for i=1:numel(p_lon)
    m_plot(p_lon{i},p_lat{i},'-k','linewidth',1.5)
m_text( mean(p_lon{i}) , mean(p_lat{i}),num2str(i),'fontweight','bold','fontsize',15)

end
 
 a1=m_plot([90,90],[90,90],'r-','linewidth',2,'color',[255, 102, 102]./255)
 a2=m_plot([90,90],[90,90],'g-','linewidth',2,'color',[153, 255, 153]./255 )
 a3=m_plot([90,90],[90,90],'b-','linewidth',2,'color',[102, 153, 255]./255 )
 a4=m_plot([90,90],[90,90],'m-','linewidth',2,'color',[255, 153, 255]./255 )
 
legend([a1,a2,a3,a4],'2014','2015','2016','2017')


% fram strait
stationlist=[540:547];
clear s a ix_transect
[a]=find(ismember(ctd2014_station,stationlist));
s(ctd2014_station(a))=a;
ix_transect=s(stationlist);
la=ctd2014_lat(ix_transect);
lo=ctd2014_lon(ix_transect);
h=m_plot(lo,la,'^r','markerfacecolor','r','markersize',6)
uistack(h,'top');

%hinlopen
stationlist=[551   552   553   554     556   557];
clear s a ix_transect
[a]=find(ismember(ctd2014_station,stationlist));
s(ctd2014_station(a))=a;
ix_transect=s(stationlist);
la=ctd2014_lat(ix_transect);
lo=ctd2014_lon(ix_transect);
h=m_plot(lo,la,'^r','markerfacecolor','r','markersize',6)
uistack(h,'top');

%fram2
stationlist=[597,596,595,594,593,591,590];
clear s a ix_transect
[a]=find(ismember(ctd2014_station,stationlist));
s(ctd2014_station(a))=a;
ix_transect=s(stationlist);
la=ctd2014_lat(ix_transect);
lo=ctd2014_lon(ix_transect);
h=m_plot(lo,la,'^r','markerfacecolor','r','markersize',6)
uistack(h,'top');

% 
% %framstrait2017
stationlist=[95,126];
clear s a ix_transect
[a]=find(ismember(ctd2017_station,stationlist));
s(ctd2017_station(a))=a;
ix_transect=s(stationlist);
la=ctd2017_lat(ix_transect);
lo=ctd2017_lon(ix_transect);
h=m_plot(lo,la,'^m','markerfacecolor','r','markersize',6)
uistack(h,'top');

% h=m_plot(6+44/60,79+40/60,'pb','markerfacecolor','b','markersize',12)
% uistack(h,'top');


m_text(15,79,'Isforden')
m_text(15,79,'Prins Karls Forlandet')
m_text(15,79,'Hinlopen strait and trench')
m_text(15,79,'Yermak Plateau')
m_text(15,79,'Spitsbergen')
m_text(15,79,'Nordauslandet')
% 
% 
%   set(gcf,'PaperPositionMode','auto')
%   print(gcf,'-dpng',['imagefolder/9boxes_map_3'],'-r600') 

%%


addpath(genpath('C:\Users\a5278\Documents\MATLAB\tidal_model\tmd_toolbox'))
Model='C:\Users\a5278\Documents\MATLAB\tidal_model\aotim5_tmd\Model_AOTIM5';

[tide.u,a]=tmd_tide_pred(Model,adcp.time,adcp.lat,adcp.lon,'u',[]);
[tide.v,a]=tmd_tide_pred(Model,adcp.time,adcp.lat,adcp.lon,'v',[]);

ix_d=adcp.depth>00 & adcp.depth<500;
u=nanmean(adcp.u_detide(ix_d,:));
v=nanmean(adcp.v_detide(ix_d,:));

adcp.u_mean_detide=u - (tide.u./100);
adcp.v_mean_detide=v - (tide.v./100);



 
 factor= deg2km(distance(80,10,81,10))/deg2km(distance(80,10,80,11)) 

latlim = [77.9 82.3];
lonlim = [2 25];
[g_lat,g_lon]=meshgrid(latlim(1):.08:latlim(2),lonlim(1):.2:lonlim(2));
errorthreshold=.4;
corr_length=km2deg(25);
dv=datevec(adcp.time);
ix_adcp=adcp.dist>0.5 & ~isnan(adcp.u_mean_detide)' &  ~isnan(adcp.v_mean_detide)' & dv(:,2)>6 &  dv(:,2)<10;
[xi,yi,ui,emu] = objmap(adcp.lat(ix_adcp),adcp.lon(ix_adcp),adcp.u_mean_detide(ix_adcp),g_lat,g_lon,[corr_length,corr_length*factor],errorthreshold);
[xi,yi,vi,emv] = objmap(adcp.lat(ix_adcp),adcp.lon(ix_adcp),adcp.v_mean_detide(ix_adcp),g_lat,g_lon,[corr_length,corr_length*factor],errorthreshold);

load('p_boxes.mat')


bottomdepth = ltln2val(Z, refvec, xi, yi);

uu=ui;
ix=emu<errorthreshold & emv<errorthreshold & bottomdepth<0;
uu(~ix)=NaN;
 [u_grid, u_grid_r] = geoloc2grid( xi,  yi, uu,.05);
 
vv=vi;
ix=emu<errorthreshold & emv<errorthreshold & bottomdepth<0;
vv(~ix)=NaN;
 [v_grid, v_grid_r] = geoloc2grid( xi,  yi, vv,.05);
 
 %% average slope profiles
 depth_vec=-2400:200:-100;

clear c
depth_intervals=[-2500:200:0];
for i_depth=1:numel(depth_intervals)
  [c{i_depth}, ~]=contour(ibcao.lat, ibcao.lon, ibcao.depth,'LevelList',depth_intervals(i_depth) );
end

 % loop trhough boxes
for i_box=1:numel(p_lat)

for i=1:numel(depth_intervals)-1
    
    lat1=c{i}(1,2:end);
    lon1=c{i}(2,2:end);
    lat2=c{i+1}(1,2:end);
    lon2=c{i+1}(2,2:end);   
    ix_box1 = find(inpolygon(lat1,lon1,p_lat{i_box},p_lon{i_box}));
    ix_box2 = find(inpolygon(lat2,lon2,p_lat{i_box},p_lon{i_box}));

%     figure(1)
%     clf
%     hold on
%     plot(lon1(ix_box1),lat1(ix_box1),'.r')
%     plot(lon2(ix_box2),lat2(ix_box2),'.b')

    
    clear d min_d
    if numel(lat1(ix_box1))>0
        
    for k=1:numel(lat1(ix_box1))
    d=distance( lat1(ix_box1(k)), lon1(ix_box1(k)) , lat2(ix_box2), lon2(ix_box2) );
    min_d(k)=deg2km(min(d));   
    end
    segment_km(i_box,i)=median(min_d(min_d<30));
    
    else
      segment_km(i_box,i)=NaN;
    end
    
end

end

%%


clear along_slope_current_median along_slope_current_mean along_slope_current_mode 
clear across_slope_current_median across_slope_current_mean across_slope_current_mode 

iy=1;
for currrent_year=2014:2017



depth_intervals=[-2500:200:0];
depth_vec=-2400:200:-100;

for i_box=1:numel(p_lon)
    
ix_box = inpolygon(adcp.lat,adcp.lon,p_lat{i_box},p_lon{i_box}) &  dv(:,2)>6 &  dv(:,2)<10 & dv(:,1)==currrent_year & adcp.dist>0.5;
box_number_of_profiles(i_box)=sum(ix_box);


  ix_box = inpolygon(adcp.lat,adcp.lon,p_lat{i_box},p_lon{i_box}) &  dv(:,2)>6 &  dv(:,2)<10 & dv(:,1)==currrent_year;
%     ix_box = inpolygon(adcp.lat,adcp.lon,p_lat{i_box},p_lon{i_box}) &  dv(:,2)>6 &  dv(:,2)<10 ;




clear u_dist_gauss_mean v_dist_gauss_mean u_dist_gauss_std v_dist_gauss_std uv_distrib_n_profiles
clear rotated_current across along perc_across perc_along perc_u perc_v



for i_depth=1:(numel(depth_intervals)-1)
    
ix=adcp.dist>0.5 & ix_box & adcp.bottomdepth>depth_intervals(i_depth) & adcp.bottomdepth<depth_intervals(i_depth+1);

asp = ltln2val(aspect, refvec, adcp.lat(ix), adcp.lon(ix));


ixx=find(ix);
clear rotated_current across along 

    
try
    
for k=1:numel(ixx)
    
 for i_deptbin=1:numel(adcp.depth)

    theta = asp(k);
R = [cosd(theta) -sind(theta); sind(theta) cosd(theta)];
rotated_current=[adcp.v_detide(i_deptbin,ixx(k)),adcp.u_detide(i_deptbin,ixx(k))]*R;
across(k,i_deptbin)=(rotated_current(1));
along(k,i_deptbin)=(rotated_current(2));
 end
 
along_slope_current_median(iy,i_box,i_depth,1:numel(adcp.depth))=nanmedian(along);
% along_slope_current_mode(i_box,i_depth,1:numel(adcp.depth))=mode(along);
% along_slope_current_mean(i_box,i_depth,1:numel(adcp.depth))=nanmean(along);

across_slope_current_median(iy,i_box,i_depth,1:numel(adcp.depth))=nanmedian(across);
% across_slope_current_mode(i_box,i_depth,1:numel(adcp.depth))=mode(across);
% across_slope_current_mean(i_box,i_depth,1:numel(adcp.depth))=nanmean(across);
% figure(1)
% clf
% subplot(211)
% imagesc(along')
% colormap(gca,cmocean('balance','pivot'))

% colorbar
% subplot(212)
% imagesc(across')
% colormap(gca,cmocean('balance','pivot')) 
% colorbar

end
end

end

end

iy=iy+1;
end

along_slope_current_median(along_slope_current_median==0)=NaN;
across_slope_current_median(across_slope_current_median==0)=NaN;

% calc transport
iy=1;
for currrent_year=2014:2017
% calc transport
for i=1:numel(p_lon)
    clear cur
    
    cur=nanmean( squeeze(along_slope_current_median(iy,i,:,:))' );
    cur(cur==0)=NaN;
    depth=depth_vec;
    depth(depth<-500)=-500;   
  area(i,:)=abs( depth .*  segment_km(i,:)*1000 );
  transport_alongslope_median(iy,i,:)= area(i,:) .* cur ./10^6;
    
end


iy=iy+1;
end

%%

%% plots

levs=-.5:.05:.5;

figure(3)
clf

ha = tight_subplot(9,4,[.01 .01],[.1 .01],[.01 .01])

axes(ha(1))
set(gca,'clim',[-.5 .5],'color','w','Xcolor','w','Ycolor','w')
colormap(gca,cmocean('balance','pivot')) 
cb=colorbar('south')
ylabel(cb,'Along slope current speed in m s^{-1}')

for iy=2:4

axes(ha(iy))

hold on
box on
grid on
contourf(depth_vec,-adcp.depth, squeeze(along_slope_current_median(iy,9,:,:))',levs)
set(gca,'xtick',depth_vec,'ytick',-700:100:0,'xticklabel',[],'yticklabel',[]) 

set(gca,'clim',[-.5 .5])
colormap(gca,cmocean('balance','pivot')) 
xlim([min(depth_vec),0])
ylim([-700,0])
patch([-700 0 0],[-700 0 -700],[.9 .9 .9])

text(-60,-530,{num2str(2013+iy),'Box 9'},'fontweight','bold','horizontalalignment','right')
end

for iy=1:4
axes(ha(iy+4))
hold on
box on
grid on
contourf(depth_vec,-adcp.depth, squeeze(along_slope_current_median(iy,8,:,:))',levs)

set(gca,'clim',[-.5 .5])
colormap(gca,cmocean('balance','pivot')) 
set(gca,'xtick',depth_vec,'ytick',-700:100:0,'xticklabel',[],'yticklabel',[]) 
xlim([min(depth_vec),0])
ylim([-700,0])
patch([-700 0 0],[-700 0 -700],[.9 .9 .9])

text(-60,-530,{num2str(2013+iy),'Box 8'},'fontweight','bold','horizontalalignment','right')

end


for iy=1:4
axes(ha(iy+8))
hold on
box on
grid on
contourf(depth_vec,-adcp.depth, squeeze(along_slope_current_median(iy,7,:,:))',levs)

set(gca,'clim',[-.5 .5])
colormap(gca,cmocean('balance','pivot')) 
set(gca,'xtick',depth_vec,'ytick',-700:100:0,'xticklabel',[],'yticklabel',[]) 
xlim([min(depth_vec),0])
ylim([-700,0])
patch([-700 0 0],[-700 0 -700],[.9 .9 .9])

text(-60,-530,{num2str(2013+iy),'Box 7'},'fontweight','bold','horizontalalignment','right')

end


for iy=1:4
axes(ha(iy+12))
hold on
box on
grid on
contourf(depth_vec,-adcp.depth, squeeze(along_slope_current_median(iy,6,:,:))',levs)

set(gca,'clim',[-.5 .5])
colormap(gca,cmocean('balance','pivot')) 
set(gca,'xtick',depth_vec,'ytick',-700:100:0,'xticklabel',[],'yticklabel',[]) 
xlim([min(depth_vec),0])
ylim([-700,0])
patch([-700 0 0],[-700 0 -700],[.9 .9 .9])
text(-60,-530,{num2str(2013+iy),'Box 6'},'fontweight','bold','horizontalalignment','right')

end


for iy=1:4
axes(ha(iy+16))
hold on
box on
grid on
contourf(depth_vec,-adcp.depth, squeeze(along_slope_current_median(iy,5,:,:))',levs)
set(gca,'clim',[-.5 .5])
colormap(gca,cmocean('balance','pivot')) 
xlim([min(depth_vec),0])
ylim([-700,0])
set(gca,'xtick',depth_vec,'ytick',-700:100:0,'xticklabel',[],'yticklabel',[]) 
patch([-700 0 0],[-700 0 -700],[.9 .9 .9])
text(-60,-530,{num2str(2013+iy),'Box 5'},'fontweight','bold','horizontalalignment','right')

end

for iy=1:4
axes(ha(iy+20))
hold on
box on
grid on
contourf(depth_vec,-adcp.depth, squeeze(along_slope_current_median(iy,4,:,:))',levs)
set(gca,'clim',[-.5 .5])
colormap(gca,cmocean('balance','pivot')) 
set(gca,'xtick',depth_vec,'ytick',-700:100:0,'xticklabel',[],'yticklabel',[]) 
xlim([min(depth_vec),0])
ylim([-700,0])
patch([-700 0 0],[-700 0 -700],[.9 .9 .9])
text(-60,-530,{num2str(2013+iy),'Box 4'},'fontweight','bold','horizontalalignment','right')

end


for iy=1:4
axes(ha(iy+24))
hold on
box on
grid on
contourf(depth_vec,-adcp.depth, squeeze(along_slope_current_median(iy,3,:,:))',levs)

set(gca,'clim',[-.5 .5])
colormap(gca,cmocean('balance','pivot')) 
set(gca,'xtick',depth_vec,'ytick',-700:100:0,'xticklabel',[],'yticklabel',[]) 
xlim([min(depth_vec),0])
ylim([-700,0])
patch([-700 0 0],[-700 0 -700],[.9 .9 .9])
text(-60,-530,{num2str(2013+iy),'Box 3'},'fontweight','bold','horizontalalignment','right')

end


for iy=1:4
axes(ha(iy+28))
hold on
box on
grid on
contourf(depth_vec,-adcp.depth, squeeze(along_slope_current_median(iy,2,:,:))',levs)
set(gca,'clim',[-.5 .5])
colormap(gca,cmocean('balance','pivot')) 
set(gca,'xtick',depth_vec,'ytick',-700:100:0,'xticklabel',[],'yticklabel',[]) 
xlim([min(depth_vec),0])
ylim([-700,0])
patch([-700 0 0],[-700 0 -700],[.9 .9 .9])
text(-60,-530,{num2str(2013+iy),'Box 2'},'fontweight','bold','horizontalalignment','right')

end


for iy=[1,3,4]
axes(ha(iy+32))
hold on
box on
grid on
contourf(depth_vec,-adcp.depth, squeeze(along_slope_current_median(iy,1,:,:))',levs)
set(gca,'clim',[-.5 .5])
colormap(gca,cmocean('balance','pivot')) 
set(gca,'xtick',depth_vec,'ytick',-700:100:0,'xticklabel',[],'yticklabel',[])

xlim([min(depth_vec),0])
ylim([-700,0])
patch([-700 0 0],[-700 0 -700],[.9 .9 .9])
text(-60,-530,{num2str(2013+iy),'Box 1'},'fontweight','bold','horizontalalignment','right')

end


set(ha(33),'xticklabel',strsplit(num2str(depth_vec),' '),'xticklabelrotation',70)
set(ha(35),'xticklabel',strsplit(num2str(depth_vec),' '),'xticklabelrotation',70)
set(ha(36),'xticklabel',strsplit(num2str(depth_vec),' '),'xticklabelrotation',70)

axes(ha(34))
box on
grid on

set(gca,'position',[0.34 0.10 0.15 0.0900],'xtick',depth_vec,'ytick',-700:100:0,...
    'xticklabel',[],'xticklabelrotation',0,...
    'yticklabelrotation',20,'yticklabel',strsplit(num2str(-700:100:0),' '))

xlim([min(depth_vec),0])
ylim([-700,0])
patch([-700 0 0],[-700 0 -700],[.9 .9 .9])
ylabel({'Water depth','in m'})
xlabel('Bottom (isobath) depth in m')

%  [left bottom width height]
% 
% mkdir('imagefolder')
%   set(gcf,'PaperPositionMode','auto')
%   print(gcf,'-dpng',['imagefolder/9boxes_mean_along_slope_current_comparison6'],'-r700') 
% 
%   

%%

%% second fran strait transect 2014


ctd_depth=ctd2014_depth;

scale=[-.3:.03:.3];

figure(4)
clf
set(gcf,'color',[1 1 1])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

stationlist=[597,596,595,594,593,591,590];
clear s a ix_transect
[a]=find(ismember(ctd2014_station,stationlist));
s(ctd2014_station(a))=a;
ix_transect=s(stationlist);

la=ctd2014_lat(ix_transect);
lo=ctd2014_lon(ix_transect);
[course,dist] = legs(la,lo)
track_dist_km=[0;nm2km(cumsum(dist))];

grid on; box on;
hold on

contourf(ctd2014_lon(ix_transect),-ctd2014_ladcp.depth,ctd2014_ladcp.along_detide(ix_transect,:)',scale,'edgecolor','none')
colormap(cmocean('balance'))
[C,H]=contour(ctd2014_lon(ix_transect),-ctd_depth,ctd2014_pottemp(ix_transect,:)',[-100,2],'color',[0 .7 0],'linewidth',1.5)

[C,H]=contour(ctd2014_lon(ix_transect),-ctd_depth,ctd2014_sigma(ix_transect,:)',[27.7:.1:28],'k')
clabel(C,H,'color','k')



% cb=colorbar

ylim([-1200 0])
% xlim([ctd2014_lon(ix_transect(1)),ctd2014_lon(ix_transect(end))])
xlim([5.5,9.5])

set(gca,'clim',[scale(1) scale(end)])
plot(ctd2014_lon(ix_transect),zeros(size(ix_transect)),'vk')
tic=get(gca,'xtick')
[course,dist] = legs(ones(size(tic))*mean(la),tic)
tickm=[0;nm2km(cumsum(dist))];
% tickm=deg2km(tic-tic(1));
for i=1:numel(tic)
tic_label{i}=[num2str(tic(i)),'^\circ \newline',num2str(tickm(i),'%2.1f'),' km'];
end
set(gca,'tickdir','out','xticklabels',tic_label)
ylabel('Depth in m')
% 
% %%%%%%%%%%%%
%  xlim([5,10])
% tic=get(gca,'xtick')
% [course,dist] = legs(ones(size(tic))*mean(la),tic)
% tickm=[0;nm2km(cumsum(dist))];
% % tickm=deg2km(tic-tic(1));
% for i=1:numel(tic)
% tic_label{i}=[num2str(tic(i)),'^\circ \newline',num2str(tickm(i),'%2.1f'),' km'];
% end
% set(gca,'tickdir','out','xticklabels',tic_label)
% ylabel('Depth in m')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
title('2014')
% 
% 
%   set(gcf,'PaperPositionMode','auto')
%   print(gcf,'-dpng',['imagefolder/fram_box2_currents_and_sigma_comparison'],'-r400') 

%%

%% currents and sigma fram strait section
ctd_depth=ctd2016_depth;

scale=[-.3:.03:.3];

figure(4)
clf
set(gcf,'color',[1 1 1])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(221)

stationlist=[540:547];
clear s a ix_transect
[a]=find(ismember(ctd2014_station,stationlist));
s(ctd2014_station(a))=a;
ix_transect=s(stationlist);

la=ctd2014_lat(ix_transect);
lo=ctd2014_lon(ix_transect);
[course,dist] = legs(la,lo)
track_dist_km=[0;nm2km(cumsum(dist))];

grid on; box on;
hold on

contourf(ctd2014_lon(ix_transect),-ctd2014_ladcp.depth,ctd2014_ladcp.along_detide(ix_transect,:)',scale,'edgecolor','none')
colormap(cmocean('balance'))
[C,H]=contour(ctd2014_lon(ix_transect),-ctd_depth,ctd2014_pottemp(ix_transect,:)',[-100,2],'color',[0 .7 0],'linewidth',1.5)

[C,H]=contour(ctd2014_lon(ix_transect),-ctd_depth,ctd2014_sigma(ix_transect,:)',[27.7:.1:28],'k')
clabel(C,H,'color','k')



% cb=colorbar

ylim([-1200 0])
set(gca,'clim',[scale(1) scale(end)])
plot(ctd2014_lon(ix_transect),zeros(size(ix_transect)),'vk')
plot(ctd2014_lon(ix_transect(6)),zeros(size(ix_transect(6))),'vk','markerfacecolor','k')

 xlim([5,10])
tic=get(gca,'xtick')
[course,dist] = legs(ones(size(tic))*mean(la),tic)
tickm=[0;nm2km(cumsum(dist))];
% tickm=deg2km(tic-tic(1));
for i=1:numel(tic)
tic_label{i}=[num2str(tic(i)),'^\circ \newline',num2str(tickm(i),'%2.1f'),' km'];
end
set(gca,'tickdir','out','xticklabels',tic_label)
ylabel('Depth in m')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
title('2014')

subplot(222)
stationlist=[86:93];
clear s a ix_transect
[a]=find(ismember(ctd2015_station,stationlist));
s(ctd2015_station(a))=a;
ix_transect=s(stationlist);

la=ctd2015_lat(ix_transect);
lo=ctd2015_lon(ix_transect);
[course,dist] = legs(la,lo)

track_dist_km=[0;nm2km(cumsum(dist))];


grid on; box on;
hold on

contourf(ctd2015_lon(ix_transect),-ctd2015_ladcp.depth,ctd2015_ladcp.along_detide(ix_transect,:)',scale,'edgecolor','none')
colormap(cmocean('balance'))
[C,H]=contour(ctd2015_lon(ix_transect),-ctd_depth,ctd2015_pottemp(ix_transect,:)',[-100,2],'color',[0 .7 0],'linewidth',1.5)


[C,H]=contour(ctd2015_lon(ix_transect),-ctd_depth,ctd2015_sigma(ix_transect,:)',[27.7:.1:28],'k')
clabel(C,H,'color','k')

% cb=colorbar

ylim([-1200 0])
set(gca,'clim',[scale(1) scale(end)])
plot(ctd2015_lon(ix_transect),zeros(size(ix_transect)),'vk')
plot(ctd2015_lon(ix_transect(3)),zeros(size(ix_transect(3))),'vk','markerfacecolor','k')

 xlim([5,10])
tic=get(gca,'xtick')
[course,dist] = legs(ones(size(tic))*mean(la),tic)
tickm=[0;nm2km(cumsum(dist))];
% tickm=deg2km(tic-tic(1));
for i=1:numel(tic)
tic_label{i}=[num2str(tic(i)),'^\circ \newline',num2str(tickm(i),'%2.1f'),' km'];
end
set(gca,'tickdir','out','xticklabels',tic_label)
ylabel('Depth in m')

title('2015')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(223)

stationlist=[58,59,60,61];
clear s a ix_transect
[a]=find(ismember(ctd2016_station,stationlist));
s(ctd2016_station(a))=a;
ix_transect=s(stationlist);

la=ctd2016_lat(ix_transect);
lo=ctd2016_lon(ix_transect);
[course,dist] = legs(la,lo)

track_dist_km=[0;nm2km(cumsum(dist))];


grid on; box on;
hold on

contourf(ctd2016_lon(ix_transect),-ctd2016_ladcp.depth,ctd2016_ladcp.along_detide(ix_transect,:)',scale,'edgecolor','none')
colormap(cmocean('balance'))

[C,H]=contour(ctd2016_lon(ix_transect),-ctd_depth,ctd2016_pottemp(ix_transect,:)',[-100,2],'color',[0 .7 0],'linewidth',1.5)

[C,H]=contour(ctd2016_lon(ix_transect),-ctd_depth,ctd2016_sigma(ix_transect,:)',[27.7:.1:28],'k')
clabel(C,H,'color','k')

% cb=colorbar

ylim([-1200 0])
set(gca,'clim',[scale(1) scale(end)])
plot(ctd2016_lon(ix_transect),zeros(size(ix_transect)),'vk')
  xlim([5,10])
tic=get(gca,'xtick')
[course,dist] = legs(ones(size(tic))*mean(la),tic)
tickm=[0;nm2km(cumsum(dist))];
% tickm=deg2km(tic-tic(1));
for i=1:numel(tic)
tic_label{i}=[num2str(tic(i)),'^\circ \newline',num2str(tickm(i),'%2.1f'),' km'];
end
set(gca,'tickdir','out','xticklabels',tic_label)
ylabel('Depth in m')

patch([5.5,7.4,7.4,5.5],[-5,-5,-800,-800],'w','edgecolor','none') 

title('2016')

subplot(224)
grid on; box on;
hold on

[a]=find(ismember(ctd2017_station,95));

scatter(ctd2017_ladcp.along_detide(a,:),-ctd2017_ladcp.depth,30,ctd2017_ladcp.along_detide(a,:),'filled')
plot(ctd2017_ladcp.along_detide(a,:),-ctd2017_ladcp.depth,'-k')
text(ctd2017_ladcp.along_detide(a,30)-.12,-100,'23.08')

[a]=find(ismember(ctd2017_station,126));
scatter(ctd2017_ladcp.along_detide(a,:),-ctd2017_ladcp.depth,30,ctd2017_ladcp.along_detide(a,:),'filled')
plot(ctd2017_ladcp.along_detide(a,:),-ctd2017_ladcp.depth,'-k')
text(ctd2017_ladcp.along_detide(a,30)+.09,-100,'05.09')

set(gca,'clim',[scale(1) scale(end)])
xlim([scale(1) scale(end)])
ylim([-1200,0])
xlabel('Along slope current speed in m s^{-1}')
title('2017')

cb=colorbar
 ylabel(cb,'Along slope current speed in m s^{-1}')
% 
%   set(gcf,'PaperPositionMode','auto')
%   print(gcf,'-dpng',['imagefolder/fram_before_yermak_currents_and_sigma_comparison2'],'-r500') 



%% hinnlopen
ctd_depth=ctd2016_depth;

scale=[-.3:.03:.3];
% figure(13)
% clf
% hold on 
stationlist=[551   552   553 ];
clear s a ix_transect
[a]=find(ismember(ctd2014_station,stationlist));
s(ctd2014_station(a))=a;
ix_transect=s(stationlist);
asp = mean( ltln2val(aspect, refvec, ctd2014_lat(ix_transect), ctd2014_lon(ix_transect)) );
% plot(ctd2014_lat(ix_transect),asp,'-r')





figure(11)
clf


stationlist=[551   552   553   554     556   557];
clear s a ix_transect
[a]=find(ismember(ctd2014_station,stationlist));
s(ctd2014_station(a))=a;
ix_transect=s(stationlist);

clear rotated_current  h_across h_along
for k=1:numel(ix_transect) 
 for i_deptbin=1:numel(ctd2014_ladcp.depth)
    theta = asp;
R = [cosd(theta) -sind(theta); sind(theta) cosd(theta)];
rotated_current=[ctd2014_ladcp.v_detide(ix_transect(k),i_deptbin),ctd2014_ladcp.u_detide(ix_transect(k),i_deptbin)]*R;
h_across(k,i_deptbin)=(rotated_current(1));
h_along(k,i_deptbin)=(rotated_current(2));
 end 
end

subplot(321)
grid on; box on;
hold on

contourf(ctd2014_lat(ix_transect),-ctd2014_ladcp.depth,h_along',scale,'edgecolor','none')
colormap(cmocean('balance'))
[C,H]=contour(ctd2014_lat(ix_transect),-ctd_depth,ctd2014_pottemp(ix_transect,:)',[-100,2],'color',[0 .7 0],'linewidth',1.5)

[C,H]=contour(ctd2014_lat(ix_transect),-ctd_depth,ctd2014_sigma(ix_transect,:)',[27.8:.1:28],'k')
clabel(C,H,'color','k')
% cb=colorbar

set(gca,'clim',[scale(1) scale(end)])
plot(ctd2014_lat(ix_transect),zeros(size(ix_transect)),'vk')

tic=get(gca,'xtick')
tickm=deg2km(tic-tic(1));
for i=1:numel(tic)
tic_label{i}=[num2str(tic(i)),'^\circ \newline',num2str(tickm(i),'%2.1f'),' km'];
end
set(gca,'tickdir','out','xticklabels',tic_label)
ylabel('Depth in m')
set(gca,'xdir','reverse')
ylim([-1000 0])
title('2014 Along slope - northeastward')

subplot(322)
grid on; box on;
hold on

contourf(ctd2014_lat(ix_transect),-ctd2014_ladcp.depth,h_across',scale,'edgecolor','none')
colormap(cmocean('balance'))
[C,H]=contour(ctd2014_lat(ix_transect),-ctd_depth,ctd2014_pottemp(ix_transect,:)',[-100,2],'color',[0 .7 0],'linewidth',1.5)

[C,H]=contour(ctd2014_lat(ix_transect),-ctd_depth,ctd2014_sigma(ix_transect,:)',[27.8:.1:28],'k')
clabel(C,H,'color','k')
% cb=colorbar

set(gca,'clim',[scale(1) scale(end)])
plot(ctd2014_lat(ix_transect),zeros(size(ix_transect)),'vk')

tic=get(gca,'xtick')
tickm=deg2km(tic-tic(1));
for i=1:numel(tic)
tic_label{i}=[num2str(tic(i)),'^\circ \newline',num2str(tickm(i),'%2.1f'),' km'];
end
set(gca,'tickdir','out','xticklabels',tic_label)
ylabel('Depth in m')
set(gca,'xdir','reverse')
ylim([-1000 0])
title('2014 Across slope - southwestward')



%%%%%%
stationlist=[42,44,46,47,49];
clear s a ix_transect
[a]=find(ismember(ctd2016_station,stationlist));
s(ctd2016_station(a))=a;
ix_transect=s(stationlist);

clear rotated_current  h_across h_along
for k=1:numel(ix_transect) 
 for i_deptbin=1:numel(ctd2016_ladcp.depth)
    theta = asp;
R = [cosd(theta) -sind(theta); sind(theta) cosd(theta)];
rotated_current=[ctd2016_ladcp.v_detide(ix_transect(k),i_deptbin),ctd2016_ladcp.u_detide(ix_transect(k),i_deptbin)]*R;
h_across(k,i_deptbin)=(rotated_current(1));
h_along(k,i_deptbin)=(rotated_current(2));
 end 
end

subplot(323)
grid on; box on;
hold on
contourf(ctd2016_lat(ix_transect),-ctd2016_ladcp.depth,h_along',scale,'edgecolor','none')
colormap(cmocean('balance'))
[C,H]=contour(ctd2016_lat(ix_transect),-ctd_depth,ctd2016_pottemp(ix_transect,:)',[-100,2],'color',[0 .7 0],'linewidth',1.5)

[C,H]=contour(ctd2016_lat(ix_transect),-ctd_depth,ctd2016_sigma(ix_transect,:)',[27.8:.1:28],'k')
clabel(C,H,'color','k')
% cb=colorbar

set(gca,'clim',[scale(1) scale(end)])
plot(ctd2016_lat(ix_transect),zeros(size(ix_transect)),'vk')

tic=get(gca,'xtick')
tickm=deg2km(tic-tic(1));
for i=1:numel(tic)
tic_label{i}=[num2str(tic(i)),'^\circ \newline',num2str(tickm(i),'%2.1f'),' km'];
end
set(gca,'tickdir','out','xticklabels',tic_label)
ylabel('Depth in m')
set(gca,'xdir','reverse')
ylim([-1000 0])

title('2016 Along slope - northeastward')

subplot(324)
grid on; box on;
hold on
contourf(ctd2016_lat(ix_transect),-ctd2016_ladcp.depth,h_across',scale,'edgecolor','none')
colormap(cmocean('balance'))
[C,H]=contour(ctd2016_lat(ix_transect),-ctd_depth,ctd2016_pottemp(ix_transect,:)',[-100,2],'color',[0 .7 0],'linewidth',1.5)

[C,H]=contour(ctd2016_lat(ix_transect),-ctd_depth,ctd2016_sigma(ix_transect,:)',[27.8:.1:28],'k')
clabel(C,H,'color','k')
% cb=colorbar

set(gca,'clim',[scale(1) scale(end)])
plot(ctd2016_lat(ix_transect),zeros(size(ix_transect)),'vk')

tic=get(gca,'xtick')
tickm=deg2km(tic-tic(1));
for i=1:numel(tic)
tic_label{i}=[num2str(tic(i)),'^\circ \newline',num2str(tickm(i),'%2.1f'),' km'];
end
set(gca,'tickdir','out','xticklabels',tic_label)
ylabel('Depth in m')
set(gca,'xdir','reverse')
ylim([-1000 0])

title('2016 Across slope - southwestward')

%%%%%%%%%%%%%%%%%%%%%%%%%


 stationlist=[99:105];
% stationlist=[99,100,101,119,102,118,103,104,105];
clear s a ix_transect
[a]=find(ismember(ctd2017_station,stationlist));
s(ctd2017_station(a))=a;
ix_transect=s(stationlist);

clear rotated_current  h_across h_along
for k=1:numel(ix_transect) 
 for i_deptbin=1:numel(ctd2017_ladcp.depth)
    theta = asp;
R = [cosd(theta) -sind(theta); sind(theta) cosd(theta)];
rotated_current=[ctd2017_ladcp.v_detide(ix_transect(k),i_deptbin),ctd2017_ladcp.u_detide(ix_transect(k),i_deptbin)]*R;
h_across(k,i_deptbin)=(rotated_current(1));
h_along(k,i_deptbin)=(rotated_current(2));
 end 
end

subplot(325)

grid on; box on;
hold on

contourf(ctd2017_lat(ix_transect),-ctd2017_ladcp.depth,h_along',scale,'edgecolor','none')
colormap(cmocean('balance'))
[C,H]=contour(ctd2017_lat(ix_transect),-ctd_depth,ctd2017_pottemp(ix_transect,:)',[-100,2],'color',[0 .7 0],'linewidth',1.5)

[C,H]=contour(ctd2017_lat(ix_transect),-ctd_depth,ctd2017_sigma(ix_transect,:)',[27.8:.1:28],'k')
clabel(C,H,'color','k')
% cb=colorbar

set(gca,'clim',[scale(1) scale(end)])
plot(ctd2017_lat(ix_transect),zeros(size(ix_transect)),'vk')

tic=get(gca,'xtick')
tickm=deg2km(tic-tic(1));
for i=1:numel(tic)
tic_label{i}=[num2str(tic(i)),'^\circ \newline',num2str(tickm(i),'%2.1f'),' km'];
end
set(gca,'tickdir','out','xticklabels',tic_label)
ylabel('Depth in m')
set(gca,'xdir','reverse')
ylim([-1000 0])
title('2017 Along slope - northeastward')

subplot(326)

grid on; box on;
hold on

h_across(h_across<scale(1))=scale(1);

contourf(ctd2017_lat(ix_transect),-ctd2017_ladcp.depth,h_across',scale,'edgecolor','none')
colormap(cmocean('balance'))
[C,H]=contour(ctd2017_lat(ix_transect),-ctd_depth,ctd2017_pottemp(ix_transect,:)',[-100,2],'color',[0 .7 0],'linewidth',1.5)

[C,H]=contour(ctd2017_lat(ix_transect),-ctd_depth,ctd2017_sigma(ix_transect,:)',[27.8:.1:28],'k')
clabel(C,H,'color','k')
% cb=colorbar

set(gca,'clim',[scale(1) scale(end)])
plot(ctd2017_lat(ix_transect),zeros(size(ix_transect)),'vk')

tic=get(gca,'xtick')
tickm=deg2km(tic-tic(1));
for i=1:numel(tic)
tic_label{i}=[num2str(tic(i)),'^\circ \newline',num2str(tickm(i),'%2.1f'),' km'];
end
set(gca,'tickdir','out','xticklabels',tic_label)
ylabel('Depth in m')
set(gca,'xdir','reverse')
ylim([-1000 0])
title('2017 Across slope - southwestward')
% 
% mkdir('imagefolder')
%   set(gcf,'PaperPositionMode','auto')
%   print(gcf,'-dpng',['imagefolder/ladcp_hinlopen_rotated_330deg_2_detided'],'-r500') 

%%

figure(2)
clf
hold on
set(gcf,'color',[1 1 1])


m_proj('lambert','long',lonlim,'lat',latlim);
m_gshhs_h('patch',[.8 .8 .8]);
m_grid('xlabeldir','end','fontsize',10);


% orig1=[80]
% orig2=[5]
% [latout,lonout] = reckon(orig1,orig2,km2deg(30),90);


 [C,h]=m_contour(ibcao.lon,ibcao.lat,ibcao.depth,[-4000,-3000,-2000,-1000,-800,-600,-400,-200],'color',[.5 .5 .5]);
 clabel(C,h,'color',[.5 .5 .5]);

 
 factor= deg2km(distance(80,10,81,10))/deg2km(distance(80,10,80,11)) 

latlim = [77.9 82.3];
lonlim = [2 25];
[g_lat,g_lon]=meshgrid(latlim(1):.08:latlim(2),lonlim(1):.2:lonlim(2));
errorthreshold=.4;
corr_length=km2deg(25);
dv=datevec(adcp.time);
ix_adcp=adcp.dist>0.5 & ~isnan(adcp.u_mean_detide)' &  ~isnan(adcp.v_mean_detide)' & dv(:,2)>6 &  dv(:,2)<10;
[xi,yi,ui,emu] = objmap(adcp.lat(ix_adcp),adcp.lon(ix_adcp),adcp.u_mean_detide(ix_adcp),g_lat,g_lon,[corr_length,corr_length*factor],errorthreshold);
[xi,yi,vi,emv] = objmap(adcp.lat(ix_adcp),adcp.lon(ix_adcp),adcp.v_mean_detide(ix_adcp),g_lat,g_lon,[corr_length,corr_length*factor],errorthreshold);

% ui=ui-gu;
% vi=vi-gv;

bottomdepth = ltln2val(Z, refvec, xi, yi);
ix=emu<errorthreshold & emv<errorthreshold & bottomdepth<0;

c=sqrt(ui(ix).^2+vi(ix).^2);
 vecs = m_vec(1, yi(ix),xi(ix),ui(ix),vi(ix),c, 'shaftwidth', .7, 'headangle', 30, 'edgeclip', 'on');
% vecs = m_vec(1, yi(ix),xi(ix),gu(ix),gv(ix),'k', 'shaftwidth', .1, 'headangle', 10);

uistack(vecs);
cb=colorbar('north')
xlabel(cb,'Interpolated current speed in m s^{-1}')

set(gca,'clim',[0.1 .4])
  colormap(gca,cmocean('thermal'))

            m_text(22,79.92,'0.1 m s^{-1}')
 vecs = m_vec(1, 22,79.8,.1,0,'shaftwidth', 0.7,  'headangle', 30, 'edgeclip', 'on');
uistack(vecs);
 vecs = m_vec(1, 22,79.7,.4,0,'shaftwidth', 0.7,  'headangle', 30, 'edgeclip', 'on');
uistack(vecs);
m_text(23,79.6,'0.4 m s^{-1}')

% 
% mkdir('imagefolder')
%   set(gcf,'PaperPositionMode','auto')
%   print(gcf,'-dpng',['objective_mapping_current_speed2'],'-r500') 

%%


%% averge interpofield in boxes

interpfield.xi=xi(:);
interpfield.yi=yi(:);
interpfield.ui=ui(:);
interpfield.vi=vi(:);

interpfield.bottomdepth = ltln2val(Z, refvec, interpfield.xi, interpfield.yi);


ma=ones(20); % smooth out variation smaller then 18.46 km!
 h = 1/numel(ma)*ma;
   z_smooth = filter2(h,Z);
   [aspect, slope, gradN, gradE] = gradientm(z_smooth, refvec);
   
load('p_boxes.mat')

clear along_slope_current_median along_slope_current_mean along_slope_current_mode 
clear across_slope_current_median across_slope_current_mean across_slope_current_mode 

depth_intervals=[-2500:200:0];
depth_vec=-2400:200:-100;

for i_box=1:numel(p_lon)
        
  ix_box = inpolygon(interpfield.xi,interpfield.yi,p_lat{i_box},p_lon{i_box});

clear u_dist_gauss_mean v_dist_gauss_mean u_dist_gauss_std v_dist_gauss_std uv_distrib_n_profiles
clear rotated_current across along perc_across perc_along perc_u perc_v



for i_depth=1:(numel(depth_intervals)-1)
    
ix= ix_box & interpfield.bottomdepth>depth_intervals(i_depth) & interpfield.bottomdepth<depth_intervals(i_depth+1);

asp = ltln2val(aspect, refvec, interpfield.xi(ix),  interpfield.yi(ix));

ixx=find(ix);
clear rotated_current across along 
if isempty(ixx)
    along_slope_current_median(i_box,i_depth)=NaN;
across_slope_current_median(i_box,i_depth)=NaN;
else
for k=1:numel(ixx)
   
    theta = asp(k);
R = [cosd(theta) -sind(theta); sind(theta) cosd(theta)];
rotated_current=[interpfield.vi(ixx(k)),interpfield.ui(ixx(k))]*R;
across(k)=(rotated_current(1));
along(k)=(rotated_current(2));
 end

interpfield.along_slope_current_median(i_box,i_depth)=nanmedian(along);
interpfield.across_slope_current_median(i_box,i_depth)=nanmedian(across);
end
end

end

%%% calc transport
for i=1:numel(p_lon)
    clear cur area
    
    cur=interpfield.along_slope_current_median(i,:);
    depth=depth_vec;
    depth(depth<-500)=-500;   
  area(i,:)=abs( depth .*  segment_km(i,:)*1000 );
  interpfield.transport_alongslope_median(i,:)= area(i,:) .* cur ./10^6;
    
end

%%


%% transport plot2
figure(6)
clf


for i=1:numel(p_lat)
subplot(3,3,i)

hold on
box on
grid on

plot(depth_vec,  squeeze(transport_alongslope_median(1,i,:)),'.r','markersize',15)
plot(depth_vec,  squeeze(transport_alongslope_median(2,i,:)),'.g','markersize',15)
plot(depth_vec,  squeeze(transport_alongslope_median(3,i,:)),'.b','markersize',15)
plot(depth_vec,  squeeze(transport_alongslope_median(4,i,:)),'.m','markersize',15)
plot(depth_vec,  squeeze(interpfield.transport_alongslope_median(i,:)),'.-k','linewidth',2,'markersize',15)

xlim([min(depth_vec),0])
ylim([-.5, 1.65])
text(-2250,1.4,['Box Nr. ',num2str(i)],'fontweight','bold')
% ylabel('Sv (10^6 m^3 s^{-1})')
p_box_interptrans_sum(i)=nansum(interpfield.transport_alongslope_median(i,:));
p_box_interptrans_sumnorth(i)=nansum(interpfield.transport_alongslope_median(i,interpfield.transport_alongslope_median(i,:)>0));
p_box_interptrans_sumsouth(i)=nansum(interpfield.transport_alongslope_median(i,interpfield.transport_alongslope_median(i,:)<0));


end

xlabel('Isobath depth in m')
ylabel('Sv (10^6 m^3 s^{-1})')

legend('2014','2015','2016','2017','Interpolation')
% 
% 
% mkdir('imagefolder')
%   set(gcf,'PaperPositionMode','auto')
%   print(gcf,'-dpng',['imagefolder/9boxes_mean_across_slope_transport5'],'-r400') 
%   savefig(gcf,'imagefolder/9boxes_mean_across_slope_transport5')
% %

%%

s_lat{1}=[78.1,78.1]; s_lon{1}=[5,10];
s_lat{2}=[78.33,78.33]; s_lon{2}=[5,10];
s_lat{3}=[78.66,78.66]; s_lon{3}=[5,10];
s_lat{4}=[79,79]; s_lon{4}=[4,10];
s_lat{5}=[79.33,79.33]; s_lon{5}=[4,10.5];
s_lat{6}=[79.66,79.66]; s_lon{6}=[4,10.5];

s_lat{7}=[81,80.3]; s_lon{7}=[11,13];

s_lat{8}=[81,80.4]; s_lon{8}=[12.7,14.7];

s_lat{9}=[81,80.4]; s_lon{9}=[14.2,16.2];

s_lat{10}=[81.3,81]; s_lon{10}=[15,19];


%branches
s_lat{11}=[80.1,79.9]; s_lon{11}=[9.5,12];
s_lat{12}=[80.3,81]; s_lon{12}=[10,10];
s_lat{13}=[80.3,79.88]; s_lon{13}=[7,3];

%%


a=deg2km(distance(xi(1,1),yi(1,1),xi(:,1),yi(:,1)) );
b=deg2km(distance(xi(1,1),yi(1,1),xi(1,:),yi(1,:)) );

xlen=mean(diff(a))
ylen=mean(diff(b))

bottomd=ltln2val(Z, refvec, xi,yi);
bottomd(bottomd<-500)=500;

area_x=abs(bottomd)*xlen*1000;
area_y=abs(bottomd)*ylen*1000;

trans_x=area_x.*vi./10^6;
trans_y=area_y.*ui./10^6;

trans_abs=sqrt(trans_x.^2 + trans_y.^2);

figure(2)
clf
hold on
set(gcf,'color',[1 1 1])


m_proj('lambert','long',lonlim,'lat',latlim);
%     m_coast('patch',[0.9 0.95 0.9]);
[C,h]=m_contour(ibcao.lon,ibcao.lat,ibcao.depth,[-6000:200:0],'color',[.5 .5 .5]);
clabel(C,h,[-4000,-2000,-1000:200:-100],'color',[.5 .5 .5]);


m_gshhs_h('patch',[.9 .9 .9]);
m_grid('xlabeldir','end','fontsize',10);

bottomdepth = ltln2val(Z, refvec, xi, yi);

ix=emu<errorthreshold & emv<errorthreshold & bottomdepth<0;
C=trans_abs;
vecs = m_vec(1, yi(ix),xi(ix),ui(ix),vi(ix),C(ix), 'shaftwidth', 1.7, 'headlength', 5);

uistack(vecs);
cmap=cmocean('speed',100);
  colormap(gca,cmap(20:100,:))
    set(gca,'clim',[0 1])


cb=colorbar('north');
ylabel(cb,'Transport in the upper 500 m in Sv (10^6 m^3 s^{-1})')

%%%%%%%%%%%%%%%%%%%%%
for i=1:10
    m_plot(s_lon{i},s_lat{i},'-b','linewidth',2)
     m_text(s_lon{i}(1), s_lat{i}(1),num2str(i),'color','b','fontweight','bold','fontsize',15)
end

for i=11:13
    m_plot(s_lon{i},s_lat{i},'-b','linewidth',2)
end
     m_text(s_lon{11}(1), s_lat{11}(1),'SB','color','b','fontweight','bold','fontsize',15)
     m_text(s_lon{12}(1), s_lat{12}(1),'YPB','color','b','fontweight','bold','fontsize',15)
     m_text(s_lon{13}(1), s_lat{13}(1),'YB','color','b','fontweight','bold','fontsize',15)

for i=1:numel(p_lon)
    m_plot(p_lon{i},p_lat{i},'-k','linewidth',1.5)
m_text( mean(p_lon{i}) , mean(p_lat{i}),num2str(i),'fontweight','bold','fontsize',15)

end

% 
% mkdir('imagefolder')
%   set(gcf,'PaperPositionMode','auto')
%   print(gcf,'-dpng',['imagefolder/transport_estimates_boxes_vs_transects_map'],'-r300') 
%   savefig(gcf,'imagefolder/transport_estimates_boxes_vs_transects_map')

%%



for i=1:numel(s_lat)

    clear section

section.origin=[s_lat{i}(1),s_lon{i}(1)];
section.end=[s_lat{i}(2),s_lon{i}(2)];

nbins=50;

section.az=azimuth(section.origin(1),section.origin(2),section.end(1),section.end(2));
theta = 90-section.az;
R = [cosd(theta) -sind(theta); sind(theta) cosd(theta)];
 

[section.v,section.degdist,section.lat,section.lon]  = mapprofile(v_grid, v_grid_r, [section.origin(1),section.end(1)], [section.origin(2),section.end(2)],nbins);
section.depth= ltln2val(Z, refvec,section.lat,section.lon);
section.depth(section.depth<-500)=-500;
section.area= [deg2km(section.degdist(2));diff(deg2km(section.degdist))] *1000.* abs(section.depth);
section.transp_v=section.area.*section.v;

[section.u,section.degdist,section.lat,section.lon]  = mapprofile(u_grid, u_grid_r, [section.origin(1),section.end(1)], [section.origin(2),section.end(2)],nbins);
section.transp_u=section.area.*section.u;

clear rotated_current
for d=1:size(section.u,1)
rotated_current(d,:)=[section.u(d),section.v(d)]*R;
end
section.transp_across=rotated_current(:,2).*section.area ;

sum_transp_across(i)=nansum(section.transp_across)./10^6;
sum_transp_u(i)=nansum(section.transp_u)./10^6;
sum_transp_v(i)=nansum(section.transp_v)./10^6;

pos_transp_across(i)=nansum(section.transp_across(section.transp_across>0))./10^6;
pos_transp_u(i)=nansum(section.transp_u(section.transp_u>0))./10^6;
pos_transp_v(i)=nansum(section.transp_v(section.transp_v>0))./10^6;

neg_transp_across(i)=nansum(section.transp_across(section.transp_across<0))./10^6;
neg_transp_u(i)=nansum(section.transp_u(section.transp_u<0))./10^6;
neg_transp_v(i)=nansum(section.transp_v(section.transp_v<0))./10^6;


end

%%


depth_intervals=[-2500:200:0];
depth_vec=-2400:200:-100;

for i_box=1:numel(p_lon)
        
  ix_box = inpolygon(interpfield.xi,interpfield.yi,p_lat{i_box},p_lon{i_box});

clear u_dist_gauss_mean v_dist_gauss_mean u_dist_gauss_std v_dist_gauss_std uv_distrib_n_profiles
clear rotated_current across along perc_across perc_along perc_u perc_v



for i_depth=1:(numel(depth_intervals)-1)
    
ix= ix_box & interpfield.bottomdepth>depth_intervals(i_depth) & interpfield.bottomdepth<depth_intervals(i_depth+1);

% asp = ltln2val(aspect, refvec, interpfield.xi(ix),  interpfield.yi(ix));

ixx=find(ix);
clear rotated_current across along 
if isempty(ixx)
    along_slope_current_median(i_box,i_depth)=NaN;
across_slope_current_median(i_box,i_depth)=NaN;
else
for k=1:numel(ixx)
   
%     theta = asp(k);
    theta = p_az(i_box);

R = [cosd(theta) -sind(theta); sind(theta) cosd(theta)];
rotated_current=[interpfield.vi(ixx(k)),interpfield.ui(ixx(k))]*R;
across(k)=(rotated_current(1));
along(k)=(rotated_current(2));
 end

interpfield.along_slope_current_median(i_box,i_depth)=nanmedian(along);
interpfield.across_slope_current_median(i_box,i_depth)=nanmedian(across);
end
end

end

%%% calc transport
for i=1:numel(p_lon)
    clear cur area
    
    cur=interpfield.across_slope_current_median(i,:);
    depth=depth_vec;
    depth(depth<-500)=-500;   
  area(i,:)=abs( depth .*  segment_km(i,:)*1000 );
  interpfield.transport_alongslope_median_constant_az(i,:)= area(i,:) .* cur ./10^6;
    
end

for i=1:numel(p_lat)
p_box_interptrans_sum_constantaz(i)=nansum(interpfield.transport_alongslope_median_constant_az(i,:));
% p_box_interptrans_sumnorth(i)=nansum(interpfield.transport_alongslope_median_constant_az(i,interpfield.transport_alongslope_median_constant_az(i,:)>0));
% p_box_interptrans_sumsouth(i)=nansum(interpfield.transport_alongslope_median_constant_az(i,interpfield.transport_alongslope_median_constant_az(i,:)<0));
% 

end

%%

figure(9)
clf
hold on
box on
grid on

% boxes
plot(1:4,p_box_interptrans_sumnorth(1:4),'.--r','markersize',20)
p1=plot(4:8,p_box_interptrans_sumnorth(4:8),'.-r','markersize',20)
plot(8:9,p_box_interptrans_sumnorth(8:9),'.--r','markersize',20)

%plot(1:9,p_box_interptrans_sumsouth,'.-b','markersize',20)
%plot(1:9,p_box_interptrans_sum,'.-k','markersize',20)
plot(1:4,p_box_interptrans_sum_constantaz(1:4),'.--b','markersize',20)
p2=plot(4:8,p_box_interptrans_sum_constantaz(4:8),'.-b','markersize',20)
plot(8:9,p_box_interptrans_sum_constantaz(8:9),'.--b','markersize',20)

x=[linspace(1,5,6),linspace(7,9,4)];
%plot(x,pos_transp_across(1:10),'^--r','markersize',4,'markerfacecolor','r')
%plot(x,neg_transp_across(1:10),'^--b','markersize',4,'markerfacecolor','b')
plot(x(1:5),sum_transp_across(1:5),'.--g','markersize',20,'color',[0 0.7 0])
p3=plot(x(5:9),sum_transp_across(5:9),'.-g','markersize',20,'color',[0 0.7 0])
plot(x(9:10),sum_transp_across(9:10),'.--g','markersize',20,'color',[0 0.7 0])

p4=plot(5,1.7,'^k','markerfacecolor','k')
plot(5,0.7,'^k','markerfacecolor','k')
text(5.15,1.7,'2015')
text(5.15,0.7,'2014')
plot(2,1.2,'^k','markerfacecolor','k')
text(2,1.2,'2014')

plot(5.5,sum_transp_across(11),'.','color',[0 0.7 0],'markersize',20)
plot(5.5,-sum_transp_across(12),'.','color',[0 0.7 0],'markersize',20)
plot(5.5,-sum_transp_across(13),'.','color',[0 0.7 0],'markersize',20)

text(5.5,sum_transp_across(11),'SB','fontweight','bold','color',[0 0.7 0])
text(5.5,-sum_transp_across(12),'YPB','fontweight','bold','color',[0 0.7 0])
text(5.5,-sum_transp_across(13),'YB','fontweight','bold','color',[0 0.7 0])

% text(5.5,-sum_transp_across(13)+sum_transp_across(11),'SB+YPB')

for i=1:9
   text(i-.2,4.1,['Box ',num2str(i)],'rotation',30) 
end

ylim([0 4])
ylabel('Sv (10^6 m^3 s^{-1})')
legend([p1,p2,p3,p4],'Obj. map. box averged net transport - local azimuth','Obj. map. box averged net transport - box averaged azimuth','Obj. map. transect net transport','L-ADCP net transport')
set(gca,'xticklabels',{'Isfjorden 78^\circ N','P.K. forlandet 78.5^\circ N','Kongsfjorden 79^\circ N','79.3^\circ N','79.6^\circ N','Yermak Plateau 80^\circ N','13^\circ E', 'Hinlopen trench 15^\circ E','18^\circ E' },'xticklabelrotation',30)

% %
% 
% mkdir('imagefolder')
%   set(gcf,'PaperPositionMode','auto')
%   print(gcf,'-dpng',['transport_estimates_boxes_vs_transects_comparison2'],'-r300') 
%   savefig(gcf,'transport_estimates_boxes_vs_transects_comparison2')

%%


%% get current field
latlim = [77.9 82.3];
lonlim = [2 25];

[g_lat,g_lon]=meshgrid(latlim(1):.08:latlim(2),lonlim(1):.2:lonlim(2));

errorthreshold=.4;

factor= deg2km(distance(80,10,81,10))/deg2km(distance(80,10,80,11)) 
corr_length=km2deg(25);

dv=datevec(adcp.time);

ix_d=adcp.depth>00 & adcp.depth<500;
u=nanmean(adcp.u_detide(ix_d,:));
v=nanmean(adcp.v_detide(ix_d,:));

ix_adcp=adcp.dist>0.5 & ~isnan(u)' &  ~isnan(v)' & dv(:,2)>6 &  dv(:,2)<10;
[xi,yi,ui,emu] = objmap(adcp.lat(ix_adcp),adcp.lon(ix_adcp),u(ix_adcp),g_lat,g_lon,[corr_length,corr_length*factor],errorthreshold);
[xi,yi,vi,emv] = objmap(adcp.lat(ix_adcp),adcp.lon(ix_adcp),v(ix_adcp),g_lat,g_lon,[corr_length,corr_length*factor],errorthreshold);

div = divergence(ui,vi);
[gu,gv]=gradient(div);

ui=ui+gu;
vi=vi+gv;

 [Z_u, refvec_u] = geoloc2grid( xi,  yi, ui,.1);
 [Z_v, refvec_v] = geoloc2grid( xi,  yi, vi,.1);
 [Z_eu, refvec_u] = geoloc2grid( xi,  yi, emu,.1);
 [Z_ev, refvec_v] = geoloc2grid( xi,  yi, emv,.1);
 
  %% calc worms for bin plot
latlim = [77.9 82.3];
lonlim = [2 25];
n=7000;
t_lat=rand([1,n])*(latlim(2)-latlim(1)) + latlim(1);
t_lon=rand([1,n])*(lonlim(2)-lonlim(1)) + lonlim(1);

bottomdepth = ltln2val(Z, refvec, t_lat,t_lon);
emu_t = ltln2val(Z_eu, refvec_u, t_lat,t_lon);
emv_t = ltln2val(Z_ev, refvec_u, t_lat,t_lon);

errorthreshold=.4

ix=bottomdepth<0 ;%& emu_t<errorthreshold & emv_t<errorthreshold;
t_lat( ~ix)=[];
t_lon(~ix)=[];

t_u=ltln2val(Z_u, refvec_u, t_lat,t_lon);
t_v=ltln2val(Z_v, refvec_v, t_lat,t_lon);

% [t_angle,t_arclen]=cart2pol(t_u,t_v);

a=atan2(t_u,t_v);
t_angle=real(asin(t_u.*sin(a)./t_v))
t_arclen=sqrt(t_u.^2 + t_v.^2);

coe=0.1;

k=3000
clear worm_lat worm_lon worm_uv
worm_lat(1,:)=t_lat;
worm_lon(1,:)=t_lon;
% worm_uv(1,:)=zeros(size(t_lon));

for i=1:k
    
t_u=ltln2val(Z_u, refvec_u, worm_lat(i,:),worm_lon(i,:));
t_v=ltln2val(Z_v, refvec_v, worm_lat(i,:),worm_lon(i,:));
a=atan2(t_u,t_v);
a= (a - pi/2);
a(a<0)=a(a<0)+2*pi;
t_angle=a + pi/2;
t_arclen=sqrt(t_u.^2 + t_v.^2);

[t_lat_new,t_lon_new] = reckon(worm_lat(i,:),worm_lon(i,:),t_arclen*coe,rad2deg(t_angle) );

worm_lat(i+1,:)=t_lat_new + randn(size(t_lat_new))*.005;
worm_lon(i+1,:)=t_lon_new + randn(size(t_lon_new))*.005;

worm_uv(i,:)=sqrt(power(t_u,2)+power(t_v,2));

end

bottomdepth = ltln2val(Z, refvec, worm_lat,worm_lon);
emu_t = ltln2val(Z_eu, refvec_u, worm_lat,worm_lon);
emv_t = ltln2val(Z_ev, refvec_u, worm_lat,worm_lon);

errorthreshold2=.5

ix=bottomdepth<0 & emu_t<errorthreshold2 & emv_t<errorthreshold2;
worm_lat( ~ix)=NaN;
worm_lon(~ix)=NaN;

% p_lat=[78,80.6,81.5,81,78,70,78];
% p_lon=[3,3,15,22,22,15,3];
% 
% % ix=bottomdepth<0 & emu_t<errorthreshold2 & emv_t<errorthreshold2;
% 
% in = inpolygon(worm_lat,worm_lon,p_lat,p_lon) & bottomdepth<-50 ;
% 
% worm_lat( ~in)=NaN;
% worm_lon(~in)=NaN;


%% svalbard branch and YPB
clear ix_svalbard_branch ix_yp_branch ix_y_branch
p_lat=[80.2,80.05];
p_lon=[10,12];

for i=1:size(worm_lat,2) 
  ixv=isnan(worm_lat(:,i));
[xi,yi] = polyxpoly(worm_lat(~ixv,i),worm_lon(~ixv,i),p_lat,p_lon);
if isempty(xi)
   ix_svalbard_branch(i)=0; 
else
   ix_svalbard_branch(i)=1;     
end
end

% p_lat=[80.5,80.9];
% p_lon=[10,10];
p_lat=[80.2,80.8];
p_lon=[8,9];
% m_plot(p_lon,p_lat,'-c')

for i=1:size(worm_lat,2) 
  ixv=isnan(worm_lat(:,i));
[xi,yi] = polyxpoly(worm_lat(~ixv,i),worm_lon(~ixv,i),p_lat,p_lon);
if isempty(xi)
   ix_yp_branch(i)=0; 
else
   ix_yp_branch(i)=1;     
end
end

% p_lat=[79.88,80.3];
% p_lon=[3,7];
p_lat=[79.88,80.3];
p_lon=[3,6.6];
%  m_plot(p_lon,p_lat,'-k')

for i=1:size(worm_lat,2) 
  ixv=isnan(worm_lat(:,i));
[xi,yi] = polyxpoly(worm_lat(~ixv,i),worm_lon(~ixv,i),p_lat,p_lon);
if isempty(xi)
   ix_y_branch(i)=0; 
else
   ix_y_branch(i)=1;     
end
end

%%

%% bin tracks 
latlim = [77.9 82.3];
lonlim = [2 25];
b_vec_lat=latlim(1):.05:latlim(2);
b_vec_lon=lonlim(1):.15:lonlim(2);

[b_lat,b_lon]=meshgrid(b_vec_lat,b_vec_lon);

ixx_sb=repmat(ix_svalbard_branch,[size(worm_lat,1),1]);
ixx_ypb=repmat(ix_yp_branch,[size(worm_lat,1),1]);
ixx_yb=repmat(ix_y_branch,[size(worm_lat,1),1]);

b_counts= nan(size(b_lat));
b_counts_sb= nan(size(b_lat));
b_counts_ypb= nan(size(b_lat));
b_counts_yb= nan(size(b_lat));

for i1=1:numel(b_vec_lat)
    for i2=1:numel(b_vec_lon)
        
 n_worms= sum( worm_lat(:) > b_vec_lat(i1)-.025  &  worm_lat(:) < b_vec_lat(i1)+.025 & worm_lon(:) > b_vec_lon(i2)-.075 &  worm_lon(:) < b_vec_lon(i2)+.075 );
 b_counts(i2,i1)=n_worms;
 
 n_worms= sum( ixx_sb(:) & worm_lat(:) > b_vec_lat(i1)-.025  &  worm_lat(:) < b_vec_lat(i1)+.025 & worm_lon(:) > b_vec_lon(i2)-.075 &  worm_lon(:) < b_vec_lon(i2)+.075 );
 b_counts_sb(i2,i1)=n_worms;
 
  n_worms= sum( ixx_ypb(:) & worm_lat(:) > b_vec_lat(i1)-.025  &  worm_lat(:) < b_vec_lat(i1)+.025 & worm_lon(:) > b_vec_lon(i2)-.075 &  worm_lon(:) < b_vec_lon(i2)+.075 );
 b_counts_ypb(i2,i1)=n_worms;
 
  n_worms= sum( ixx_yb(:) & worm_lat(:) > b_vec_lat(i1)-.025  &  worm_lat(:) < b_vec_lat(i1)+.025 & worm_lon(:) > b_vec_lon(i2)-.075 &  worm_lon(:) < b_vec_lon(i2)+.075 );
 b_counts_yb(i2,i1)=n_worms;
    end
    i1./numel(b_vec_lat)
end

%%



b_counts_yb_n=b_counts_yb ./ max(b_counts_yb(:));
b_counts_ypb_n=b_counts_ypb ./ max(b_counts_ypb(:));
b_counts_sb_n=b_counts_sb ./ max(b_counts_sb(:));

plotfield=nan(size(b_counts_yb_n));
cbins=10;
for i1=1:numel(b_vec_lat)
    for i2=1:numel(b_vec_lon)
       
  [a,b]=max([b_counts_sb_n(i2,i1),b_counts_ypb_n(i2,i1),b_counts_yb_n(i2,i1)]);
  
  if b==2 & b_counts_yb_n(i2,i1)>b_counts_sb_n(i2,i1)
      b=3;
  end
        if a==0
plotfield(i2,i1)=NaN;
        else
            switch  b             
                case 1
[c]=hist(b_counts_sb_n(i2,i1),[linspace(0,.002,cbins-1),1] );
plotfield(i2,i1)=find(c);
                  case 2
[c]=hist(b_counts_ypb_n(i2,i1),[linspace(0,.01,cbins-1),1] );
plotfield(i2,i1)=find(c)+cbins;
                  case 3
[c]=hist(b_counts_yb_n(i2,i1),[linspace(0,.01,cbins-1),1]   );
plotfield(i2,i1)=find(c)+cbins*2;
            end
        end
      
    end
end

c1=cmocean('amp',cbins)
c2=cmocean('speed',cbins)
c3=cmocean('dense',cbins)

cmap=[c1;c2;c3];

% cut off
plotfield(b_lon>18)=NaN;
plotfield(b_lat>81.3)=NaN;

figure(4)
clf
hold on
set(gcf,'color',[1 1 1])
latlim = [77.9 81.5];
lonlim = [2 20];
m_proj('lambert','long',lonlim,'lat',latlim);
%     m_coast('patch',[0.9 0.95 0.9]);
 [C,h]=m_contour(ibcao.lon,ibcao.lat,ibcao.depth,[-5000:200:0],'color',[.5 .5 .5]);
clabel(C,h,[-4000,-2000,-1000:50:-100],'color',[.5 .5 .5]);

 m_pcolor(b_lon,b_lat,plotfield);
 shading flat
colormap(gca,cmap)
  set(gca,'clim',[0 cbins*3])
 m_gshhs_h('patch',[.9 .9 .9]);
m_grid('xlabeldir','end','fontsize',10);

p_lat=[80.2,80.05];
p_lon=[10,12];
m_plot(p_lon,p_lat,'-r','linewidth',2)

p_lat=[80.2,80.8];
p_lon=[8,9];
m_plot(p_lon,p_lat,'-g','linewidth',2)

p_lat=[79.88,80.3];
p_lon=[3,6.6];
m_plot(p_lon,p_lat,'-b','linewidth',2)


% % mkdir('imagefolder')
%   set(gcf,'PaperPositionMode','auto')
%   print(gcf,'-dpng',['imagefolder/particel_tracking_pathways6'],'-r500') 
%   savefig(gcf,'imagefolder/particel_tracking_pathways6')

