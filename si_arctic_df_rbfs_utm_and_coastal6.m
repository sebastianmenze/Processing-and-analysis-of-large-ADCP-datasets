
 clear all
addpath(genpath('C:\Users\a5278\Documents\MATLAB\matlab_functions'))


%svalbard
% latlim = [78 81.5];
% lonlim = [0 25];

%yermak
latlim = [79.5 81];
lonlim = [5 13];

%hinlopen
% latlim = [79.5 81.2];
% lonlim = [12 20];

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
dv=datevec(adcp.time);

%%
clear lon_coast lat_coast
 ma=ones(30); % smooth out variation smaller then 18.46 km!
 h = 1/numel(ma)*ma;
   z_smooth = filter2(h,ibcao.depth);
   
   counter=1;
   for cval=[0];
%    cval=0
[C,h]=contour(ibcao.lat,ibcao.lon,z_smooth,[cval,cval]);

% figure(1)
% clf
% hold on

ix0=find(C(1,:)==cval);
ix0len=C(2,ix0)

ixdel=ix0len<100;
ix0(ixdel)=[];
ix0len(ixdel)=[];


for i=1:numel(ix0)
    
       ix= ix0(i)+1 : 10 : ix0(i)+ix0len(i);
 
    lat_coast{counter}=C(1, ix ) ;
    lon_coast{counter}=C(2, ix ) ;

% plot(C(1, ix ),C(2,ix),'.-')
counter=counter+1;
end
   end




%%%% utm conversions
 
 z1 = utmzone( [mean(lat_coast{1}),mean(lon_coast{1})  ]);
 [ellipsoid,estr] = utmgeoid(z1);
 utmstruct = defaultm('utm');
utmstruct.zone = '18T';
utmstruct.geoid = ellipsoid;
utmstruct = defaultm(utmstruct);


    for i=1:numel(lat_coast)
[x_utm_coast{i},y_utm_coast{i}] = mfwdtran(utmstruct,lat_coast{i},lon_coast{i});
        
    end
    
    
%%

i=1;
for y=2014:2017
    for m=1:12
 nobs(i)= sum( dv(:,1)==y & dv(:,2)==m & adcp.dist>0.3 & adcp.lat>latlim(1) & adcp.lat<latlim(2) & adcp.lon>lonlim(1) & adcp.lon<lonlim(2) )  ;
 datenums(i)=datenum(y,m,1);
  i=i+1;
end
end

figure(8)
clf
bar(datenums,nobs)

datamonth=datenums(nobs>100);
dvdatamonth=datevec(datamonth);



for itime=1:numel(datamonth)
% the observations
 %ix=adcp.dist>0.5  & dv(:,2)>6 &  dv(:,2)<10 & ~isnan(adcp.u_mean)' & adcp.lat>latlim(1) & adcp.lat<latlim(2) & adcp.lon>lonlim(1) & adcp.lon<lonlim(2);
  ix=adcp.dist>0.3 & dv(:,1)==dvdatamonth(itime,1) &  dv(:,2)==dvdatamonth(itime,2) & ~isnan(adcp.u_mean)' & adcp.lat>latlim(1) & adcp.lat<latlim(2) & adcp.lon>lonlim(1) & adcp.lon<lonlim(2);
a_lat=adcp.lat(ix);
a_lon=adcp.lon(ix);
ix_d=adcp.depth>100 & adcp.depth<200;
a_u=nanmean(adcp.u_detide(ix_d,ix));
a_v=nanmean(adcp.v_detide(ix_d,ix));

ix_nan=isnan(a_u)  | isnan(a_v);

a_v(ix_nan)=[];
a_u(ix_nan)=[];
a_lat(ix_nan)=[];
a_lon(ix_nan)=[];

% a_v=[a_v(1:10:end)];
% a_u=[a_u(1:10:end)];
% a_lat=[a_lat(1:10:end)];
% a_lon=[a_lon(1:10:end)];

%   [a_utm_x,a_utm_y,f]=ll2utm([a_lat,a_lon]);


% the centers
% n_centers=600;
n_centers=round(numel(a_u)*.07 );

ixr=round(linspace(1,numel(a_lat),n_centers)) ;
c_lat=a_lat(ixr);
c_lon=a_lon(ixr);

% [idx,C] = kmeans([a_lon,a_lat],n_centers);
% c_lat=C(:,2);
% c_lon=C(:,1);

for ic=1:numel(lat_coast)
c_lat=[c_lat ; lat_coast{ic}(1:10:end)' ];
c_lon=[c_lon ; lon_coast{ic}(1:10:end)' ];
end
% ixr=randi(numel(a_lat),[n_centers,1])
% c_lat=a_lat(ixr);
% c_lon=a_lon(ixr);
%   [x_c_utm,y_c_utm,f]=ll2utm([c_lat,c_lon]);

% [c_lat,c_lon] = utm2ll(C(:,1),C(:,2),32);
% the grid
% latlim = [77.9 82.3];
% lonlim = [2 25];
[g_lat,g_lon]=meshgrid(latlim(1):.08:latlim(2),lonlim(1):.2:lonlim(2));

x_c=c_lon;
y_c=c_lat;

%plot 

figure(1)
clf
hold on

plot(a_lon,a_lat,'.k')
plot(c_lon,c_lat,'or')
 plot(g_lon,g_lat,'.b')
 
 %%%% utm conversions
 
 z1 = utmzone( [mean(a_lat),mean(a_lon)  ]);
 [ellipsoid,estr] = utmgeoid(z1);
 utmstruct = defaultm('utm');
utmstruct.zone = '18T';
utmstruct.geoid = ellipsoid;
utmstruct = defaultm(utmstruct);

[x_utm_a,y_utm_a] = mfwdtran(utmstruct,a_lat,a_lon);
[x_c_utm,y_c_utm] = mfwdtran(utmstruct,c_lat,c_lon);

% x_c_utm=[x_c_utm ; x_utm_coast{1}(1:10:end)' ; x_utm_coast{2}(1:10:end)' ]
% y_c_utm=[y_c_utm ; y_utm_coast{1}(1:10:end)' ; y_utm_coast{2}(1:10:end)' ]

figure(1)
clf
hold on

plot(x_utm_a,y_utm_a,'.k')
plot(x_c_utm,y_c_utm,'or')
     for i=1:numel(lat_coast)
plot(x_utm_coast{i},y_utm_coast{i},'.-');
        
    end
 %% train RBFs with observations UTM
% 
% a=[1,1,1;2,2,2]
% a=a(:)
% reshape(a,[2,3])

data=[a_v;a_u];
data_adcp=data(:);
x_obs=x_utm_a;
y_obs=y_utm_a;

% determine radiuses/differences
clear r rbf_adcp

dataindex=1;
for i=1:numel(x_obs)
    
        %distance to centers
%         r=sqrt( (x_obs(i)-x_c).^2 + (y_obs(i)-y_c).^2 )';
r=deg2km(distance(a_lat(i),a_lon(i),c_lat,c_lon))'.*1000;
        logr=log(r);
        logr(r==0)=0;
             
        % u part
   rbf_adcp(dataindex,:)  =   [ 1,0,x_obs(i),y_obs(i),0 , [r.^2.*(12*logr+7) , zeros(size(r))]  + [ -(8*logr+6) .* (x_obs(i)-x_c_utm)' .* (x_obs(i)-x_c_utm)' ,  -(8*logr+6) .* (x_obs(i)-x_c_utm)' .* (y_obs(i)-y_c_utm)'   ]  ] ;
        %v part
   rbf_adcp(dataindex+1,:)= [0,1,-y_obs(i),0,x_obs(i) , [zeros(size(r)), r.^2.*(12*logr+7)]   + [ -(8*logr+6) .* (y_obs(i)-y_c_utm)' .* (x_obs(i)-x_c_utm)' ,  -(8*logr+6) .* (y_obs(i)-y_c_utm)' .* (y_obs(i)-y_c_utm)'   ]  ] ;
    dataindex=dataindex+2;    
end

%% costal segments

coastcounter=1;
clear rbf_coast

for icoast=1:numel(x_utm_coast)
    
    r2=deg2km(distance( lat_coast{icoast}(1) , lon_coast{icoast}(1) ,c_lat,c_lon))'.*1000;
        logr2=log(r2);
        logr2(r2==0)=0;
        
for i=2:numel(lat_coast{icoast})
    
        %distance to centers
        r1= deg2km(distance( lat_coast{icoast}(i) , lon_coast{icoast}(i) ,c_lat,c_lon))'.*1000; ;
        logr1=log(r1);
        logr1(r1==0)=0;

        alphapart= [ y_utm_coast{icoast}(i)-y_utm_coast{icoast}(1) , -x_utm_coast{icoast}(i)+x_utm_coast{icoast}(1) , x_utm_coast{icoast}(i)*y_utm_coast{icoast}(i)- x_utm_coast{icoast}(1)*y_utm_coast{icoast}(1) , .5*y_utm_coast{icoast}(i).^2-.5*y_utm_coast{icoast}(1).^2 , -.5*x_utm_coast{icoast}(i).^2+.5*x_utm_coast{icoast}(1).^2 ] ;     
   
%         centerpart=[ r1.^2.*(1+4*logr1) .* (x_utm_coast{icoast}(i)-x_c)' - r2.^2.*(1+4*logr2) .* (x_utm_coast{icoast}(1)-x_c)' , -(  r1.^2.*(1+4*logr1) .* (y_utm_coast{icoast}(i)-y_c)' - r2.^2.*(1+4*logr2) .* (y_utm_coast{icoast}(1)-y_c)'  ) ];
        centerpart=[ r1.^2.*(1+4*logr1) .* (y_utm_coast{icoast}(i)-y_c_utm)' - r2.^2.*(1+4*logr2) .* (y_utm_coast{icoast}(1)-y_c_utm)'  ,  - (  r1.^2.*(1+4*logr1) .* (x_utm_coast{icoast}(i)-x_c_utm)' - r2.^2.*(1+4*logr2) .* (x_utm_coast{icoast}(1)-x_c_utm)' ) ];

        rbf_coast(coastcounter,:)  = [alphapart,centerpart];   
        coastcounter=coastcounter+1;
end
end

data_coast=zeros(coastcounter-1,1)

 %%  ectra constraint

n_params=size(rbf_adcp,2);
dataindex=1;
clear a b rbf_side
    a=zeros([1,n_params]);
    b=zeros([1,n_params]);
for i=1:numel(x_c_utm)
    a(5+dataindex)=sum( [1,0,x_c_utm(i),y_c_utm(i),0] );
    b(5+dataindex+1)=sum( [0,1,-y_c_utm(i),0,x_c_utm(i)] );
        dataindex=dataindex+2;    
end
rbf_side(1,:)  =  a ;
rbf_side(2,:)  =  b ;

data_side=[0;0];

%%
%   params= [data_adcp',data_coast' ]/ [ rbf_adcp',rbf_coast']  ;%use mrdivide to solve system of equations. For large systems it may 
 
params= [data_adcp',data_coast',data_side' ]/ [ rbf_adcp',rbf_coast',rbf_side']  ;%use mrdivide to solve system of equations. For large systems it may 

% params2 = linsolve(rbf,data);

% sum(params1==params1)/numel(params1)

%% evaluate rbfs


% eval rbfs for grid

% latlim = [77.9 82.3];
% lonlim = [2 25];
 [g_lat,g_lon]=meshgrid(latlim(1):.04:latlim(2),lonlim(1):.15:lonlim(2));
% [g_lat,g_lon]=meshgrid(latlim(1):.06:latlim(2),lonlim(1):.3:lonlim(2));
latgridv=g_lat(:);
longridv=g_lon(:);
% xgridv=g_lon(:);
% ygridv=g_lat(:);

[xgridv,ygridv] = mfwdtran(utmstruct,latgridv,longridv);

clear r rbf dist2centers
dataindex=1;
for i=1:size(xgridv,1)
    
%         r=sqrt( (xgridv(i)-x_c).^2 + (ygridv(i)-y_c).^2 )';
r=deg2km(distance(latgridv(i),longridv(i),c_lat,c_lon))'.*1000;

dist2centers(i,:)=r;
       logr=log(r);
        logr(r==0)=0;
            % u part
   rbf(dataindex,:)  =   [ 1,0,xgridv(i),ygridv(i),0 , [r.^2.*(12*logr+7) , zeros(size(r))]  + [ -(8*logr+6) .* (xgridv(i)-x_c_utm)' .* (xgridv(i)-x_c_utm)' ,  -(8*logr+6) .* (xgridv(i)-x_c_utm)' .* (ygridv(i)-y_c_utm)'   ]  ] ;
        %v part
   rbf(dataindex+1,:)= [0,1,-ygridv(i),0,xgridv(i) , [zeros(size(r)), r.^2.*(12*logr+7)]   + [ -(8*logr+6) .* (ygridv(i)-y_c_utm)' .* (xgridv(i)-x_c_utm)' ,  -(8*logr+6) .* (ygridv(i)-y_c_utm)' .* (ygridv(i)-y_c_utm)'   ]  ] ;
    dataindex=dataindex+2   ;
end

estimate=params*rbf';%apply weights
estimate=reshape(estimate,[2,size(xgridv,1)]);

%%

ugrid=reshape( estimate(2,:) , size(g_lat) );
vgrid=reshape( estimate(1,:) , size(g_lat) );
dgrid=reshape( sum( dist2centers(:,1:n_centers)<10000 ,2) , size(g_lat) );
ixclose=dgrid>0;
% 
% figure(9)
% clf
% dgrid=reshape( sum(dist2centers<20000 ,2) , size(g_lat) );
% imagesc(dgrid)
% colorbar

m_ugrid(itime,:,:)=ugrid;
m_vgrid(itime,:,:)=vgrid;
m_ixclose(itime,:,:)=ixclose;


% end
% 
% figure(2)
% clf
% hold on
% % quiver(x,y,u,v,'k')
% %  quiver(x_obs,y_obs,u_obs,v_obs,'r')
% % plot(x_c,y_c,'or')
% quiver(xgridv,ygridv,ugrid(:),vgrid(:),'k')
% 
% figure(3)
% clf
% hold on
% div=-divergence(vgrid,ugrid);
% pcolor(g_lon,g_lat,div)
% shading flat
% colormap(cmocean('balance'))
% plot(c_lon,c_lat,'or')
%      for i=1:numel(lat_coast)
% plot(lon_coast{i},lat_coast{i},'.-');        
%     end
% colorbar


%%

% for itime=1:numel(datamonth)

%     
%    ugrid=squeeze( m_ugrid(itime,:,:));
% vgrid=squeeze( m_vgrid(itime,:,:));
% ixclose=squeeze( m_ixclose(itime,:,:) );

figure(4)
clf
hold on
set(gcf,'color',[1 1 1])


m_proj('lambert','long',lonlim,'lat',latlim);
m_grid('xlabeldir','end','fontsize',10);

% m_pcolor(g_lon,g_lat,div)
% shading flat
% colormap(cmocean('balance'))
set(gca,'clim',[-0.2 .2])

m_gshhs_h('patch',[.8 .8 .8]);

[C,h]=m_contour(ibcao.lon,ibcao.lat,ibcao.depth,[-4000,-3000,-2000,-1000,-800,-600,-400,-200],'color',[.5 .5 .5]);
clabel(C,h,'color',[.5 .5 .5]);

bottomdepth = ltln2val(Z, refvec, g_lat,g_lon);
ui=ugrid;
vi=vgrid;



c=sqrt(ui.^2+vi.^2);
ix= bottomdepth<0 & c<1 & ixclose==1 ;


 vecs = m_vec(1, g_lon(ix),g_lat(ix),ui(ix),vi(ix),'k', 'shaftwidth', .7, 'headangle', 30, 'edgeclip', 'on');
  vecs = m_vec(1, a_lon(1:5:end),a_lat(1:5:end),a_u(1:5:end),a_v(1:5:end),'r','shaftwidth', .7, 'headangle', 30, 'edgeclip', 'on');

uistack(vecs);
%  cb=colorbar('north')
 % xlabel(cb,'Interpolated current speed in m s^{-1}')

%   colormap(gca,cmocean('thermal'))

            m_text(16,79.8,'0.1 m s^{-1}')
 vecs = m_vec(1, 16,79.78,.1,0,'shaftwidth', 0.7,  'headangle', 30, 'edgeclip', 'on');
uistack(vecs);
 vecs = m_vec(1, 16,79.79,.4,0,'shaftwidth', 0.7,  'headangle', 30, 'edgeclip', 'on');
uistack(vecs);
m_text(23,79.4,'0.4 m s^{-1}')
title(datestr(datamonth(itime),'yyyy-mm'))

  set(gcf,'PaperPositionMode','auto')
  print(gcf,'-dpng',[datestr(datamonth(itime),'yyyy-mm')],'-r300') 

end