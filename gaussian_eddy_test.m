
clc;clear;
 
% Set coefficients
eta0=0.01; % highest height of eddy
R=300; % eddy's radius
f=0.01; % coriolis coefficient
g=9.81; % gravity constant
 
% Set grid
x=-500:10:500;
y=-500:10:500;
[x y]=meshgrid(x,y);
 
% Set sea surface height
eta=eta0*exp(-(x.^2+y.^2)/R^2);
 
% Calcualte geostrophic velocity
v=(g*eta0)/f*(-2*x/R^2).*exp(-(x.^2+y.^2)/R^2);
u=-(g*eta0)/f*(-2*y/R^2).*exp(-(x.^2+y.^2)/R^2);


% Visualization
figure(1)
clf
set(gcf,'color','w')
subplot(2,1,1)
surf(x,y,eta)
view([-45 80])
title('\eta, sea surface height','fontweight','bold')
subplot(2,1,2)
z0=zeros(size(x));
quiver3(x,y,z0,u,v,z0)
view([-45 80])
axis([-500 500 -500 500 0 0.01])
title('u and v, velocities','fontweight','bold')

%% get observations

n_obs=200;

x_vec=x(:);
y_vec=y(:);
u_vec=u(:);
v_vec=v(:);

ixr=randi(numel(x_vec),[n_obs,1]);
x_obs=x_vec(ixr);
y_obs=y_vec(ixr);
u_obs=u_vec(ixr);
v_obs=v_vec(ixr);

figure(2)
clf
hold on
quiver(x,y,u,v,'k')
quiver(x_obs,y_obs,u_obs,v_obs,'r')

figure(3)
clf
div=divergence(u,v);
surf(div)
shading flat
%% determine centers

% 50 centers
n_centers=50;

ixr=randi(numel(x_obs),[n_centers,1]);
x_c=x_obs(ixr);
y_c=y_obs(ixr);
u_c=u_obs(ixr);
v_c=v_obs(ixr);


figure(2)
clf
hold on
quiver(x,y,u,v,'k')
% quiver(x_obs,y_obs,u_obs,v_obs,'r')
plot(x_c,y_c,'or')


%% train RBFs

data=[u_obs,v_obs]';
data=data(:);

% determine radiuses/differences
clear r rbf

dataindex=1;
for i=1:size(x_obs,1)
    
        %distance to centers
        r=sqrt( (x_obs(i)-x_c).^2 + (y_obs(i)-y_c).^2 )';
        logr=log(r);
        logr(r==0)=0;
             
        % u part
   rbf(dataindex,:)  =   [ 1,0,x_obs(i),y_obs(i),0 , [r.^2.*(12*logr+7) , zeros(size(r))]  + [ -(8*logr+6) .* (x_obs(i)-x_c)' .* (x_obs(i)-x_c)' ,  -(8*logr+6) .* (x_obs(i)-x_c)' .* (y_obs(i)-y_c)'   ]  ] ;
        %v part
   rbf(dataindex+1,:)= [0,1,-y_obs(i),0,x_obs(i) , [zeros(size(r)), r.^2.*(12*logr+7)]   + [ -(8*logr+6) .* (y_obs(i)-y_c)' .* (x_obs(i)-x_c)' ,  -(8*logr+6) .* (y_obs(i)-y_c)' .* (y_obs(i)-y_c)'   ]  ] ;
    dataindex=dataindex+2    
end


params=data'/rbf'  ;%use mrdivide to solve system of equations. For large systems it may 

%% evaluate rbfs


% eval rbfs for grid

[xgrid,ygrid]=meshgrid(-500:50:500,-500:50:500);
xgridv=xgrid(:);
ygridv=ygrid(:);

clear r rbf
dataindex=1;
for i=1:size(xgridv,1)
    
        r=sqrt( (xgridv(i)-x_c).^2 + (ygridv(i)-y_c).^2 )';
        logr=log(r);
        logr(r==0)=0;
            % u part
   rbf(dataindex,:)  =   [ 1,0,xgridv(i),ygridv(i),0 , [r.^2.*(12*logr+7) , zeros(size(r))]  + [ -(8*logr+6) .* (xgridv(i)-x_c)' .* (xgridv(i)-x_c)' ,  -(8*logr+6) .* (xgridv(i)-x_c)' .* (ygridv(i)-y_c)'   ]  ] ;
        %v part
   rbf(dataindex+1,:)= [0,1,-ygridv(i),0,xgridv(i) , [zeros(size(r)), r.^2.*(12*logr+7)]   + [ -(8*logr+6) .* (ygridv(i)-y_c)' .* (xgridv(i)-x_c)' ,  -(8*logr+6) .* (ygridv(i)-y_c)' .* (ygridv(i)-y_c)'   ]  ] ;
    dataindex=dataindex+2   ;
end

estimate=params*rbf';%apply weights
estimate=reshape(estimate,[2,size(xgridv,1)]);
ugrid=reshape( estimate(1,:) , size(xgrid) );
vgrid=reshape( estimate(2,:) , size(xgrid) );

% a=[1,1,1;2,2,2]
% a=a(:)
% reshape(a,[2,3])

figure(2)
clf
hold on
% quiver(x,y,u,v,'k')
 quiver(x_obs,y_obs,u_obs,v_obs,'r')
plot(x_c,y_c,'or')
quiver(xgrid,ygrid,ugrid,vgrid,'k')

figure(3)
clf
div=divergence(ugrid,vgrid);
surf(div)

