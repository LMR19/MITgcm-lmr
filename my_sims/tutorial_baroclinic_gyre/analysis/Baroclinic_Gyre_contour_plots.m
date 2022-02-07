% tutorial baroclinic gyre  
% '/Users/Lina/Data/MITgcm_project/MITgcm-lmr/my_sims/tutorial_baroclinic_gyre/'

addpath( '/Users/Lina/Data/MITgcm_project/MITgcm-lmr/utils/matlab')

dat2d=rdmnc('outs_2D.0000000000*',8640); % 8640 for run < 06 233280 for run06
dat3d=rdmnc('outs_3D.0000000000*',8640);
% dat=rdmnc('surfDiag.0000000000*',8640);
% dat=rdmnc('dynDiag.0000000000*',8640);

figure

X=dat2d.X;
Y=dat2d.Y;
Eta=dat2d.ETAN;
Xp1=dat3d.Xp1;
Yp1=dat3d.Yp1;
U=dat3d.UVEL;
V=dat3d.VVEL;
W=dat3d.WVEL;

subplot(221)
pcolor(X,Y,Eta(:,:)')
hold on
title('Eta [m]')
xlabel('X')
ylabel('Y')
shading flat
colorbar
ylim([10 80])
xlim([0 61])
set(gca,'fontsize',16,'fontname','Times New Roman');
caxis([-round(max(max(Eta)),1) round(max(max(Eta)),1)])
T=text(-29.5000, 97.0153,'Matlab code: Baroclinic\_Gyre\_contour\_plots.m')

subplot(222)
pcolor(Xp1,Y,U(:,:,1)')
hold on
title('U [m/s]')
xlabel('X')
ylabel('Yp1')
shading flat
colorbar
ylim([10 80])
xlim([0 61])
set(gca,'fontsize',16,'fontname','Times New Roman');
caxis([-round(max(max(max(U))),1) round(max(max(max(U))),1)])

subplot(223)
pcolor(X,Yp1,V(:,:,1)')
hold on
title('V [m/s]')
xlabel('Xp1')
ylabel('Y')
shading flat
colorbar
ylim([10 80])
xlim([0 61])
set(gca,'fontsize',16,'fontname','Times New Roman');
caxis([-round(max(max(max(V))),1) round(max(max(max(V))),1)])

subplot(224)
pcolor(X,Y,W(:,:,1)')
hold on
title('W [m/s]')
xlabel('X')
ylabel('Y')
shading flat
colorbar
ylim([10 80])
xlim([0 61])
set(gca,'fontsize',16,'fontname','Times New Roman');
caxis([-round(max(max(max(W))),1) round(max(max(max(W))),1)])

f=pwd;
% saveas(gcf,['Eta_UVW_' f(end-4:end)],'jpg')
%%
dat2d_full=rdmnc('outs_2D.0000000000*'); % 8640 for run < 06 233280 for run06
dat3d_full=rdmnc('outs_3D.0000000000*');

T=dat2d_full.T;
X=dat2d_full.X;
Y=dat2d_full.Y;
Eta=dat2d_full.ETAN;
Xp1=dat3d_full.Xp1;
Yp1=dat3d_full.Yp1;
U=dat3d_full.UVEL;
V=dat3d_full.VVEL;
W=dat3d_full.WVEL;

n=10;

for ts=1:length(T)
    
    Eta_hm(ts,:)=Eta(:,n,1,ts);
    U_hm(ts,:)=U(:,n,1,ts);
    V_hm(ts,:)=V(:,n,1,ts);
    W_hm(ts,:)=W(:,n,1,ts);
    
end


figure 

% sgtitle('Subplot Grid Title','location','northeast')

subplot(221)
pcolor(T/(3600*24),Y,Eta_hm')
hold on
title('Eta [m]')
ylabel('Y')
xlabel('Time [days]')
shading flat
colorbar
xlim([0 T(end)/(3600*24)])
ylim([10 80])
set(gca,'fontsize',16,'fontname','Times New Roman');
caxis([-round(max(max(Eta_hm)),1) round(max(max(Eta_hm)),1)])
text(-599, 97.0153,'Matlab code: Baroclinic\_Gyre\_contour\_plots.m')

subplot(222)
pcolor(T/(3600*24),Yp1,U_hm')
hold on
title('U [m/s]')
ylabel('Yp1')
xlabel('Time [days]')
shading flat
colorbar
xlim([0 T(end)/(3600*24)])
ylim([10 80])
set(gca,'fontsize',16,'fontname','Times New Roman');
caxis([-round(max(max(max(U_hm))),1) round(max(max(max(U_hm))),1)])

subplot(223)
pcolor(T/(3600*24),Y,V_hm')
hold on
title('V [m/s]')
ylabel('Y')
xlabel('Time [days]')
shading flat
colorbar
xlim([0 T(end)/(3600*24)])
ylim([10 80])
set(gca,'fontsize',16,'fontname','Times New Roman');
caxis([-round(max(max(max(V_hm))),1) round(max(max(max(V_hm))),1)])

subplot(224)
pcolor(T/(3600*24),Y,W_hm')
hold on
title('W [m/s]')
ylabel('Y')
xlabel('Time [days]')
shading flat
colorbar
xlim([0 T(end)/(3600*24)])
ylim([10 80])
set(gca,'fontsize',16,'fontname','Times New Roman');
caxis([-round(max(max(max(W_hm))),1) round(max(max(max(W_hm))),1)])

f=pwd;
saveas(gcf,['Eta_UVW_hm_' f(end-4:end)],'jpg')


%% 

% reading in data
% ncdisp('grid.t001.nc')

% Z=ncread('grid.t001.nc','Z');               % 'vertical coordinate of cell center'
% RC=ncread('grid.t001.nc','RC');             % 'R coordinate of cell center'
% Zp1=ncread('grid.t001.nc','Zp1');           % 'vertical coordinate of cell interface'
% RF=ncread('grid.t001.nc','RF');             % 'R coordinate of cell interface'
% Zu=ncread('grid.t001.nc','Zu');             % 'vertical coordinate of lower cell interface'
% RU=ncread('grid.t001.nc','RU');             % 'R coordinate of upper interface'
% Zl=ncread('grid.t001.nc','Zl');             % 'vertical coordinate of upper cell interface'
% RL=ncread('grid.t001.nc','RF');             % 'R coordinate of lower interface'
% drC=ncread('grid.t001.nc','drC');           % 'r cell center separation'
% drF=ncread('grid.t001.nc','drF');           % 'r cell face separation'
% X=ncread('grid.t001.nc','X');               % 'longitude of cell center'
% Y=ncread('grid.t001.nc','Y');               % 'latitude of cell center'
% XC=ncread('grid.t001.nc','XC');             % 'X coordinate of cell center (T-P point)'
% YC=ncread('grid.t001.nc','YC');             % 'Y coordinate of cell center (T-P point)'
% Xp1=ncread('grid.t001.nc','Xp1');           % 'longitude of cell corner'
% Yp1=ncread('grid.t001.nc','Yp1');           % 'latitude of cell corner'
% XG=ncread('grid.t001.nc','XG');             % 'X coordinate of cell corner (Vorticity point)'
% YG=ncread('grid.t001.nc','YG');             % 'Y coordinate of cell corner (Vorticity point)'
% dxC=ncread('grid.t001.nc','dxC');           % 'x cell center separation'
% dyC=ncread('grid.t001.nc','dyC');           % 'y cell center separation'
% dxF=ncread('grid.t001.nc','dxF');           % 'x cell face separation'
% dyF=ncread('grid.t001.nc','dyF');           % 'y cell face separation'
% dxG=ncread('grid.t001.nc','dxG');           % 'x cell corner separation'
% dyG=ncread('grid.t001.nc','dyG');           % 'y cell corner separation'
% dxV=ncread('grid.t001.nc','dxV');           % 'x v-velocity separation'
% dyU=ncread('grid.t001.nc','dyU');           % 'y u-velocity separation'
% rA=ncread('grid.t001.nc','rA');             % 'r-face area at cell center'
% rAw=ncread('grid.t001.nc','rAw');           % 'r-face area at U point'
% rAs=ncread('grid.t001.nc','rAs');           % 'r-face area at V point'
% rAz=ncread('grid.t001.nc','rAz');           % 'r-face area at cell corner'
% fCori=ncread('grid.t001.nc','fCori');       % 'Coriolis f at cell center'
% fCoriG=ncread('grid.t001.nc','fCoriG');     % 'Coriolis f at cell corner'
% R_low=ncread('grid.t001.nc','R_low');       % 'base of fluid in r-units'
% Ro_surf=ncread('grid.t001.nc','Ro_surf');   % 'surface reference (at rest) position'
% Depth=ncread('grid.t001.nc','Depth');       % 'fluid thickness in r coordinates (at rest)'
% HFacC=ncread('grid.t001.nc','HFacC');       % 'vertical fraction of open cell at cell center'
% HFacW=ncread('grid.t001.nc','HFacW');       % 'vertical fraction of open cell at West face'
% HFacS=ncread('grid.t001.nc','HFacS');       % 'vertical fraction of open cell at South face'

% T=ncread('phiHyd.0000000000.t001.nc','T');       % 'model_time'  units     = 's'
% iter=ncread('phiHyd.0000000000.t001.nc','iter');       % 'iteration_count'
% X=ncread('phiHyd.0000000000.t001.nc','X');       % 'longitude of cell center'  units     = 'degrees_east'
% Y=ncread('phiHyd.0000000000.t001.nc','Y');       % 'latitude of cell center' units     = 'degrees_north'
% Z=ncread('phiHyd.0000000000.t001.nc','Z');       % 'vertical coordinate of cell center' units     = 'meters' positive  = 'up'
% phiHyd=ncread('phiHyd.0000000000.t001.nc','phiHyd');       % ????
% phiHydLow=ncread('phiHyd.0000000000.t001.nc','phiHydLow');       % ????

% T=ncread('state.0000000000.t001.nc','T');       % 'model_time'  units     = 's'
% iter=ncread('state.0000000000.t001.nc','iter');       % 'iteration_count'
% Xp1=ncread('state.0000000000.t001.nc','Xp1');           % 'longitude of cell corner'
% Yp1=ncread('state.0000000000.t001.nc','Yp1');           % 'latitude of cell corner'
% X=ncread('state.0000000000.t001.nc','X');               % 'longitude of cell center'
% Y=ncread('state.0000000000.t001.nc','Y');       % 'latitude of cell center' units     = 'degrees_north'
% Z=ncread('state.0000000000.t001.nc','Z');       % 'vertical coordinate of cell center' units     = 'meters' positive  = 'up'
% U=ncread('state.0000000000.t001.nc','U');       % 'm/s' 'XU YU RC iter'
% V=ncread('state.0000000000.t001.nc','V');       % 'm/s' 'XV YV RC iter'
% W=ncread('state.0000000000.t001.nc','W');       % 'm/s''XC YC RC iter'
% Temp=ncread('state.0000000000.t001.nc','Temp');       % 'degC' 'potential_temperature'  'XC YC RC iter'
% S=ncread('state.0000000000.t001.nc','S');       % 'salinity' 'XC YC RC iter'
% Eta=ncread('state.0000000000.t001.nc','Eta');   % 'free-surface_r-anomaly' 'm' 'XC YC RC iter'
% Zl=ncread('state.0000000000.t001.nc','Zl');       % 'vertical coordinate of upper cell interface' 'm' 'up'                       


% file=dir('state*.nc');
% 
% figure
% 
% for fl=1:length(file);
% 
%     file_name=file(fl).name;
% 
%     T=ncread(file_name,'T');       % 'model_time'  units     = 's'
%     iter=ncread(file_name,'iter');       % 'iteration_count'
%     Xp1=ncread(file_name,'Xp1');           % 'longitude of cell corner'
%     Yp1=ncread(file_name,'Yp1');           % 'latitude of cell corner'
%     X=ncread(file_name,'X');               % 'longitude of cell center'
%     Y=ncread(file_name,'Y');       % 'latitude of cell center' units     = 'degrees_north'
%     Z=ncread(file_name,'Z');       % 'vertical coordinate of cell center' units     = 'meters' positive  = 'up'
%     U=ncread(file_name,'U');       % 'm/s' 'XU YU RC iter'
%     V=ncread(file_name,'V');       % 'm/s' 'XV YV RC iter'
%     W=ncread(file_name,'W');       % 'm/s''XC YC RC iter'
%     Temp=ncread(file_name,'Temp');       % 'degC' 'potential_temperature'  'XC YC RC iter'
%     S=ncread(file_name,'S');       % 'salinity' 'XC YC RC iter'
%     Eta=ncread(file_name,'Eta');   % 'free-surface_r-anomaly' 'm' 'XC YC RC iter'
%     Zl=ncread(file_name,'Zl');       % 'vertical coordinate of upper cell interface' 'm' 'up'                       
% 
%     t=2
% 
%     subplot(221)
%     pcolor(X,Y,Eta(:,:,t)')
%     hold on
%     title('Eta [m]')
%     xlabel('X')
%     ylabel('Y')
%     shading flat
%     colorbar
%     ylim([10 80])
%     xlim([0 61])
% 
%     subplot(222)
%     pcolor(Xp1,Y,U(:,:,1,t)')
%     hold on
%     title('U [m/s]')
%     xlabel('X')
%     ylabel('Yp1')
%     shading flat
%     colorbar
%     ylim([10 80])
%     xlim([0 61])
% 
%     subplot(223)
%     pcolor(X,Yp1,V(:,:,1,t)')
%     hold on
%     title('V [m/s]')
%     xlabel('Xp1')
%     ylabel('Y')
%     shading flat
%     colorbar
%     ylim([10 80])
%     xlim([0 61])
% 
%     subplot(224)
%     pcolor(X,Y,W(:,:,1,t)')
%     hold on
%     title('W [m/s]')
%     xlabel('X')
%     ylabel('Y')
%     shading flat
%     colorbar
%     ylim([10 80])
%     xlim([0 61])
% 
% end