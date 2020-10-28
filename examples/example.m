clear
close all
%% Input Parameters
% Inclination angle (unit:degree)
a=0; 
% Burgers vector (unit: nm)
b=[0.2556,0,0]; 
% Poisson's ratio
v=0.3;
% Grid step size (unit: nm)
grid_step=3; % nm
% Number of grid points (even number to avoid singularity) (unit:)
% Total size of the grid=grid_step*grid_size
grid_size=334;
% Young's Modulus (unit:GPa)
E=200;
% step size in depth (unit: nm)
z_step=1;
% maximum depth to consider (unit: nm)
max_z=10;
%%
% convert inclination angle it to radians
a=deg2rad(a);
% stiffness matrix
la=E*v/((1+v)*(1-2*v));
G=E/(2*(1+v));
C=[la+2*G,     la,         la,  0,  0,  0;
    la,     la+2*G,         la,  0,  0,  0;
    la,           la,    la+2*G,  0,  0,  0;
    0,           0,          0, G,  0,  0;
    0,           0,          0,  0, G,  0;
    0,           0,          0,  0,  0, G];
% grid with origin at the (0,0) 
bound=(grid_size*grid_step-grid_step)/2;
x0=[-bound:grid_step:bound];
y0=[-bound:grid_step:bound];
% generate the x y grid mesh
[x_grid,y_grid]=meshgrid(x0,y0);
y_grid=flipud(y_grid);
% the depth grid
depth_z=-[0:z_step:10];

%% plots for strain and rotation field (z=0)
ca=[-0.001,0.001];
% pure edge dislocation with Burgers vector in x direction
[beta,defTensor]=YSH(b,a,G,v,x_grid,y_grid,depth_z,'Edge-xz');
% surface rotation and strain at depth = 0 [beta(:,:,:,:,1)]
strainrot=YSHStrainRotation(beta(:,:,:,:,1));
plotStrainRotation(strainrot,ca);
%% plots for stress field (z=0)
stress=YSHStress(C,strainrot);
ca=[-0.08,0.08];
plotStress(stress,ca)
%% functions

function [strainrot]=YSHStrainRotation(beta)
strainrot=zeros(size(beta));
strainrot(1,1,:,:)=reshape(beta(1,1,:,:),[size(beta,3),size(beta,4)]);
strainrot(2,2,:,:)=reshape(beta(2,2,:,:),[size(beta,3),size(beta,4)]);
strainrot(3,3,:,:)=reshape(beta(3,3,:,:),[size(beta,3),size(beta,4)]);
strainrot(1,2,:,:)=0.5*(reshape(beta(1,2,:,:),[size(beta,3),size(beta,4)])+reshape(beta(2,1,:,:),[size(beta,3),size(beta,4)]));
strainrot(1,3,:,:)=0.5*(reshape(beta(1,3,:,:),[size(beta,3),size(beta,4)])+reshape(beta(3,1,:,:),[size(beta,3),size(beta,4)]));
strainrot(2,3,:,:)=0.5*(reshape(beta(2,3,:,:),[size(beta,3),size(beta,4)])+reshape(beta(3,2,:,:),[size(beta,3),size(beta,4)]));
strainrot(2,1,:,:)=0.5*(reshape(beta(1,2,:,:),[size(beta,3),size(beta,4)])-reshape(beta(2,1,:,:),[size(beta,3),size(beta,4)]));
strainrot(3,2,:,:)=0.5*(reshape(beta(2,3,:,:),[size(beta,3),size(beta,4)])-reshape(beta(3,2,:,:),[size(beta,3),size(beta,4)]));
strainrot(3,1,:,:)=0.5*(reshape(beta(3,1,:,:),[size(beta,3),size(beta,4)])-reshape(beta(1,3,:,:),[size(beta,3),size(beta,4)]));

end

function plotStrainRotation(beta,ca,print_title)
figure
subplot(3,3,1);
imagesc(reshape(beta(1,1,:,:),[size(beta,3),size(beta,4)]));
colormap jet;caxis(ca);
axis image; axis ij; axis off;
title('\epsilon_{11}','Fontsize',14);

subplot(3,3,2);
imagesc(reshape(beta(1,2,:,:),[size(beta,3),size(beta,4)]));colormap jet;caxis(ca);
axis image; axis ij;  axis off;
title('\epsilon_{12}','Fontsize',14);

subplot(3,3,3);
imagesc(reshape(beta(1,3,:,:),[size(beta,3),size(beta,4)]));colormap jet;caxis(ca);
axis image; axis ij;  axis off;
title('\epsilon_{13}','Fontsize',14);

subplot(3,3,4);
imagesc(reshape(beta(2,1,:,:),[size(beta,3),size(beta,4)]));colormap jet;caxis(ca);
axis image; axis ij; axis off;
title('\omega_{12}','Fontsize',14);

subplot(3,3,5);
imagesc(reshape(beta(2,2,:,:),[size(beta,3),size(beta,4)]));colormap jet;caxis(ca);
axis image; axis ij; axis off;
title('\epsilon_{22}','Fontsize',14);

subplot(3,3,6);
imagesc(reshape(beta(2,3,:,:),[size(beta,3),size(beta,4)]));colormap jet;caxis(ca);
axis image; axis ij; axis off;
title('\epsilon_{23}','Fontsize',14);

subplot(3,3,7);
imagesc(reshape(beta(3,1,:,:),[size(beta,3),size(beta,4)]));colormap jet;caxis(ca);
axis image; axis ij;  axis off;
title('\omega_{31}','Fontsize',14);

subplot(3,3,8);
imagesc(reshape(beta(3,2,:,:),[size(beta,3),size(beta,4)]));colormap jet;caxis(ca);
axis image; axis ij;  axis off;
title('\omega_{23}','Fontsize',14);

subplot(3,3,9);
imagesc(reshape(beta(3,3,:,:),[size(beta,3),size(beta,4)]));colormap jet;caxis(ca);
axis image; axis ij;  axis off;
title('\epsilon_{33}','Fontsize',14);
set(gcf,'position',[375 87.4000 910.4000 678]);

hp4 = get(subplot(3,3,9),'Position');
colorbar('Position', [hp4(1)+hp4(3)-0.85 hp4(2)-0.01+0.14  0.04 hp4(2)+hp4(3)*2.1])
figure_handle = gcf;
cbar_handle = findobj(figure_handle,'tag','Colorbar');
set(cbar_handle, 'YAxisLocation','right','YTick',ca, 'FontSize', 10)
switch nargin
    case 3
        print(print_title,'-dtiff','-r300');
end
end
function [stress]=YSHStress(C,strain)
stress=zeros(size(strain));
for i=1:size(strain,3)
    for j=1:size(strain,4)
        tempstrain=[strain(1,1,i,j);strain(2,2,i,j); strain(3,3,i,j);2*strain(2,3,i,j);2*strain(1,3,i,j); 2*strain(1,2,i,j)];
        tempstress=C*tempstrain;
        stress(1,1,i,j)=tempstress(1);
        stress(2,2,i,j)=tempstress(2);
        stress(3,3,i,j)=tempstress(3);
        stress(2,3,i,j)=tempstress(4);
        stress(1,3,i,j)=tempstress(5);
        stress(1,2,i,j)=tempstress(6);
        clear tempstress tempstrain
    end
end
end

function plotStress(beta,ca,print_title)
figure
subplot(3,3,1);
imagesc(reshape(beta(1,1,:,:),[size(beta,3),size(beta,4)]));
colormap jet;caxis(ca);
axis image; axis ij; axis off;
title('\sigma_{11}','Fontsize',14);

subplot(3,3,2);
imagesc(reshape(beta(1,2,:,:),[size(beta,3),size(beta,4)]));colormap jet;caxis(ca);
axis image; axis ij;  axis off;
title('\sigma_{12}','Fontsize',14);

subplot(3,3,3);
imagesc(reshape(beta(1,3,:,:),[size(beta,3),size(beta,4)]));colormap jet;caxis(ca);
axis image; axis ij;  axis off;
title('\sigma_{13}','Fontsize',14);



subplot(3,3,5);
imagesc(reshape(beta(2,2,:,:),[size(beta,3),size(beta,4)]));colormap jet;caxis(ca);
axis image; axis ij; axis off;
title('\sigma_{22}','Fontsize',14);

subplot(3,3,6);
imagesc(reshape(beta(2,3,:,:),[size(beta,3),size(beta,4)]));colormap jet;caxis(ca);
axis image; axis ij; axis off;
title('\sigma_{23}','Fontsize',14);

subplot(3,3,9);
imagesc(reshape(beta(3,3,:,:),[size(beta,3),size(beta,4)]));colormap jet;caxis(ca);
axis image; axis ij;  axis off;
title('\sigma_{33}','Fontsize',14);
set(gcf,'position',[375 87.4000 910.4000 678]);

hp4 = get(subplot(3,3,9),'Position');
colorbar('Position', [hp4(1)+hp4(3)-0.85 hp4(2)-0.01+0.14  0.04 hp4(2)+hp4(3)*2.1])
figure_handle = gcf;
cbar_handle = findobj(figure_handle,'tag','Colorbar');
set(cbar_handle, 'YAxisLocation','right','YTick',ca, 'FontSize', 10)
switch nargin
    case 3
        print(print_title,'-dtiff','-r300');
end
end