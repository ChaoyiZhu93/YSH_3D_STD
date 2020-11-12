clear
close all
%% Input Parameters
% lattice parameter (cubic, nm)
a=0.3615;
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
z_step=5;
% maximum depth to consider (unit: nm)
max_z=100;
%%

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

new_x_grid=zeros(size(y_grid));
new_y_grid=zeros(size(y_grid));
% the depth grid
depth_z=-[0:z_step:max_z];

%% edge type dislocations (12)
burger_slip_edge=0.5*a*[-1 0 1;-1 1 0;0 -1 1;0 1 1;-1 1 0;1 0 1;
    1 0 1;1 1 0;0 -1 1;0 1 1;1 1 0;-1 0 1].';

plane_slip_edge=[1 1 1;1 1 1; 1 1 1;-1 -1 1;-1 -1 1;-1 -1 1;
    -1 1 1;-1 1 1;-1 1 1;1 -1 1;1 -1 1;1 -1 1].';

line_slip_screw=[1 1 0;1 0 1;0 1 1;-1 1 0;1 0 -1;0 -1 1];
%%
e1=[1,0,0];
e2=[0,1,0];
e3=[0,0,1];

for i=1
    % dislocation line of edge dislocation
    line_slip_edge(:,i)=cross(plane_slip_edge(:,i),burger_slip_edge(:,i));
    % line direction projected on the z plane (normalzed)
    line_slip_edge_zp(:,i)=[line_slip_edge(1,i),line_slip_edge(2,i),0]./norm([line_slip_edge(1,i),line_slip_edge(2,i),0]);
    % normal to the plane containing z and dislocation line (yz plane)
    x_slip_edge_zp(:,i)=cross(line_slip_edge_zp(:,i),e3)/norm(cross(line_slip_edge_zp(:,i),e3));
    
    p_1=x_slip_edge_zp(1:2,i)';
    p_2=line_slip_edge_zp(1:2,i)';
    p1=[1,0];
    p2=[0,1];
    R(:,:,i)=[dot(p1,p_1),dot(p1,p_2);
        dot(p2,p_1),dot(p2,p_2)];
    
    for j=1:numel(x_grid)
        temp_position=R(:,:,i)'*[x_grid(j);y_grid(j)];
        new_x_grid(j)=temp_position(1);
        new_y_grid(j)=temp_position(2);
    end
    
    
    
    e_2=line_slip_edge_zp(:,i);
    e_1=x_slip_edge_zp(:,i);
    e_3=e3;
    Q(:,:,i)=[dot(e_1,e1),dot(e_1,e2),dot(e_1,e3);
        dot(e_2,e1),dot(e_2,e2),dot(e_2,e3);
        dot(e_3,e1),dot(e_3,e2),dot(e_3,e3)];
    edge_yz(:,i)=burger_slip_edge(:,i)-dot(burger_slip_edge(:,i),x_slip_edge_zp(:,i))*x_slip_edge_zp(:,i);
    edge_xz(:,i)=dot(burger_slip_edge(:,i),x_slip_edge_zp(:,i))*x_slip_edge_zp(:,i);
    
 
    if line_slip_edge(3,i)>0
        inclination(i)=-acosd(dot(line_slip_edge(:,i),e3)/norm(line_slip_edge(:,i)));
    elseif line_slip_edge(3,i)<0
        inclination(i)=180-acosd(dot(line_slip_edge(:,i),e3)/norm(line_slip_edge(:,i)));
    end
    
    
    
    [beta_yz,~]=YSH(Q(:,:,i)*edge_yz(:,i),deg2rad(inclination(i)),G,v,new_x_grid,new_y_grid,depth_z,'Edge-yz');
    [beta_xz,~]=YSH(Q(:,:,i)*edge_xz(:,i),deg2rad(inclination(i)),G,v,new_x_grid,new_y_grid,depth_z,'Edge-xz');
    
    for k=1
        
        beta=beta_xz(:,:,:,:,k)+beta_yz(:,:,:,:,k);
        beta=replaceNaN(beta);
        
        ca=[-0.001,0.001];
        
        strainrot=YSHStrainRotation(beta);
        plotStrainRotation(strainrot,ca);
        
        stress=YSHStress(C,strainrot);
        
        ca=[-0.08,0.08];
        plotStress(stress,ca)
        
    end
    clear beta beat_yz beta_xz
    

end

function [beta]=replaceNaN(beta)

for i=1:size(beta,1)
    for j=1:size(beta,2)
        for k=1:size(beta,3)
            for l=1:size(beta,4)
                if isnan(beta(i,j,k,l))
                    if l<size(beta,4)
                        beta(i,j,k,l)=beta(i,j,k,l+1);
                    elseif l==size(beta,4)
                        beta(i,j,k,l)=beta(i,j,k,l-1);
                    end
                end
            end
        end
    end
end
end

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

function DefTensorVolume(name,np,step,eu,pc,defTensor)
warning('off')
h5create(name,'/DeformationFieldInfo/npix',[1,1],'Datatype','int32');
h5create(name,'/DeformationFieldInfo/npiy',[1,1],'Datatype','int32');
h5create(name,'/DeformationFieldInfo/npiz',[1,1],'Datatype','int32');
h5create(name,'/DeformationFieldInfo/stepx',[1,1],'Datatype','double');
h5create(name,'/DeformationFieldInfo/stepy',[1,1],'Datatype','double');
h5create(name,'/DeformationFieldInfo/stepz',[1,1],'Datatype','double');
h5create(name,'/DeformationFieldInfo/eu',[1,3],'Datatype','double');
h5create(name,'/DeformationFieldInfo/pcx',[np(1),np(2)],'Datatype','double');
h5create(name,'/DeformationFieldInfo/pcy',[np(1),np(2)],'Datatype','double');
h5create(name,'/DeformationFieldInfo/L',[1,1],'Datatype','double');
h5create(name,'/DeformationFieldInfo/deftensor',[9,np(3),np(1),np(2)],'Datatype','double');
h5write(name, '/DeformationFieldInfo/npix', np(1));
h5write(name, '/DeformationFieldInfo/npiy', np(2));
h5write(name, '/DeformationFieldInfo/npiz', np(3));
h5write(name, '/DeformationFieldInfo/stepx', step(1));
h5write(name, '/DeformationFieldInfo/stepy', step(2));
h5write(name, '/DeformationFieldInfo/stepz', step(3));
h5write(name, '/DeformationFieldInfo/eu', eu);
h5write(name, '/DeformationFieldInfo/pcx', pc(:,:,1));
h5write(name, '/DeformationFieldInfo/pcy', pc(:,:,2));
h5write(name, '/DeformationFieldInfo/L', pc(1,1,3));
h5write(name, '/DeformationFieldInfo/deftensor', defTensor);
end