
function [beta]=YSH(b,angle,mu,p,x_grid,y_grid,z_grid,type)
% Compute the distortion field around surface threading dislocations
% Reference: Shaibani, S.J. and Hazzledine, P.M., 1981. 

% Input parameters:
% b: Burgers vector (values for substitution)
% angle: inclination angle relative to -z (values for substitution)
% p: Poisson's ratio (values for substitution)
% mu: shear modulus (values for substitution)
% x_grid,y_grid,z_grid: 3D grid for the dislocation
% type: type of YSH dislocation (screw, edge_yz, edge_xz)

% Output parameter: 
% beta: numerical values of the distortion field (gradient in the
% displacement)

% define symbolic variables
% x y z coordinates, inclination angle
% Burgers vector and the Poissons ratio
syms x y z a 
syms bx by bz v

% coordinate system: x, eta, zeta (m)
eta=y*cos(a)-z*sin(a);
zeta=z*cos(a)+y*sin(a);

% coordinate system: x, zeta, eta (m)
eta_=-y*cos(a)-z*sin(a);
zeta_=-z*cos(a)+y*sin(a);

% sqrt of sum of squared (m)
R=sqrt(x^2+y^2+z^2);
R1=sqrt(x^2+eta^2+zeta^2);
R2=sqrt(x^2+eta_^2+zeta_^2);

% Burgers vector of screw (along zeta_) (m)
b_s=bz*cos(a)-by*sin(a);
% Burgers vector of edge (along eta_) (m)
b_e=by*cos(a)+bz*sin(a);

% parameters (rad)
omega=atan(y/x)-atan(eta/x)+atan(x*R*sin(a)/(y*eta+x^2*cos(a)));
omega_=atan(y/x)-atan(eta_/x)+atan(x*R*sin(a)/(y*eta_-x^2*cos(a)));

% (m)
A=R-z;B=R-zeta;B_=R-zeta_;
depth=z_grid*ones(size(x_grid));

% screw parameters
switch type
    case 'Screw' % y-z plane
        m_s=x*sin(2*a)/(R*B); % (m^-1)
        S=(mu*b_s/(2*pi))^(-1);% (m^-1)
        
        % inclined screw dislocation displacement (u_x, u_y , u_z)
        u_x=(x*m_s+(2*eta*(cos(a))^2)/B+2*(1-2*v)*cot(a)*(-1+cos(a)+cos(a)*log(A)-...
            y*sin(a)/A-log(B))-sin(2*a))/(2*S*mu);
        
        u_y=(y*m_s-(2*x*cos(a))/B-sin(a)*(omega_-omega)+2*(1-2*v)*cot(a)*(x*sin(a)/A-...
            omega*cos(a)))/(2*S*mu);
        
        u_z=(z*m_s+cos(a)*(omega_-omega)-2*(1-2*v)*omega*cos(a))/(2*S*mu);
        
        % displacement field with inclination angle substitution
        [u_x,u_y,u_z]=DisplacementSym(u_x,u_y,u_z,angle);
        % distortion field
        Beta=DistortionSym(b,p,u_x,u_y,u_z);
        
        beta=YSHDistortion(Beta,x_grid,y_grid,depth);
        

    case 'Edge-yz'
        % inclined edge dislocation displacement (b in the y-z plane)
        % parameters
        q_e=x*((1/B_)-(1/B)+2*z*cos(a)/B^2);
        m_e=-q_e/R-4*(1-v)*x*(cos(a))^2/(R*B);
        D_e=(mu*b_e/(4*pi*(1-v)))^(-1);%different to defectmodule
        lamda=(1-2*v)*log(B_/B);
        theta=2*(1-v)*(omega_-omega);
        k=4*(1-v)*(1-2*v)*(cot(a))^2;
        %displacement
        u_x=(x*m_e+lamda+2*cos(a)/B*(z+2*(1-v)*eta*sin(a))-4*(1-v)*(sin(a))^2+...
            k*(1-cos(a)-cos(a)*log(A)+y*sin(a)/A+log(B)))/(D_e*2*mu);
        
        u_y=(y*m_e+q_e*sin(a)+theta*cos(a)+k*(-x*sin(a)/A+omega*cos(a)))/(D_e*2*mu);
        
        u_z=(z*m_e+q_e*cos(a)+theta*sin(a)-2*x*cos(a)*(1/B_+(1-2*v)/B)+k*omega*sin(a))/(D_e*2*mu);
        
        [u_x,u_y,u_z]=DisplacementSym(u_x,u_y,u_z,angle);
        Beta=DistortionSym(b,p,u_x,u_y,u_z);
        beta=YSHDistortion(Beta,x_grid,y_grid,depth);
        
        
    case 'Edge-xz' % x-z plane
        % inclined edge dislocation displacement (b in the x-z plane)
        
        % parameters
        q_x=eta_/B_-eta/B-2*z*eta*cos(a)/(B^2);
        m_x=-q_x/R+2*(1-2*v)*y*cos(a)/(R*B);
        D_x=(mu*bx/(4*pi*(1-v)))^(-1);
        theta=2*(1-v)*(omega_-omega);
        k=4*(1-v)*(1-2*v)*(cot(a))^2;
        lamda=(1-2*v)*log(B_/B);
        
        u_x=(x*m_x+theta+k*(x*tan(a)/A-omega))/(D_x*2*mu);
        u_y=(y*m_x+q_x*sin(a)-lamda*cos(a)-2*cos(a)/B*(z*cos(a)+(1-2*v)*y*sin(a))+...
            k*(-1+cos(a)-log(A)+y*tan(a)/A+cos(a)*log(B)))/(D_x*2*mu);
        u_z=(z*m_x+q_x*cos(a)-lamda*sin(a)-2*eta_*cos(a)/B_+4*cos(a)/B*((1-v)*y*cos(a)...
            -(1-2*v)*z*sin(a))+k*tan(a)*(cos(a)-log(A)+cos(a)*log(B))+4*(1-v)*cos(a)*cot(a))/(D_x*2*mu);

        [u_x,u_y,u_z]=DisplacementSym(u_x,u_y,u_z,angle);
        Beta=DistortionSym(b,p,u_x,u_y,u_z);
        beta=YSHDistortion(Beta,x_grid,y_grid,depth);
        
end
end
function [u_x,u_y,u_z]=DisplacementSym(u_x,u_y,u_z,angle)
% DisplacementSym deals with inclination angle substitution and the special 
% case when inclination angle is zero.
%tic
syms x y z a
syms bx by bz v
if angle==0 
    u_x=limit(u_x,a,0);
    u_y=limit(u_y,a,0);
    u_z=limit(u_z,a,0); 
else
    u_x=subs(u_x,a,angle);
    u_y=subs(u_y,a,angle);
    u_z=subs(u_z,a,angle);
end
end

function [Beta]=DistortionSym(b,p,u,v1,w)
% DistortionSym obtains the distortion field by taking derivarives of the 
% displacement field and substitute Burgers vector and Poisson's ratio

syms x y z 
syms bx by bz v
% betaij = dui/dxj (deformation gradient)
B(1,1) = diff(u,x);
B(1,2) = diff(u,y);
B(1,3) = diff(u,z);
B(2,1) = diff(v1,x);
B(2,2) = diff(v1,y);
B(2,3) = diff(v1,z);
B(3,1) = diff(w,x);
B(3,2) = diff(w,y);
B(3,3) = diff(w,z);
% substitue values for Burgers vector and Poisson's ratio
Beta=subs(B,{bx,by,bz,v},{b(1),b(2),b(3),p});
end

function [beta]=YSHDistortion(Beta,x_grid,y_grid,depth)
% YSHDistortion substitute the spatial grids and obtain numerical 
% distortion field 
beta=zeros(3,3,size(x_grid,1),size(x_grid,2));
syms x y z
b11(x,y,z)=Beta(1,1);
b12(x,y,z)=Beta(1,2);
b13(x,y,z)=Beta(1,3);
b21(x,y,z)=Beta(2,1);
b22(x,y,z)=Beta(2,2);
b23(x,y,z)=Beta(2,3);
b31(x,y,z)=Beta(3,1);
b32(x,y,z)=Beta(3,2);
b33(x,y,z)=Beta(3,3);

hb11=matlabFunction(b11);
hb12=matlabFunction(b12);
hb13=matlabFunction(b13);
hb21=matlabFunction(b21);
hb22=matlabFunction(b22);
hb23=matlabFunction(b23);
hb31=matlabFunction(b31);
hb32=matlabFunction(b32);
hb33=matlabFunction(b33);

beta(1,1,:,:)=hb11(x_grid,y_grid,depth);
beta(1,2,:,:)=hb12(x_grid,y_grid,depth);
beta(1,3,:,:)=hb13(x_grid,y_grid,depth);
beta(2,1,:,:)=hb21(x_grid,y_grid,depth);
beta(2,2,:,:)=hb22(x_grid,y_grid,depth);
beta(2,3,:,:)=hb23(x_grid,y_grid,depth);
beta(3,1,:,:)=hb31(x_grid,y_grid,depth);
beta(3,2,:,:)=hb32(x_grid,y_grid,depth);
beta(3,3,:,:)=hb33(x_grid,y_grid,depth);
end