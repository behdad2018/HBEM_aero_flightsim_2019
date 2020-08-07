%     Copyright (C) 12/12/2018  Regents of the University of Michigan
%     Aerospace Engineering Department
%     Computational Aeroscience Lab, written by Behdad Davoudi

%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <https://www.gnu.org/licenses/>.

function varargout=try2(T,varargin)
%#codegen
global geometry

if nargin==3
    rpm=varargin{1};
    V_rel_B=varargin{2};
end

% prcoedures:

% the geometry is a cell array, should already be created in the workspace
% and defined as global variable
% find the uniform lamb0
% used linear inflow model for forward model
% calculate thrust, roll and moment coefficiencts
% using blade element theory

%% rotor characteristics

% list={R,nb,A,rho,nr,npsi,M,alp0,aLeq0,th,c,cla,r,psi}

R=geometry(1,1);
nb=geometry(2,1);
A=geometry(3,1);
rho=geometry(4,1);
nr=geometry(5,1);
npsi=geometry(6,1);
M=geometry(7,1);
alp0=geometry(8,1);
aLeq0=geometry(9,1);
th=geometry(10,1:nr);
c=geometry(11,1:nr);
cla=geometry(12,1:nr);
r=geometry(13,1:nr);
psi=geometry(14,:);

sig=nb*mean(c)/(pi*R);        % solidity
%%

Vx=V_rel_B(1);Vy=V_rel_B(2);Vz=V_rel_B(3);
Vinf=[Vx,Vy,Vz];                            % relative incoming velocity

om=rpm*2*pi/60;                             % rad per second
vt=om*R ;                                   % tip velocity
mu= sqrt(Vx^2+Vy^2) / vt ;                  % advance ratio
lamb_tot = Vinf(3)/vt ;                     % climb inflow ratio

% the pitch angle of the quad-copter is already embded becuse the Vinf is
% in body frame!

alpha = atan(Vinf(3)/norm(Vinf,2));     % copter angle of attack

ct=T/(rho*pi*R^2*vt^2);

% climb ratio, prependicular to free stream
lambc=lamb_tot/cos(alpha);
%% finding lamb0

% lambda= mu*tan(alpha) + ct/(2*sqrt(mu^2+lambda(n)^2)) + lamc * cos(alpha)

er1=1 ;

% inflow for hover
lah = sqrt(0.5*ct);
% using exact solution for a more simplified case to behoove the initial guess
%lamb0 = lah * ((0.25*(mu/lah)^4+1)^0.5-0.5*(mu/lah)^2)^0.5;
lamb0=lah;
%   lamb0 = mu * tan (alpha) + lambc * cos(alpha) + 0.5 *  ct * (mu^2 + lamb0 ^2) ^ (-0.5) ;
%   it is noted that lamb_tot = mu * tan (alpha) + lambc * cos(alpha)


while er1>0.00001
    
%    f = lamb0 - mu * tan (alpha) - lambc * cos(alpha) - 0.5 *  ct * (mu^2 + lamb0 ^2) ^ (-0.5) ;
   
    f= lamb0 -  lamb_tot - 0.5 *  ct * (mu^2 + lamb0 ^2) ^ (-0.5) ;
    
    fp = 1 + 0.5 * ct * (mu^2 + lamb0 ^2) ^ (-1.5) * lamb0 ;
    
    lambnew = lamb0 - f/fp;
    
    er1 = abs((lambnew - lamb0)/lamb0) ;
    
    lamb0 = lambnew ;
    
end

%lambda_one=lamb0;

%dct=2*lamb0^2-ct;

%% page 158 - 160

% muz is the advance ratio prependicular to the rotor
% mux is the advance ratio parallel to the rotor

% muz = lambc * cos(alpha) ;
% mux = mu;

mux = mu ;
muz = lamb_tot;
lambi=lamb0;

% wake skew angle
x = atan (mux/(lamb0));

% Drees:
% kx = 4/3 * (1 - cos(x) - 1.8 * mu^2) / sin(x) ;
% 
% if x==0
%     kx=0;
% end
% 
% ky = -2*mu;

% Pitt and Peters
kx = (15*pi/23)*tan(x/2) ;
ky = 0;

%% fixing psi given a Vy velcoity, psi is defined when x is alinged with free stream

psi_new=psi-atan(-Vy/Vx);

% constructing inflow ratio based on a linear model
% lambda is a matrix, rows r, and columns azimuth, (nr * npsi)

lam= lamb0 * (1 + kx * r' * cos(psi_new) + ky * r' * sin(psi_new));

%% intergration to find roll and pitch for a rotor

% the incoming velocity seen by the blade
% for a quad-copter (no flpapping)

beta=0 ;
beta_dot=0 ;

r2=repmat(r',1,npsi);        % nr*npsi
psi2=repmat(psi_new,nr,1);       % nr*npsi

ut = r2 + mu*sin(psi2);

% up = lamb0 + lambc + mu * beta * cos (psi2) + r2 * beta_dot / om(1) ;

%%%%%%%% make sure to change om to om (t)

up = lam + mu * beta * cos (psi2) + r2 * beta_dot / om(1) ;

%ur = mu *cos(psi2) ;

u = (ut.^2 +up.^2) .^0.5 ;

% matrix form of the variables

%phi=atan(lam./r') ;                % nr*npsi
%phi=atan(bsxfun(@rdivide, lam, r'));
phi=atan(up./ut);
th2=repmat(th',1,npsi);      % nr*npsi
cla2=repmat(cla',1,npsi);     % nr*npsi
c2=repmat(c',1,npsi);        % nr*npsi

%diff2=diff(r2);
%dr2=[diff2(1,:);diff(r2)];
%da2=2*pi.*r2.*dr2;

% cl=cl(r,psi)
cl1 = cla2 .* (th2 - phi);

% based on Beard-McLain book page 47 to model stall condition
alp=(th2 - phi)-aLeq0;
sigma=(1+exp(-M*(alp-alp0)) + exp(M*(alp+alp0)))./...
      ((1+exp(-M*(alp-alp0))).*(1+exp(M*(alp+alp0))));
cl = (1-sigma).*cl1 + sigma .* (2* sign(alp) .* sin(alp).^2 .* cos(alp));
aoaeff=(th2 - phi)*180/pi;

%plot(cl);hold on;

% lift generated by every annulus
% intergral int_0^{2pi} cl * 1/2 * rho * U^2 * c2 / (2pi)
l = nb * trapz(psi, cl * 0.5 .* rho .* (vt*u).^2 .*c2.*cos(atan(up./ut)) , 2) / (2*pi);

%ct = l ./ (rho * 2 * pi * r' .* dr2(:,1) .* (r'*R*om).^2 * R^2) ;

%cT=trapz(r*R,ct)

% note that vt wil be cancel out in ct calculations

ct_two=trapz(r*R,l')/(rho*pi*R^2*vt^2);

% the error parameter defined between momentum and blade element theories
er_ct=abs(ct_two-ct)/ct;

% x point to the nose of the vehicle, y to the right wing (pilot's right),
% z pointing downward - right wing dowm, nose up, right wing back are
% positive roll, pitch, yaw moments. psi = 0, starts farther from the
% inflow, x<0

x=-r2.*cos(psi2)*R;
y=r2.*sin(psi2)*R;

%note, I already cancel out Vt
cmroll = trapz(r*R,...
    nb/(2*pi)*trapz(psi,-y.*cl*0.5*rho.*u.^2.*c2.*cos(atan(up./ut)),2))/(rho*pi*R^3);

cmpitch = trapz(r*R,...
    nb/(2*pi)*trapz(psi,x.*cl*0.5*rho.*u.^2.*c2.*cos(atan(up./ut)),2))/(rho*pi*R^3);


% blade profile drag
% Cd0=0.0081-0.0216*aoa+0.4*aoa^2; from helicopter literature, note aoa is the angle of attack - page 14, chapter 2 of Friedmann notes
Cd0=0.008;

% f is the equivalent flat plate area that copter frame occupies in space
% -- note that a quarter of the area should be used since this is model for
% for a helicopter, f/A is between 0.004 to 0.025

f=0.005*A * 0.25;

cmyaw=sig * Cd0 * 0.125 * (1 + 4.6* mu^2) + ...
    1.15 * 0.5 * ct_two^2 / sqrt(mu^2+lamb0^2) + ...
    f * mu^3 /(2*pi*R^2) + ...
    ct_two*lambc ;

% cmroll= (trapz(psi(1:end/2),r2(:,[1:end/2]).*cl(:,[1:end/2])*0.5*rho.*u(:,[1:end/2]).^2.*c2(:,[1:end/2]),2) -...
%               trapz(psi(end/2+1:end),r2(:,[end/2+1:end]).*cl(:,[end/2+1:end])*0.5*rho.*u(:,[end/2+1:end]).^2.*c2(:,[end/2+1:end]),2)) / (rho*pi*R^3*vt*2);

% cmpitch=0;
%
% mroll= trapz(r,trapz(psi(1:end/2),r2(:,[1:end/2]).*cl(:,[1:end/2]).*u(:,[1:end/2]).^2,2) -...
%                 trapz(psi(end/2+1:end),r2(:,[end/2+1:end]).*cl(:,[end/2+1:end]).*u(:,[1:end/2]).^2,2));

%outputs:
% er_ct    the difference between ct by momentum and blade element theories
% Mr        The roll moment
% Mp        The pitch moment
% My        The yaw moment 
% T         Thrust 
% cT        The thrust coefficient
% cr        The roll moment coefficient
% cp        The pitch moment coefficient

varargout{1}=er_ct;
if nargout>1
    varargout{2}=cmroll*(rho*pi*R^3*vt^2);
    varargout{3}=cmpitch*(rho*pi*R^3*vt^2);
    varargout{4}=cmyaw*(rho*pi*R^3*vt^2);
    varargout{5}=ct_two*(rho*pi*R^2*vt^2);
    varargout{6}=ct_two;
    varargout{7}=cmroll;
    varargout{8}=cmpitch;
end
