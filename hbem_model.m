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

function varargout=hbem_model(T,V_rel_B)
%% Geometry of the blade

global geometry2

R=4*0.0254;                          % propeller radius  [m]
nb=2;                                % number of blade
A=pi*R^2;                            % rotor disk area
M=50;                                  % stall transition rate
alp0=20.6*pi/180;                      % stall cut-off aoa
aLeq0=2*pi/180;                        % absolute values of angle of attack where lift is zero

% blade characteristics

load('MR8x45');
c=MR8x45(:,2)*0.0254;
th=[MR8x45(:,8)]*pi/180+aLeq0;         % aLeq0 must be added here!
r=MR8x45(:,1)/4;                       % normolized radial locations, R=4in

rho=1.2;
nr=size(MR8x45,1);                    % number of data points in radial location
npsi=60;                               % number of data point in azimuth

% plotting chord and twist
% reset(gcf);reset(gca)
% set(0,'defaultLineLineWidth',2)
% set(0,'defaultAxesFontSize',12)
% set(0,'defaultTextFontWeight','bold')
% set(0,'DefaultAxesFontWeight','bold')
%
% h=subplot(2,1,1)
% set(h, 'Position', [.12 .55 .8 .4]);
% plot(r,c*100)
% set(gca,'xticklabel',{[]})
% ylabel('Chord [cm]')
% %%legend('Ascend-straight-descend path','Location','best')
% h=subplot(2,1,2)
% set(h, 'Position', [.12 .12 .8 .4]);
% plot(r,MR8x45(:,8))
% %%legend('Circular path','Location','best')
% ylabel('Twist angle [deg]')
% xlabel('Normalized span radial location (r/R) [cm]')

cla=repmat(1.8*pi,nr,1)';             % 2-D lift curve slope

psi=linspace(0,2*pi,npsi);             % azimuth angle

% constructing all geometry and aerodynamic data as one variable
% "geometry2"
maxsize=max(nr,npsi);
numvar=14;
geometry2=zeros(numvar,maxsize);

list={R,nb,A,rho,nr,npsi,M,alp0,aLeq0,th,c,cla,r,psi};

for  i=1:numvar   
    geometry2(i,1:length(list{i}))=[list{i}];    
end

subfcn2 = @(rpm) try11(rpm,V_rel_B,T);
rpm=golden_search(100,16000,subfcn2);

mu=sqrt(V_rel_B(1)^2+V_rel_B(2)^2)/(rpm*2*pi/60*R);
lambda=-V_rel_B(3)/(rpm*2*pi/60*R);

[~,~,~,~,~,Re,aoaeff,cl]=try1(rpm,V_rel_B,T);

varargout{1}=rpm;
if nargout>1
    varargout{2}=mu;
    varargout{3}=lambda;
    varargout{4}=Re;
    varargout{5}=aoaeff;
    varargout{6}=cl;
end
end