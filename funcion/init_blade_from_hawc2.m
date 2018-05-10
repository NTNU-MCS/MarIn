function [ bladePar ] = init_blade_from_hawc2(data_Blade_ae,data_Blade_clcdcm,data_Blade_st,data_Blade_posCenterLine_i)
% To generate the parameters for a blade by the inputs.
%         o  o  -------------------------------------------`
%      o  /  /  o    \        \        ......   \           \
%     o /     /  o    \        \       ......    \           \
%     o  /    /  o -------------------------------------------)
%      o  / /   o     /        /       ......    /          /
%         o  o  ------------------------------------------'
%           ?           ?      ?       ......       ?     ?
%  node#1:Root(Cicle)  node#2 node#3   ......      node#50 node#51:Tips(NACA0025)
%                       (FFA-W3-xxx)   ......   (FFA-W3-xxx)
% 
% Reference frame
%   | y
%   |   /  x
%   | /
%   o------> z


%% Import data
% n_node  =  200;

% r: curved length distance from main_body node 1 [m]
% bladePar.r_i    =   linspace(1,data_Blade_st(end,1),n_node)';    

bladePar.n_i    =  length(data_Blade_st(:,1));

bladePar.r_i    =  data_Blade_st(:,1);

% m: mass per unit length, kg/m
bladePar.rho_i  =   interp1(data_Blade_st(:,1),data_Blade_st(:,2),bladePar.r_i);    

% x_m,y_m: coordinate from C1/2 to mass center
bladePar.x_m_i =   interp1(data_Blade_st(:,1),data_Blade_st(:,3),bladePar.r_i);  
bladePar.y_m_i	=   interp1(data_Blade_st(:,1),data_Blade_st(:,4),bladePar.r_i); 

% r_ix,r_iy: radius of gyration related to elastic center [m].
bladePar.r_ix_i	=   interp1(data_Blade_st(:,1),data_Blade_st(:,5),bladePar.r_i); 
bladePar.r_iy_i	=   interp1(data_Blade_st(:,1),data_Blade_st(:,6),bladePar.r_i); 

% x_s,y_s: coordinate from C1/2 to shear center [m]
bladePar.x_s_i	=   interp1(data_Blade_st(:,1),data_Blade_st(:,7),bladePar.r_i); 
bladePar.y_s_i	=   interp1(data_Blade_st(:,1),data_Blade_st(:,8),bladePar.r_i); 

% E: modulus of elasticity [N/m2]
bladePar.E_i	=   interp1(data_Blade_st(:,1),data_Blade_st(:,9),bladePar.r_i); 

% G: shear modulus of elasticity
bladePar.G_i	=   interp1(data_Blade_st(:,1),data_Blade_st(:,10),bladePar.r_i); 

% I_x,I_y: area moment of inertia with respect to principal bending xe axis [m4]
bladePar.I_x_i =    interp1(data_Blade_st(:,1),data_Blade_st(:,11),bladePar.r_i);
bladePar.I_y_i	=   interp1(data_Blade_st(:,1),data_Blade_st(:,12),bladePar.r_i);

% K: torsional stiffness constant with respect to ze axis at the shear center [m4/rad]. 
bladePar.K_i	=   interp1(data_Blade_st(:,1),data_Blade_st(:,13),bladePar.r_i);   

% k_x,k_y: shear factor for force in principal bending xe/ye direction.
bladePar.k_x_i	=   interp1(data_Blade_st(:,1),data_Blade_st(:,14),bladePar.r_i);
bladePar.k_y_i	=   interp1(data_Blade_st(:,1),data_Blade_st(:,15),bladePar.r_i);

% A: cross sectional area [m2]
bladePar.A_i	=   interp1(data_Blade_st(:,1),data_Blade_st(:,16),bladePar.r_i);

% theta_s: structural pitch about zc2 axis. 
bladePar.theta_s_i	=   interp1(data_Blade_st(:,1),data_Blade_st(:,17),bladePar.r_i);

% x_e,y_e: coordinate from C1/2 to center of elasticity [m]. 
bladePar.x_e_i	=   interp1(data_Blade_st(:,1),data_Blade_st(:,18),bladePar.r_i);
bladePar.y_e_i	=   interp1(data_Blade_st(:,1),data_Blade_st(:,19),bladePar.r_i);

% chord length
bladePar.chordLength_i  = interp1(data_Blade_ae.r,data_Blade_ae.chordLength,bladePar.r_i);

% thickness ratio 
bladePar.thicknessRatio_i  = interp1(data_Blade_ae.r,data_Blade_ae.thicknessRatio,bladePar.r_i);

% wind attack angle (deg)
bladePar.alpha_i = data_Blade_clcdcm.alpha;

bladePar.data_Blade_clcdcm.thickness = data_Blade_clcdcm.thickness;
bladePar.data_Blade_clcdcm.alpha = data_Blade_clcdcm.alpha;
bladePar.data_Blade_clcdcm.Cl = data_Blade_clcdcm.Cl;
bladePar.data_Blade_clcdcm.Cd = data_Blade_clcdcm.Cd;
bladePar.data_Blade_clcdcm.Cm = data_Blade_clcdcm.Cm;

%% Blade aerodynamic coefficients
% Air density
bladePar.rho_air = 1.225;   

% lifitng coefficient CL
bladePar.Cl_i = interp2(data_Blade_clcdcm.thickness,data_Blade_clcdcm.alpha,data_Blade_clcdcm.Cl,bladePar.thicknessRatio_i,bladePar.alpha_i);

% drag coefficient CD
bladePar.Cd_i = interp2(data_Blade_clcdcm.thickness,data_Blade_clcdcm.alpha,data_Blade_clcdcm.Cd,bladePar.thicknessRatio_i,bladePar.alpha_i);

% Moment coefficient CM
bladePar.Cm_i = interp2(data_Blade_clcdcm.thickness,data_Blade_clcdcm.alpha,data_Blade_clcdcm.Cm,bladePar.thicknessRatio_i,bladePar.alpha_i);

% position of blade center line
bladePar.posCenterLine_i(:,1) = interp1(data_Blade_posCenterLine_i(:,3),data_Blade_posCenterLine_i(:,1),bladePar.r_i);
bladePar.posCenterLine_i(:,3) = -interp1(data_Blade_posCenterLine_i(:,3),data_Blade_posCenterLine_i(:,2),bladePar.r_i);
bladePar.posCenterLine_i(:,2) = bladePar.r_i;
bladePar.posCenterLine_i_thetay = interp1(data_Blade_posCenterLine_i(:,3),data_Blade_posCenterLine_i(:,4),bladePar.r_i);

%% Blade geometry infomation
% Blade curved length
bladePar.L = bladePar.r_i(end);         
                                                                                                                                                                                      
% Number of nodes on the blade
bladePar.no_node = length(bladePar.r_i);

% Distance between a node and its next node [no_node-1,1]
delta_l_i = bladePar.r_i(2:end) - bladePar.r_i(1:end-1);

% Distance of acting mass of the node
bladePar.delta_r_i =  0.5*[delta_l_i;0] + 0.5*[0;delta_l_i];

% Mass of each node
bladePar.m_i = bladePar.delta_r_i.*bladePar.rho_i;

% total mass
bladePar.totalMass = sum(bladePar.m_i);     

% Caluclate the centerline position in the blade coordinate
bladePar.posMassCenter_i(:,1) = bladePar.posCenterLine_i(:,1) + bladePar.x_m_i.*cosd(bladePar.posCenterLine_i_thetay) - bladePar.y_m_i.*sind(bladePar.posCenterLine_i_thetay);
bladePar.posMassCenter_i(:,3) = bladePar.posCenterLine_i(:,3) - bladePar.x_m_i.*sind(bladePar.posCenterLine_i_thetay) - bladePar.y_m_i.*cosd(bladePar.posCenterLine_i_thetay);
bladePar.posMassCenter_i(:,2) = bladePar.posCenterLine_i(:,2);

% Calculate the position of COG
bladePar.posCOG = zeros(3,1);
bladePar.posCOG(1) = sum(bladePar.posMassCenter_i(:,1).*bladePar.m_i)/bladePar.totalMass;
bladePar.posCOG(2) = sum(bladePar.posMassCenter_i(:,2).*bladePar.m_i)/bladePar.totalMass;
bladePar.posCOG(3) = sum(bladePar.posMassCenter_i(:,3).*bladePar.m_i)/bladePar.totalMass;

% Caluclate the elasticity center position in the blade coordinate
bladePar.posCE_E(:,1) = bladePar.posCenterLine_i(:,1) + bladePar.x_e_i.*cosd(bladePar.posCenterLine_i_thetay) - bladePar.y_e_i.*sind(bladePar.posCenterLine_i_thetay);
bladePar.posCE_E(:,3) = bladePar.posCenterLine_i(:,3) - bladePar.x_e_i.*sind(bladePar.posCenterLine_i_thetay) - bladePar.y_e_i.*cosd(bladePar.posCenterLine_i_thetay);
bladePar.posCE_E(:,2) = bladePar.posCenterLine_i(:,2);


% Caluclate the 1/4 center position in the blade coordinate
bladePar.posC_14_i(:,1) = bladePar.posCenterLine_i(:,1) + bladePar.chordLength_i*0.25.*cosd(bladePar.posCenterLine_i_thetay) - bladePar.chordLength_i*0.25.*sind(bladePar.posCenterLine_i_thetay);
bladePar.posC_14_i(:,3) = bladePar.posCenterLine_i(:,3) - bladePar.chordLength_i*0.25.*sind(bladePar.posCenterLine_i_thetay) - bladePar.chordLength_i*0.25.*cosd(bladePar.posCenterLine_i_thetay);
bladePar.posC_14_i(:,2) = bladePar.posCenterLine_i(:,2);

% Caluclate the 3/4 center position in the blade coordinate
bladePar.posC_34_i(:,1) = bladePar.posCenterLine_i(:,1) - bladePar.chordLength_i*0.25.*cosd(bladePar.posCenterLine_i_thetay) + bladePar.chordLength_i*0.25.*sind(bladePar.posCenterLine_i_thetay);
bladePar.posC_34_i(:,3) = bladePar.posCenterLine_i(:,3) + bladePar.chordLength_i*0.25.*sind(bladePar.posCenterLine_i_thetay) + bladePar.chordLength_i*0.25.*cosd(bladePar.posCenterLine_i_thetay);
bladePar.posC_34_i(:,2) = bladePar.posCenterLine_i(:,2);


for i = 1:length(bladePar.r_i)
    bladePar.posC_14_i_rltv_COG(i,:) =  bladePar.posC_14_i(i,:) - bladePar.posCOG';
    bladePar.posC_34_i_rltv_COG(i,:) =  bladePar.posC_34_i(i,:) - bladePar.posCOG';
end

Ixx_o = sum(bladePar.m_i.* ((bladePar.posMassCenter_i(:,2) ).^2 + ( bladePar.posMassCenter_i(:,3) ).^2 ) ); 
Iyy_o = sum(bladePar.m_i.* (bladePar.r_ix_i.^2 + bladePar.r_iy_i.^2) );
Izz_o = sum( (bladePar.posCE_E(:,1).^2 + bladePar.posCE_E(:,2).^2).*bladePar.m_i );

Ixy_o = sum( bladePar.m_i.* (bladePar.posMassCenter_i(:,1)).* (bladePar.posMassCenter_i(:,2) ) );
Iyz_o = sum( bladePar.m_i.* bladePar.posMassCenter_i(:,2).* bladePar.posMassCenter_i(:,3) );

% second area moment of enertia
I_xx =   bladePar.r_ix_i.^2.*bladePar.A_i;
I_zz =   bladePar.r_iy_i.^2.*bladePar.A_i;
Ixz_o =  2*sum((I_zz-I_xx).*tand(2*bladePar.posCenterLine_i_thetay).*bladePar.m_i);


% Moment of inertia at the center of the root circle [0;0;0]
bladePar.I_o = [Ixx_o -Ixy_o -Ixz_o; -Ixy_o Iyy_o -Iyz_o; -Ixz_o -Iyz_o Izz_o];

% Moment of inertia at COG
bladePar.I_COG = bladePar.I_o - bladePar.totalMass*(bladePar.posCOG'*bladePar.posCOG*eye(3) - bladePar.posCOG*bladePar.posCOG');

bladePar.M = [bladePar.totalMass*eye(3),zeros(3,3);zeros(3,3),bladePar.I_COG];

bladePar.invM = inv(bladePar.M);

bladePar.posMassCenter_i = bladePar.posMassCenter_i';
bladePar.posCE_E = bladePar.posCE_E';
bladePar.posC_14_i = bladePar.posC_14_i';
bladePar.posC_34_i = bladePar.posC_34_i';
bladePar.posC_14_i_rltv_COG = bladePar.posC_14_i_rltv_COG';
bladePar.posC_34_i_rltv_COG = bladePar.posC_34_i_rltv_COG';

%%
% % Moment of inertia in x-axis at COG
% Ixx = sum(bladePar.m_i.*(...
%                       (bladePar.posMassCenter_i(:,2) - bladePar.posCOG(2)).^2 + ( bladePar.posMassCenter_i(:,3) - bladePar.posCOG(3)).^2  ...
%                         )  );
% 
% % Moment of inertia in y-axis at COG
% Iyy = sum(bladePar.m_i.*(...
%             (bladePar.posMassCenter_i(:,1) - bladePar.posCOG(1)).^2 + (bladePar.posMassCenter_i(:,3) - bladePar.posCOG(3)).^2 ...
%                         )   );
% 
% % Moment of inertia in z-axis at COG
% Izz = sum(bladePar.m_i.*((bladePar.posMassCenter_i(:,1) - bladePar.posCOG(1)).^2+(bladePar.posMassCenter_i(:,2) - bladePar.posCOG(1)).^2));
% 
% % I_xy
% Ixy = - sum( bladePar.m_i.* (bladePar.posMassCenter_i(:,1) - bladePar.posCOG(1)).* (bladePar.posMassCenter_i(:,2) - bladePar.posCOG(2)));
% 
% % I_xz
% Ixz = - sum( bladePar.m_i.* (bladePar.posMassCenter_i(:,1) - bladePar.posCOG(1)).* (bladePar.posMassCenter_i(:,3) - bladePar.posCOG(3)));
% 
% % I_yz
% Iyz = - sum( bladePar.m_i.* (bladePar.posMassCenter_i(:,2) - bladePar.posCOG(2)).* (bladePar.posMassCenter_i(:,3) - bladePar.posCOG(3)));
% bladePar.I_COG = [Ixx Ixy Ixz; Ixy Iyy Iyz; Ixz Iyz Izz];


% if plotBlade_flag == 1
%     subplot(2,1,1)
%     plot(x,z1); hold on; plot(x,z2)
%     axis equal;
%     subplot(2,1,2)
%     
% end