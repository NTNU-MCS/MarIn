function R = Rzxy(phi,theta,psi)
% R = Rzxy(phi,theta,psi) computes the Euler angle
% rotation matrix R in SO(3) using the zyx convention
%
% Author:   Zhengru Ren
% Date:     18.9.2017



R = [   cos(psi)*cos(theta) - sin(phi)*sin(psi)*sin(theta),     -cos(phi)*sin(psi),     cos(psi)*sin(theta) + cos(theta)*sin(phi)*sin(psi);
        cos(theta)*sin(psi) + cos(psi)*sin(phi)*sin(theta),     cos(phi)*cos(psi),    	sin(psi)*sin(theta) - cos(psi)*cos(theta)*sin(phi);
    	-cos(phi)*sin(theta),                                   sin(phi),           	cos(phi)*cos(theta)];
