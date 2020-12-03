% linktrans
%
%   linktrans(a,d,alpha,theta)
%   takes in the dh parameters and returns a homogenous matrix 0T1 that
%   maps from 0 to 1.
%
%   a:float | distance along xi from zi-1 to zi
%   d:float | distance along zi-1 from xi-1 to xi
%   alpha:float | radian angle around xi from zi-1 to zi
%   theta:float | radian angle around zi-1 from xi-1 to xi
%
% linktrans
%
%   linktrans(a,d,alpha,theta)
%   takes in the dh parameters and returns a homogenous matrix 0T1 that
%   maps from 0 to 1.
%
%   a:float | distance along xi from zi-1 to zi
%   d:float | distance along zi-1 from xi-1 to xi
%   alpha:float | radian angle around xi from zi-1 to zi
%   theta:float | radian angle around zi-1 from xi-1 to xi
%
% ikinelbow
%
%   ikinelbow(ais, dis, Tw_tool)
%   takes in a T matrix (4x4) and returns the angles necessary to move
%   the arm into that position
%
%   ais:array | distances along xi from zi-1 to zi
%   dis:array | distances along zi-1 from xi-1 to xi
%   Tw_tool:float | T matrix designating end position of the tool in world view
%
%   Iain Lee
%   Robotics
%   10-22-20

a = [69 370.82 0 0 0 0 0];
d = [270.35 0 0 374.29 0 368.3 0];

Tw_tool = [0.984  0.177  0.024 727.089;
           0.177 -0.984  0.020 418.581;
           0.027 -0.015 -0.999 732.233;
             0      0      0      1   ];

theta = ikinelbow(a, d, Tw_tool, 1, 0, 0);

% tests = load('test_cases.mat')
% fileID = fopen('joint_angles.txt','w');
% 
% theta = ikinelbow(a, d, tests.T1, 1, 0, 0);
% fprintf(fileID,'%8.3f,',theta);
% fprintf(fileID,'\n');
% 
% theta = ikinelbow(a, d, tests.T2, 1, 0, 0);
% fprintf(fileID,'%8.3f,',theta);
% fprintf(fileID,'\n');

poses = load('poses.mat');
fileID = fopen('solutions.txt','w');

theta = ikinelbow(a, d, poses.above_goal, 1, 0, 0);
fprintf(fileID,'%8.3f,',theta);
fprintf(fileID,'\n');

theta = ikinelbow(a, d, poses.above_object, 1, 0, 0);
fprintf(fileID,'%8.3f,',theta);
fprintf(fileID,'\n');

theta = ikinelbow(a, d, poses.goal, 1, 0, 0);
fprintf(fileID,'%8.3f,',theta);
fprintf(fileID,'\n');

theta = ikinelbow(a, d, poses.midpoint, 1, 0, 0);
fprintf(fileID,'%8.3f,',theta);
fprintf(fileID,'\n');

theta = ikinelbow(a, d, poses.object, 1, 0, 0);
fprintf(fileID,'%8.3f,',theta);
fprintf(fileID,'\n');


function R = rotmatx(alpha)
    R = [1          0           0;
         0 cos(alpha) -sin(alpha);
         0 sin(alpha)  cos(alpha)];
end

function R = rotmatzx(theta, alpha)
    R = [cos(theta) -sin(theta)*cos(alpha)  sin(theta)*sin(alpha);
         sin(theta)  cos(theta)*cos(alpha) -cos(theta)*sin(alpha);
                  0             sin(alpha)             cos(alpha)];
end

function thetas = ikinelbow(ais, dis, Tw_tool,LR,UD,NF)
    alphas = [-pi/2 0 -pi/2 pi/2 -pi/2 0 0];
    thetas = zeros(6,0);
    Tw_0 = [ sqrt(2)/2 sqrt(2)/2 0  221;
            -sqrt(2)/2 sqrt(2)/2 0   22;
                 0         0     1 1104;
                 0         0     0    1];
    T0_w = inv(Tw_0);

%   remove world view
    T0_tool = mtimes(T0_w, Tw_tool);
    R0_tool = T0_tool(1:3,1:3);
    d0_tool = T0_tool(1:3,4);

%   find R0_6 and d0_4 (ie remove tool view)
    R0_6 = mtimes(R0_tool, rotmatx(alphas(7)).');
    z0_6 = R0_6(1:3,3);
    x0_6 = R0_6(1:3,1);
    d4_tool = dis(6)*z0_6 + dis(7)*z0_6 + ais(7)*x0_6;
    d0_4 = d0_tool - d4_tool;

%   solve for joint angles 1-3
    x1_1 = [1 0 0].';
    z0_0 = [0 0 1].';
    if LR == 1
      thetas(1) = atan2(d0_4(2), d0_4(1));
    else
      thetas(1) = atan(d0_4(2) / d0_4(1));
    end
    R0_1 = rotmatzx(thetas(1), alphas(1));
    d1_4 = mtimes( R0_1.', (d0_4 - dis(1)*z0_0) ) - ais(1)*x1_1;

    x = d1_4(1);
    y = d1_4(2);
    a1h = ais(2);
    a2h = dis(4);
    a12p = (a1h + a2h)^2;
    a12m = (a1h - a2h)^2;
    xy2 = x^2 + y^2;
    top = a12p - xy2;
    bot = xy2 - a12m;
    thetas(3) = 2 * atan( sqrt( top / bot ) );
    if UD == 1
      thetas(3) = -thetas(3);
    end

    phi = atan2(y, x);
    yh = a2h * sin(thetas(3));
    xh = a1h + a2h*cos(thetas(3));
    psi = atan2(yh, xh);
    thetas(2) = phi - psi;

%   subtract pi/2 to account for 0 position
    thetas(3) = thetas(3) - pi/2;

%   solve for joint angles 4-6
    R1_2 = rotmatzx(thetas(2), alphas(2));
    R2_3 = rotmatzx(thetas(3), alphas(3));
    R0_3 = mtimes(R0_1, mtimes(R1_2, R2_3));
    R3_6 = mtimes(R0_3.', R0_6);

    thetas(4) = atan(-R3_6(2,3) / -R3_6(1,3));
    if NF == 0
      if thetas(4) > (pi/2)
          thetas(4) = thetas(4) - pi;
      elseif thetas(4) < (-pi/2)
          thetas(4) = thetas(4) + pi;
      end
    else
      if thetas(4) < (pi/2) & thetas(4) > (-pi/2)
        thetas(4) = thetas(4) + pi;
      end
    end

    yh = -R3_6(1,3)*cos(thetas(4)) - R3_6(2,3)*sin(thetas(4));
    xh = R3_6(3,3);
    thetas(5) = atan2(yh, xh);

    yh = -R3_6(1,1)*sin(thetas(4)) + R3_6(2,1)*cos(thetas(4));
    xh = -R3_6(1,2)*sin(thetas(4)) + R3_6(2,2)*cos(thetas(4));
    thetas(6) = atan2(yh, xh);
end
