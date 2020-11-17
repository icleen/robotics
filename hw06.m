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
%   Iain Lee
%   Robotics
%   10-22-20


% test
% temp = linktrans(1, 1, 0, pi/2)
% temp2 = linktrans(0, 2, pi/2, 0)
% comb = mtimes(temp, temp2)
% comb2 = fkine([1 0], [1 2], [0 pi/2], [pi/2 0])
a = [69 370.82 0 0 0 0];
d = [270.35 0 0 374.29 0 368.2];
alpha = [-pi/2 0 -pi/2 pi/2 -pi/2 0];
theta = [pi/2 -pi/4 pi/3 pi pi/6 pi/2]

T = fkine(a, d, alpha, theta)

% calculated answer:
% T = [-1.0000   -0.0000   -0.0000   -0.0000;
%      -0.0000    0.9659    0.2588  329.6331;
%       0.0000    0.2588   -0.9659 -184.6309;
%            0         0         0    1.0000];

theta = ikinebaxter(a(1), a(2), d(1), d(4), d(6), T, 1, 0, 1)


function T = linktrans(a,d,alpha,theta)

    ctheta = cos(theta);
    stheta = sin(theta);

    calpha = cos(alpha);
    salpha = sin(alpha);

    T = [ctheta -stheta*calpha   stheta*salpha a*ctheta;
         stheta  ctheta*calpha  -ctheta*salpha a*stheta;
         0       salpha          calpha        d;
         0       0               0             1];

end


function T = fkine(a,d,alpha,theta)

    i = 1;
    oT = linktrans(a(i), d(i), alpha(i), theta(i));
    for i=2:length(a)
        nT = linktrans(a(i), d(i), alpha(i), theta(i));
        oT = mtimes(oT, nT);
    end
    T = oT;

end

function R = rotmatzx(theta, alpha)
    R = [cos(theta) -sin(theta)*cos(alpha)  sin(theta)*sin(alpha);
         sin(theta)  cos(theta)*cos(alpha) -cos(theta)*sin(alpha);
                  0             sin(alpha)             cos(alpha)];
end

function thetas = ikinebaxter(a1,a2,d1,d4,d6,T,LR,UD,NF)
    alphas = [-pi/2 0 -pi/2 pi/2 -pi/2 0];
    thetas = zeros(6,0);
    R0_6 = T(1:3,1:3);
    d0_6 = T(1:3,4);

%   check workspace
    maxdist = a2 + d4 + d6;
    xydist = sqrt(d0_6(1)^2 + d0_6(2)^2);
    if xydist > (maxdist+a1)
        "return None"
    end
    if d0_6(3) > (maxdist+d1)
        "return None"
    end

%   find d0_4
    z0_6 = R0_6(1:3,3);
    d0_4 = d0_6 - d6*z0_6;

%   find thetas 1-3
    x1_1 = [1 0 0].';
    z0_0 = [0 0 1].';
%   LR decides lefty or righty
    if LR == 1
      thetas(1) = atan2(d0_4(2), d0_4(1));
    else
      thetas(1) = atan(d0_4(2) / d0_4(1));
    end
    R0_1 = rotmatzx(thetas(1), alphas(1));
    d1_4 = mtimes( R0_1.', (d0_4 - d1*z0_0) ) - a1*x1_1;

    x = d1_4(1);
    y = d1_4(2);
    a1h = a2;
    a2h = d4;
    a12p = (a1h + a2h)^2;
    a12m = (a1h - a2h)^2;
    xy2 = x^2 + y^2;
    top = a12p - xy2;
    bot = xy2 - a12m;
    thetas(3) = 2 * atan( sqrt( top / bot ) );
%   whether you choose + or minus needs to depend on UD (elbow up or down)
    if UD == 1
      thetas(3) = -thetas(3);
    end
%     thetas(3) = pi - thetas(3);

    phi = atan2(y, x);
    yh = a2h * sin(thetas(3));
    xh = a1h + a2h*cos(thetas(3));
    psi = atan2(yh, xh);
    thetas(2) = phi - psi;

    thetas(3) = thetas(3) - pi/2;

%   find thetas 4-6
    R1_2 = rotmatzx(thetas(2), alphas(2));
    R2_3 = rotmatzx(thetas(3), alphas(3));
    R0_3 = mtimes(R0_1, mtimes(R1_2, R2_3));
    R3_6 = mtimes(R0_3.', R0_6);

    y = R3_6(2,3);
    x = R3_6(1,3);
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
