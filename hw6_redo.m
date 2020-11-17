

function R = rotmatzx(theta, alpha)
    R = [cos(theta) -sin(theta)*cos(alpha)  sin(theta)*sin(alpha);
         sin(theta)  cos(theta)*cos(alpha) -cos(theta)*sin(alpha);
                  0             sin(alpha)             cos(alpha)];
end

function thetas = ikinelbow(ais, dis, Tw_tool):
    alpha = [-pi/2 0 -pi/2 pi/2 -pi/2 0 0];
    theta = zeros(6,0);
    Tw_0 = [ sqrt(2)/2 sqrt(2)/2 0  221;
            -sqrt(2)/2 sqrt(2)/2 0   22;
                 0         0     1 1104;
                 0         0     0    1];
    T0_w = inv(Tw_0);

    # revome world view
    T0_tool = mtimes(T0_w, Tw_tool);
    R0_tool = T0_tool(1:3,1:3)
    d0_tool = T0_tool(1:3,4)

    # find R0_6 and d0_4 (ie remove tool view)
    R0_6 = mtimes(R0_tool, rotmatx(alphas(7)).');
    z0_6 = R0_6(:,4);
    x0_6 = R0_6(:,1);
    d4_tool = dis(6)*z0_6 + dis(7)*z0_6 + ais(7)*x0_6;
    d0_4 = d0_tool - d4_tool;

    # solve for joint angles 1-3
    x1_1 = [1 0 0].';
    z0_0 = [0 0 1].';
    thetas(1) = atan2(d0_4(2), d0_4(1));
    R0_1 = rotmatzx(thetas(1), alphas(1));
    d1_4 = mtimes( R0_1.', (d0_4 - dis(1)*z0_0) ) - ais(1)*x1_1;

    x = d1_4(1);
    y = d1_4(2);
    a1h = a2;
    a2h = d4;
    a12p = (a1h + a2h)^2;
    a12m = (a1h - a2h)^2;
    xy2 = x^2 + y^2;
    top = a12p - xy2;
    bot = xy2 - a12m;
    theta(3) = 2 * atan( sqrt( top / bot ) );

    phi = atan2(y, x);
    yh = a2h * sin(thetas(3));
    xh = a1h + a2h*cos(thetas(3));
    psi = atan2(yh, xh);
    thetas(2) = phi - psi;

    # subtract pi/2 to account for 0 position
    thetas(3) = thetas(3) - pi/2;

    # solve for joint angles 4-6
    R1_2 = rotmatzx(thetas(2), alphas(2));
    R2_3 = rotmatzx(thetas(3), alphas(3));
    R0_3 = mtimes(R0_1, mtimes(R1_2, R2_3));
    R3_6 = mtimes(R0_3.', R0_6);

    thetas(4) = atan(-R3_6(2,3) / -R3_6(1,3));
    if thetas(4) > (pi/2)
        thetas(4) = thetas(4) - pi
    elseif thetas(4) < (-pi/2)
        thetas(4) = thetas(4) + pi
    end

    yh = -R3_6(1,3)*cos(thetas(4)) - R3_6(2,3)*sin(thetas(4));
    xh = R3_6(3,3);
    thetas(5) = atan2(yh, xh);

    yh = -R3_6(1,1)*sin(thetas(4)) + R3_6(2,1)*cos(thetas(4));
    xh = -R3_6(1,2)*sin(thetas(4)) + R3_6(2,2)*cos(thetas(4));
    thetas(6) = atan2(yh, xh);
end
