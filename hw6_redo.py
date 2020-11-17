
from math import cos, sin, sqrt, atan, atan2, pi
import numpy as np

def linktrans(a, d, alpha, theta):
    ctheta = cos(theta)
    stheta = sin(theta)
    calpha = cos(alpha)
    salpha = sin(alpha)
    return np.array(
      [[ctheta, -stheta*calpha,   stheta*salpha, a*ctheta],
       [stheta,  ctheta*calpha,  -ctheta*salpha, a*stheta],
       [0,       salpha,          calpha,        d],
       [0,       0,               0,             1]] )

def fkine(ai, di, alpha, theta):
    mat = linktrans(ai[0], di[0], alpha[0], theta[0])
    for i in range(1,len(ai)):
        mat = np.matmul(mat, linktrans(ai[i], di[i], alpha[i], theta[i]))
    return mat


def revmat(theta, alpha):
    ctheta = cos(theta)
    stheta = sin(theta)
    calpha = cos(alpha)
    salpha = sin(alpha)
    return np.array(
      [[ctheta, -stheta*calpha,  stheta*salpha],
       [stheta,  ctheta*calpha, -ctheta*salpha],
       [     0,         salpha,         calpha]]
    )

def rotmatz(theta):
    return np.array([
      [cos(theta), -sin(theta), 0],
      [sin(theta),  cos(theta), 0],
      [    0,           0,      1]
    ])

def rotmatx(alpha):
    return np.array([
      [1,          0,           0],
      [0, cos(alpha), -sin(alpha)],
      [0, sin(alpha),  cos(alpha)]
    ])

def rotmatzx(theta, alpha):
    ctheta = cos(theta)
    stheta = sin(theta)
    calpha = cos(alpha)
    salpha = sin(alpha)
    return np.array(
      [[ctheta, -stheta*calpha,  stheta*salpha],
       [stheta,  ctheta*calpha, -ctheta*salpha],
       [     0,         salpha,         calpha]]
    )

def ikinebaxter(a1,a2,d1,d4,d6,T,LR,UD,NF,alphas=None):
    if alphas is None:
        alphas = np.array([-pi/2, 0, -pi/2, pi/2, -pi/2, 0])
    thetas = np.zeros(6)

    R0_6 = T[:3,:3]
    d0_6 = T[:3,-1]

    # check workspace
    maxdist = a2 + d4 + d6
    xydist = sqrt(d0_6[0]**2 + d0_6[1]**2)
    if xydist > (maxdist+a1):
        return None
    if d0_6[2] > (maxdist+d1):
        return None

    x0_6 = R0_6[:,0]
    z0_6 = R0_6[:,-1]
    x0_6 = R0_6[:,0]
    d0_4 = d0_6 - d4*z0_6

    # solve for joint angles 1-3
    x1_1 = np.array([1, 0, 0])
    z0_0 = np.array([0, 0, 1])
    thetas[0] = atan(d0_4[1] / d0_4[0])
    if LR==1:
        thetas[0] = atan2(d0_4[1], d0_4[0])
    R0_1 = rotmatzx(thetas[0], alphas[0])
    d1_4 = np.matmul( R0_1.T, (d0_4 - d1*z0_0) ) - a1*x1_1

    x, y = d1_4[:2]
    a1h, a2h = a2, d4
    a12p = (a1h + a2h)**2
    a12m = (a1h - a2h)**2
    xy2  = x**2 + y**2
    top = a12p - xy2
    bot = xy2 - a12m
    thetas[2] = 2*atan(sqrt( top / bot ))
    elbowup = UD==1
    if elbowup:
        thetas[2] *= -1
    # thetas[2] = pi - thetas[2]

    phi = atan2(y, x)
    yh = a2h * sin(thetas[2])
    xh = a1h + a2h*cos(thetas[2])
    psi = atan2(yh, xh)
    thetas[1] = phi - psi

    # subtract pi/2 to account for 0 position
    thetas[2] = thetas[2] - pi/2

    # solve for joint angles 4-6
    R1_2 = rotmatzx(thetas[1], alphas[1])
    R2_3 = rotmatzx(thetas[2], alphas[2])
    R0_3 = np.matmul(R0_1, np.matmul(R1_2, R2_3))
    R3_6 = np.matmul(R0_3.T, R0_6)

    thetas[3] = atan(-R3_6[1,2] / -R3_6[0,2])
    if NF == 0:
        if thetas[3] > (pi/2):
            thetas[3] -= pi
        elif thetas[3] < (-pi/2):
            thetas[3] += pi
    else:
        if thetas[3] < (pi/2) and thetas[3] > (-pi/2):
            # flip = 1 if thetas[3] < 0 else -1
            # thetas[3] += pi*flip
            thetas[3] += pi

    yh = -R3_6[0,2]*cos(thetas[3]) - R3_6[1,2]*sin(thetas[3])
    xh = R3_6[2,2]
    thetas[4] = atan2(yh, xh)

    yh = -R3_6[0,0]*sin(thetas[3]) + R3_6[1,0]*cos(thetas[3])
    xh = -R3_6[0,1]*sin(thetas[3]) + R3_6[1,1]*cos(thetas[3])
    thetas[5] = atan2(yh, xh)

    return thetas

if __name__ == '__main__':
    # pi = np.pi
    # acc = 0
    #
    # ais = [0, 1, 0, 0, 0, 0]
    # dis = [1, 0, 0, 1, 0, 0]
    # alphas = [-pi/2, 0, -pi/2, pi/2, -pi/2, 0]

    print('\nactual')
    ais = [69, 370.82, 0, 0, 0, 0]
    dis = [270.35, 0, 0, 374.29, 0, 368.2]
    alphas = [-pi/2, 0, -pi/2, pi/2, -pi/2, 0]
    thetas = [pi/2, -pi/4, pi/3, pi, pi/6, pi/2]
    mat = fkine( ais, dis, alphas, thetas )
    thets = ikinebaxter(
      ais[0], ais[1], dis[0], dis[3], dis[5], mat, 1, 0, 1, alphas )
    print('thetas:', np.round(thetas, 2))
    print('fkine:', np.round(mat[:-1,-1], 2))
    print('ikine:', np.round(thets, 2))
