
from math import cos, sin, sqrt, atan, atan2, pi
import numpy as np


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

def ikinelbow(ais, dis, Tw_tool):
    thetas = np.zeros(6)
    alphas = np.array([-pi/2, 0, -pi/2, pi/2, -pi/2, 0, 0])
    Tw_0 = np.array([
      [ sqrt(2)/2, sqrt(2)/2, 0,  221],
      [-sqrt(2)/2, sqrt(2)/2, 0,   22],
      [     0,         0,     1, 1104],
      [     0,         0,     0,    1]
    ])
    T0_w = np.linalg.inv(Tw_0)
    # T0_w2 = Tw_0.copy()
    # T0_w2[:3,:3] = T0_w2[:3,:3].T
    # T0_w2[:3, -1] = -1*np.matmul(T0_w2[:3, :3], T0_w2[:3, -1])
    # import pdb; pdb.set_trace()

    # revome world view
    T0_tool = np.matmul(T0_w, Tw_tool)
    R0_tool = T0_tool[:3,:3]
    d0_tool = T0_tool[:3,-1]

    # find R0_6 and d0_4 (ie remove tool view)
    R0_6 = np.matmul(R0_tool, rotmatx(alphas[6]).T)
    z0_6 = R0_6[:,-1]
    x0_6 = R0_6[:,0]
    d4_tool = dis[5]*z0_6 + dis[6]*z0_6 + ais[6]*x0_6
    d0_4 = d0_tool - d4_tool

    # solve for joint angles 1-3
    x1_1 = np.array([1, 0, 0])
    z0_0 = np.array([0, 0, 1])
    thetas[0] = atan(d0_4[1] / d0_4[0])
    R0_1 = rotmatzx(thetas[0], alphas[0])
    d1_4 = np.matmul( R0_1.T, (d0_4 - dis[0]*z0_0) ) - ais[0]*x1_1

    x, y = d1_4[:2]
    a1, a2 = ais[1], dis[3]
    a12p = (a1 + a2)**2
    a12m = (a1 - a2)**2
    xy2  = x**2 + y**2
    top = a12p - xy2
    bot = xy2 - a12m
    thetas[2] = 2*atan(sqrt( top / bot ))
    elbowup = False
    if elbowup:
        thetas[2] *= -1

    phi = atan2(y, x)
    yh = a2 * sin(thetas[2])
    xh = a1 + a2*cos(thetas[2])
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
    if thetas[3] > (pi/2):
        thetas[3] -= pi
    elif thetas[3] < (-pi/2):
        thetas[3] += pi

    yh = -R3_6[0,2]*cos(thetas[3]) - R3_6[1,2]*sin(thetas[3])
    xh = R3_6[2,2]
    thetas[4] = atan2(yh, xh)

    yh = -R3_6[0,0]*sin(thetas[3]) + R3_6[1,0]*cos(thetas[3])
    xh = -R3_6[0,1]*sin(thetas[3]) + R3_6[1,1]*cos(thetas[3])
    thetas[5] = atan2(yh, xh)

    return thetas


def main():
    ais = [69, 370.82, 0, 0, 0, 0, 0]
    dis = [270.35, 0, 0, 374.29, 0, 368.3, 0]
    Tw_tool = np.array([
      [0.984,  0.177,  0.024, 727.089],
      [0.177, -0.984,  0.020, 418.581],
      [0.027, -0.015, -0.999, 732.233],
      [  0,      0,      0,      1   ]
    ])

    above_goal = np.array([
      [1,    0,    0,  150],
      [0,   -1,    0, -800],
      [0,    0,   -1, 1100],
      [0,    0,    0,    1]
    ])
    above_object = np.array([
      [1,    0,    0,  150],
      [0,   -1,    0, -500],
      [0,    0,   -1, 1100],
      [0,    0,    0,    1]
    ])
    goal = np.array([
      [1,    0,    0,  150],
      [0,   -1,    0, -800],
      [0,    0,   -1, 1000],
      [0,    0,    0,    1]
    ])
    midpoint = np.array([
      [1,    0,    0,  150],
      [0,   -1,    0, -650],
      [0,    0,   -1, 1200],
      [0,    0,    0,    1]
    ])
    object = np.array([
      [1,    0,    0,  150],
      [0,   -1,    0, -500],
      [0,    0,   -1, 1000],
      [0,    0,    0,    1]
    ])

    thetas = ikinelbow(ais, dis, Tw_tool)
    print(thetas)


if __name__ == '__main__':
    main()
