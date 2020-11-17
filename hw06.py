from math import cos, sin, sqrt, atan, atan2
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

def ikinebaxter(a1,a2,d1,d4,d6,T,LR,UD,NF,alpha=[-np.pi/2, 0, -np.pi/2, np.pi/2, -np.pi/2, 0]):
    pi = np.pi
    theta = np.zeros(6)
    R06 = T[:-1,:-1]
    d06 = T[:-1,-1]

# %   find d04
    d56 = np.array([0, 0, d6])
    d56 = np.matmul(R06, d56)
    d04 = d06 - d56

# %   find theta 1-3

    x = d04[0]
    y = d04[1]
# %   LR decides lefty or righty
    if LR == 1:
        theta[0] = atan2(y, x)
    else:
        theta[0] = atan(y / x)

    R01 = revmat(theta[0], alpha[0])

    d01 = np.array([a1*cos(theta[0]), a1*sin(theta[0]), d1])
    d14 = np.matmul(R01.T, (d04 - d01))
    x = d14[0]
    y = d14[1]

    xy = x**2 + y**2
    ad = a2**2 + d4**2
    ad2 = 2*a2*d4
    temp = (ad2 - ad + xy) / (ad2 + ad - xy)
    theta[2] = 2 * atan( sqrt( temp ) )
# %   whether you choose + or minus needs to depend on UD (elbow up or down)
    if UD == 1:
        theta[2] = -theta[2]
    theta[2] = pi - theta[2]

    phi = atan2(y, x)
    omega = atan2( (d4*sin(theta[2])), (a2 + d4*cos(theta[2])) )
    theta[1] = phi - omega

    theta[2] = (theta[2] - pi/2)

# %   find theta 4-6

    R12 = revmat(theta[1], alpha[1])
    R23 = revmat(theta[2], alpha[2])
    R03 = np.matmul(np.matmul(R01, R12), R23)
    R36 = np.matmul(R03.T, R06)

    y = R36[1,2]
    x = R36[0,2]
    if (y == 0) and (x == 0):
        if True:
            theta[3] = atan2(R36[1,0], R36[0,0])
        else:
            theta[3] = atan2(-R36[1,0], -R36[0,0])

    else:
        theta[3] = atan2(-y, -x)
        if NF == 1:
            theta[3] = -theta[3]


        y = -R36[0,2]*cos(theta[3]) - R36[1,2]*sin(theta[3])
        x = R36[2,2]
        theta[4] = atan2(y, x)

        y = -R36[0,0]*sin(theta[3]) + R36[1,0]*cos(theta[3])
        x = -R36[1,0]*sin(theta[3]) + R36[1,1]*cos(theta[3])
        theta[5] = atan2(y, x)

    return theta

def get_d06(a1,a2,d1,d4,thetas):
    xy = (
      a1 + a2*cos(thetas[1]) + d4*(
        cos(thetas[2])*sin(thetas[1]) - sin(thetas[2])*cos(thetas[1])
      )
    )

    dx = cos(thetas[0])*( xy )

    dy = sin(thetas[0])*( xy )

    dz = (
      d1 + a2*sin(thetas[1]) + d4*(
        -cos(thetas[1])*cos(thetas[2]) - sin(thetas[1])
      )
    )

    return np.round([dx, dy, dz], 2)

if __name__ == '__main__':
    pi = np.pi
    acc = 0

    ais = [0, 1, 0, 0, 0, 0]
    dis = [1, 0, 0, 1, 0, 0]
    alphas = [-pi/2, 0, -pi/2, pi/2, -pi/2, 0]

    print('test 1')
    thetas = np.zeros(6)
    mat = fkine(ais, dis, alphas, thetas)
    thets = ikinebaxter(ais[0], ais[1], dis[0], dis[3], dis[5], mat, 1, 0, 0, alphas)
    print('thetas:', np.round(thetas, 2))
    print('fkine:', np.round(mat[:-1,-1], 2))
    print('ikine:', np.round(thets, 2))
    print('d06:', [1, 0, 0])
    acc += (thets == thetas).sum() / 6
    # print('d06c:', get_d06(ais[0], ais[1], dis[0], dis[3], thetas) )


    print('\ntest 2')
    thetas = [pi/2, -pi/2, 0, 0, 0, 0]
    mat = fkine(ais, dis, alphas, thetas)
    thets = ikinebaxter(ais[0], ais[1], dis[0], dis[3], dis[5], mat, 1, 0, 0, alphas)
    print('thetas:', np.round(thetas, 2))
    print('fkine:', np.round(mat[:-1,-1], 2))
    print('ikine:', np.round(thets, 2))
    print('d06:', [0, 1, 2])
    acc += (thets == thetas).sum() / 6
    # print('d06c:', get_d06(ais[0], ais[1], dis[0], dis[3], thetas) )

    print('\ntest 3')
    thetas = [pi/4, -pi/4, 0, 0, 0, 0]
    mat = fkine(ais, dis, alphas, thetas)
    thets = ikinebaxter(ais[0], ais[1], dis[0], dis[3], dis[5], mat, 1, 0, 0, alphas)
    print('thetas:', np.round(thetas, 2))
    print('fkine:', np.round(mat[:-1,-1], 2))
    print('ikine:', np.round(thets, 2))
    print('d06:', [1, 1, 1])
    acc += (thets == thetas).sum() / 6
    # print('d06c:', get_d06(ais[0], ais[1], dis[0], dis[3], thetas) )

    print('\ntest 4')
    dis = [2, 0, 0, 1, 0, 1]
    thetas = [pi/4, -pi/4, 0, -pi/4, 0, 0]
    mat = fkine(ais, dis, alphas, thetas)
    thets = ikinebaxter(ais[0], ais[1], dis[0], dis[3], dis[5], mat, 1, 0, 0, alphas)
    print('thetas:', np.round(thetas, 2))
    print('fkine:', np.round(mat[:-1,-1], 2))
    print('ikine:', np.round(thets, 2))
    print('d06:', np.round([1.5, 1.5, 2-(sqrt(2)/2)], 2))
    acc += (thets == thetas).sum() / 6
    # print('d06c:', get_d06(ais[0], ais[1], dis[0], dis[3], thetas) )

    print('\ntest 5')
    ais = [1, 1, 0, 0, 0, 0]
    mat = fkine(ais, dis, alphas, thetas)
    thets = ikinebaxter(ais[0], ais[1], dis[0], dis[3], dis[5], mat, 1, 0, 0, alphas)
    print('thetas:', np.round(thetas, 2))
    print('fkine:', np.round(mat[:-1,-1], 2))
    print('ikine:', np.round(thets, 2))
    print('d06:', np.round([1.5+(sqrt(2)/2), 1.5+(sqrt(2)/2), 2-(sqrt(2)/2)], 2))
    acc += (thets == thetas).sum() / 6
    # print('d06c:', get_d06(ais[0], ais[1], dis[0], dis[3], thetas) )

    print('\ntest 6')
    thetas = [pi/4, -pi/4, 0, 0, pi/4, 0]
    mat = fkine(ais, dis, alphas, thetas)
    thets = ikinebaxter(ais[0], ais[1], dis[0], dis[3], dis[5], mat, 1, 0, 0, alphas)
    print('thetas:', np.round(thetas, 2))
    print('fkine:', np.round(mat[:-1,-1], 2))
    print('ikine:', np.round(thets, 2))
    print('d06:', np.round([1+(sqrt(2)/2), 1+(sqrt(2)/2), 1], 2))
    acc += (thets == thetas).sum() / 6
    # print('d06c:', get_d06(ais[0], ais[1], dis[0], dis[3], thetas) )



    print('\ntest 7')
    ais = [69, 370.82, 0, 0, 0, 0]
    dis = [270.35, 0, 0, 374.29, 0, 0]
    mat = fkine(ais, dis, alphas, thetas)
    thets = ikinebaxter(ais[0], ais[1], dis[0], dis[3], dis[5], mat, 1, 0, 0, alphas)
    print('thetas:', np.round(thetas, 2))
    print('fkine:', np.round(mat[:-1,-1], 2))
    print('ikine:', np.round(thets, 2))
    acc += (thets == thetas).sum() / 6
    # print('d06c:', get_d06(ais[0], ais[1], dis[0], dis[3], thetas) )

    print('\nactual')
    ais = [69, 370.82, 0, 0, 0, 0]
    dis = [270.35, 0, 0, 374.29, 0, 368.2]
    alphas = [-pi/2, 0, -pi/2, pi/2, -pi/2, 0]
    thetas = [pi/2, -pi/4, pi/3, pi, pi/6, pi/2]
    mat = fkine(ais, dis, alphas, thetas)
    thets = ikinebaxter(ais[0], ais[1], dis[0], dis[3], dis[5], mat, 1, 0, 1, alphas)
    print('thetas:', np.round(thetas, 2))
    print('fkine:', np.round(mat[:-1,-1], 2))
    print('ikine:', np.round(thets, 2))
    acc += (thets == thetas).sum() / 6

    print('\naccuracy:', (acc / 8))
