import math
import numpy as np

def spatialWindow(alpha,delta,spx,spy,alpha_axis,delta_axis):
    """ Frequency window.
       """
    Dlta_x=alpha_axis[1] - alpha_axis[0]
    Dlta_y=delta_axis[1] - delta_axis[0]
    xbord=[int(round((alpha - spx - alpha_axis[0])/Dlta_x)), int(round((alpha + spx - alpha_axis[0])/Dlta_x))]
    ybord=[int(round((delta - spy - delta_axis[0])/Dlta_y)),int(round((delta + spy - delta_axis[0])/Dlta_y))]
    if xbord[0] < 0:
       xbord[0]=0
    if ybord[0] < 0:
       ybord[0]=0
    if xbord[1] > len(alpha_axis):
       xbord[1]=len(alpha_axis)
    if ybord[1] > len(delta_axis):
       ybord[1]=len(delta_axis)
    return (xbord,ybord)


def genSurface(form,alpha,delta,alpha_axis,delta_axis):
    "Create a gaussian surface over a mesh created by x and y axes"
    sx = form[1]
    sy = form[2]
    theta = form[3]
    spx= abs(3*sx*math.cos(theta)) + abs(3*sy*math.sin(theta))
    spy= abs(3*sx*math.sin(theta)) + abs(3*sy*math.cos(theta))
    xbord,ybord=spatialWindow(alpha,delta,spx,spy,alpha_axis,delta_axis)
    alpha_axis=alpha_axis[xbord[0]:xbord[1]]
    delta_axis=delta_axis[ybord[0]:ybord[1]]
#    r = 3 * math.sqrt(sx ** 2 + sy ** 2)
    alpha_mesh, delta_mesh = np.meshgrid(alpha_axis, delta_axis, sparse = False, indexing = 'xy')
    Xc = alpha_mesh.flatten() - alpha * np.ones(len(alpha_axis) * len(delta_axis))
    Yc = delta_mesh.flatten() - delta * np.ones(len(alpha_axis) * len(delta_axis))
    XX = (Xc) * math.cos(theta) - (Yc) * math.sin(theta);
    YY = (Xc) * math.sin(theta) + (Yc) * math.cos(theta);
    u = (XX / sx) ** 2 + (YY / sy) ** 2;
    sol = sx * sy * np.exp(-u / 2) / (2 * math.pi);
    #res = np.transpose(np.reshape(sol, (len(alpha_axis), len(delta_axis))))
    res = np.transpose(np.reshape(sol, (len(delta_axis), len(alpha_axis))))
    res = res / res.max()
    return res,[xbord,ybord]


