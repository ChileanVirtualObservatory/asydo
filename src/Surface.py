import math
import numpy as np

def genSurface(form,alpha,delta,alpha_axis,delta_axis):
        "Create a gaussian surface over a mesh created by x and y axes"
        sx = form[1]
        sy = form[2]
        theta = form[3]
        r = 3 * math.sqrt(sx ** 2 + sy ** 2)
        alpha_mesh, delta_mesh = np.meshgrid(alpha_axis, delta_axis, sparse = False, indexing = 'xy')
        Xc = alpha_mesh.flatten() - alpha * np.ones(len(alpha_axis) * len(delta_axis))
        Yc = delta_mesh.flatten() - delta * np.ones(len(alpha_axis) * len(delta_axis))
        XX = (Xc) * math.cos(theta) - (Yc) * math.sin(theta);
        YY = (Xc) * math.sin(theta) + (Yc) * math.cos(theta);
        u = (XX / sx) ** 2 + (YY / sy) ** 2;
        sol = sx * sy * np.exp(-u / 2) / (2 * math.pi);
        res = np.transpose(np.reshape(sol, (len(delta_axis), len(delta_axis))))
        res = res / res.max()
        return res,[[0, len(alpha_axis)],[0,len(delta_axis)]]


