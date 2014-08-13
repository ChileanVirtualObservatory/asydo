import math
import numpy as np

def genGradient(form,alpha,delta,alpha_axis,delta_axis,Tbord):
        "Create a gradient matrix"
        gtype = form[0]
        theta = form[1]
        km_sarcs = form[2]
        xbord=Tbord[0]
        ybord=Tbord[1]
        alpha_axis=alpha_axis[xbord[0]:xbord[1]]
        delta_axis=delta_axis[ybord[0]:ybord[1]]
        alpha_mesh, delta_mesh = np.meshgrid(alpha_axis, delta_axis, sparse = False, indexing = 'xy')
        Xc = alpha_mesh.flatten() - alpha * np.ones(len(alpha_axis) * len(delta_axis))
        Yc = delta_mesh.flatten() - delta * np.ones(len(alpha_axis) * len(delta_axis))
        XX = (Xc) * math.cos(theta) - (Yc) * math.sin(theta);
        YY = (Xc) * math.sin(theta) + (Yc) * math.cos(theta);
        u = km_sarcs*(XX + YY)
        res = np.transpose(np.reshape(u, (len(alpha_axis), len(delta_axis))))
        return res


