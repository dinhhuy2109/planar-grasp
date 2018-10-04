import numpy as np

class SE2Algebra(object):

    def  __init__(self):

    def twist_matrix(vec):
        return np.array([[0     ,-vec[2],vec[0]],
                         [vec[2],0      ,vec[1]],
                         [0     ,      0,1     ]])

    def exponetial_map(vec):
        theta =  vec[2]
        if theta > 1e-8:
            V = (1./theta)*np.array([[np.sin(theta), -(1-np.cos(theta))],
                                     [1-np.cos(theta), np.sin(theta)]])
            R = np.array([[np.cos(theta),-np.sin(theta)],
                          [np.sin(theta), np.cos(theta)]])
            T = np.hstack((R,np.dot(V,vec[:2])))
            T = np.vstack((T,np.array([0,0,1])))
        else:
            T = np.eye(3)
            T[:2,2] = vec[:2]
        return T

    
