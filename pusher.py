import numpy as np
from shapely.geometry import Polygon, Point, LineString
import IPython

class Pusher(object):
    """
    Pushed object
    """
    def __init__(self,poly,pose,v_push):
        self.poly = poly
        self.points = np.asarray(self.poly.coords)[:2]
        self.pose = pose
        self.v_push = v_push
        vec = self.points[0]-self.points[1]
        self.normal = np.array([vec[1],-vec[0]])
        if np.dot(self.normal,v_push)>0:
            self.normal = self.normal
        else:
            self.normal = -self.normal
        self.normal = self.normal/np.linalg.norm(self.normal)

    def eval_pusher(self,time):
        delta = self.v_push*time
        new_pts = []
        for pt in self.points:
            new_pts.append(pt + delta)
        return LineString(new_pts)

    def update_pusher(self,pusher):
        self.poly = pusher
        self.points = np.asarray(self.poly.coords)[:2]
        return True
    
