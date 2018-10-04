import numpy as np
from shapely.geometry import Polygon, Point, LineString
from shapely import affinity

import matplotlib.pyplot as plt
from descartes.patch import PolygonPatch
from scipy import spatial

import motionmodel
import pushedobject
import pusher
import IPython

poly = Polygon([(63./2.,26./2.),(-63./2.,26./2),(-63./2.,-26./2.),(63./2.,-26./2.)])
pose = [0.,0.,-np.pi/3.] #x,y,theta in world frame
pose = np.array([[np.cos(pose[2]),-np.sin(pose[2]),pose[0]],
              [np.sin(pose[2]), np.cos(pose[2]),pose[1]],
              [0,0,1]])
num_support_points = 10
pressure_option = 'uniform'

obj = pushedobject.PushedObject('polygon',poly,pose,num_support_points,pressure_option)

contact_mu = 0.5

motion = motionmodel.MotionModel()

# #test 1 translation
# A = np.eye(3)
# v_pushes = np.array([[1,0.2],[1,0.2]])
# contact_points = np.array([[-1,0.5],[-1,-0.5]])
# contact_normals = np.array([[-1,0],[-1,0]])

#test 2 translate and rotate
A = np.eye(3)
v_pushes = np.array([[1,0.],[0.5,0.]])
contact_points = np.array([[-1,0.5],[-1,-0.5]])
contact_normals = np.array([[-1,0],[-1,0]])

# #test3 jamming when increase mu to tan(pi/6)
# A = np.eye(3)
# v_pushes = np.array([[0.,-1],[np.sqrt(3)/2,0.5]])
# contact_points = np.array([[0.,1],[-np.sqrt(3)/2,-0.5]])
# contact_normals = np.array([[0.,1.],[-np.sqrt(3)/2,-0.5]])
# # contact_mu = np.tan(np.pi/6)
print motion.compute_vel_multi_contacts(v_pushes, contact_points,contact_normals,contact_mu,A)

