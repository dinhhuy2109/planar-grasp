import numpy as np
from scipy.integrate import ode
from shapely.geometry import Polygon, Point, LineString
from shapely import affinity

import matplotlib.pyplot as plt
from descartes.patch import PolygonPatch
from scipy import spatial

import motionmodel
import pushedobject
import pusher
import simengine
import IPython

def vec_to_pose(vec):
    return np.array([[np.cos(vec[2]),-np.sin(vec[2]),vec[0]],
                        [np.sin(vec[2]), np.cos(vec[2]),vec[1]],
                        [0,0,1]])

IN_CONTACT = 3
V_PUSHES = 2
NORMALS = 1
POINTS = 0

poly = Polygon([(63./2.,26./2.),(-63./2.,26./2),(-63./2.,-26./2.),(63./2.,-26./2.)])
pose = [0.,0.,0.] # x,y,theta in world frame
pose = np.array([[np.cos(pose[2]),-np.sin(pose[2]),pose[0]],
              [np.sin(pose[2]), np.cos(pose[2]),pose[1]],
              [0,0,1]])
num_support_points = 10
pressure_option = 'uniform'

obj = pushedobject.PushedObject('polygon',poly,pose,num_support_points,pressure_option)

contact_mu = 0.5

motion = motionmodel.MotionModel()

v_push1 = np.array([0.,10.]) # mm/s (in world)
pusher1 = pusher.Pusher(poly=LineString([(-50,-40),(50,-40)]),pose=np.eye(3),v_push=v_push1)

v_push2 = np.array([0.,-10.]) # mm/s (in world)
pusher2 = pusher.Pusher(poly=LineString([(-50,40),(50,40)]),pose=np.eye(3),v_push=v_push2)

sim_engine = simengine.SimulationEngine(obj,pusher1,pusher2,motion,contact_mu)

# pose = [3.,10.,-np.pi/6.] # x,y,theta in world frame
# pose = np.array([[np.cos(pose[2]),-np.sin(pose[2]),pose[0]],
#               [np.sin(pose[2]), np.cos(pose[2]),pose[1]],
#               [0,0,1]])
# sim_engine.roll_out(pose, show_final_only = False, show_init = True, discrtimestep= 2e-1, t_max = 10)
# sim_engine.show(True)
              
cov = np.diag([16,9,0.01])
mean_pose = [0.,0.,0]
# corrupted_pose = np.random.multivariate_normal(mean_pose,cov)
# corrupted_pose = vec_to_pose(corrupted_pose)
# print sim_engine.roll_out(corrupted_pose,show_final_only = True, show_init = True, discrtimestep= 1e-2, t_max = 5)
    
for i in range(50):
    corrupted_pose = np.random.multivariate_normal(mean_pose,cov)
    corrupted_pose = vec_to_pose(corrupted_pose)
    sim_engine.roll_out(corrupted_pose,show_final_only = True, show_init = False, discrtimestep= 1e-2, t_max = 5)

sim_engine.show(True)
# IPython.embed()
