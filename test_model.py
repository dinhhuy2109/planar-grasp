import numpy as np
from shapely.geometry import Polygon, Point
from shapely import affinity

import matplotlib.pyplot as plt
from descartes.patch import PolygonPatch
from scipy import spatial

import motionmodel
import pushedobject


poly = Polygon([(63./2.,26./2.),(-63./2.,26./2),(-63./2.,-26./2.),(63./2.,-26./2.)])
pose = [15.,0.,0.] # [x,y,theta]
num_support_points = 100
pressure_option = 'uniform'
obj = pushedobject.PushedObject('polygon',poly,pose,num_support_points,pressure_option)

contact_mu = 0.5

motion = motionmodel.MotionModel(obj, contact_mu)
v_push = np.array([0.,40.]) # mm/s
contact_point = np.array([0.,-26./2.])
contact_normal = np.array([0.,-1.])
v_o, omega, v_slip, contact_mode = motion.compute_vel_single_contact(v_push,contact_point,contact_normal,contact_mu)

fig = plt.figure(1)
ax = fig.add_subplot(111)

new_poly = affinity.rotate(obj.poly,obj.pose[2],origin='centroid',use_radians=True)
new_poly = affinity.translate(obj.poly,obj.pose[0],obj.pose[1])
patch = PolygonPatch(new_poly, alpha=0.25, zorder=2)
ax.add_patch(patch)

# ax.plot(contact_point[0],contact_point[1], 'o', color='#999999')

print v_o, omega, v_slip, contact_mode


# discrtimestep= 1e-2

# while contact_mode is 'sticking':
#     new_poly = affinity.rotate(obj.poly,obj.pose[2],origin='centroid',use_radians=True)
#     new_poly = affinity.translate(obj.poly,obj.pose[0],obj.pose[1])
#     patch = PolygonPatch(new_poly, alpha=0.25, zorder=2)
#     ax.add_patch(patch)
#     ax.plot(contact_point[0],contact_point[1], 'o', color='#999999')
#     # obj_coords = list(obj.poly.exterior.coords)
#     # tree = spatial.KDTree(points)
#     # nearest = tree.query([contact_point])
#     # print Point(contact_point).within(obj.poly)
#     print contact_mode
    
#     v_o, omega, v_slip, contact_mode = motion.compute_vel_single_contact(v_push,contact_point,contact_normal,contact_mu)
#     obj.update_pose(discrtimestep,v_o,omega)
#     contact_point, contact_normal = obj.update_contact(contact_point, contact_normal,discrtimestep,v_o,omega,v_slip)
    


# ax.set_xlim([-50,50])
# ax.set_ylim([-20,50])
# ax.set_title('Visualization')
# plt.show(False)
# raw_input()
