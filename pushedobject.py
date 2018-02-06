import numpy as np
from shapely.geometry import Polygon, Point


class PushedObject(object):
    """
    Pushed object
    """
    def __init__(self,object_type,poly,pose,num_support_points,pressure_option):
        self.object_type = object_type
        if self.object_type is 'polygon':
            self.pose = pose
            self.poly = poly
            self.support_points = []
            self.num_support_points = num_support_points
            self.grid_support_points(num_support_points,'polygon')
            self.pressure_option = pressure_option
            self.presure_weights = self.assign_pressure(self.support_points)
            self.ls_coeff = self.compute_limit_surface_coeffs()
        
        print 'Object init DONE'

    def update_pose(self,  discrtimestep, v_o, omega):
        displacements = discrtimestep*np.asarray([v_o[0],v_o[1],omega])
        self.pose = self.pose + displacements

    def update_contact(self, contact_point, contact_normal, discrtimestep, v_o, omega, v_slip):
        new_contact_point = contact_point + discrtimestep*np.asarray(v_slip)
        theta = self.pose[2]
        T = np.array([[np.cos(theta),-np.sin(theta)],
                      [np.sin(theta), np.cos(theta)]])
        new_contact_normal = np.dot(T,contact_normal)
        return new_contact_point, new_contact_normal
    
    def assign_pressure(self,support_points):
        if self.pressure_option is 'uniform':
            pressure_weights = np.ones(len(support_points))
            pressure_weights = pressure_weights/sum(pressure_weights)
            return pressure_weights        
    
    def grid_support_points(self,num_points,option):
        if option is 'polygon':
            object_area = self.poly.area
            
            min_x,min_y,max_x,max_y = self.poly.bounds # enclosing rectangle
            rec_area = (max_x - min_x)*(max_y-min_y)

            grid_size = np.sqrt(object_area/num_points)
            rec_edge_x = np.arange(min_x,max_x,grid_size)
            rec_edge_y = np.arange(min_y,max_y,grid_size)
            for x in rec_edge_x:
                for y in rec_edge_y:
                    grid_center = Point(x,y)
                    if grid_center.within(self.poly):
                        self.support_points.append(np.asarray(grid_center.xy).reshape(2))
            grid = [option,[self.support_points,grid_size]]
            return grid
        # if option == 'rim':
        #     grid = [option]
        #     return 

    def compute_limit_surface_coeffs(self):
        c = 0 # tau_max/f_max
        num_support_points = len(self.support_points)
        for pt in self.support_points:
            c += np.linalg.norm(pt)
        c = c/num_support_points
        return c
