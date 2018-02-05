import numpy as np
from shapely.geometry import Polygon, Point

class PushedObject(object):
    """
    Pushed object
    """
    def __init__(self,object_type,geometry):
        self.object_type = object_type
        if self.object_type == 'polygon':
            self.poly = geometry #poly
        self.support_points = []
        num_support_points = 100
        grid_support_points(num_support_points,'polygon')
        pressure_option = 'uniform'
        self.presure_weights = assign_pressure(support_points,pressure_option)
        self.ls_constant = compute_limit_surface_coeffs()
        
        print 'init DONE'

    def assign_pressure(self,support_points,pressure_option):
        if pressure_option == 'uniform':
            pressure_weights = np.ones(len(support_points))
            pressure_weights = pressure_weights/sum(pressure_weights)
            return pressure_weights
        
    
    def grid_support_points(self,num_points,option):
        if option == 'polygon':
            obj_area = self.poly.area
            
            min_x,min_y,max_x,max_y = self.poly.bounds # enclosing rectangle
            rec_area = (max_x - min_x)*(max_y-min_y)

            grid_size = np.sqrt(object_area/num_points)
            rec_edge_x = np.arange(min_x,max_x,grid_size)
            rec_edge_y = np.arange(min_y,max_y,grid_size)
            for i in rec_edge_x:
                for j in rec_edge_y:
                    grid_center = Point(rec_edge_x[i],rec_edge_y[j])
                    if grid_center.within(self.poly):
                        self.support_points.append(np.asarray(pt.xy).reshape(2))
            grid = [option,[self.support_points,grid_size]]
            return grid
        # if option == 'rim':
        #     grid = [option]
        #     return 

    def compute_limit_surface_coeffs(self):
        c = 0# tau_max/f_max
        num_support_points = len(self.support_points)
        for pt in self.support_points:
            c += np.linalg.norm(pt)
        c = c/num_support_points
        return c
