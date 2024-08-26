import numpy as np
from scipy.optimize import root
import matplotlib.pyplot as plt
def slope_circle(x,y,h,k):
  return -(x-h)/(y-k)
    
def solve_system_circle(x1, y1, x2, y2, dy_dx_1, dy_dx_2,spline):
    """  
    Calculate the slope of the tangent line to a circle at a given point.  

    Parameters:  
    x (float): X-coordinate of the point on the circle.  
    y (float): Y-coordinate of the point on the circle.  
    h (float): X-coordinate of the circle's center.  
    k (float): Y-coordinate of the circle's center.  

    Returns:  
    float: The slope of the tangent line to the circle at the specified point.  
    """  
    def equations(vars):
      h, k, r = vars

      eq1 = (x1-h)**2 + (y1-k)**2 - r**2
      eq2 =slope_circle(x1,y1,h,k) - (dy_dx_1)
      eq3 = (x2-h)**2 + (y2-k)**2 - r**2

      return [eq1, eq2, eq3]
  
    h_guess = (x1 + x2) / 2
    k_guess = (y1 + y2) / 2
    r_guess = ((x1-x2)**2 + (y1-y2)**2)**0.5
    initial_guess = [h_guess, k_guess, r_guess]
    solution = root(equations, initial_guess)
    return solution
  
def solve_system_TE_Circle(x_points_s,x_points_p,spline_s, spline_p):
    """  
    Solve the system to determine the parameters for a circle at the trailing edge  
    based on the provided points and their splines.  

    Parameters:  
    x_points_s (list): X-coordinates of the suction side spline control points.  
    x_points_p (list): X-coordinates of the pressure side spline control points.  
    spline_s: Spline function for the suction side.  
    spline_p: Spline function for the pressure side.  

    Returns:  
    OptimizeResult: The solution to the optimization problem including parameters  
    for the center (h, k) and the radius (R).  
    """ 
    def slope_circle(x,y,h,k):
        return -(x-h)/(y-k)
  
    x_TE = (x_points_s[-1]+x_points_p[-1])/2
    y_TE_p = spline_p(x_points_p[-1])
    y_TE_s = spline_s(x_points_s[-1])
    y_TE = (y_TE_p+y_TE_s)/2  
    TE_thickness = ((x_points_s[-1]-x_points_p[-1])**2 + (y_TE_p - y_TE_s)**2)**0.5

    def equations(vars):
        xs, xp, h, k, R = vars
        ys = spline_s(xs)
        yp = spline_p(xp)
        slope_s = spline_s.derivative()(xs)
        slope_p = spline_p.derivative()(xp)

        eq1 =(x_TE - h)**2 + (y_TE - k)**2 - R**2 
        eq2 =(xs - h)**2 + (ys - k)**2 - R**2+ 100 * ((R<0))
        eq3 =(xp - h)**2 + (yp - k)**2 - R**2
        eq4 = slope_s - slope_circle(xs,ys,h,k)
        eq5 = slope_p - slope_circle(xp,yp,h,k)
        return [eq1, eq2, eq3, eq4, eq5]
    initial_guess = [x_points_s[-1]-TE_thickness/2, x_points_p[-1]-TE_thickness/2,x_TE,y_TE,TE_thickness/2 ]
    solution = root(equations, initial_guess)
    print(f"circle_TE_solution: {solution.success}")
    if not solution.success:
        print(solution)
    return solution


    
def plot_circle(h, k, r, alpha_start=0, alpha_end=360, label="", n=15):
    """  
    Plot a circle based on its center and radius.  

    Parameters:  
    h (float): X-coordinate of the circle's center.  
    k (float): Y-coordinate of the circle's center.  
    r (float): Radius of the circle.  
    alpha_start (float): Starting angle in degrees for the plot (default is 0).  
    alpha_end (float): Ending angle in degrees for the plot (default is 360).  
    label (str): Label for the plot legend (default is empty string).  
    n (int): Number of points to plot the circle (default is 15).  
    """ 
    alpha_vec = np.linspace(alpha_start, alpha_end, n)
    x_circle = r * np.cos(alpha_vec*np.pi/180) + h
    y_circle = r * np.sin(alpha_vec*np.pi/180) + k
    plt.plot(x_circle,y_circle,c="k",label=label)
    if label != "":
        plt.legend()

def find_circle_start_end(x_s, x_p, spline_s, spline_p, h, k, r):
    """  
    Find the starting and ending angles of the circle based on the given points.  

    Parameters:  
    x_s (float): X-coordinate of a point on the suction side.  
    x_p (float): X-coordinate of a point on the pressure side.  
    spline_s: Spline function for the suction side.  
    spline_p: Spline function for the pressure side.  
    h (float): X-coordinate of the circle's center.  
    k (float): Y-coordinate of the circle's center.  
    r (float): Radius of the circle.  

    Returns:  
    tuple: Minimum and maximum angles (in degrees) corresponding to the start and end of the circle.  
    """ 
    y_s = spline_s(x_s)
    y_p = spline_p(x_p)
    alpha_start_c = ( np.arctan2((y_p - k), (x_p - h) )+2*np.pi ) *180/np.pi
    alpha_end_c = ( np.arctan2((y_s - k), (x_s - h )  ) +2*np.pi) *180/np.pi
    return min(alpha_start_c, alpha_end_c), max(alpha_start_c, alpha_end_c)



