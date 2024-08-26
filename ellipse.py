import numpy as np
from scipy.optimize import root


def slope_ellipse(x, y, h, k, a, b, theta):
    """  
    Calculate the slope of the tangent line to an ellipse at a given point.  

    Parameters:  
    x (float): X-coordinate of the point on the ellipse.  
    y (float): Y-coordinate of the point on the ellipse.  
    h (float): X-coordinate of the ellipse's center.  
    k (float): Y-coordinate of the ellipse's center.  
    a (float): Semi-major axis length of the ellipse.  
    b (float): Semi-minor axis length of the ellipse.  
    theta (float): Rotation angle of the ellipse in radians.  

    Returns:  
    float: The slope of the tangent line to the ellipse at the specified point.  
           If the denominator is zero, it returns ±inf based on the sign of the numerator.  
    """  
    import numpy as np
    u = (x - h) * np.cos(theta) + (y - k) * np.sin(theta)
    v = - (x - h) * np.sin(theta) + (y - k) * np.cos(theta)
    numerator = - b**2 * u * np.cos(theta) + a**2 * v * np.sin(theta)
    denominator = b**2 * u * np.sin(theta) + a**2 * v * np.cos(theta)
    if denominator != 0:
        return numerator/denominator
    elif numerator > 0:
        return np.inf
    else:
        return -np.inf


def solve_system(x1, y1, x2, y2, dy_dx_1, dy_dx_2, er):
    def equations(vars):
        h, k, a, theta = vars
        # er = a/b # Ellipse ratio
        a_limit = (((x1-x2)**2 + (y1-y2)**2)**0.5)/2
        b = a/er

        u = (x1-h)*np.cos(theta) + (y1-k)*np.sin(theta)
        v = -(x1-h)*np.sin(theta) + (y1-k)*np.cos(theta)

        eq1 = ((u/a)**2 + (v/b)**2 - 1)*2 + (a > 5*a_limit) * 100
        u = (x2-h)*np.cos(theta) + (y2-k)*np.sin(theta)
        v = -(x2-h)*np.sin(theta) + (y2-k)*np.cos(theta)

        eq2 = ((u/a)**2 + (v/b)**2 - 1)*2

        eq3 = dy_dx_1 - slope_ellipse(x1, y1, h, k, a, b, theta)
        eq4 = dy_dx_2 - slope_ellipse(x2, y2, h, k, a, b, theta)

        return [eq1, eq2, eq3, eq4]
        
    a_guess = (((x1-x2)**2 + (y1-y2)**2)**0.5)/2 * er
    h_guess = (x1 + x2) / 2 * 0.9
    k_guess = (y1 + y2) / 2
    theta_guess = (dy_dx_1 + dy_dx_2)/2
    initial_guess = [h_guess, k_guess, a_guess, theta_guess]
    solution = root(equations, initial_guess)
    # print(solution.success)
    return solution


def plot_ellipse(h, k, a, b, theta, alpha_start=0, alpha_end=360, title="", n=15):
    """  
    Plot an ellipse based on its parameters.  

    Parameters:  
    h (float): X-coordinate of the ellipse's center.  
    k (float): Y-coordinate of the ellipse's center.  
    a (float): Semi-major axis length of the ellipse.  
    b (float): Semi-minor axis length of the ellipse.  
    theta (float): Rotation angle of the ellipse in radians.  
    alpha_start (float): Starting angle (in degrees) for the plot (default is 0).  
    alpha_end (float): Ending angle (in degrees) for the plot (default is 360).  
    title (str): Title for the plot (default is an empty string).  
    n (int): Number of points used to plot the ellipse (default is 15).  

    Returns:  
    None  
    """ 
    import numpy as np
    import matplotlib.pyplot as plt
    alpha = np.deg2rad(np.linspace(alpha_start, alpha_end, n))
    x_vec = h + a * np.cos(alpha) * np.cos(theta) - b * np.sin(alpha) * np.sin(theta)
    y_vec = k + a * np.cos(alpha) * np.sin(theta) + b * np.sin(alpha) * np.cos(theta)
    plt.plot(x_vec, y_vec, label=title)
    plt.axis("equal")
    # plt.axhline(k, xmin=0.1 ,xmax=0.9,ls='--')
    # plt.axvline(h, ymin=0.05 ,ymax=0.95,ls='--')

    if title != "":
        plt.legend()


def solve_system_LE_ellipse(x_points_s, x_points_p, spline_s, spline_p, theta, er):
    """  
    Solve for the parameters of an ellipse at the leading edge based on the   
    provided points and their splines.  

    Parameters:  
    x_points_s (list): X-coordinates of the suction side points.  
    x_points_p (list): X-coordinates of the pressure side points.  
    spline_s: Spline function for the suction side.  
    spline_p: Spline function for the pressure side.  
    theta (float): Angle of the ellipse in degrees.  
    er (float): Ellipse ratio (a/b).  

    Returns:  
    OptimizeResult: The solution of the optimization problem including parameters  
    for the center (h, k) and the semi-major axis (a).  
    """  
    ## theta in degree (input)
    theta = (np.deg2rad(theta)) # convert to radian
    # slope_normal_to_LE = -1/slope_LE

    x_LE_s = x_points_s[0]
    y_LE_s = spline_s(x_LE_s)
    x_LE_p = x_points_p[0]
    y_LE_p = spline_p(x_LE_p)
    xLE = (x_LE_s+x_LE_p)/2
    yLE = (y_LE_s+y_LE_p)/2
    LE_thickness = ( (x_LE_s - x_LE_p)**2 + (y_LE_s - y_LE_p)**2 )**0.5

  
    def equations(vars):
        xs, xp, h, k, a = vars
        # er = a/b # Ellipse ratio
        b = a/er
        ys = spline_s(xs)
        yp = spline_p(xp)
        dy_dx_s = spline_s.derivative()(xs)
        dy_dx_p = spline_p.derivative()(xp)
        slope_LE = np.tan(theta)
        slope_normal_to_LE = -1/slope_LE

        u = (xs-h)*np.cos(theta) + (ys-k)*np.sin(theta)
        v = -(xs-h)*np.sin(theta) + (ys-k)*np.cos(theta)
        eq1 = ((u/a)**2 + (v/b)**2 - 1)

        u = (xp-h)*np.cos(theta) + (yp-k)*np.sin(theta)
        v = -(xp-h)*np.sin(theta) + (yp-k)*np.cos(theta)
        eq2 = ((u/a)**2 + (v/b)**2 - 1)

        u = (xLE-h)*np.cos(theta) + (yLE-k)*np.sin(theta)
        v = -(xLE-h)*np.sin(theta) + (yLE-k)*np.cos(theta)
        eq3 = ((u/a)**2 + (v/b)**2 - 1)

        eq4 = dy_dx_s - slope_ellipse(xs, ys, h, k, a, b, theta)
        eq5 = dy_dx_p - slope_ellipse(xp, yp, h, k, a, b, theta)

        return [eq1, eq2, eq3, eq4, eq5]
    
    xs_guess = x_points_s[0] + LE_thickness / 2 * np.cos(theta) * 2 * np.random.rand()
    xp_guess = x_points_p[0] + LE_thickness / 2 * np.cos(theta) * 2 * np.random.rand()
    a_guess = LE_thickness/2*(1 + er * 2 * np.random.rand())
    h_guess = xLE + LE_thickness/2 * er * np.cos(theta)
    k_guess = yLE + LE_thickness/2 * np.sin(theta)
    initial_guess = [xs_guess, xp_guess, h_guess, k_guess, a_guess]
    solution = root(equations, initial_guess)
    # print(solution.success)
    return solution

def find_ellipse_start_end(h, k, a, b, theta, xs, ys, xp, yp, geo="LE"):
    """  
    Find the starting and ending angles of the ellipse based on the given points.  

    Parameters:  
    h (float): X-coordinate of the ellipse's center.  
    k (float): Y-coordinate of the ellipse's center.  
    a (float): Semi-major axis length of the ellipse.  
    b (float): Semi-minor axis length of the ellipse.  
    theta (float): Rotation angle of the ellipse in radians.  
    xs (float): X-coordinate of a point on the suction side.  
    ys (float): Y-coordinate of a point on the suction side.  
    xp (float): X-coordinate of a point on the pressure side.  
    yp (float): Y-coordinate of a point on the pressure side.  
    geo (str): Geometry type ("LE" or "TE", default is "LE").  

    Returns:  
    tuple: Minimum and maximum angles (in degrees) corresponding to the start and end of the ellipse.  
    """

    a11 = a * np.cos(theta)
    a12 = - b * np.sin(theta)
    a21 = a * np.sin(theta)
    a22 = b * np.cos(theta)
    A = np.array([[a11, a12],
              [a21, a22]])
    inv_A = np.linalg.inv(A)
    B = np.array([[xs - h],
              [ys - k]])
    cos_alpha_s, sin_alpha_s = np.dot(inv_A, B)
#    match geo:
    if geo == "LE":     #case "LE":
            alpha_s = (np.arctan2(sin_alpha_s, cos_alpha_s) + 2 * np.pi * (np.arctan2(sin_alpha_s, cos_alpha_s) < 0)) * 180/np.pi
    else:
  #      case "TE":
            alpha_s = (np.arctan2(sin_alpha_s, cos_alpha_s) + 2 * np.pi) * 180/np.pi
    B = np.array([[xp - h],
              [yp - k]])
    cos_alpha_p, sin_alpha_p = np.dot(inv_A, B)
 #   match geo:
  #      case "LE":
    if geo == "LE":
            alpha_p = (np.arctan2(sin_alpha_p, cos_alpha_p) + 2 * np.pi * (np.arctan2(sin_alpha_p, cos_alpha_p)<0)) *180/np.pi
            return alpha_s, alpha_p
  #      case "TE":
    else:
            alpha_p = (np.arctan2(sin_alpha_p, cos_alpha_p) + 2 * np.pi ) *180/np.pi
            return min(alpha_s, alpha_p), max(alpha_s, alpha_p)
  



def solve_system_TE_ellipse(x_points_s, x_points_p, spline_s, spline_p, theta, er):
    """  
    Solve for the parameters of an ellipse at the trailing edge based on the   
    provided points and their splines.  

    Parameters:  
    x_points_s (list): X-coordinates of the suction side points.  
    x_points_p (list): X-coordinates of the pressure side points.  
    spline_s: Spline function for the suction side.  
    spline_p: Spline function for the pressure side.  
    theta (float): Angle of the ellipse in degrees.  
    er (float): Ellipse ratio (a/b).  

    Returns:  
    OptimizeResult: The solution of the optimization problem including parameters for the ellipse.  
    """
    ## theta in degree (input)
    theta = (np.deg2rad(theta)) # convert to radian
    # slope_normal_to_LE = -1/slope_LE

    x_TE_s = x_points_s[-1]
    y_TE_s = spline_s(x_TE_s)
    x_TE_p = x_points_p[-1]
    y_TE_p = spline_p(x_TE_p)
    xTE = (x_TE_s+x_TE_p)/2
    yTE = (y_TE_s+y_TE_p)/2
    TE_thickness = ( (x_TE_s - x_TE_p)**2 + (y_TE_s - y_TE_p)**2 )**0.5

  
    def equations(vars):
        xs, xp, h, k, a = vars
        # er = a/b # Ellipse ratio
        b = a/er
        ys = spline_s(xs)
        yp = spline_p(xp)
        dy_dx_s = spline_s.derivative()(xs)
        dy_dx_p = spline_p.derivative()(xp)
        slope_TE = np.tan(theta)
        slope_normal_to_TE = -1/slope_TE

        u =  (xs-h)*np.cos(theta) + (ys-k)*np.sin(theta)
        v = -(xs-h)*np.sin(theta) + (ys-k)*np.cos(theta)
        eq1 = ((u/a)**2 + (v/b)**2 - 1) + (a < 0) * 100

        u =  (xp-h)*np.cos(theta) + (yp-k)*np.sin(theta)
        v = -(xp-h)*np.sin(theta) + (yp-k)*np.cos(theta)
        eq2 = ((u/a)**2 + (v/b)**2 - 1)

        u =  (xTE-h)*np.cos(theta) + (yTE-k)*np.sin(theta)
        v = -(xTE-h)*np.sin(theta) + (yTE-k)*np.cos(theta)
        eq3 = ((u/a)**2 + (v/b)**2 - 1)

        eq4 = dy_dx_s - slope_ellipse(xs, ys, h, k, a, b, theta)
        eq5 = dy_dx_p - slope_ellipse(xp, yp, h, k, a, b, theta)

        return [eq1, eq2, eq3, eq4, eq5]
    
    xs_guess = x_points_s[-1] - TE_thickness/2  * np.cos(theta) * 2 * np.random.rand()
    xp_guess = x_points_p[-1] - TE_thickness/2  * np.cos(theta)  * 2 * np.random.rand()
    a_guess = TE_thickness/2*(1+ er * 2 *np.random.rand())
    h_guess = xTE - TE_thickness/2 * er * np.cos(theta)
    k_guess = yTE - TE_thickness/2 * np.sin(theta)
    initial_guess = [xs_guess,xp_guess,h_guess, k_guess, a_guess]
    solution = root(equations, initial_guess)
    # print(solution.success)
    return solution


def ellipse_points(h, k, a, b, theta, alpha_start=0, alpha_end=360, n=15):
    """  
    Generate points on an ellipse defined by its parameters.  

    Parameters:  
    h (float): X-coordinate of the ellipse's center.  
    k (float): Y-coordinate of the ellipse's center.  
    a (float): Semi-major axis length of the ellipse.  
    b (float): Semi-minor axis length of the ellipse.  
    theta (float): Rotation angle of the ellipse in radians.  
    alpha_start (float): Starting angle (in degrees) for generating points (default is 0).  
    alpha_end (float): Ending angle (in degrees) for generating points (default is 360).  
    n (int): Number of points to generate (default is 15).  

    Returns:  
    tuple: X-coordinates and Y-coordinates of the points on the ellipse.  
    """ 
    import numpy as np
    alpha = np.deg2rad(np.linspace(alpha_start, alpha_end, n))
    x_vec = h + a * np.cos(alpha) * np.cos(theta) - b * np.sin(alpha) * np.sin(theta)
    y_vec = k + a * np.cos(alpha) * np.sin(theta) + b * np.sin(alpha) * np.cos(theta)
    return x_vec, y_vec