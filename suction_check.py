
def suction_check(x_points, spline):
    
    """
    check the number of change in  derivitation of suction side spline sign
    args: 
        spline
    return:
        change direction of curve derivative
    """
    import numpy as np
    from scipy.interpolate import CubicSpline
    spline_derivative = spline.derivative()  # This returns a new CubicSpline object for the derivative
    
    x_LE = x_points[0]
    x_TE = x_points[-1]
    
    x_vec = np.linspace(x_LE, x_TE, 201)
    y_derivative = spline_derivative(x_vec)
    dy_diff = np.diff(y_derivative)
    signs = np.sign(dy_diff)  
    sign_changes = np.sum(np.abs(np.diff(signs)) > 0) 
    return sign_changes
