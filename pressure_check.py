def pressure_check(x_points, spline):
    
    """
    args: 
        spline
    return:
        change direction of curve derivative
    """
    import numpy as np
    from scipy.interpolate import CubicSpline
    spline_derivative = spline.derivative()  # This returns a new CubicSpline object for the derivative
    

    
    x_vec = np.linspace(x_points.min(), x_points.max(), 101)
    y_derivative = spline_derivative(x_vec)
    
    dy_diff = np.diff(y_derivative)
    signs = np.sign(dy_diff)  
    sign_changes = np.sum(np.abs(np.diff(signs)) > 0) 
    return sign_changes
