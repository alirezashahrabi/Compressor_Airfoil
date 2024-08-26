
def pressure_side(kappa1,
                  kappa2,
                  wedge_angle_LE,
                  wedge_angle_TE,
                  LE_th,
                  TE_th,
                  chord,
                  stagger,
                  x0,
                  y0,
                  slope_tmax,
                  x_p,
                  y_p):
    """
    returns two splines: 
    first from LE to maximum thickness location,
    second : from maximum thickness location to TE
    
    """
    import numpy as np
    from scipy.interpolate import CubicHermiteSpline
    
    kappa1 = np.deg2rad(kappa1)
    kappa2 = np.deg2rad(kappa2)
    wedge_angle_LE = np.deg2rad(wedge_angle_LE)
    wedge_angle_TE = np.deg2rad(wedge_angle_TE)
    stagger = np.deg2rad(stagger)
    
    x1 = x0 + chord * np.cos(stagger)
    y1 = y0 + chord * np.sin(stagger)
    
    x_LE = x0 + LE_th / 2 * np.sin(kappa1)
    y_LE = y0 - LE_th / 2 * np.cos(kappa1)
    
    slope_LE = np.tan(kappa1 - wedge_angle_LE)
    
    x_TE = x1 + TE_th /2 * np.sin(kappa2) 
    y_TE = y1 - TE_th /2 * np.cos(kappa2)  
    slope_TE = np.tan(kappa2 + wedge_angle_TE)
    
      
    x_points = np.array([x_LE, x_p, x_TE])
    y_points = np.array([y_LE, y_p, y_TE])

    slopes = np.array([slope_LE, slope_tmax, slope_TE])
        
    spline = CubicHermiteSpline(x_points, y_points, slopes) 
    
    return spline, x_points, y_points 
    
    
