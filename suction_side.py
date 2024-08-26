def suction_side(kappa1 : float,
                 kappa2: float,
                 wedge_angle_LE,
                 wedge_angle_TE,
                 LE_th,
                 TE_th,
                 chord,
                 stagger,
                 x_sp,
                 y_sp,
                 x0,
                 y0):
    """
    returns suction side spline and its control points
    args:
    kappa1 : inlet_metal_angle (degree)
    kappa2 : outlet_metal_angle (degree)
    wedge_angle_LE: LE Wedge angle (degree)
    wedge_angle_TE: TE Wedge angle (degree)
    LE_th: Leading Edge thickness
    TE_th: Trailing Edge thickness
    chord: chord_length
    stagger: stagger_angle (degree )
    x_sp : x position of spline points
    y_sp : y position of spline points
    x0 : x-coordinate of the origin
    y0 : y-coordinate of the origin
    return:

    spline : spline fitted to first, last and mid-points and considering wedge angles as boundary conditions for the spline.
    x_points: x-coordinate of construction point [x_LE , x_control_points, x_TE]
    y_points: y-coordinate of construction point [y_LE , y_control_points, y_TE]
    """
    import numpy as np
    from scipy.interpolate import CubicSpline
    kappa1 = np.deg2rad(kappa1)
    kappa2 = np.deg2rad(kappa2)
    wedge_angle_LE = np.deg2rad(wedge_angle_LE)
    wedge_angle_TE = np.deg2rad(wedge_angle_TE)
    stagger = np.deg2rad(stagger)
    
    x1 = x0 + chord * np.cos(stagger)
    y1 = y0 + chord * np.sin(stagger)
    
    x_LE = x0 - LE_th / 2 * np.sin(kappa1)
    y_LE = y0 + LE_th / 2 * np.cos(kappa1)
    slope_LE = np.tan(kappa1 + wedge_angle_LE)
    
    x_TE = x1 - TE_th /2 * np.sin(kappa2) 
    y_TE = y1 + TE_th /2 * np.cos(kappa2)  
    slope_TE = np.tan(kappa2 - wedge_angle_TE)
    x_points = np.concatenate(([x_LE], x_sp, [x_TE]))
    y_points = np.concatenate(([y_LE], y_sp, [y_TE]))

    spline = CubicSpline(x_points, y_points, bc_type=((1, slope_LE), (1, slope_TE)))  

    return spline, x_points, y_points 
    
