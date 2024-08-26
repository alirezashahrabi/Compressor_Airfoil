def max_thick_loc(spline, chord, stagger, tmax_percent, tmax_loc, x0, y0):
    print( chord, stagger, tmax_percent, tmax_loc, x0, y0)
    import matplotlib.pyplot as plt
    import numpy as np
    """
    find the location and slope of suction and pressure side @ maximum thickness location  
    """
    tmax = tmax_percent * chord / 100    
    def max_thick(spline, x_suc):
        import numpy as np
        from scipy.interpolate import CubicSpline
        
        y_suc = spline(x_suc)
        
        spline_derivative = spline.derivative()
        m1 = spline_derivative(x_suc) # slope of tangent line @ x_suc
        m2 = -1/m1
        
        y_p = y_suc + np.sin(np.arctan(m2)) * tmax 
        x_p = x_suc + np.cos(np.arctan(m2)) * tmax 
        
        x_t_max = x_suc + np.cos(np.arctan(m2)) * tmax /2
        y_t_max = y_suc + np.sin(np.arctan(m2)) * tmax /2
        
        return x_t_max, y_t_max, x_p, y_p, m1
  
    axial_chord = chord * np.cos(np.deg2rad(stagger))
    print(f"axial_chord={axial_chord}")
    x_suc = x0 + tmax_loc * axial_chord / 100

    error = 100 
    x_t_max, y_t_max, x_p, y_p, slope = max_thick(spline, x_suc)
    print(x_t_max, y_t_max, x_p, y_p, slope)
    while error > 5e-1:
        x_t_max, y_t_max, x_p, y_p, slope = max_thick(spline, x_suc)
        x_suc +=  ( (x_t_max-x0)/axial_chord *100 - tmax_loc ) 
        error = abs( (x_t_max - x0) / axial_chord * 100 - tmax_loc)
        print(error)
    return x_suc, slope, x_p, y_p, x_t_max, y_t_max


def max_thick_loc_v02(spline, chord, stagger, tmax_percent, tmax_loc, x0, y0):
    """  
    Alternative method to find the location and slope of the suction and pressure   
    sides at the maximum thickness location of an airfoil profile using root finding.  

    Parameters:  
    spline (CubicSpline): Spline representing the airfoil profile.  
    chord (float): Chord length of the airfoil.  
    stagger (float): Stagger angle of the airfoil in degrees.  
    tmax_percent (float): Maximum thickness percentage of the chord (0-100).  
    tmax_loc (float): Target maximum thickness location along the chord (0-100).  
    x0 (float): X-coordinate of the starting point for the search.  
    y0 (float): Y-coordinate of the starting point for the search.  

    Returns:  
    tuple: Contains the x-location of the maximum thickness point, the slope  
           at that point, the x and y coordinates of the pressure side, and the   
           y-coordinate of the maximum thickness point.  
    """ 
    from scipy.optimize import root
    import numpy as np
    axial_chord = chord * np.cos(np.deg2rad(stagger))
    tmax = tmax_percent * chord / 100    
  
    def equations(vars):
        x_s = vars
        m1 = spline.derivative()(x_s)
        y_s = spline(x_s)
        m2 = -1/m1
        y_p = y_s + np.sin(np.arctan(m2)) * tmax 
        x_p = x_s + np.cos(np.arctan(m2)) * tmax 
            
        x_t_max = x_s + np.cos(np.arctan(m2)) * tmax /2
        y_t_max = y_s + np.sin(np.arctan(m2)) * tmax /2

        eq1 = tmax_loc - ((x_t_max-x0)/axial_chord*100)
        return eq1
    initial_guess = x0 + chord * np.cos(np.deg2rad(stagger)) * tmax_loc /100
    solution = root(equations, initial_guess)
    x_s = solution.x[0]
    m1 = spline.derivative()(x_s)
    y_s = spline(x_s)
    m2 = -1/m1
    y_p = y_s + np.sin(np.arctan(m2)) * tmax 
    x_p = x_s + np.cos(np.arctan(m2)) * tmax 
        
    x_t_max = x_s + np.cos(np.arctan(m2)) * tmax /2
    y_t_max = y_s + np.sin(np.arctan(m2)) * tmax /2
    return x_s, m1, x_p, y_p, x_t_max, y_t_max


def solve_system_camber_thickness(xs,x_points_p,spline_s, spline_p):
  
  import numpy as np
  from scipy.optimize import root

  dys = spline_s.derivative()(xs)
  ys = spline_s(xs)
  bs = -dys * xs + ys
  def equations(vars):
    xp = vars
    yp = spline_p(xp)
    dyp = spline_p.derivative()(xp)
    bp = -dyp * xp + yp
    A = np.array([[-dys.item(), 1],
                   [-dyp.item(), 1]])

    B = np.array([[bs.item()],
     [bp.item()]])
    inv_A = np.linalg.inv(A)
    x0, y0 = np.dot(inv_A, B)

    d1 = (x0-xs)**2 + (y0-ys)**2
    d2 = (x0-xp)**2 + (y0-yp)**2
    penalty = 100000 if xp < x_points_p[0] or xp > x_points_p[-1] else 0 #or abs(xp-xs)> tmax else 0


    eq1 = abs(d1 - d2) + penalty
    return eq1

  initial_guess = xs
  solution = root(equations, initial_guess)
  return solution

def cam_thick_dist(x_start,x_end,x_points_s,x_points_p,spline_s, spline_p,n=51):
  import numpy as np
  import matplotlib.pyplot as plt
  from scipy.interpolate import CubicSpline
  from circle import plot_circle
  x_LE = (x_points_s[0]+ x_points_p[0])/2
  x_TE = (x_points_s[-1]+ x_points_p[-1])/2

  t_vec=[]
  x_t_vec=[]
  x_camber = []
  y_camber = []
  x_suction = []
  x_pressure = []
  y_suction = []
  y_pressure = []
  for xs in np.linspace(x_start, x_end, n):
    solution = solve_system_camber_thickness(xs,x_points_p,spline_s, spline_p)
    if solution.success:
      xp = solution.x[0]
      yp = spline_p(xp)
      ys = spline_s(xs)
      diameter = ( (xs - xp)**2 + (ys - yp)**2 )**0.5
      t_vec.append(diameter)
      x_t_vec.append((xs+xp)/2)
      x_camber.append((xs+xp)/2)
      y_camber.append((ys+yp)/2)
      x_suction.append(xs)
      x_pressure.append(xp)
      y_suction.append(ys)
      y_pressure.append(yp)
  camber_spline = CubicSpline(np.array(x_camber), np.array(y_camber))
  camber_spline_derivative = camber_spline.derivative()
  theta_camber = np.arctan(camber_spline_derivative(np.array(x_camber)))*180/np.pi

  return np.array(t_vec),np.array(x_t_vec),np.array(x_camber),np.array(y_camber),np.array(x_suction) \
    ,np.array(x_pressure),np.array(y_suction),np.array(y_pressure), camber_spline, theta_camber
