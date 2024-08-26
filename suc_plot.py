def suc_plot(x_points_s,
             spline_s,
             x_points_p,
             spline_p,
             ):
    
    import matplotlib.pyplot as plt
    import numpy as np
    x_LE_s = x_points_s[0]
    x_TE_s = x_points_s[-1]
    x_vec_s = np.linspace(x_LE_s, x_TE_s, 200)
    y_vec_s = spline_s(x_vec_s)

    x_LE_p = x_points_p[0]
    x_TE_p = x_points_p[-1]
    x_vec_p = np.linspace(x_LE_p, x_TE_p, 200)
    y_vec_p = spline_p(x_vec_p)
    plt.plot(x_vec_s, y_vec_s, label="Suction_side")
    plt.plot(x_vec_p, y_vec_p, label="Pressure Side")
    plt.axis("equal")
    plt.grid()
    plt.legend()
    # plt.show()

def check_plots(x_points_s, spline_s, x_points_p, spline_p,plot=1):
    import matplotlib.pyplot as plt
    import numpy as np
    x_vec_s = np.linspace(x_points_s.min(), x_points_s.max(), 300)
    spline_derv_s = spline_s.derivative()
    second_derv_s = spline_derv_s.derivative()
    y_vec_s = spline_derv_s(x_vec_s)
    y2_vec_s = second_derv_s(x_vec_s)
    R_curvature_s = ( (1 + y_vec_s**2)**1.5) / abs( y2_vec_s)

    x_vec_p = np.linspace(x_points_p.min(), x_points_p.max(), 300)
    spline_derv_p = spline_p.derivative()
    second_derv_p = spline_derv_p.derivative()
    y_vec_p = spline_derv_p(x_vec_p)
    y2_vec_p = second_derv_p(x_vec_p)
    R_curvature_p = ( (1 + y_vec_p**2)**1.5) / abs(y2_vec_p)

    plt.subplots(2,3, figsize=(15,10))
    plt.subplot(2,3,1)
    plt.plot(x_vec_s, y_vec_s)
    plt.title("Derivative of suction side")
    plt.xlabel("x",fontsize=15)
    plt.ylabel("y'",fontsize=15)
    plt.subplot(2,3,2)
    plt.plot(x_vec_s, y2_vec_s)
    plt.title("Second derivative of suction side")
    plt.xlabel("x",fontsize=15)
    plt.ylabel("y''",fontsize=15)
    plt.tight_layout()
    plt.subplot(2,3,3)
    plt.plot(x_vec_s, 1/R_curvature_s)
    plt.title("Curvature of suction side")
    plt.xlabel("x",fontsize=15)
    plt.ylabel("1/R($\kappa$)",fontsize=15)

    plt.subplot(2,3,4)

    plt.title("Derivative of Pressure side")
    plt.plot(x_vec_p, y_vec_p)
    plt.xlabel("x",fontsize=15)
    plt.ylabel("y'",fontsize=15)

    plt.subplot(2,3,5)
    plt.plot(x_vec_p, y2_vec_p)
    plt.title("Second derivative of pressure side")
    plt.xlabel("x",fontsize=15)
    plt.ylabel("y''",fontsize=15)
    plt.tight_layout()

    plt.subplot(2,3,6)
    plt.plot(x_vec_p, 1/R_curvature_p)
    plt.title("Curvature of pressure side")
    plt.xlabel("x",fontsize=15)
    plt.ylabel("1/R($\kappa$)",fontsize=15)

    plt.show()
    return x_vec_s, y_vec_s, y2_vec_s, R_curvature_s