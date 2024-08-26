"""  
BladeGenClass.py  

This module contains the definition of the CompressorAirfoil class, which is designed to model and analyze  
compressor airfoil profiles. The class utilizes various mathematical and scientific libraries, including  
NumPy for numerical operations, Pandas for data manipulation, and SciPy for optimization and interpolation.  

Key functionalities include:  
- Calculation of airfoil geometries based on input parameters.  
- Implementation of aerodynamic properties and analyses.  
- Visualization of airfoil shapes and characteristics through plotting functions.  

Dependencies:  
- NumPy  
- Pandas  
- SciPy  
- Matplotlib  
- Shapely.geometry  

To use this module, import the CompressorAirfoil class and initialize it with required parameters such as  
kappa values, chord length, stagger angle, and leading edge wedge angle. The class methods will then provide  
tools to analyze the airfoil's aerodynamic performance and geometrical properties.  

Author: Alireza Shahrabi  
Date: August 25, 2024  
"""
import os
import sys
import numpy as np
import pandas as pd
from shapely.geometry import Polygon
from scipy.optimize import root
from scipy.interpolate import CubicSpline, CubicHermiteSpline
import matplotlib.pyplot as plt
from suction_side import suction_side
from suction_check import suction_check
from max_thick_loc import max_thick_loc, solve_system_camber_thickness, cam_thick_dist, max_thick_loc_v02
from pressure_side import pressure_side
from suc_plot import suc_plot
from ellipse import plot_ellipse, slope_ellipse, solve_system, solve_system_LE_ellipse,find_ellipse_start_end,solve_system_TE_ellipse, ellipse_points
from pressure_check import pressure_check
from circle import plot_circle, solve_system_circle,solve_system_TE_Circle, find_circle_start_end


class CompressorAirfoil:
    def __init__(
            self, kappa1, kappa2, chord, stagger, wedge_angle_LE,
            wedge_angle_TE, tmax, tmax_loc, LE_th, TE_th, er_LE,
            er_TE, x_s_sp, y_s_sp, x0, y0):
        self.kappa1 = kappa1
        self.kappa2 = kappa2
        self.chord = chord
        self.stagger = stagger
        self.wedge_angle_LE = wedge_angle_LE
        self.wedge_angle_TE = wedge_angle_TE
        self.tmax = tmax
        self.tmax_loc = tmax_loc
        self.LE_th = LE_th
        self.TE_th = TE_th
        self.x_s_sp = x_s_sp
        self.y_s_sp = y_s_sp
        self.er_LE = er_LE
        self.er_TE = er_TE
        self.x0 = x0
        self.y0 = y0
        self.LE_thickness = LE_th * chord / 100
        self.TE_thickness = TE_th * chord / 100
        self.chord_axial = self.chord * np.cos(np.deg2rad(self.stagger))
        self.parameters = np.array([self.kappa1,
                                    self.kappa2,
                                    self.chord,
                                    self.stagger,
                                    self.wedge_angle_LE,
                                    self.wedge_angle_TE,
                                    self.tmax,
                                    self.tmax_loc,
                                    self.LE_th,
                                    self.TE_th,
                                    self.er_LE,
                                    self.er_TE,
                                    self.x_s_sp[0],
                                    self.x_s_sp[1],
                                    self.y_s_sp[0],
                                    self.y_s_sp[1],
                                    self.x0,
                                    self.y0])

    def BladeGen(self):
        self.spline_s, self.x_points_s, self.y_points_s = suction_side(
            kappa1=self.kappa1,
            kappa2=self.kappa2,
            wedge_angle_LE=self.wedge_angle_LE,
            wedge_angle_TE=self.wedge_angle_TE,
            LE_th=self.LE_thickness,
            TE_th=self.TE_thickness,
            chord=self.chord,
            stagger=self.stagger,
            x_sp=self.x_s_sp,
            y_sp=self.y_s_sp,
            x0=self.x0,
            y0=self.y0)
        self.suction_sign_change = suction_check(
            self.x_points_s,
            self.spline_s)
        self.x_suc_max_th, self.slope_max_th, self.x_p_max_th, self.y_p_max_th, self.x_t_max, self.y_t_max = max_thick_loc_v02(
            spline=self.spline_s,
            chord=self.chord,
            stagger=self.stagger,
            tmax_percent=self.tmax,
            tmax_loc=self.tmax_loc,
            x0=self.x0,
            y0=self.y0)
        self.spline_p, self.x_points_p, self.y_points_p = pressure_side(
            kappa1=self.kappa1,
            kappa2=self.kappa2,
            wedge_angle_LE=self.wedge_angle_LE,
            wedge_angle_TE=self.wedge_angle_TE,
            LE_th=self.LE_thickness,
            TE_th=self.TE_thickness,
            chord=self.chord,
            stagger=self.stagger,
            x0=self.x0,
            y0=self.y0,
            slope_tmax=self.slope_max_th,
            x_p=self.x_p_max_th,
            y_p=self.y_p_max_th)
        
        self.pressure_sign_change = pressure_check(
            self.x_points_p, self.spline_p)
    
    def calculate_camber_thick_distribution(self):
        self.t_vec,self.x_t_vec,self.x_camber,self.y_camber,self.x_suction, \
            self.x_pressure,self.y_suction,self.y_pressure,self.camber_spline, self.theta_camber  = cam_thick_dist(
            self.x_points_s[0],
            self.x_points_s[-1],
            self.x_points_s,
            self.x_points_p,
            self.spline_s,
            self.spline_p,
            n=31)
        plt.subplots(2,2, figsize=(10,8))
        plt.subplot(2,2,1)
        # plt.figure(num="Thickness Distribution")
        plt.plot(self.x_t_vec/(self.chord_axial) * 100, self.t_vec/self.chord * 100, "k-")
        plt.title("Thickness Distribution")
        plt.ylabel("Thickness % chord")
        plt.xlabel("% chord")
        plt.grid()
        plt.subplot(2,2,2)
        # plt.figure(num="Camber_Distribution")
        plt.plot(self.x_camber/(self.chord_axial) * 100, self.theta_camber, "k-")
        plt.title("Camber Distribution")
        plt.xlabel("% chord")
        plt.ylabel(r"$\theta$")
        plt.grid()
        plt.subplot(2,2,4).remove()
        plt.subplot(2,2,3).remove()
        plt.subplot(2, 1, 2)
        for i in range(self.t_vec.shape[0]):
            plt.plot(
            [self.x_suction[i], self.x_pressure[i]],
            [self.y_suction[i], self.y_pressure[i]], 'k-', label = "Thickness" if i == 0 else "")
        plt.plot(self.x_suction, self.y_suction, label="Suction Side", lw=2.5)
        plt.plot(self.x_pressure, self.y_pressure, label="Pressure Side", lw=2.5)
        plt.plot(self.x_camber, self.y_camber, 'k--', label="Camber", lw=1.5)
        plt.legend()
        plt.axis('equal')
        plt.tight_layout()
        
        plt.show()

    def plot_suction_pressure(self):
        x_vec = np.linspace(self.x_points_s[0], self.x_points_s[-1], 51)
        y_vec = self.spline_s(x_vec)
        plt.figure(num="suction_pressure",figsize=(6,5))
        plt.plot(x_vec, y_vec, label="suction side")
        plt.axis('equal')
        plt.plot(self.x_points_s, self.y_points_s, 'ko', label='suction side control points')
        plot_circle(
            self.x_t_max,
            self.y_t_max,
            self.tmax*self.chord/100/2,
            n=51)
        plt.text(self.x_t_max,self.y_t_max,r"$t_{max}$",horizontalalignment='center')
        x_vec = np.linspace(self.x_points_p[0], self.x_points_p[-1], 51)
        y_vec = self.spline_p(x_vec)
        plt.plot(x_vec, y_vec, label="pressure side")
        plt.plot(self.x_points_p, self.y_points_p, 'ro', label='pressure side control points')        
        plt.xlabel("x [mm]")
        plt.ylabel("y [mm]")
        plt.legend()
        plt.show()
    
    def calculate_LE_TE(self):
        # fit an ellipse to the L.E.
        self.success = False
        iteration = 0
        while iteration<=10:
            Ellipse_LE = solve_system_LE_ellipse(
                self.x_points_s,
                self.x_points_p,
                self.spline_s,
                self.spline_p,
                self.kappa1,
                self.er_LE)
            self.success = Ellipse_LE.success
            if self.success:
                print("Leading Edge is fitted successfully!")
                break
            else:
                iteration += 1
        if not self.success:
            print("Leading Edge is not fitted!\nTry new input parameters")
            
        self.success = Ellipse_LE.success
        self.theta_ellipse = np.deg2rad(self.kappa1)
        self.x_s_ellipse = Ellipse_LE.x[0]
        self.x_p_ellipse = Ellipse_LE.x[1]
        self.y_s_ellipse = self.spline_s(self.x_s_ellipse)
        self.y_p_ellipse = self.spline_p(self.x_p_ellipse)
        self.h_ellipse = Ellipse_LE.x[2]
        self.k_ellipse = Ellipse_LE.x[3]
        self.a_ellipse = Ellipse_LE.x[4]
        self.b_ellipse = self.a_ellipse/self.er_LE
    # fit an ellipse (er =1 ~ circle) to the T.E.
        self.success = False
        iteration = 0
        while iteration<=5:
            Ellipse_TE = solve_system_TE_ellipse(
            self.x_points_s,
            self.x_points_p,
            self.spline_s,
            self.spline_p,
            self.kappa2,
            self.er_TE)
            self.success = Ellipse_TE.success
            if self.success:
                print("Trailing Edge is fitted successfully!")
                break
            else:
                iteration += 1
        if not self.success:
            print("Trailing-Edge is not fitted!\nTry new input parameters")
            
        self.theta_ellipse_TE = np.deg2rad(self.kappa2)
        self.x_s_ellipse_TE = Ellipse_TE.x[0]
        self.x_p_ellipse_TE = Ellipse_TE.x[1]
        self.y_s_ellipse_TE = self.spline_s(self.x_s_ellipse_TE)
        self.y_p_ellipse_TE = self.spline_p(self.x_p_ellipse_TE)
        self.h_ellipse_TE = Ellipse_TE.x[2]
        self.k_ellipse_TE = Ellipse_TE.x[3]
        self.a_ellipse_TE = Ellipse_TE.x[4]
        self.b_ellipse_TE = self.a_ellipse_TE/self.er_TE


    ### find alpha start and end of LE_ellipse
        self.alpha_start_ellipse, self.alpha_end_ellipse = find_ellipse_start_end(
            self.h_ellipse,
            self.k_ellipse,
            self.a_ellipse,
            self.b_ellipse,
            self.theta_ellipse,
            self.x_s_ellipse,
            self.y_s_ellipse,
            self.x_p_ellipse,
            self.y_p_ellipse,
            "LE")
### find alpha start and end of TE_ellipse
        self.alpha_start_ellipse_TE, self.alpha_end_ellipse_TE = find_ellipse_start_end(
            self.h_ellipse_TE,
            self.k_ellipse_TE,
            self.a_ellipse_TE,
            self.b_ellipse_TE,
            self.theta_ellipse_TE,
            self.x_s_ellipse_TE,
            self.y_s_ellipse_TE,
            self.x_p_ellipse_TE,
            self.y_p_ellipse_TE,
            "TE")
        
    def plot_blade(self):
        self.x_points_s_final = self.x_points_s.copy()
        self.x_points_s_final[0] =self.x_s_ellipse
        self.x_points_s_final[-1]=self.x_s_ellipse_TE
        self.x_points_p_final = self.x_points_p.copy()
        self.x_points_p_final[0] =self.x_p_ellipse
        self.x_points_p_final[-1]=self.x_p_ellipse_TE
        x_vec = np.linspace(self.x_points_s_final[0], self.x_points_s_final[-1], 51)
        y_vec = self.spline_s(x_vec)

        plt.figure(num="final_Blade",figsize=(6,5))
        plt.plot(x_vec, y_vec, label="suction side")
        plt.axis('equal')
        # plt.plot(self.x_points_s, self.y_points_s,'ko', label='suction side control points')

        x_vec = np.linspace(self.x_points_p_final[0], self.x_points_p_final[-1], 51)
        y_vec = self.spline_p(x_vec)
        plt.plot(x_vec, y_vec, label="pressure side")
        # plt.plot(self.x_points_p, self.y_points_p,'ro', label='pressure side control points')
        plot_ellipse(
            self.h_ellipse,
            self.k_ellipse,
            self.a_ellipse,
            self.b_ellipse,
            self.theta_ellipse,
            self.alpha_start_ellipse,
            self.alpha_end_ellipse,
            title="L.E.")
        plot_ellipse(
            self.h_ellipse_TE,
            self.k_ellipse_TE,
            self.a_ellipse_TE,
            self.b_ellipse_TE,
            self.theta_ellipse_TE,
            self.alpha_start_ellipse_TE,
            self.alpha_end_ellipse_TE,
            title="T.E.")        
        
        plt.xlabel("x [mm]")
        plt.ylabel("y [mm]")
        plt.legend()
        plt.show()

    def post_process(self,plot):
        self.x_vec_s = np.linspace(self.x_points_s_final.min(), self.x_points_s_final.max(), 300)
        self.spline_derv_s = self.spline_s.derivative()
        self.second_derv_s = self.spline_derv_s.derivative()
        self.y_vec_s = self.spline_derv_s(self.x_vec_s)
        self.y2_vec_s = self.second_derv_s(self.x_vec_s)
        self.R_curvature_s = ( (1 + self.y_vec_s**2)**1.5) / abs( self.y2_vec_s)

        self.x_vec_p = np.linspace(self.x_points_p_final.min(), self.x_points_p_final.max(), 300)
        self.spline_derv_p = self.spline_p.derivative()
        self.second_derv_p = self.spline_derv_p.derivative()
        self.y_vec_p = self.spline_derv_p(self.x_vec_p)
        self.y2_vec_p = self.second_derv_p(self.x_vec_p)
        self.R_curvature_p = ( (1 + self.y_vec_p**2)**1.5) / abs(self.y2_vec_p)
        if plot==1:

            plt.subplots(2,3, figsize=(15,10))
            plt.subplot(2,3,1)
            plt.plot(self.x_vec_s, self.y_vec_s)
            plt.title("Derivative of suction side")
            plt.xlabel("x",fontsize=15)
            plt.ylabel("y'",fontsize=15)
            plt.subplot(2,3,2)
            plt.plot(self.x_vec_s, self.y2_vec_s)
            plt.title("Second derivative of suction side")
            plt.xlabel("x",fontsize=15)
            plt.ylabel("y''",fontsize=15)
            plt.tight_layout()
            plt.subplot(2,3,3)
            plt.plot(self.x_vec_s, 1/self.R_curvature_s)
            plt.title("Curvature of suction side")
            plt.xlabel("x",fontsize=15)
            plt.ylabel("1/R($\kappa$)",fontsize=15)

            plt.subplot(2,3,4)

            plt.title("Derivative of Pressure side")
            plt.plot(self.x_vec_p, self.y_vec_p)
            plt.xlabel("x",fontsize=15)
            plt.ylabel("y'",fontsize=15)

            plt.subplot(2,3,5)
            plt.plot(self.x_vec_p, self.y2_vec_p)
            plt.title("Second derivative of pressure side")
            plt.xlabel("x",fontsize=15)
            plt.ylabel("y''",fontsize=15)
            plt.tight_layout()

            plt.subplot(2,3,6)
            plt.plot(self.x_vec_p, 1/self.R_curvature_p)
            plt.title("Curvature of pressure side")
            plt.xlabel("x",fontsize=15)
            plt.ylabel("1/R($\kappa$)",fontsize=15)

            plt.show()
        return self.x_vec_s, self.y_vec_s, self.y2_vec_s, self.R_curvature_s, self.x_vec_p, self.y_vec_p, self.y2_vec_p, self.R_curvature_p
    
    def save_parameters(self, file_name):
        self.parameters = np.array([self.kappa1,
                                         self.kappa2,
                                         self.chord,
                                         self.stagger,
                                         self.wedge_angle_LE,
                                         self.wedge_angle_TE,
                                         self.tmax,
                                         self.tmax_loc,
                                         self.LE_th,
                                         self.TE_th,
                                         self.er_LE,
                                         self.er_TE,
                                         self.x_points_s[1],
                                         self.x_points_s[2],
                                         self.y_points_s[1],
                                         self.y_points_s[2],
                                         self.x0,
                                         self.y0])
        print(self.parameters)
        self.new_df = pd.DataFrame([self.parameters], columns=[
            'kappa1', 'kappa2', 'chord', 'stagger', 'wedge_angle_LE',
            'wedge_angle_TE', 'tmax', 'tmax_loc', 'LE_th', 'TE_th',
            'er_LE', 'er_TE', 'x_s_sp1','x_s_sp2','y_s_sp1','y_s_sp2','x0','y0'])
        self.new_df.to_csv(file_name,
                           mode='a',
                           header=False,
                           index=False)

    def calculate_centroid(self):
        x_vec_s = np.linspace(self.x_points_s_final[0], self.x_points_s_final[-1], 101)
        y_vec_s = self.spline_s(x_vec_s)
        self.suction_points = np.column_stack((x_vec_s, y_vec_s))
        x_vec_p = np.linspace(self.x_points_p_final[-1], self.x_points_p_final[0], 101)
        y_vec_p = self.spline_p(x_vec_p)
        self.pressure_points = np.column_stack((x_vec_p, y_vec_p))
        plt.figure(figsize=(6,5))
        plt.plot(self.suction_points[:, 0], self.suction_points[:, 1], 'k-')
        plt.plot(self.pressure_points[:, 0], self.pressure_points[:, 1], 'k-')
        self.x_LE, self.y_LE = ellipse_points(
            self.h_ellipse,
            self.k_ellipse,
            self.a_ellipse,
            self.b_ellipse,
            self.theta_ellipse,
            self.alpha_start_ellipse,
            self.alpha_end_ellipse)
        self.ellipse_LE_points = np.column_stack((self.x_LE, self.y_LE))
        plt.plot(self.ellipse_LE_points[:, 0], self.ellipse_LE_points[:, 1], 'k-')
        self.x_TE, self.y_TE = ellipse_points(
            self.h_ellipse_TE,
            self.k_ellipse_TE,
            self.a_ellipse_TE,
            self.b_ellipse_TE,
            self.theta_ellipse_TE,
            self.alpha_start_ellipse_TE,
            self.alpha_end_ellipse_TE)
        self.ellipse_TE_points = np.column_stack((self.x_TE, self.y_TE))
        plt.plot(self.ellipse_TE_points[:, 0], self.ellipse_TE_points[:, 1], 'k-')
        plt.grid('minor')
        plt.axis('equal')
        
        airfoil_points = np.concatenate([self.suction_points, self.ellipse_TE_points , self.pressure_points, self.ellipse_LE_points])

        # Create a polygon from these points
        airfoil_polygon = Polygon(airfoil_points)
        # Calculate the centroid
        self.centroid = airfoil_polygon.centroid
        print(f"The centroid of the airfoil is at: ({round(self.centroid.x,3)}, {round(self.centroid.y,3)})")
        plt.plot(self.centroid.x, self.centroid.y, 'rx', label="Centroid")
        plt.legend()
        plt.show()
        return self.centroid.x, self.centroid.y


    def export_points(self, n_suction = 51, n_pressure= 51, n_LE = 13, n_TE = 9):
        x_vec_s = np.linspace(self.x_points_s_final[0], self.x_points_s_final[-1], n_suction)
        y_vec_s = self.spline_s(x_vec_s)
        self.suction_points = np.column_stack((x_vec_s, y_vec_s))
        x_vec_p = np.linspace(self.x_points_p_final[-1], self.x_points_p_final[0], n_pressure)
        y_vec_p = self.spline_p(x_vec_p)
        self.pressure_points = np.column_stack((x_vec_p, y_vec_p))
        self.x_LE, self.y_LE = ellipse_points(
            self.h_ellipse,
            self.k_ellipse,
            self.a_ellipse,
            self.b_ellipse,
            self.theta_ellipse,
            self.alpha_start_ellipse,
            self.alpha_end_ellipse,
            n=n_LE)
        self.ellipse_LE_points = np.column_stack((self.x_LE, self.y_LE))
        self.x_TE, self.y_TE = ellipse_points(
            self.h_ellipse_TE,
            self.k_ellipse_TE,
            self.a_ellipse_TE,
            self.b_ellipse_TE,
            self.theta_ellipse_TE,
            self.alpha_start_ellipse_TE,
            self.alpha_end_ellipse_TE,
            n=n_TE)
        self.ellipse_TE_points = np.column_stack((self.x_TE, self.y_TE))
        np.savetxt('TE.txt', self.ellipse_TE_points, delimiter=',',
                   fmt='%.4f', header='x y', comments='')
        np.savetxt('LE.txt', self.ellipse_LE_points, delimiter=',',
                   fmt='%.4f', header='x y', comments='')
        np.savetxt('Suction.txt', self.suction_points, delimiter=',',
                   fmt='%.4f', header='x y', comments='')
        np.savetxt('Pressure.txt', self.pressure_points, delimiter=',',
                   fmt='%.4f', header='x y', comments='')   
        
        airfoil_points = np.concatenate([self.suction_points, self.ellipse_TE_points , self.pressure_points, self.ellipse_LE_points])



