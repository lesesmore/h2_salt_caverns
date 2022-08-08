#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 20 14:27:04 2022

@author: lesarmstrong
"""

# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import seaborn as sbn

from eos import create_eos
import preos


##############################################################################################################





def pressure_to_moles(p, V, T):
    
    # Peng Robinson coefficients for H2

    # Critical temperature and pressure
    Tc = 33 # [K]
    Pc = 1.297 * 10  # [MPa] --> [bar]
    # accentric factor
    omega = -0.215
    
    # Change MPa --> bar
    p = p*10
    
    hydrogen = preos.Molecule("hydrogen", Tc, Pc, omega)    
    
    props = preos.preos(hydrogen, T, p, plotcubic=False, printresults=False)
    density = props['density(mol/m3)']
    
    moles = density * V
    
    return(moles)


def pressure_from_PR(n, T, V):
    
    # Peng Robinson coefficients for H2

    # Critical temperature and pressure
    Tc = 33 # K
    Pc = 1.297 * 1000000 # [MPa] --> [Pa]
    # accentric factor
    omega = -0.215

    
    # Calculates specific volume from fixed Volume and variable mass
    Vn = V/n
    
    # Peng-Robinson Equation of State
    pr_eos = create_eos('PR', Pc, Tc, omega=omega)
    
   # Calculate pressure p from given temperature and Vm (which is related to mass since V is constant)
    p = pr_eos.pressure(T, Vn) #[Pa]   
    p = p / 1000000 # [Pa] -- > [MPa]
    
   # P = ((R*T)/ (Vm - b)) - ((a * alpha)/(Vm**2 + 2 * b * Vm - b**2))
    
    return(p)
    



################################################################################################################################



class Cavern:    
    
    def __init__(self, depth, shape, height):
        
        
        # Varies between 2200 and 2300
        row_rock = 2250


        # Hydrogen viscosity (Pa m)
        mu = 9.37 * 10**(-6)

        # Molarity (kg/mol)
        M = 1.00794 * 10**(-3)

        # Gas constant (m^3 Pa / K mol) or (J / K mol)
        R = 8.3144

        # Gravity (m/s^2)
        g = 9.81
        
        # MA SAYS 10^-6, but 10^-7 makes more sense??
        self.D = 6 * 10**(-7) #MPa^-n year^-1 
        
        # Convergence number. Has to do with salt properties / Ref: ???? Nye et al
        self.n = 4.089
        
    
        self.depth = depth
        self.height = height
        self.shape = shape


        # Derived values
        # Volume of cylinder [m^3]
        

        # Overburden pressure[MPa]
        #self.P_over_top = (row_rock * g * self.depth) / 1000000
        #self.P_over_bot = (row_rock * g * (self.depth + self.height)) / 1000000



        
        # Temperature (K). Assumes average surface temp of 15C
        self.temperature = 288 + 0.0025 * ((self.depth - self.height)/2)
        
        


    def __str__(self):
        return (f'r={self.radius}, d={self.depth} h={self.height}, shape={self.shape}')


    '''
    
    def pressure_from_PR(self):
        
        # Peng Robinson coefficients for H2
    
        # Critical temperature and pressure
        Tc = 33 # K
        Pc = 1.297 * 1000000 # [MPa] --> [Pa]
        # accentric factor
        omega = -0.215
    
        
        # Calculates specific volume from fixed Volume and variable mass
        Vn = self.volume/self.moles
        
        # Peng-Robinson Equation of State
        pr_eos = create_eos('PR', Pc, Tc, omega=omega)
        
       # Calculate pressure p from given temperature and Vm (which is related to mass since V is constant)
        p = pr_eos.pressure(self.Temperature, Vn) #[Pa]   
        p = p / 1000000 # [Pa] -- > [MPa]
        
        self.pressure = p
       # P = ((R*T)/ (Vm - b)) - ((a * alpha)/(Vm**2 + 2 * b * Vm - b**2))
        
        return(p)
    '''



# The legal maximum diameter for caverns in the Netherlands is 90 m


class horizontal_cylinder(Cavern):
    

        def __init__(self, radius, depth, lenght, shape='horizontal_cylinder'):
        

            Cavern.__init__(self, depth, shape, height=radius*2,)    
        
            self.radius = radius
            self.lenght = lenght 
     
        
            T = self.temperature
            self.volume = (np.pi) * (self.radius**2) * self.lenght
            
            # Steady state creep (Lankof 2022)
            n = 4.089
            A_optimistic = 0.3423 #[%]
            A_pessimistic = 0.6846 #[%]
            
            # Varies between 2200 and 2300
            row_rock = 2250


            # Hydrogen viscosity (Pa m)
            mu = 9.37 * 10**(-6)

            # Molarity (kg/mol)
            M = 1.00794 * 10**(-3)

            # Gas constant (m^3 Pa / K mol) or (J / K mol)
            R = 8.3144

            # Gravity (m/s^2)
            g = 9.81
            
            
            
            
            
           # Overburden pressures [MPa] / Ref: Caglayan
            self.P_over_top = (row_rock * g * self.depth) / 1000000
            self.P_over_bot = (row_rock * g * (self.depth + self.radius * 2)) / 1000000
    
        
            # Max pressure [MPa] / Can change safety constant 0.8 to something else
            self.P_max = self.P_over_top * 0.8 
    
            ## ESPECIALLY MIN PRESSURE
            # Min pressure [MPa]. Two kinds of pressures that compete with each other: Stress & Creep. We take the largest value of the two as the base of the working gas pressure.
            
            # pmin1 is related to stress. Minumum pressure that is needed to avert cavern collapse / Ref: ?????
            p_min1 = self.P_over_bot * 0.24
            self.p_min_stress = p_min1
            
           # pmin_2 is related to creep.           
           
            # Overburden pressure at the bottom is the geostatic pressure 
            p_geostatic = self.P_over_bot
            
            # Equation is set up to cutoff so that delta_v of loss of cavern over t years is avoided.
            
            # t years of convergence. delta_v volume lossed
            t = 30
            delta_v = 0.3
           
            
            # Algebra eq (29) from Ma et al paper where p_internal is Pc. / Ref: Xuqiang Ma et al., “Creep Deformation Analysis of Gas Storage in Salt Caverns,” International Journal of Rock Mechanics and Mining Sciences 139 (March 1, 2021): 104635, https://doi.org/10.1016/j.ijrmms.2021.104635.
            p_min2 = p_geostatic - ( np.log(1 - delta_v) / ( t * (-np.sqrt(3)) * self.D * (np.sqrt(3)/(self.n))**self.n )  )**(1/self.n)
            self.p_min_convergence = p_min2
            
            # There are 2 parameters affecting min pressure. One is the convergence and loss of cavern volume given by creep established by pmin2 (30% loss if )
            if p_min2 > p_min1:
                
                self.P_min = p_min2
                
            else:
                self.P_min = p_min1
            
            # If P_min > P_max it means that convergence is happening at a rate greather than 1% per year (30% lose of volume in 30 years)
            # THIS SETS THE PRESSURES AS EQUAL TO EACH OTHER SO WE WILL GET 0 MASS
            if self.P_min >= self.P_max:
                self.P_min = self.P_max
            
        
            # Calculates max and min possible mass given pressure contraints
            # Goes onto calculating working gas capacity and cushion gas
            self.Mass_min = pressure_to_moles(self.P_min, self.volume, self.temperature) * 2.016 / 1000000 #[kg]
            self.Mass_max = pressure_to_moles(self.P_max, self.volume, self.temperature) * 2.016 / 1000000 #[kg]
            
            # Max and min mass give us the paramount Working Gas Capacity
            self.working_gas_capacity = self.Mass_max - self.Mass_min #[kg]
            
            ### ????? ###
            #self.cushion_gas =  self.working_gas_capacity / self.Mass_min
            
            # Cushion gas is the amount of hydrogen that always needs to be present in the well
            self.cushion_gas = self.Mass_min



            ## ASSUMING THAT THE INITIAL PRESSURE IS PMIN
            self.pressure = self.P_min 


            # Assigning variables to cavern after assuming that the initial pressure is pmi
            self.moles = pressure_to_moles(self.pressure, self.volume, self.temperature)
            self.mass = self.moles * 2.016 / 1000000 #[Kg]
            self.energy = self.mass * 33.33 / 1000 #[MW/h] 
            
            
            ## FOLLOWING FOR PLOTTING POSSIBLE VALUES, NOT ACTUAL STATE OF THE SYSTEM.
            ## TRUE VALUES OF STATE OF SYSTEM WILL USE WHOLE NAME, NOT ABBREVIATIONS
            
            ###### CAN CHANGE
            ##### Number of intervals of moles calculated between max and min pressure
            n_intervals = 10
            n_min = pressure_to_moles(self.P_min, self.volume, self.temperature)
            n_max = pressure_to_moles(self.P_max, self.volume, self.temperature)
            n_array = np.linspace(n_min,n_max,n_intervals)
            
            self.n_array = n_array
            
            PR_pressures = np.array([])
            for i in n_array:
                PR_pressures = np.append(pressure_from_PR(n=i, T=self.temperature, V=self.volume),PR_pressures)
            PR_pressures = np.flip(PR_pressures)
            self.PR_pressures = PR_pressures
            

            self.m_array =  n_array * 2.016 / 1000000 #[metric tons]
            
            self.e_array = self.m_array * 33.33 / 1000 #[MW/h] 
            
      
class vertical_cylinder(Cavern):
    

        def __init__(self, radius, depth, height, shape='vertical cylinder'):
        

            Cavern.__init__(self, depth, shape, height)    
        
            self.radius = radius
     
        
            T = self.temperature
            self.volume = (np.pi) * (self.radius**2) * self.height
            
            # Steady state creep (Lankof 2022)
            n = 4.089
            A_optimistic = 0.3423 #[%]
            A_pessimistic = 0.6846 #[%]
            
            # Varies between 2200 and 2300
            row_rock = 2250


            # Hydrogen viscosity (Pa m)
            mu = 9.37 * 10**(-6)

            # Molarity (kg/mol)
            M = 1.00794 * 10**(-3)

            # Gas constant (m^3 Pa / K mol) or (J / K mol)
            R = 8.3144

            # Gravity (m/s^2)
            g = 9.81
            
            
            
            
           # Overburden pressures [MPa] / Ref: Caglayan
            self.P_over_top = (row_rock * g * self.depth) / 1000000
            self.P_over_bot = (row_rock * g * (self.depth + self.height)) / 1000000
    
        
            # Max pressure [MPa] / Can change safety constant 0.8 to something else
            self.P_max = self.P_over_top * 0.8 
    
            ## ESPECIALLY MIN PRESSURE
            # Min pressure [MPa]. Two kinds of pressures that compete with each other: Stress & Creep. We take the largest value of the two as the base of the working gas pressure.
            
            # pmin1 is related to stress. Minumum pressure that is needed to avert cavern collapse / Ref: ?????
            p_min1 = self.P_over_bot * 0.24
            self.p_min_stress = p_min1
            
           # pmin_2 is related to creep.           
           
            # Overburden pressure at the bottom is the geostatic pressure 
            p_geostatic = self.P_over_bot
            
            # Equation is set up to cutoff so that delta_v of loss of cavern over t years is avoided.
            
            # t years of convergence. delta_v volume lossed
            t = 30
            delta_v = 0.3
           
            
            # Algebra eq (29) from Ma et al paper where p_internal is Pc. / Ref: Xuqiang Ma et al., “Creep Deformation Analysis of Gas Storage in Salt Caverns,” International Journal of Rock Mechanics and Mining Sciences 139 (March 1, 2021): 104635, https://doi.org/10.1016/j.ijrmms.2021.104635.
            p_min2 = p_geostatic - ( np.log(1 - delta_v) / ( t * (-np.sqrt(3)) * self.D * (np.sqrt(3)/(self.n))**self.n )  )**(1/self.n)
            self.p_min_convergence = p_min2
            
            # There are 2 parameters affecting min pressure. One is the convergence and loss of cavern volume given by creep established by pmin2 (30% loss if )
            if p_min2 > p_min1:
                
                self.P_min = p_min2
                
            else:
                self.P_min = p_min1
            
            # If P_min > P_max it means that convergence is happening at a rate greather than 1% per year (30% lose of volume in 30 years)
            # THIS SETS THE PRESSURES AS EQUAL TO EACH OTHER SO WE WILL GET 0 MASS
            if self.P_min >= self.P_max:
                self.P_min = self.P_max
            
        
            # Calculates max and min possible mass given pressure contraints
            # Goes onto calculating working gas capacity and cushion gas
            self.Mass_min = pressure_to_moles(self.P_min, self.volume, self.temperature) * 2.016 / 1000000 #[kg]
            self.Mass_max = pressure_to_moles(self.P_max, self.volume, self.temperature) * 2.016 / 1000000 #[kg]
            
            # Max and min mass give us the paramount Working Gas Capacity
            self.working_gas_capacity = self.Mass_max - self.Mass_min #[kg]
            
            ### ????? ###
            #self.cushion_gas =  self.working_gas_capacity / self.Mass_min
            
            # Cushion gas is the amount of hydrogen that always needs to be present in the well
            self.cushion_gas = self.Mass_min



            ## ASSUMING THAT THE INITIAL PRESSURE IS PMIN
            self.pressure = self.P_min 


            # Assigning variables to cavern after assuming that the initial pressure is pmi
            self.moles = pressure_to_moles(self.pressure, self.volume, self.temperature)
            self.mass = self.moles * 2.016 / 1000000 #[Kg]
            self.energy = self.mass * 33.33 / 1000 #[MW/h] 
            
            
            ## FOLLOWING FOR PLOTTING POSSIBLE VALUES, NOT ACTUAL STATE OF THE SYSTEM.
            ## TRUE VALUES OF STATE OF SYSTEM WILL USE WHOLE NAME, NOT ABBREVIATIONS
            
            ###### CAN CHANGE
            ##### Number of intervals of moles calculated between max and min pressure
            n_intervals = 10
            n_min = pressure_to_moles(self.P_min, self.volume, self.temperature)
            n_max = pressure_to_moles(self.P_max, self.volume, self.temperature)
            n_array = np.linspace(n_min,n_max,n_intervals)
            
            self.n_array = n_array
            
            PR_pressures = np.array([])
            for i in n_array:
                PR_pressures = np.append(pressure_from_PR(n=i, T=self.temperature, V=self.volume),PR_pressures)
            PR_pressures = np.flip(PR_pressures)
            self.PR_pressures = PR_pressures
            

            self.m_array =  n_array * 2.016 / 1000000 #[metric tons]
            
            self.e_array = self.m_array * 33.33 / 1000 #[MW/h] 
            
          
             
      
class sphere(Cavern):
    

        def __init__(self, radius, depth, shape='sphere'):
        

            Cavern.__init__(self, depth, shape, height=radius*2)    
        
            self.radius = radius
     
            self.volume = 4  * (np.pi) * (self.radius**3) / 3
            T = self.temperature
          
            # Steady state creep (Lankof 2022)
            n = 4.089
            A_optimistic = 0.3423 #[%]
            A_pessimistic = 0.6846 #[%]
            
            
            # Varies between 2200 and 2300
            row_rock = 2250


            # Hydrogen viscosity (Pa m)
            mu = 9.37 * 10**(-6)

            # Molarity (kg/mol)
            M = 1.00794 * 10**(-3)

            # Gas constant (m^3 Pa / K mol) or (J / K mol)
            R = 8.3144

            # Gravity (m/s^2)
            g = 9.81
            
            
            
           # Overburden pressures [MPa] / Ref: Caglayan
            self.P_over_top = (row_rock * g * self.depth) / 1000000
            self.P_over_bot = (row_rock * g * (self.depth + self.radius*2)) / 1000000
    
        
            # Max pressure [MPa] / Can change safety constant 0.8 to something else
            self.P_max = self.P_over_top * 0.8 
    
            ## ESPECIALLY MIN PRESSURE
            # Min pressure [MPa]. Two kinds of pressures that compete with each other: Stress & Creep. We take the largest value of the two as the base of the working gas pressure.
            
            # pmin1 is related to stress. Minumum pressure that is needed to avert cavern collapse / Ref: ?????
            p_min1 = self.P_over_bot * 0.24
            self.p_min_stress = p_min1
            
           # pmin_2 is related to creep.           
           
            # Overburden pressure at the bottom is the geostatic pressure 
            p_geostatic = self.P_over_bot
            
            # Equation is set up to cutoff so that delta_v of loss of cavern over t years is avoided.
            
            # t years of convergence. delta_v volume lossed
            t = 30
            delta_v = 0.3
           
            
            # Algebra eq (29) from Ma et al paper where p_internal is Pc. / Ref: Xuqiang Ma et al., “Creep Deformation Analysis of Gas Storage in Salt Caverns,” International Journal of Rock Mechanics and Mining Sciences 139 (March 1, 2021): 104635, https://doi.org/10.1016/j.ijrmms.2021.104635.
            p_min2 = p_geostatic - ( np.log(1 - delta_v) / ( t * (-3/2) * self.D * (3/(2*self.n))**self.n )  )**(1/self.n)
            
            self.p_min_convergence = p_min2
            
            # There are 2 parameters affecting min pressure. One is the convergence and loss of cavern volume given by creep established by pmin2 (30% loss if )
            if p_min2 > p_min1:
                
                self.P_min = p_min2
                
            else:
                self.P_min = p_min1
            
            # If P_min > P_max it means that convergence is happening at a rate greather than 1% per year (30% lose of volume in 30 years)
            # THIS SETS THE PRESSURES AS EQUAL TO EACH OTHER SO WE WILL GET 0 MASS
            if self.P_min >= self.P_max:
                self.P_min = self.P_max
            
        
            # Calculates max and min possible mass given pressure contraints
            # Goes onto calculating working gas capacity and cushion gas
            self.Mass_min = pressure_to_moles(self.P_min, self.volume, self.temperature) * 2.016 / 1000000 #[kg]
            self.Mass_max = pressure_to_moles(self.P_max, self.volume, self.temperature) * 2.016 / 1000000 #[kg]
            
            # Max and min mass give us the paramount Working Gas Capacity
            self.working_gas_capacity = self.Mass_max - self.Mass_min #[kg]
            
            ### ????? ###
            #self.cushion_gas =  self.working_gas_capacity / self.Mass_min
            
            # Cushion gas is the amount of hydrogen that always needs to be present in the well
            self.cushion_gas = self.Mass_min



            ## ASSUMING THAT THE INITIAL PRESSURE IS PMIN
            self.pressure = self.P_min 


            # Assigning variables to cavern after assuming that the initial pressure is pmi
            self.moles = pressure_to_moles(self.pressure, self.volume, self.temperature)
            self.mass = self.moles * 2.016 / 1000000 #[Kg]
            self.energy = self.mass * 33.33 / 1000 #[MW/h] 
            
            
            ## FOLLOWING FOR PLOTTING POSSIBLE VALUES, NOT ACTUAL STATE OF THE SYSTEM.
            ## TRUE VALUES OF STATE OF SYSTEM WILL USE WHOLE NAME, NOT ABBREVIATIONS
            
            ###### CAN CHANGE
            ##### Number of intervals of moles calculated between max and min pressure for plotting
            n_intervals = 10
            n_min = pressure_to_moles(self.P_min, self.volume, self.temperature)
            n_max = pressure_to_moles(self.P_max, self.volume, self.temperature)
            n_array = np.linspace(n_min,n_max,n_intervals)
            
            self.n_array = n_array
            
            PR_pressures = np.array([])
            for i in n_array:
                PR_pressures = np.append(pressure_from_PR(n=i, T=self.temperature, V=self.volume),PR_pressures)
            PR_pressures = np.flip(PR_pressures)
            self.PR_pressures = PR_pressures
            

            self.m_array =  n_array * 2.016 / 1000000 #[metric tons]
            
            self.e_array = self.m_array * 33.33 / 1000 #[MW/h] 
            
          
             
          





########################################################################################################################################################################################################################################
########################################################################################################################################################################################################################################



## SIZE OF GRID from arcgis?? 

a = 540
A = a*a


def capacity_calculation_vertical_cavern_from_csv(path,a=540,range_outputs=False):
    
    A = a*a
    
    df = pd.read_csv(path)

    # List of working gas capacity of single caverns
   
    working_gas_sphere_single_cavern_list = []
    working_gas_sphere_list = []
    
    working_gas_horizontal_cylinder_single_cavern_list = []
    working_gas_horizontal_cylinder_list = []
    
    working_gas_vertical_cylinder_single_cavern_list = []
    working_gas_vertical_cylinder_list = []


    for i in range(len(df)):
        
        # Buffer in salt that cavern needs [m]
        buffer = 15
        
        
        
        ## For Spheres
        radius = (df['thickness'][i] - buffer*2) / 2
        cavern = sphere(radius, df['depth'][i])
       
        masskg = cavern.working_gas_capacity #[kg]
        massmetric_ton = masskg / 1000  # [kg] --> [metric ton]
        working_gas_sphere_single_cavern_list.append(massmetric_ton)
        
        
        Ncavernsthatfitinasinglecell_sphere = A / (np.pi * ((cavern.radius*8)**2))
        total_working_gas_sphere_capacity =  massmetric_ton * Ncavernsthatfitinasinglecell_sphere
        working_gas_sphere_list.append(total_working_gas_sphere_capacity)
        
        
        
        ## For horizontal caverns
        radius = (df['thickness'][i] - buffer*2) / 2
        cavern = horizontal_cylinder(radius, df['depth'][i], a - buffer*2)
        
        masskg = cavern.working_gas_capacity #[kg]
        massmetric_ton = masskg / 1000  # [kg] --> [metric ton]
        working_gas_horizontal_cylinder_single_cavern_list.append(massmetric_ton)
        
        
        Ncavernsthatfitinasinglecell_horizantalcylinder = A / (cavern.lenght * (cavern.radius*8))
        total_working_gas_horizontal_cylinder_capacity =  massmetric_ton * Ncavernsthatfitinasinglecell_horizantalcylinder
        working_gas_horizontal_cylinder_list.append(total_working_gas_horizontal_cylinder_capacity)
        
        
        ## For vertical caverns
        

        height = (df['thickness'][i] - buffer*2) / 2
       
       ## HOW DO WE KNOW THE RADIUS OF VERTICAL CAVERN??  THIS IS A PLACEHOLDER BASED ON HEURISTICS
        radius = height / 4
        
        
        cavern = vertical_cylinder(radius,df['depth'][i],height)
        
        masskg = cavern.working_gas_capacity #[kg]
        massmetric_ton = masskg / 1000  # [kg] --> [metric ton]
        working_gas_vertical_cylinder_single_cavern_list.append(massmetric_ton)
        
        
        Ncavernsthatfitinasinglecell_verticalcylinder = A / (np.pi * ((cavern.radius*8)**2))
        total_working_gas_vertical_cylinder_capacity =  massmetric_ton * Ncavernsthatfitinasinglecell_verticalcylinder
        working_gas_vertical_cylinder_list.append(total_working_gas_vertical_cylinder_capacity)
        
      
        
    df['sphere_working_gas_total']  = working_gas_sphere_list 
    df['sphere_working_gas_single_cavern']  = working_gas_sphere_single_cavern_list 
    
    df['horizontal_cylinder_working_gas_total']  = working_gas_horizontal_cylinder_list 
    df['horizontal_cylinder_working_gas_single_cavern']  = working_gas_horizontal_cylinder_single_cavern_list 
        
    df['vertical_cylinder_working_gas_total']  = working_gas_vertical_cylinder_list 
    df['vertical_cylinder_working_gas_single_cavern']  = working_gas_vertical_cylinder_single_cavern_list  
    
    processed_path = path[:-4] + '_processed.csv'
    
    df.to_csv(processed_path)
        
        
        
        
path='/Users/lesarmstrong/Documents/UHS/SYR_dummy.csv'       
        
        
    
    
    





















