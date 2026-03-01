import numpy as np
import timeit
import math

import pandas as pd

from functions import interpolation
from functions import math_functions
from functions import linear_algebra as la

import matplotlib.pyplot as plt

from functions.linear_algebra import forward_sub, backward_sub


def import_data():
    url = 'https://home.strw.leidenuniv.nl/~daalen/Handin_files/Vandermonde.txt'
    data=pd.DataFrame(np.genfromtxt(url,comments='#',dtype=np.float64),columns=['x','y'])
    return data



def question_1() -> pd.DataFrame:
    question_1_values = np.array([[1,0],[5,10],[3,21],[2.6,40],[100,5],[101,200]],dtype=np.float32)

    results = []
    for values in question_1_values:

        #defining how the data is stored in 32 bits
        lam_raw = np.float32(values[0])
        if lam_raw.is_integer():
            lam = np.int32(lam_raw)
        else:
            lam = lam_raw
        k = np.int32(values[1])

        result = math_functions.poisson_distribution(lam,k)
        results.append([lam,k,result])
    results = pd.DataFrame(results,columns=['lambda','k','Poisson result'])
    return results

def question_2() -> tuple[pd.DataFrame]:
    data = import_data()
    x,y = np.array(data['x']),np.array(data['y'])

    x_points_to_interpolate = np.linspace(x[0],x[-1],1001)

    vandermonde = math_functions.get_vandermonde_matrix(x)

    L,U = la.crout(vandermonde)
    A = la.crouts_algorithm(vandermonde)
    L1 = L.copy()
    U1 = U.copy()
    y1 = y.copy()

    out = la.forward_sub(A,y1)
    c = la.backward_sub(A,out)
    for ans in c:
        print(ans)

    ys = math_functions.get_y(x,c,1001)

    y1 = y.copy()
    L1 = L.copy()
    U1 = U.copy()
    for n in range(0,10):
        z = forward_sub(L1,y1)
        c_1 = backward_sub(U1,z)
    yyy = math_functions.get_y(x,c_1,1001)
    print(yyy)
    M=20 #change to 3/5/8 for smooth function
    nevilles    = []
    lagrange    = []
    linear      = []
    cubic       = []
    akima       = []
    x_points = []
    x_y = np.array([x, y])

    #need to split this up for timing purposes
    for point_index in range(0,len(x_points_to_interpolate)):

        point_request = x_points_to_interpolate[point_index]

        interpolated_values_linear  = interpolation.linear_interpolation(point_request,x_y)
        interpolated_values_nev     = interpolation.nevilles_algorithm(point_request,M,x_y)
        interpolated_values_lag     = interpolation.lagrange_polynomial(point_request,19,x_y)
        interpolated_cubic_spline   = interpolation.cubic_spline(point_request,x_y)
        # interpolated_akima_spline   = interpolation.

        linear.append(interpolated_values_linear[1])
        nevilles.append(interpolated_values_nev[1])
        lagrange.append(interpolated_values_lag[1])
        # cubic.append(interpolated_cubic_spline[1])
        # akima.append(interpolated_akima_spline[1])
        x_points.append(x_points_to_interpolate[point_index])

    fig, axes = plt.subplots(2,1)
    axes = axes.flatten()
    axes[0].scatter(x_points,yyy, marker = 'o', alpha = 0.4, label = 'yyy')
    axes[0].scatter(x_points,ys, marker = '^', alpha = 0.4, label = 'ys')
    # axes[0].scatter(x_points,linear, marker = 'o', alpha = 0.4, label = 'linear')
    # axes[0].scatter(x_points,nevilles, marker = '^', alpha = 0.4, label = 'nevilles')
    # axes[0].scatter(x_points,lagrange, marker = 'o', alpha = 0.2, label = 'lagrange')
    axes[0].set_ylim(-200,200)
    axes[0].scatter(x,y,label = 'true values')


    # axes[1].plt()

    fig.legend()
    plt.show()


question_1_results = question_1()
print('\n')
print('Question 1 values and results')
print('------------------------------------------------------')
print(question_1_results)
print('------------------------------------------------------')
print('\n')

question_2_results = question_2()
print('\n')
print('Question 2 values and results')
print('------------------------------------------------------')
print(question_2_results)
print('------------------------------------------------------')
print('\n')