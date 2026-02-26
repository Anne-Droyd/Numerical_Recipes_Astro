import numpy as np
import timeit
import math

import pandas as pd

from functions import interpolation
from functions import linear_algebra
from functions import math_functions

import matplotlib.pyplot as plt

from functions.linear_algebra import crouts_algorithm


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

    vandermonde = math_functions.get_vandermonde_matrix(x,y)

    A = crouts_algorithm(vandermonde)
    M=20
    nevilles    = []
    lagrange    = []
    linear      = []
    x_points = []
    x_y = np.array([x, y])

    #need to split this up for timing purposes
    for point_index in range(0,len(x_points_to_interpolate)):

        point_request = x_points_to_interpolate[point_index]

        interpolated_values_nev     = interpolation.nevilles_algorithm(point_request,M,x_y)
        interpolated_values_lag     = interpolation.lagrange_polynomial(point_request,19,x_y)
        interpolated_cubic_spline   = interpolation.
        interpolated_akima_spline   = interpolation.

        lagrange.append(interpolated_values_lag[1])
        nevilles.append(interpolated_values_nev[1])
        x_points.append(x_points_to_interpolate[point_index])

    plt.scatter(x_points,nevilles, marker = '^', alpha = 0.4, label = 'nevilles')
    plt.scatter(x_points,lagrange, marker = 'o', alpha = 0.2, label = 'lagrange')
    plt.ylim(-200,200)
    plt.scatter(x,y,label = 'true values')
    plt.legend()
    plt.show()


question_1_results = question_1()
print('\n')
print('Question 1 values and results')
print('------------------------------------------------------')
print(question_1_results)
print('------------------------------------------------------')
print('\n')

question_2()
print('\n')
print('Question 2 values and results')
print('------------------------------------------------------')

print('------------------------------------------------------')
print('\n')