
import os
import time
import numpy as np
import pandas as pd

from functions import interpolation
from functions import math_functions
from functions import linear_algebra as la

import matplotlib.pyplot as plt

# get the folder where the script is located
script_dir = os.path.dirname(__file__)
data_dir = os.path.join(script_dir, "../data")
plots_dir = os.path.join(script_dir, "../plots")

os.makedirs(plots_dir, exist_ok=True)  # just in case

def import_data():

    data_path = os.path.join(data_dir, "Vandermonde.txt")
    data=pd.DataFrame(np.genfromtxt(data_path,comments='#',dtype=np.float64),columns=['x','y'])
    return data

def LU_1_iter(vandermonde,x,y,points = 1001):
    A = la.crouts_algorithm(vandermonde)
    y1 = y.copy()
    out = la.forward_sub(A, y1)
    c = la.backward_sub(A, out)
    y = math_functions.get_y(x,c,points)
    return c, y

def LU_10_iter(vandermonde,x,y,points = 1001):
    y1 = y.copy()
    A = la.crouts_algorithm(vandermonde)
    z = la.forward_sub(A, y1)
    c = la.backward_sub(A, z)

    for _ in range(10):
        #todo ask for clarification.
        # this was applying to my x? i don't understand how, i thought this was y, I need to ask for clarification

        residual = y1 - vandermonde @ c
        z = la.forward_sub(A, residual)
        delta = la.backward_sub(A, z)
        c += delta
    y = math_functions.get_y(x,c,points)
    return c, y

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
    results.to_csv(os.path.join(data_dir, 'poisson.csv'))
    return results

def question_2() -> tuple[pd.DataFrame]:
    data = import_data()
    x,y = np.array(data['x']),np.array(data['y'])
    x_points_to_interpolate = np.linspace(x[0],x[-1],1001)
    vandermonde = math_functions.get_vandermonde_matrix(x)

    start = time.perf_counter()
    c, ys = LU_1_iter(vandermonde, x, y)
    end = time.perf_counter()
    elapsed_1_iter = end - start

    # math_functions.get_y_check(vandermonde,c)

    start = time.perf_counter()
    c_1, yyy = LU_10_iter(vandermonde, x, y)
    end = time.perf_counter()
    elapsed_10_iter = end - start


    M=20 #change to 3/5/8 for smooth function
    nevilles    = []
    lagrange    = []
    linear      = []
    x_points = []
    x_y = np.array([x, y])


    start = time.perf_counter()
    for point_index in range(0,len(x_points_to_interpolate)):

        point_request = x_points_to_interpolate[point_index]
        interpolated_values_nev     = interpolation.nevilles_algorithm(point_request,M,x_y)
    end = time.perf_counter()
    elapsed_nev = end - start

    #need to split this up for timing purposes
    for point_index in range(0,len(x_points_to_interpolate)):

        point_request = x_points_to_interpolate[point_index]

        interpolated_values_linear  = interpolation.linear_interpolation(point_request,x_y)
        interpolated_values_nev     = interpolation.nevilles_algorithm(point_request,M,x_y)
        interpolated_values_lag     = interpolation.lagrange_polynomial(point_request,19,x_y)
        interpolated_cubic_spline   = interpolation.cubic_spline(point_request,x_y)

        linear.append(interpolated_values_linear[1])
        nevilles.append(interpolated_values_nev[1])
        lagrange.append(interpolated_values_lag[1])
        x_points.append(x_points_to_interpolate[point_index])


    lagrange = np.array(lagrange)
    linear = np.array(linear)
    nevilles = np.array(nevilles)

    y_check=[]
    for point_request in x:
        interpolated_values_nev = interpolation.nevilles_algorithm(point_request,M,x_y)
        y_check.append(interpolated_values_nev[1])
    diff_nev = abs(y_check - y)

    y_check = []
    c_1, y_check = LU_1_iter(vandermonde, x, y,20)
    diff_1st = abs(y_check - y)

    y_check = []
    c_1, y_check = LU_10_iter(vandermonde, x, y,20)
    diff_10th = abs(y_check - y)

    # plotting
    fig, axes = plt.subplots(2,1)
    axes = axes.flatten()
    axes[0].scatter(x_points,yyy, marker = 'o', alpha = 0.4, label = '10th')
    axes[0].scatter(x_points,ys, marker = '^', alpha = 0.4, label = '1st')
    axes[0].scatter(x_points,linear, marker = 'o', alpha = 0.4, label = 'linear')
    axes[0].scatter(x_points,nevilles, marker = '^', alpha = 0.4, label = 'nevilles')
    axes[0].scatter(x_points,lagrange, marker = 'o', alpha = 0.2, label = 'lagrange')
    axes[0].set_ylim(-200,200)
    axes[0].set_title('Interpolated values')
    axes[0].scatter(x,y,label = 'true values')
    axes[0].legend()

    axes[1].plot(x,diff_nev, marker = '^', alpha = 0.4, label = 'nevilles')
    axes[1].plot(x,diff_1st ,marker = '^', alpha = 0.4, label = '1st')
    axes[1].plot(x,diff_10th, marker = 'o', alpha = 0.4, label = '10th')
    axes[1].set_yscale('log')
    axes[1].set_title('Error in log scale')
    axes[1].set_ylabel('|$y_i$ - y|')
    axes[1].legend()
    fig.suptitle('Results')
    plt.tight_layout()
    plt.savefig(os.path.join(plots_dir, "results.png"))

    #saving results
    np.savetxt(os.path.join(data_dir, 'c_1_iter.csv'), c, delimiter=',')
    np.savetxt(os.path.join(data_dir, 'c_10_iter.csv'), c_1, delimiter=',')
    np.savetxt(os.path.join(data_dir, '1_iter_time.txt'), [elapsed_1_iter], delimiter=',')
    np.savetxt(os.path.join(data_dir, '10_iter_time.txt'), [elapsed_10_iter], delimiter=',')
    np.savetxt(os.path.join(data_dir, 'nev_time.txt'), [elapsed_nev], delimiter=',')

    return c, c_1, elapsed_1_iter, elapsed_10_iter, elapsed_nev

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
print('C from 1 iteration of crouts algorithm\n',question_2_results[0])
print('\n')
print('C from 10 iterations of crouts algorithm\n',question_2_results[1])
print('\n')
print('Elapsed time 1 iteration of crouts algorithm (sec)\n',question_2_results[2])
print('\n')
print('Elapsed time 10 iterations of crouts algorithm (sec)\n',question_2_results[3])
print('\n')
print('Elapsed time for Neville\'s algorithm (sec)\n',question_2_results[4])
print('------------------------------------------------------')
print('\n')