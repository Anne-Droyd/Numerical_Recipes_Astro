import numpy as np
import math
from functions import linear_algebra as la

#handles 1d or 2d
def bisection(point_request:float,
              data_input:np.ndarray,
              points_around:int|None = None #this is value left/right ie 2 would be n-2 n-1 n n+1 n+2
              ) -> tuple[np.array,bool]:

    data = data_input.copy()
    edge_missing = False
    try:
        left_index = 0
        right_index = len(data[0]) - 1
        left_edge = data[0,left_index]
        right_edge = data[0,right_index]
    except:

        left_index = 0
        right_index = len(data) - 1

        left_edge = data[left_index]
        right_edge = data[right_index]

    # assume point is within range
    assert left_edge <= point_request <= right_edge, ("point out of range")

    #assume left dominant
    middle_index = math.floor((right_index-left_index) / 2)
    #continue until edges found
    while right_index - left_index > 1:

        # assume left dominant, update edges
        try:
            if point_request <= data[0,middle_index]:
                right_index = middle_index
            elif point_request > data[0,middle_index]:
                left_index = middle_index
        except:
            if point_request <= data[middle_index]:
                right_index = middle_index
            elif point_request > data[middle_index]:
                left_index = middle_index
        middle_index = math.floor((left_index+right_index)/2)

        #check stopping condition
        if right_index == left_index+1:

            #check if point falls on one of the edges, if so both points the same
            try:
                if point_request == data[0,left_index]:
                    right_index = left_index
                elif point_request == data[0,right_index]:
                    left_index = right_index
            except:
                if point_request == data[left_index]:
                    right_index = left_index
                elif point_request == data[right_index]:
                    left_index = right_index

    #if no boundary condition return edges


    if points_around is None:
        return np.array([data[0][left_index:right_index+1],data[1][left_index:right_index+1]]), edge_missing

    #else return list with boundary indices included
    else:

        assert points_around > 0, ('Points around is negative, required to be stritcly postitive')
        left_bound = left_index - points_around
        right_bound = right_index + points_around
        if left_bound < 0:
            left_bound = 0
            edge_missing = True
        if right_bound > len(data[0])-1:
            right_bound = len(data[0])-1
            edge_missing = True

        return np.array([data[0][left_bound:right_bound+1],data[1][left_bound:right_bound+1]]), edge_missing

def linear_interpolation(point_request:float,
                         data = np.ndarray,
                         ) -> np.ndarray[float,float]:

    points, egde = bisection(point_request,data)
    if points.shape[1] == 1:
        return np.array([points[0,0],points[1,0]])
    else:
        right_point = points[:,0]
        left_point = points[:,1]
    fraction = (right_point[1] - left_point[1])/(right_point[0] - left_point[0])
    interpolated_value = fraction*(point_request - left_point[0]) + left_point[1]
    return np.array([point_request,interpolated_value])

def lagrange_polynomial(point_request:float,
                        degree_of_polynomial:int,
                        full_data:np.ndarray
                        ) -> np.ndarray[float,float]:
    # assert degree_of_polynomial <= 5, ('Degree of polynomial should be less than or equal to 5')
    assert degree_of_polynomial > 0, ('Degree of polynomial must be strictly postitive and greater than 0')

    data = full_data.copy()

    nearby_samples, edge_case = bisection(point_request,data,40)
    # if edge_case:
    #     return np.array([point_request, None])
    nearby_data = nearby_samples.copy()
    interpolated_value = 0

    for i in range(degree_of_polynomial+1):
        products = nearby_data[1,i]
        for j in range(degree_of_polynomial+1):
            if j == i:
                continue
            product = (point_request - nearby_data[0,j])/(nearby_data[0,i] - nearby_data[0,j])
            products *=product

        interpolated_value += products

    return np.array([point_request, interpolated_value])

def nevilles_algorithm(point_request:float,
                       M_points_around:int,
                       data:np.ndarray
                       ) -> np.ndarray[float,float]:

    nearby_samples, edge_case = bisection(point_request,data,M_points_around)
    #save only ys
    p_i = nearby_samples[1,:].copy()

    # tried to follow the slides perfectly, used https://github.com/gisalgs/geom/blob/master/neville.py
    # to adjust my indices and sums
    for k in range(1,len(p_i)):
        for i in range(len(p_i)-k):
            j= i + k
            top_1 = (point_request - nearby_samples[0,j])*p_i[i]
            top_2 = (nearby_samples[0,i] - point_request)*p_i[i+1]
            bottom = 1/(nearby_samples[0,i]-nearby_samples[0,j])
            p_i[i] = (top_1+top_2)*bottom

    return np.array([point_request,p_i[0]])

    # heres what I had before not super different just wasn't summing correct
    # def nevilles_algorithm(point_request:float,
    #                        M_points_around:int,
    #                        data:np.ndarray
    #                        ) -> tuple[float,float]:
    #     nearby_samples, edge_case = bisection(point_request,data,M_points_around)
    #     print(nearby_samples)
    #     #save only ys
    #     p_i = nearby_samples[:,1].copy()
    #     #find the closest pair or single x value
    #     nearest_pair, _ = bisection(point_request,data)
    #     if len(nearest_pair[:,0]) > 1:
    #         closest = check_closest(nearest_pair,point_request)
    #         initial_solution = closest[1]
    #     else:
    #         initial_solution = nearest_pair[0,1]
    #
    #     for k in range(1,len(p_i)):
    #
    #         print(p_i)
    #         for i in range(0,len(p_i)-k):
    #             j= i + k
    #             top_1 = (nearby_samples[j,0]-point_request)*p_i[i]
    #             top_2 = (point_request-nearby_samples[i,0])*p_i[j]
    #             bottom = 1/(nearby_samples[j,0]-nearby_samples[i,0])
    #             p_i[i] = (top_1+top_2)*bottom
    #     print(p_i)
    #     print(point_request)
    #     print(nearby_samples)
    #     print(p_i)
    #

def cubic_spline(point_request:float,
                 data:np.ndarray
                 ) -> np.ndarray[float,float]:
    h=0
    # requires LA, do later

#todo
# def cubic_spline():
# def akima_spline():
# def bilinear_interpolation_1():
# def bilinear_interpolation_2():
# def bicubic_interpolation():

