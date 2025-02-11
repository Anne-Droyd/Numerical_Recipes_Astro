import numpy as np
import matplotlib.pyplot as plt
import timeit

mean=10**6
sigma=10**5
sample_size=1000000
s=np.random.normal(mean,sigma,sample_size)

def schwarzchild_radius_mult(data):
    result=[]
    inverse_c_square=1/(299792458**2)
    G=6.6743*(10**(-11))
    for n in data:
        radius=2*G*n*inverse_c_square
        result.append(radius)
    return result

def schwarzchild_radius_div(data):
    result=[]
    c_square=(299792458**2)
    G=6.6743*(10**(-11))
    for n in data:
        radius=2*G*n/c_square
        result.append(radius)
    return result

# Function to time
def code_to_time_mult():
    schwarzchild_radius_mult(s)

# Function to time
def code_to_time_div():
    schwarzchild_radius_div(s)

# Use timeit to time the function
execution_time_mult = timeit.timeit(code_to_time_mult, number=1)
execution_time_div = timeit.timeit(code_to_time_div, number=1)
print(f"Multiplication execution time: {execution_time_mult} seconds")
print(f"Division execution time: {execution_time_div} seconds")