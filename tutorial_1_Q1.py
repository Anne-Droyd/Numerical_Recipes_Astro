#a)  Write a piece of code that implements the sinc(x) function in two ways:
# one using the power series expansion,
# and one using a library function, e.g.  np.sin(x).
# Clearly, the accuracy of the answer given by the power series expansion depends
# on the number of terms included. Which kind of error are we dealing with in such an example?


#make plot to show the error changing as terms change

import numpy as np
from scipy.special import factorial
import matplotlib.pyplot as plt

def factorial_mine(x):
    fact=1
    for n in range(2,x+1):
        fact *= n
    return fact

def sinc_function(data,terms):

    y_data=[]
    for x in data:
        sinc=0
        x=x*3.14
        # Handle x = 0 case directly to avoid division by zero
        if x == 0:
            sinc = 1
        else:
            for n in range(0,terms):
                bottom=(2*n)+1
                bottom_fact=factorial_mine(bottom)
                bottom_inverse=1/bottom_fact
                top_1=(-1)**n
                top_2=x**(2*n)
                top=top_1*top_2
                sinc = sinc + top*bottom_inverse
        y_data.append(sinc)

    return y_data


# plt.plot(np_data, np.sinc(np_data))

data=[]
data = np.linspace(-4, 4, 100)
all_data=[]
max_errors=[]
for m in range(1,50):
    sinc_data=sinc_function(data,m)
    all_data.append(sinc_data)
    sinc_through_np = np.sinc(data)
    max_errors.append(np.max(np.abs(sinc_through_np-sinc_data)))

plt.subplot(3,1,1)
plt.plot(data,sinc_data)
plt.title("My function")
plt.subplot(3,1,2)
plt.plot(data,sinc_through_np)
plt.title("Numpy function")
plt.subplot(3,1,3)
plt.plot(range(1,50),max_errors)
plt.title("Max error per # of terms")
plt.tight_layout()
plt.show()