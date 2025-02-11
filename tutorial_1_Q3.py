import requests
import numpy as np
import matplotlib.pyplot as plt
from PIL import Image
from io import BytesIO

def grab_array():
    # Fetch the image from the URL
    url = "https://home.strw.leidenuniv.nl/~daalen/NURfiles/M42_128.jpg"
    response = requests.get(url)

    # Ensure the request was successful
    if response.status_code == 200:
        # Convert the response content to a BytesIO object
        image_stream = BytesIO(response.content)

        # Open the image with PIL
        with Image.open(image_stream) as img:
            # Convert the image into a numpy array
            image_array = np.array(img)

        # Print the shape of the image array to confirm it was read correctly
        print("Image shape:", image_array.shape)

    else:
        print(f"Failed to retrieve image. Status code: {response.status_code}")

    return image_array

image_array = grab_array()

# Extract the first row (row index 0) from the image for the linear interpolation
first_row = image_array[0]

# Define the linear interpolation function
def linear_interpolator(x_0, y_0, x_1, y_1, x):
    # Handle division by zero issue
    if x_1 == x_0:
        return y_0

    #issue here if the values are not ints, get overflow not sure why, maybe stored as 8bit
    top_1 = int(y_1) - int(y_0)
    top_2 = x - x_0
    bottom = 1 / (x_1 - x_0)

    return (top_1 * top_2 * bottom) + y_1

#make a function for Neville's algorithm, pass the data through instead of data points as need to iterate over neighbours
#this is for the linear interpolation
def nevilles_algorithm(x_values,y_values,x):

    n = len(x_values)
    # Initialize a 2D table to hold the results of Neville's algorithm
    P = np.zeros((n, n))

    # Base case: Set the first row to the given y-values
    for i in range(n):
        P[i, 0] = y_values[i]

    # Recursively calculate P[i,j] for all i < j that have been passed
    for j in range(1, n):
        for i in range(n - j):
            top_1 = (x-x_values[i+j])*P[i,j-1]
            top_2 = (x-x_values[i])*P[i+1,j-1]
            inverse_bottom = 1/(x_values[i] - x_values[i + j])
            P[i, j] = (top_1-top_2)*inverse_bottom

    # The result is the value at P[0, n-1] format[col,row]
    return P[0,n-1]

#function for bisecting
def bisection(a,b,x):
    #while the difference between the upper and lower bound is greater than one continue
    while ((int(b)-int(a)) >= 1):
        #calculate the center point
        c = (a+b)*0.5
        #if the difference is one or less we've found the bounds
        if (int(b)-int(a) <= 1):
            break

        #decide which half to focus on based on where the center point lies
        if c >= x:
            b = int(c)
        else:
            a = int(c)
    return a,b

#this was a bit confusing but essentially you create a linspace with 201 steps for 128 pixels
def parse_201_values():

    # Create new equally spaced points
    new_points = np.linspace(0, 127, 201)

    # Initialize an array to store the interpolated values
    interpolated_values = []
    interpolated_values_neville = []
    # Define a small window for interpolation (e.g., 5 points)
    window_size = 10

    # Perform linear interpolation for each new point
    for new_x in new_points:

        # Find the two enclosing grid points using the bisection method
        a=0
        b=128
        x_0, x_1 = bisection(a,b,new_x)

        #not using the bisection method
        # x_0 = int(np.floor(new_x))  # Left grid point
        # x_1 = x_0 + 1  # Right grid point

        # Ensure x_1 does not exceed the bounds
        if x_1 >= len(first_row):
            x_1 = len(first_row) - 1

        # Get the pixel value from grid points after bisection
        y_0 = first_row[x_0]  # Value at the left grid point
        y_1 = first_row[x_1]  # Value at the right grid point

        # Get the integer indices for the points to interpolate
        x_values = np.arange(128)  # Original x-values (pixel indices) just making these new because CBA figuring out shape stuff
        y_values = first_row  # The corresponding pixel values for the first row

        # Determine a local window of points to use for Neville's algorithm
        center = int(new_x)
        start = max(0, center - window_size // 2)
        end = min(128, center + window_size // 2)

        # Interpolate using Neville's algorithm over the window otherwise get crazy values
        interpolated_value_nev = nevilles_algorithm(x_values[start:end], y_values[start:end], new_x)
        interpolated_values_neville.append((interpolated_value_nev))

        # Use the linear interpolater to find the value at new_x
        interpolated_value = linear_interpolator(x_0, y_0, x_1, y_1, new_x)
        interpolated_values.append(interpolated_value)

    # Convert to numpy array for easier plotting
    interpolated_values = np.array(interpolated_values)
    interpolated_values_neville = np.array(interpolated_values_neville)

    return new_points, interpolated_values, interpolated_values_neville

#call function to clean up code
new_points, interpolated_values, interpolated_values_neville = parse_201_values()

#make new interpolater for 2D function
def true_bilinear_interpolation(data,num_points):


def two_d_interpolation_approx(data,num_points):


# Plot the original row and the interpolated result
plt.subplot(2,1,1)
plt.plot(np.arange(128), first_row, label='Original row', marker='o', linestyle='-', color='b')
plt.plot(new_points, interpolated_values, label='Interpolated row (201 points)', linestyle='-', color='r')
plt.plot(new_points, interpolated_values_neville, label='Nevilles algorithm',linestyle='-',color='g')
plt.xlabel('Pixel index')
plt.ylabel('Pixel value')
plt.legend()
plt.title('Image Row Interpolation')
plt.subplot(2,1,2)
plt.plot()
plt.show()
