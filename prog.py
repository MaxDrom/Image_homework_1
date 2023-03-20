from astropy.io import fits
import numpy as np
import math
import matplotlib.pyplot as plt

def is_in_bounds(i, j, data):
    result = (i>=0) and (i<np.shape(data)[0]) 
    result = result and (j>=0) and (j<np.shape(data)[1]) 
    return result

def I(x, y, data):
    i = math.floor(x)
    j = math.floor(y)
    f_x = x-i
    f_y = y-j
    result = 0
    areas = {}
    areas[0] = lambda x: 1-x
    areas[1] = lambda x: x
    for k in range(i, i+2):
        for l in range(j, j+2):
            if is_in_bounds(k, l, data):
                result+= areas[k-i](f_x)*areas[l-j](f_y)*data[k, l]
    return result

def do_plot_data(N, x_0, x_1, f):
    result_y = []
    result_x = []
    for i in range(N):
        x = x_0+(x_1-x_0)/(N-1)*i
        result_x.append(x)
        result_y.append(f(x))
    return (result_x, result_y)

def circle_I(R, N_r, N_phi, data):
    result = 0
    dr = R/N_r
    dphi = 2*math.pi/N_phi
    phi = 0
    for _ in range(N_phi):
        x_0 = math.cos(phi)
        y_0 = math.sin(phi)
        phi += dphi
        r = 0
        for __ in range(1,N_r+1):
            r += dr
            x = r*x_0+100
            y = r*y_0+100
            result+= r*I(x, y, data)*dr*dphi
    return result

with fits.open("data.fits") as hdu_list:
    data = [hdu.data for hdu in hdu_list[1:]]

mean_data = np.mean(data, axis = 0)
median_data = np.median(data, axis = 0)

data_dict = {}
data_dict["Mean"] = mean_data
data_dict["Median"] = median_data
count = 1
for title in data_dict.keys():
    plt.subplot(2, 2, count)
    plt.title(title+" slice")
    plt.xlabel("$x$")
    plt.ylabel("$I$")
    plt.plot(*do_plot_data(1000,0,199,lambda x: I(x,100, data_dict[title])))

    plt.subplot(2, 2, count+2)
    plt.title(title+" growth curve")
    plt.xlabel("$r$")
    plt.ylabel("$I$")
    plt.plot(*do_plot_data(100,0,99,lambda x: circle_I(x,100, 5, data_dict[title])))
    count+=1

plt.show()

