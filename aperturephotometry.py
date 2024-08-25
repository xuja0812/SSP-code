import numpy as np
import vpython as vp
from math import *
import matplotlib.pyplot as plt
from astropy.io import fits

def findCentroid(fits_file, target_x, target_y, radius, sky_radius):
    data = fits.getdata("D:\\Benedicts_07032023\\better_2002KL6_0703_001light120.fits")
    # plt.imshow(data[target_y-radius:target_y+radius+1, target_x-radius:target_x+radius+1])
    # plt.show()
    totalRadius = radius+sky_radius
    pixelCountSum = 0
    sky_sum = 0
    sky_num = 0
    ysum = 0
    xsum = 0
    num = 0
    for i in range(target_y-totalRadius,target_y+totalRadius + 1):
        for j in range(target_x-totalRadius,target_x+totalRadius + 1):
            if(sqrt((i-target_y)**2 + (j-target_x)**2) <= radius):
                pixelCountSum += data[i, j]
                num+=1
            elif(sqrt((i-target_y)**2 + (j-target_x)**2) <= totalRadius):
                sky_sum += data[i, j]
                sky_num += 1
    sky_avg = sky_sum/sky_num
    pixelCountSum -= (num * sky_avg)

    for i in range(target_x-totalRadius,target_x+totalRadius + 1):
        for j in range(target_y-totalRadius,target_y+totalRadius + 1):
            if(sqrt((j-target_y)**2 + (i-target_x)**2) <= radius):
                xsum += (data[j, i] - sky_avg) * i

    for i in range(target_y-totalRadius,target_y+totalRadius + 1):
        for j in range(target_x-totalRadius,target_x+totalRadius + 1):
            if(sqrt((i-target_y)**2 + (j-target_x)**2) <= radius):
                ysum += (data[i, j] - sky_avg) * i

    x_centroid = (xsum)/pixelCountSum
    y_centroid = (ysum)/pixelCountSum

    return x_centroid, y_centroid

def aperture_photometry(img_name, x, y, r_ap, r_out, r_in, R, D):
    img = fits.getdata("D:\\Benedicts_07032023\\better_2002KL6_0703_001light120.fits")
    x, y = findCentroid(img_name,x,y,r_ap, r_out)
    X, Y = x, y
    width = len(img[0])
    length = len(img)
    temp_arr = np.array([[0]*width]*length)
    x //= 1
    y //= 1
    x = int(x)
    y = int(y)
    print("CENTROID X:",X)
    print("CENTROID Y:",Y)
    ADU_ap = 0
    n_ap= 0
    n_an = 0
    an_values = []
    for i in range(y-r_out, y+r_out+1):
        for j in range(x-r_out, x+r_out+1):
            if(sqrt((i-Y)**2 + (j-X)**2) <= r_ap):
                ADU_ap += img[i,j]
                temp_arr[i,j]+=1
                n_ap+=1
            elif(sqrt((i-Y)**2 + (j-X)**2) > r_in and sqrt((i-Y)**2 + (j-X)**2) <= r_out):
                an_values.append(img[i,j])
                temp_arr[i,j]+=1
                n_an +=1
    an_values.sort()
    median = np.median(an_values)

    signal = ADU_ap - median*n_ap
    m_inst = -2.5 * log10(signal)
    exp = sqrt(signal + n_ap*(1 + (n_ap/n_an)) * (median + D + R**2))
    SNR = signal/exp

    uncertainty = 1.0875/SNR
    signal_uncertainty = signal/SNR

    print("SIGNAL:",signal)
    print("SIGNAL UNCERTAINTY:",signal_uncertainty)
    print("SNR:",SNR)
    print("M_inst:",m_inst)
    print("INST MAG UNCERTAINTY:",uncertainty)
    return signal
