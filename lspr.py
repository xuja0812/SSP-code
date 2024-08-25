import numpy as np
import astropy.io
import matplotlib.pyplot as plt
import math

def lspr(x, y):
    data = np.loadtxt("D:\\Xu_LSPR_Input - Copy (2).txt", dtype=str)
    xs = np.array(data[:,0])
    xs = xs.astype(float)
    ys = np.array(data[:,1])
    ys = ys.astype(float)
    RAs = np.array(data[:,2])
    decs = np.array(data[:,3])
    RA_sum = RA_value_sum(RAs)
    RAx_sum = RA_value_sum(RAs,xs)
    RAy_sum = RA_value_sum(RAs,[],ys)
    dec_sum = dec_value_sum(decs)
    decx_sum = dec_value_sum(decs,xs)
    decy_sum = dec_value_sum(decs,[],ys)
    XandYs = xs * ys
    XSquareds = xs * xs
    YSquareds = ys * ys
    XYsum = np.sum(XandYs)
    Xsum = np.sum(xs)
    Ysum = np.sum(ys)
    X2sum = np.sum(XSquareds)
    Y2sum = np.sum(YSquareds)
    
    # finally actually doing and solving the matrix multiplication
    
    mainMat = np.array([[10,Xsum,Ysum],
                        [Xsum,X2sum,XYsum],
                        [Ysum,XYsum,Y2sum]])
    # print(mainMat)
    
    mainMat = np.linalg.inv(mainMat)
    
    RAMat = np.array([[RA_sum],[RAx_sum],[RAy_sum]])
    # print(RAMat)
    
    RA_constants = np.dot(mainMat, RAMat)
    # print(RA_constants)
    
    decMat = np.array([[dec_sum],[decx_sum],[decy_sum]])
    
    dec_constants = np.dot(mainMat, decMat)
    # print(dec_constants)

    b1 = RA_constants[0,0]
    a11 = RA_constants[1,0]
    a12 = RA_constants[2,0]
    b2 = dec_constants[0,0]
    a21 = dec_constants[1,0]
    a22 = dec_constants[2,0]
    ra = b1 + a11*x + a12*y
    dec = b2 + a21*x + a22*y

    chi2RA = 0
    chi2dec = 0
    for i in range(len(xs)):
        chi2RA += (RA_value(RAs[i]) - b1 - a11*xs[i] - a12*ys[i]) ** 2
        chi2dec += (dec_value(decs[i]) - b2 - a21*xs[i] - a22*ys[i]) ** 2

    errorRA = math.sqrt(chi2RA/7) * 3600
    errordec = math.sqrt(chi2dec/7) * 3600


    return hours_minutes_seconds(ra), degrees_arcmins_arcsecs(dec), round(errorRA,2), round(errordec,2)
    
    # RA = RA_constants[0] + 
    
def RA_value(RA):
    list = RA.split(":")
    list = np.array(list)
    list = list.astype(float)
    return (list[0] + (list[1]/60) + (list[2]/3600)) * 15

def dec_value(dec):
    list = dec.split(":")
    list = np.array(list)
    list = list.astype(float)
    return list[0] + (list[1]/60) + (list[2]/3600)
    
    
def RA_value_sum(RA_list, x_list=[], y_list=[]):
    sum = 0
    for i in range(len(RA_list)):
        list = RA_list[i].split(":")
        list = np.array(list)
        list = list.astype(float)
        if(len(x_list)>0):
            sum+=(list[0] + (list[1]/60) + (list[2]/3600)) * x_list[i] * 15
        elif(len(y_list)>0):
            sum+=(list[0] + (list[1]/60) + (list[2]/3600)) * y_list[i] * 15
        else:
            sum+=(list[0] + (list[1]/60) + (list[2]/3600)) * 15
    return sum

def dec_value_sum(dec_list, x_list=[], y_list=[]):
    sum = 0
    for i in range(len(dec_list)):
        list = dec_list[i].split(":")
        list = np.array(list)
        list = list.astype(float)
        if(len(x_list)>0):
            sum+=(list[0] + (list[1]/60) + (list[2]/3600)) * x_list[i]
        elif(len(y_list)>0):
            sum+=(list[0] + (list[1]/60) + (list[2]/3600)) * y_list[i]
        else:
            sum+=list[0] + (list[1]/60) + (list[2]/3600)
    return sum

def hours_minutes_seconds(RA):
    RA /= 15
    hours = RA//1
    minutes = ((RA - hours) * 60)//1
    seconds = (((RA - hours) * 60) - minutes) * 60
    return hours, minutes, round(seconds,2)

def degrees_arcmins_arcsecs(dec):
    degrees = dec//1
    arcminutes = ((dec - degrees) * 60)//1
    arcseconds = (((dec - degrees) * 60) - arcminutes) * 60
    return degrees, arcminutes, round(arcseconds,2)

RA, dec, RA_uncertainty, dec_uncertainty = lspr(942.63423,967.02236)
print("RA:",RA)
print("dec:",dec)
print("RA uncertainty:",RA_uncertainty,"arcseconds")
print("dec uncertainty:",dec_uncertainty,"arcseconds")
