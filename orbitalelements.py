import numpy as np
import vpython 
from math import *

def findAngle(cosine,sine):
    exp = cosine ** 2 + sine ** 2
    if(abs(exp ** 0.5 - 1) > 0.01):
        return "User error: invalid input"
    angleCos = acos(cosine)
    if (sine < 0):
        return 2*pi - angleCos
    return angleCos

def angMomentum(r2, r2dot): # r2 and r2dot
    r2 = np.array([r2.x, r2.y, r2.z])
    r2dot = np.array([r2dot.x, r2dot.y, r2dot.z])
    h = np.cross(r2, r2dot)
    h = vpython.vector(h[0], h[1], h[2])
    return h

def true_anomaly(E, A, r2, r2dot): # e, a, r2, r2dot
    exp = ((A/r2.mag)*(1-E**2)-1)/E
    return acos(exp)

def incline(h): # angular momentum
    exp = ((h.x**2 + h.y**2) ** 0.5)/h.z
    return atan(exp)

def longitude_of_ascending_node(h, i): # angular momentum, incline
    i = radians(i)
    exp = (h.x)/(h.mag * sin(i))
    cosine = -1 * (h.y)/(h.mag * sin(i))
    sine = exp
    angleCos = acos(cosine)
    if (sine < 0):
        return 2*pi - angleCos
    return angleCos

def a(r2, r2dot): # r2 and r2dot
    return 1/((2/r2.mag) - r2dot.mag**2)

def e(A, h): # a, angular momentum
    return (1 - (h.mag**2/A))**0.5

def argument_of_perihelion(omega, i, r2, TA, A, E, h, r2dot): # omega, incline, r2.z, r2, true anomaly
    omega = radians(omega)
    i = radians(i)
    sine = (r2.z)/(r2.mag * sin(i))
    cosine = (r2.x * cos(omega) + r2.y * sin(omega))/r2.mag
    angleCos = acos(cosine)
    if (sine < 0):
        U=  2*pi - angleCos
    else:
        U = angleCos
    r2mag = r2.mag
    r2 =  np.array([r2.x, r2.y, r2.z])
    r2dot = np.array([r2dot.x, r2dot.y, r2dot.z])
    sine = (A * (1-E**2)/(E * h.mag)) * ((np.dot(r2, r2dot))/r2mag)
    cosine = (1/E)*(A*(1-E**2)/r2mag - 1)
    angleCos = acos(cosine)
    if (sine < 0):
        TA=  2*pi - angleCos
    else:
        TA = angleCos
    if(U - TA < 0):
        return U - TA + 2*pi
    return U - TA

def mean_anomaly(E, Eccen): # e, eccentric anomaly
    return Eccen - E*sin(Eccen) 

def eccentric_anomaly(E, r2, A, TA): # e, r2, r2dot, a
    exp = (1/E)*(1 - (r2.mag/A))
    if(TA >= 0 and TA <= pi):
        return 2*pi - acos(exp)
    return 2*pi - acos(exp)

def percent_error(calculated, expected):
    return str((abs(calculated - expected)/expected) * 100) + "%"

def do_rotation(coords): # coords must be a NUMPY ARRAY
    coords = np.array([coords.x, coords.y, coords.z])
    epsilon = radians(23.4374)
    toEcliptic = np.linalg.inv([[1,0,0],
                    [0,cos(epsilon),-sin(epsilon)],
                    [0,sin(epsilon),cos(epsilon)]])
    toEcliptic = np.array(toEcliptic)
    coords = np.matmul(toEcliptic, coords)
    coords = vpython.vector(coords[0], coords[1], coords[2])

    return coords

def find_JD(year, month, day, time):
    JulianDay = (367 * year) - (7 * (year + (month+9)//12)//4) + (275 * month)//9 + day + 1721013.5
    time = time.split(":")
    decimal_time = 0
    if(float(time[0]) < 0):
        decimal_time += float(time[0]) - float(time[1])/60 - float(time[2])/3600
    else:
        decimal_time += float(time[0]) + float(time[1])/60 + float(time[2])/3600
    decimal_time/=24
    JulianDay += decimal_time
    return JulianDay

def M(t2, t1, M2, A):
    k = 0.01720209894
    n = sqrt(k**2/(A**3))
    return n*(t1 - t2) + M2
