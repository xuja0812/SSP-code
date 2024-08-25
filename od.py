import numpy as np
from math import *
import vpython as vp
import scipy

def process_inputs(year, month, day, time, RA, dec, Rx, Ry, Rz):
    JD = find_JD(year, month, day, time)
    RA = find_RA(RA)
    dec = find_dec(dec)
    R = vp.vector(Rx, Ry, Rz)
    return [JD, RA, dec, R]

def find_JD(year, month, day, time):
    JulianDay = (367 * year) - (7 * (year + (month+9)//12)//4) + (275 * month)//9 + day + 1721013.5
    time = time.split(":")
    decimal_time = 0
    if(float(time[0]) < 0):
        decimal_time += float(time[0]) - float(time[1])/60 - float(time[2])/3600
    else:
        decimal_time += float(time[0]) + float(time[1])/60 + float(time[2])/3600
    decimal_time/=24
    # JulianDay += decimal_time
    return JulianDay

def find_RA(RA):
    RA = RA.split(":")
    deg = float(RA[0]) 
    min = float(RA[1]) 
    sec = float(RA[2]) 
    min = copysign(min, deg)
    sec = copysign(sec, deg)
    return radians((deg + min/60 + sec/3600) * 15)

def find_dec(dec):
    dec = dec.split(":")
    deg = float(dec[0]) 
    min = float(dec[1]) 
    sec = float(dec[2]) 
    min = copysign(min, deg)
    sec = copysign(sec, deg)
    return radians(deg + min/60 + sec/3600)

def rho_hat(RA, dec):
    phat = vp.vector(cos(RA)*cos(dec), sin(RA)*cos(dec), sin(dec))
    return phat

def find_taus(t1, t2, t3):
    k = 0.0172020989484
    tau3 = k*(t3 - t2)
    tau1 = k*(t1 - t2)
    tau = tau3 - tau1
    return [tau1, tau3, tau]

def find_Ds(rhohat1, R1, R2, R3, rhohat3, rhohat2):
    rhohat1 = np.array([rhohat1.x, rhohat1.y, rhohat1.z])
    rhohat2 = np.array([rhohat2.x, rhohat2.y, rhohat2.z])
    rhohat3 = np.array([rhohat3.x, rhohat3.y, rhohat3.z])
    R1 = np.array([R1.x, R1.y, R1.z])
    R2 = np.array([R2.x, R2.y, R2.z])
    R3 = np.array([R3.x, R3.y, R3.z])
    D0 = np.dot(rhohat1, np.cross(rhohat2, rhohat3))
    D21 = np.dot((np.cross(rhohat1, R1)), rhohat3)
    D22 = np.dot((np.cross(rhohat1, R2)), rhohat3)
    D23 = np.dot((np.cross(rhohat1, R3)), rhohat3)
    return [D0, D21, D22, D23]

def SEL(taus, R2, rhohat2, Ds):
    D0, D21, D22, D23 = Ds
    tau1, tau3, tau = taus
    A1 = tau3/tau
    B1 = (A1/6) * (tau**2 - tau3**2)
    A3 = (-1*tau1)/tau
    B3 = (A3/6)*(tau**2 - tau1**2)
    A = (A1*D21 - D22 + A3*D23)/(-D0)
    B = (B1*D21 + B3*D23)/(-D0)

    R2 = np.array([R2.x,R2.y,R2.z])
    rhohat2 = np.array([rhohat2.x,rhohat2.y,rhohat2.z])

    E = -2*(np.dot(R2, rhohat2))

    rhohat2 = vp.vector(rhohat2[0], rhohat2[1], rhohat2[2])
    R2 = vp.vector(R2[0], R2[1], R2[2])

    F = R2.mag ** 2

    a = -(A**2 + A*E + F)
    b = -(2*A*B + B*E)
    c = -(B**2)

    coeffs = (c, 0, 0, b, 0, 0, a, 0, 1)
    res = np.polynomial.polynomial.polyroots(coeffs)
    res = res.real

    temp = []
    for num in res:
        if num > 0:
            temp.append(num)

    res = np.array(temp)

    u2 = 1/(res[1]**3)
    c1 = (tau3/tau)*(1 + (u2/6)*(tau**2 - tau3**2))
    
    rhos = []
    for root in res:
        rhos.append(A + (B/(root**3)))
    
    return res, np.array(rhos)

def initial_f_g(tau1, tau3, r2): # second order; 
    # convert r2 into a numpy array 
    f1 = 1 - 1/(2*r2**3)*tau1**2
    g1 = tau1 - 1/(6*r2**3) * tau1**3
    f3 = 1 - 1/(2*r2**3)*tau3**2
    g3 = tau3 - 1/(6*r2**3) * tau3**3
    return f1, g1, f3, g3

# find c1 c2 c3
# find rhomag1 and rhomag3 --> combine with rhohat1 and rhohat3 to get rho1 and rho3
# find r1 and r3 using vector subtraction
# use equation 106 to find r2dot
# keep iterating using f and g to update rho2 ????

def initial_r2dot(tau1, tau3, r2, R1, R2, R3, rhohat1, rhohat2, rhohat3): # r2 is a vector
    rhohat1 = np.array([rhohat1.x, rhohat1.y, rhohat1.z])
    rhohat2 = np.array([rhohat2.x, rhohat2.y, rhohat2.z])
    rhohat3 = np.array([rhohat3.x, rhohat3.y, rhohat3.z])
    R1 = np.array([R1.x, R1.y, R1.z])
    R2 = np.array([R2.x, R2.y, R2.z])
    R3 = np.array([R3.x, R3.y, R3.z])

    f1, g1, f3, g3 = initial_f_g(tau1, tau3, r2)
    c1 = g3/(f1*g3 - g1*f3)
    c2 = -1
    c3 = -g1/(f1*g3 - g1*f3)

    # convert R vectors to numpy arrays

    D0 = np.dot(rhohat1, np.cross(rhohat2, rhohat3))
    D11 = np.dot(np.cross(R1, rhohat2),rhohat3)
    D12 = np.dot(np.cross(R2, rhohat2),rhohat3)
    D13 = np.dot(np.cross(R3, rhohat2),rhohat3)

    D31 = np.dot(rhohat1, np.cross(rhohat2, R1))
    D32 = np.dot(rhohat1, np.cross(rhohat2, R2))
    D33 = np.dot(rhohat1, np.cross(rhohat2, R3))

    Ds = [D0, D11, D12, D13, D31, D32, D33]

    rho1mag = (c1*D11 + c2*D12 + c3*D13)/(c1*D0)
    rho3mag = (c1*D31 + c2*D32 + c3*D33)/(c3*D0)

    rho1 = rhohat1 * rho1mag
    rho3 = rhohat3 * rho3mag

    d1 = -f3/(f1*g3 -  f3*g1)
    d3 = f1/(f1*g3 - f3*g1)

    # rhos and Rs are numpy arrays

    r1 = rho1 - R1
    r3 = rho3 - R3

    initial_r2dot = d1*r1 + d3*r3
    r2 = c1*r1 + c3*r3
    initial_r2dot = vp.vector(initial_r2dot[0], initial_r2dot[1], initial_r2dot[2])
    return initial_r2dot, r2, rho1mag, rho3mag

def iterate(r2dot, r2, tau1, tau3, R2, rho2, R3, R1, rhohat1, rhohat2, rhohat3, t10, t20, t30, rho1mag2, rho2mag2, rho3mag2): # all parameters are vectors except for the taus
    past_rho2mag = rho2.mag - 1
    k = 0.0172020989484
    current_rho2mag = rho2.mag
    count = 0
    rhohat1 = np.array([rhohat1.x, rhohat1.y, rhohat1.z])
    rhohat2 = np.array([rhohat2.x, rhohat2.y, rhohat2.z])
    rhohat3 = np.array([rhohat3.x, rhohat3.y, rhohat3.z])
    R1 = np.array([R1.x, R1.y, R1.z])
    R2 = np.array([R2.x, R2.y, R2.z])
    R3 = np.array([R3.x, R3.y, R3.z])

    cAU = 173.144643267
    while(abs(past_rho2mag - current_rho2mag) > 1E-12):
        t1 = t10 - rho1mag2/cAU
        t2 = t20 - rho2mag2/cAU
        t3 = t30 - rho3mag2/cAU
        tau3 = k*(t3 - t2)
        tau1 = k*(t1 - t2)
        count+=1
        # r2 needs to be a VECTOR 
        r2mag = r2.mag
        A = a(r2, r2dot)
        n = sqrt(1/(A**3))
        delta_E1 = delta_e(r2, A, r2dot, n, tau1)
        delta_E3 = delta_e(r2, A, r2dot, n, tau3)
        f1 = 1 - (A/r2mag)*(1 - cos(delta_E1))
        g1 = tau1 + (1/n)*(sin(delta_E1) - delta_E1)
        f3 = 1 - (A/r2mag)*(1 - cos(delta_E3))
        g3 = tau3 + (1/n)*(sin(delta_E3) - delta_E3)
        c1 = g3/(f1*g3 - g1*f3)
        c2 = -1
        c3 = -g1/(f1*g3 - g1*f3)

        D0 = np.dot(rhohat1, np.cross(rhohat2, rhohat3))
        D11 = np.dot(np.cross(R1, rhohat2),rhohat3)
        D12 = np.dot(np.cross(R2, rhohat2),rhohat3)
        D13 = np.dot(np.cross(R3, rhohat2),rhohat3)

        D21 = np.dot((np.cross(rhohat1, R1)), rhohat3)
        D22 = np.dot((np.cross(rhohat1, R2)), rhohat3)
        D23 = np.dot((np.cross(rhohat1, R3)), rhohat3)

        D31 = np.dot(rhohat1, np.cross(rhohat2, R1))
        D32 = np.dot(rhohat1, np.cross(rhohat2, R2))
        D33 = np.dot(rhohat1, np.cross(rhohat2, R3))

        rho1mag = (c1*D11 + c2*D12 + c3*D13)/(c1*D0)
        rho3mag = (c1*D31 + c2*D32 + c3*D33)/(c3*D0)
        rho2mag = (c1*D21 + c2*D22 + c3*D23)/(c2*D0)

        rho1 = rhohat1 * rho1mag
        rho3 = rhohat3 * rho3mag
        rho2 = rhohat2 * rho2mag

        d1 = -f3/(f1*g3 -  f3*g1)
        d3 = f1/(f1*g3 - f3*g1)
        r1 = rho1 - R1
        r3 = rho3 - R3
        r2dot = d1*r1 + d3*r3
        r2dot = vp.vector(r2dot[0], r2dot[1], r2dot[2])
        r2 = c1*r1 + c3*r3
        r2 = vp.vector(r2[0], r2[1], r2[2])
        past_rho2mag = current_rho2mag
        current_rho2mag = rho2mag

    return r2dot, r2, current_rho2mag

def a(r2, r2dot):
    position = r2.mag
    velocity = r2dot.mag
    return 1/((2/position) - velocity**2)

def delta_e(r2, a, r2dot, n, tau):
    x0 = n*tau
    r2mag = r2.mag
    r2 = np.array([r2.x, r2.y, r2.z])
    r2dot = np.array([r2dot.x, r2dot.y, r2dot.z])
    func = lambda x : x - (1 - (r2mag)/a)*sin(x) + (np.dot(r2, r2dot)/(n*a**2))*(1 - cos(x)) - n*tau
    newt = scipy.optimize.newton(func, x0, tol=1E-12)
    return newt
