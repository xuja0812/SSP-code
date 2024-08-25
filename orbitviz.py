import txaio
txaio.use_asyncio()
from vpython import *
import numpy as np
from numpy import *
import math

a = 2.773017979589484
e = 0.1750074901308245
M = radians(336.0050001501443)
Oprime = radians(108.032597191534)
iprime = radians(16.34548466739393)
wprime = radians(74.95130563682554)


def solvekep(M):
    Eguess = M
    Mguess = Eguess - e*sin(Eguess)
    Mtrue = M
    while abs(Mguess - Mtrue) > 1e-004:
        Mguess = Eguess - e*sin(Eguess)
        Eguess = Eguess - (Mtrue - (Eguess - e*sin(Eguess))) / (e*cos(Eguess) - 1)
        # Eguess = Eguess - (Eguess - e*sin(Eguess) - Mtrue) / (1 - e*cos(Eguess))
    return Eguess

sqrtmu = 0.01720209895
mu = sqrtmu**2
time = 0
dt = .05
period = sqrt(4*pi**2*a**3/mu)
# r1ecliptic = vector(0, 0, 0)
Mtrue = 2*pi/period*(time) + M
Etrue = solvekep(Mtrue)
#
r1ecliptic = vector(a*cos(Etrue) - a*e, a*sqrt(1-e**2) * sin(Etrue), 0)
r1eclipticList = [[r1ecliptic.x],[r1ecliptic.y],[r1ecliptic.z]]
r1ecliptic = array(r1eclipticList)
print(r1ecliptic)
r1ecliptic = vector(r1ecliptic[0,0], r1ecliptic[1,0], r1ecliptic[2,0])
print(r1ecliptic)
r1eclipticList = [[r1ecliptic.x],[r1ecliptic.y],[r1ecliptic.z]]
r1ecliptic = array(r1eclipticList)
#
epsilon = radians(23.4)

# spinning the r vector

firstRotation = [[cos(Oprime),-sin(Oprime), 0],[sin(Oprime),cos(Oprime),0],[0,0,1]]
firstRotation = array(firstRotation)
print(firstRotation)
secondRotation = [[1, 0, 0],[0,cos(iprime),-sin(iprime)],[0,sin(iprime),cos(iprime)]]
secondRotation = array(secondRotation)
print(secondRotation)
thirdRotation = [[cos(wprime),-sin(wprime), 0],[sin(wprime),cos(wprime),0],[0,0,1]]
thirdRotation = array(thirdRotation)
print(thirdRotation)

# vec = [[r1ecliptic.x],[r1ecliptic.y],[r1ecliptic.z]]
# vec = array(vec)

# r1ecliptic.x = firstRotation * secondRotation * thirdRotation * r1ecliptic.x
# r1ecliptic.y = firstRotation * secondRotation * thirdRotation * r1ecliptic.y
# r1ecliptic.z = firstRotation * secondRotation * thirdRotation * r1ecliptic.z
# r1ecliptic = matmul(r1ecliptic, firstRotation)
# r1ecliptic = matmul(r1ecliptic, secondRotation)
# r1ecliptic = matmul(r1ecliptic, thirdRotation)

firstRotation = matmul(firstRotation, secondRotation)
print(firstRotation)
firstRotation = matmul(firstRotation, thirdRotation)
print(firstRotation)
r1ecliptic = matmul(firstRotation, r1ecliptic)
# r1ecliptic.x = vec[0,0]
# r1ecliptic.y = vec[1,0]
# r1ecliptic.z = vec[2,0]
print(r1ecliptic)

# r1ecliptic.x = (cos(wprime)*cos(Oprime) - sin(wprime)*cos(iprime)*sin(Oprime))*(a*cos(Etrue)-a*e) - (cos(wprime)*cos(iprime)*sin(Oprime) + sin(wprime)*cos(Oprime))*(a*sqrt(1-e**2)*sin(Etrue))
# r1ecliptic.y = (cos(wprime)*sin(Oprime) + sin(wprime)*cos(iprime)*cos(Oprime))*(a*cos(Etrue)-a*e) + (cos(wprime)*cos(iprime)*cos(Oprime) - sin(wprime)*sin(Oprime))*(a*sqrt(1-e**2)*sin(Etrue))
# r1ecliptic.z = sin(wprime)*sin(iprime)*(a*cos(Etrue)-a*e) + cos(wprime)*sin(iprime)*(a*sqrt(1-e**2)*sin(Etrue))

# convert r to equatorial 

# vec = [[r1ecliptic.x],[r1ecliptic.y],[r1ecliptic.z]]
# vec = array(vec)
toEquatorial = [[1,0,0],[0,cos(epsilon),-sin(epsilon)],[0,sin(epsilon),cos(epsilon)]]
toEquatorial = array(toEquatorial)
r1ecliptic = matmul(toEquatorial, r1ecliptic)
# r1ecliptic.x = vec[0,0]
# r1ecliptic.y = vec[1,0]
# r1ecliptic.z = vec[2,0]
print(r1ecliptic)

r1ecliptic = vector(r1ecliptic[0,0], r1ecliptic[1,0], r1ecliptic[2,0])
asteroid = sphere(pos=r1ecliptic*150, radius=(15), color=color.white)
asteroid.trail = curve(color=color.white)
sun = sphere(pos=vector(0,0,0), radius=(50), color=color.yellow)
# r1eclipticList = [[r1ecliptic.x],[r1ecliptic.y],[r1ecliptic.z]]
# r1ecliptic = array(r1eclipticList)

# finding RA and dec
# June 18 2023 4:00 UTC

R = [-1.090645644786802 * 10 ** 7, -1.518916634843708 * 10 ** 8, 3.977208618997782 * 10 ** 4]
R = array(R)
px = math.sqrt(r1ecliptic.x ** 2 + R[0] ** 2)
py = math.sqrt(r1ecliptic.y ** 2 + R[1] ** 2)
pz = math.sqrt(r1ecliptic.z ** 2 + R[2] ** 2)
p = vector(px, py, pz)
pmag = math.sqrt(px ** 2 + py ** 2 + pz ** 2)
phat = p/pmag
dec = degrees(asin(phat.z))
RA = degrees(acos(phat.x / cos(dec)))
print("Declination:",dec)
print("RA:",RA)

while (1==1):
    rate(200)
    time = time + 1
    Mtrue = 2*pi/period*(time) + M
    Etrue = solvekep(Mtrue)

#     # r1ecliptic.x = (cos(wprime)*cos(Oprime) - sin(wprime)*cos(iprime)*sin(Oprime))*(a*cos(Etrue)-a*e) - (cos(wprime)*cos(iprime)*sin(Oprime) + sin(wprime)*cos(Oprime))*(a*sqrt(1-e**2)*sin(Etrue))
#     # r1ecliptic.y = ((cos(wprime)*sin(Oprime) + sin(wprime)*cos(iprime)*cos(Oprime))*(a*cos(Etrue)-a*e)) + (cos(wprime)*cos(iprime)*cos(Oprime) - sin(wprime)*sin(Oprime))*(a*sqrt(1-e**2)*sin(Etrue))
#     # r1ecliptic.z = sin(wprime)*sin(iprime)*(a*cos(Etrue)-a*e) + cos(wprime)*sin(iprime)*(a*sqrt(1-e**2)*sin(Etrue))

    # firstRotation = [[cos(Oprime),-sin(Oprime), 0],[sin(Oprime),cos(Oprime),0],[0,0,1]]
    # firstRotation = array(firstRotation)
    # secondRotation = [[1, 0, 0],[0,cos(iprime),-sin(iprime)],[0,sin(iprime),cos(iprime)]]
    # secondRotation = array(secondRotation)
    # thirdRotation = [[cos(wprime),-sin(wprime), 0],[sin(wprime),cos(wprime),0],[0,0,1]]
    # thirdRotation = array(thirdRotation)

#     # vec = [[r1ecliptic.x],[r1ecliptic.y],[r1ecliptic.z]]
#     # vec = array(vec)

#     # r1ecliptic.x = firstRotation * secondRotation * thirdRotation * r1ecliptic.x
#     # r1ecliptic.y = firstRotation * secondRotation * thirdRotation * r1ecliptic.y
#     # r1ecliptic.z = firstRotation * secondRotation * thirdRotation * r1ecliptic.z
#     # vec = firstRotation @ secondRotation @ thirdRotation @ vec

    # firstRotation = matmul(firstRotation, secondRotation)
    # firstRotation = matmul(firstRotation, thirdRotation)
    r1ecliptic = vector(a*cos(Etrue) - a*e, a*sqrt(1-e**2) * sin(Etrue), 0)
    r1eclipticList = [[r1ecliptic.x],[r1ecliptic.y],[r1ecliptic.z]]
    r1ecliptic = array(r1eclipticList)
    r1ecliptic = matmul(firstRotation, r1ecliptic)
#     # r1ecliptic.x = vec[0,0]
#     # r1ecliptic.y = vec[1,0]
#     # r1ecliptic.z = vec[2,0]

    r1ecliptic = vector(r1ecliptic[0,0], r1ecliptic[1,0], r1ecliptic[2,0])
    asteroid.pos = r1ecliptic*150
    asteroid.trail.append(pos=asteroid.pos)  
    r1eclipticList = [[r1ecliptic.x],[r1ecliptic.y],[r1ecliptic.z]]
    r1ecliptic = array(r1eclipticList)
