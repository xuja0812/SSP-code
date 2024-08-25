import numpy as np
from math import *
from XuOD import *
from astropy.io import fits
from Orbital_elements_Xu import *
import matplotlib.pyplot as plt

def find_error(file_name):
    tbl = fits.open(file_name)[1].data
    rms_ra = 3600*(np.mean((tbl.field_ra - tbl.index_ra)**2.))**0.5
    rms_dec = 3600*(np.mean((tbl.field_dec - tbl.index_dec)**2.))**0.5

    return rms_ra, rms_dec

observation3 = "D:\\corr_0703_001light120_2002KL6.fits"
observation2 = "D:\\corr (2).fits"
observation1 = "D:\\corr (3).fits"

ra_error3, dec_error3 = find_error(observation3)
ra_error2, dec_error2 = find_error(observation2)
ra_error1, dec_error1 = find_error(observation1)

ra_error3 /= 360000
dec_error3 /= 360000
ra_error2 /= 360000
dec_error2 /= 360000
ra_error1 /= 360000
dec_error1 /= 360000
# ra_error1, dec_error1 = find_error(observation1)
print(ra_error3, dec_error3)
print(ra_error2, dec_error2)
print(ra_error1, dec_error1)

# RA1 = find_RA("15:49:41.38")
# dec1 = find_dec("00:32:12.50")
# RA2 = find_RA("15:52:32.32")
# dec2 = find_dec("03:06:24.80")
# RA3 = find_RA("16:04:31.18")
# dec3 = find_dec("09:32:21.60")

a_list = []
e_list = []
i_list = []
w_list = []
om_list = []
ma_list = []

data = np.loadtxt("D:\\OD_Xu\\2002KL6_input_Xu.txt", dtype=str)
years = data[:,0].astype(float)
months = data[:,1].astype(float)
days = data[:,2].astype(float)
times = data[:,3]
RAs = data[:,4]
decs = data[:,5]
Rxs = data[:,6].astype(float)
Rys = data[:,7].astype(float)
Rzs = data[:,8].astype(float)

data1 = process_inputs(years[0], months[0], days[0], times[0], RAs[0], decs[0], Rxs[0], Rys[0], Rzs[0])
data2 = process_inputs(years[1], months[1], days[1], times[1], RAs[1], decs[1], Rxs[1], Rys[1], Rzs[1])
data3 = process_inputs(years[2], months[2], days[2], times[2], RAs[2], decs[2], Rxs[2], Rys[2], Rzs[2])

RA1 = data1[1]
RA2 = data2[1]
RA3 = data3[1]

# print("RA1",RA1)
# print("RA2",RA2)
# print("RA3",RA3)

dec1 = data1[2]
dec2 = data2[2]
dec3 = data3[2]

# print("dec1",dec1)
# print("dec2",dec2)
# print("dec3",dec3)

for i in range(10 ** 2):
    RA1 = np.random.normal(RA1,ra_error1)
    dec1 = np.random.normal(dec1,dec_error1)
    RA2 = np.random.normal(RA2,ra_error2)
    dec2 = np.random.normal(dec2,dec_error2)
    RA3 = np.random.normal(RA3,ra_error3)
    dec3 = np.random.normal(dec3,dec_error3)

    # print("RA1",RA1)
    # print("RA2",RA2)
    # print("RA3",RA3)

    # print("dec1",dec1)
    # print("dec2",dec2)
    # print("dec3",dec3)

    data = np.loadtxt("D:\\OD_Xu\\2002KL6_input_Xu.txt", dtype=str)
    years = data[:,0].astype(float)
    months = data[:,1].astype(float)
    days = data[:,2].astype(float)
    times = data[:,3]
    RAs = data[:,4]
    decs = data[:,5]
    Rxs = data[:,6].astype(float)
    Rys = data[:,7].astype(float)
    Rzs = data[:,8].astype(float)

    data1 = process_inputs(years[0], months[0], days[0], times[0], RAs[0], decs[0], Rxs[0], Rys[0], Rzs[0])
    data2 = process_inputs(years[1], months[1], days[1], times[1], RAs[1], decs[1], Rxs[1], Rys[1], Rzs[1])
    data3 = process_inputs(years[2], months[2], days[2], times[2], RAs[2], decs[2], Rxs[2], Rys[2], Rzs[2])

    t1 = data1[0]
    t2 = data2[0]
    t3 = data3[0]

    taus = find_taus(t1, t2, t3)

    # RA1 = data1[1]
    # RA2 = data2[1]
    # RA3 = data3[1]

    # dec1 = data1[2]
    # dec2 = data2[2]
    # dec3 = data3[2]

    R1 = data1[3]
    R2 = data2[3]
    R3 = data3[3]

    rhohat1 = rho_hat(RA1, dec1)
    rhohat2 = rho_hat(RA2, dec2)
    rhohat3 = rho_hat(RA3, dec3)

    Ds = find_Ds(rhohat1, R1, R2, R3, rhohat3, rhohat2)
    r2, rho2 = SEL(taus, R2, rhohat2, Ds)

    # UPDATE BELOW TO CHOOSE THE RIGHT R2 AND RHO2 VALUES

    # if(len(r2) > 1):
    #     print("ROOTS:",r2)
    #     print("RHOS:",rho2)
    #     index = int(input("Which root to use (index)? "))
    #     r2 = r2[index]
    #     rho2 = rho2[index]
    # else:
    r2 = r2[-1]
    rho2 = rho2[-1]
    r2dot, r2vec, rho1, rho3 = initial_r2dot(taus[0], taus[1], r2, R1, R2, R3, rhohat1, rhohat2, rhohat3)
    rho2vec = rho2 * rhohat2
    r2vec = vp.vector(r2vec[0], r2vec[1], r2vec[2]) 

    A = a(r2vec, r2dot)
    n = sqrt(1/(A**3))

    r2dot, r2, rho2mag2 = iterate(r2dot, r2vec, taus[0], taus[1], R2, rho2vec, R3, R1, rhohat1, rhohat2, rhohat3, t1, t2, t3, rho1, rho2vec.mag,rho3)

    # converts vectors to ecliptic 
    r2dot = do_rotation(r2dot)
    r2 = do_rotation(r2)

    h = angMomentum(r2, r2dot)
    A = a(r2, r2dot)
    E = e(A, h)
    TA = true_anomaly(E, A, r2, r2dot)
    i = degrees(incline(h))
    omega = degrees(longitude_of_ascending_node(h, i))
    w = degrees(argument_of_perihelion(omega, i, r2, TA, A, E, h, r2dot))
    EA = eccentric_anomaly(E, r2, A, TA)
    MA = degrees(mean_anomaly(E, EA))
    t0 = find_JD(2021, 7, 24, "07:00:00")
    MA2 = M(t2, t0, MA, A)

    a_list.append(A)
    e_list.append(E)
    i_list.append(i)
    w_list.append(w)
    om_list.append(omega)
    ma_list.append(MA)

a_avg = np.mean(np.array(a_list))
a_std = np.std(np.array(a_list))
a_SDOM = a_std/len(a_list)
print("A AVG:",a_avg)
print("A STD:",a_std)
e_avg = np.mean(np.array(e_list))
e_std = np.std(np.array(e_list))
print("E AVG:",e_avg)
i_avg = np.mean(np.array(i_list))
i_std = np.std(np.array(i_list))
print("I AVG:",i_avg)
w_avg = np.mean(np.array(w_list))
w_std = np.std(np.array(w_list))
print("w AVG:",w_avg)
om_avg = np.mean(np.array(om_list))
om_std = np.std(np.array(om_list))
print("om AVG:",om_avg)
ma_avg = np.mean(np.array(ma_list))
ma_std = np.std(np.array(ma_list))
print("MA AVG:",ma_avg)

ExpA = 2.306817004734406 # 2.720870
ExpE = .5483638958524116  # 0.536660
ExpIncline = 3.24614366112847 # 3.896806
ExpAoP = 98.04416849015669 # 108.222281
ExpLoAN = 213.2862496588689 # 238.484915
ExpMA = 3.508990579885306E+02 # 339.757416
ExpH = vp.vector(-0.08063796, 0.0494636, 1.38863187)
ExpMA2 = 343.994717
print()
print("INCLINE (i)")
print("Expected incline:",ExpIncline)
print("Calculated incline:",i_avg)
print("Incline percent error:",percent_error(i_avg,ExpIncline))
print()
print("SEMI-MAJOR AXIS (a)")
print("Expected a:",ExpA)
print("Calculated a:",a_avg)
print("a percent error:",percent_error(a_avg,ExpA))
print()
print("ECCENTRICITY (e)")
print("Expected e:",ExpE)
print("Calculated e:",e_avg)
print("e percent error:",percent_error(e_avg,ExpE))
print()
print("LONGITUDE OF ASCENDING NODE (OM)")
print("Expected OM:",ExpLoAN)
print("Calculated OM:",om_avg)
print("OM percent error:",percent_error(om_avg,ExpLoAN))
print()
print("ARGUMENT OF PERIHELION (w)")
print("Expected w:",ExpAoP)
print("Calculated w:",w_avg)
print("w percent error:",percent_error(w_avg,ExpAoP))
print()
print("MEAN ANOMALY (MA)")
print("Expected MA:",ExpMA)
print("Calculated MA:",ma_avg)
print("MA percent error:",percent_error(ma_avg,ExpMA))
print()


def gaussian(x, B, mean, std):
    return 1./(np.sqrt(np.pi * 2) * std)*np.exp(-(x-mean)**2/(2*std**2))

fig, ax = plt.subplots(1,1)

# hist_a = np.random.normal(a_avg, a_std, 500)

# ax.hist(hist_a, alpha=0.75, bins = 25, zorder = 2, histtype="step", linewidth=3, density=True)
# # ax.axvline(np.mean(hist_a), color="black", linestyle="--", label="Mean = " + "%.2F" % np.mean(hist_a))
# # ax.axvline(a_avg, color="black", linestyle="--", label="Mean = " + "%.2F" % a_avg)
# span = np.linspace(a_avg - a_std, a_avg + a_std,50)
# # xrange = np.linspace(a_avg - 3*a_std, a_avg + 3*a_std, 50)
# xrange = np.linspace(a_avg - 5*a_std, a_avg + 5*a_std, 50)
# ax.plot(xrange, gaussian(xrange, None, a_avg, a_std), color="black")
# plt.fill_between(span, np.zeros(50), gaussian(span, None, a_avg, a_std), alpha = 0.50, color="red", zorder=3, label=r"$\sigma$= %.4F" % a_std)
# plt.xlabel("Semi-major axis (a) [AU]")
# plt.ylabel("Frequency")

hist_e = np.random.normal(e_avg, e_std, 500)

ax.hist(hist_e, alpha=0.75, bins = 25, zorder = 2, histtype="step", linewidth=3, density=True)
ax.axvline(np.mean(hist_e), color="black", linestyle="--", label="Mean = " + "%.2F" % np.mean(hist_e))
span = np.linspace(e_avg - e_std, e_avg + e_std,50)
xrange = np.linspace(e_avg - 3*e_std, e_avg + 3*e_std, 50)
ax.plot(xrange, gaussian(xrange, None, e_avg, e_std), color="black")
plt.fill_between(span, np.zeros(50), gaussian(span, None, e_avg, e_std), alpha = 0.50, color="red", zorder=3, label=r"$\sigma$= %.4F" % e_std)
plt.xlabel("Eccentricity (e)")
plt.ylabel("Frequency")

# hist_i = np.random.normal(i_avg, i_std, 500)

# ax.hist(hist_i, alpha=0.75, bins = 25, zorder = 2, histtype="step", linewidth=3, density=True)
# ax.axvline(np.mean(hist_i), color="black", linestyle="--", label="Mean = " + "%.2F" % np.mean(hist_i))
# span = np.linspace(i_avg - i_std, i_avg + i_std,50)
# xrange = np.linspace(i_avg - 3*i_std, i_avg + 3*i_std, 50)
# ax.plot(xrange, gaussian(xrange, None, i_avg, i_std), color="black")
# plt.fill_between(span, np.zeros(50), gaussian(span, None, i_avg, i_std), alpha = 0.50, color="red", zorder=3, label=r"$\sigma$= %.4F" % i_std)
# plt.xlabel("Incline (i) [degrees]")
# plt.ylabel("Frequency")

# hist_w = np.random.normal(w_avg, w_std, 500)

# ax.hist(hist_w, alpha=0.75, bins = 25, zorder = 2, histtype="step", linewidth=3, density=True)
# ax.axvline(np.mean(hist_w), color="black", linestyle="--", label="Mean = " + "%.2F" % np.mean(hist_w))
# span = np.linspace(w_avg - w_std, w_avg + w_std,50)
# xrange = np.linspace(w_avg - 3*w_std, w_avg + 3*w_std, 50)
# ax.plot(xrange, gaussian(xrange, None, w_avg, w_std), color="black")
# plt.fill_between(span, np.zeros(50), gaussian(span, None, w_avg, w_std), alpha = 0.50, color="red", zorder=3, label=r"$\sigma$= %.4F" % w_std)
# plt.xlabel("Argument of Perihelion (w) [degrees]")
# plt.ylabel("Frequency")

# hist_om = np.random.normal(om_avg, om_std, 500)

# ax.hist(hist_om, alpha=0.75, bins = 25, zorder = 2, histtype="step", linewidth=3, density=True)
# ax.axvline(np.mean(hist_om), color="black", linestyle="--", label="Mean = " + "%.2F" % np.mean(hist_om))
# span = np.linspace(om_avg - om_std, om_avg + om_std,50)
# xrange = np.linspace(om_avg - 3*om_std, om_avg + 3*om_std, 50)
# ax.plot(xrange, gaussian(xrange, None, om_avg, om_std), color="black")
# plt.fill_between(span, np.zeros(50), gaussian(span, None, om_avg, om_std), alpha = 0.50, color="red", zorder=3, label=r"$\sigma$= %.4F" % om_std)
# plt.xlabel("Longitude of Ascending Node (OM) [degrees]")
# plt.ylabel("Frequency")

# hist_ma = np.random.normal(ma_avg, ma_std, 500)

# ax.hist(hist_ma, alpha=0.75, bins = 25, zorder = 2, histtype="step", linewidth=3, density=True)
# ax.axvline(np.mean(hist_ma), color="black", linestyle="--", label="Mean = " + "%.2F" % np.mean(hist_ma))
# span = np.linspace(ma_avg - ma_std, ma_avg + ma_std,50)
# xrange = np.linspace(ma_avg - 3*ma_std, ma_avg + 3*ma_std, 50)
# ax.plot(xrange, gaussian(xrange, None, ma_avg, ma_std), color="black")
# plt.fill_between(span, np.zeros(50), gaussian(span, None, ma_avg, ma_std), alpha = 0.50, color="red", zorder=3, label=r"$\sigma$= %.4F" % ma_std)
# plt.xlabel("Mean Anomaly (MA) [degrees]")
# plt.ylabel("Frequency")

plt.legend()
plt.show()
