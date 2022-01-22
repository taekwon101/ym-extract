from math import floor
from matplotlib import pyplot as plt
import pandas as pd
import numpy as np
from scipy import stats

# TO-DO: automate finding csvs

# create temporary dataframe
df = pd.read_csv('DAQ- Crosshead, … - (Timed).txt')
# extract filename and convert to case-normalized string
filenamex = df.loc[1].to_string()
filenamex = filenamex.lower()  
# find string "csa" and read the proceeding 5 digits into crossSection as float (be sure the format is csaXX.XX in mm^2)
index_csa = (filenamex.find('csa'))
tempstring = filenamex[index_csa+3:index_csa+8]
crossSection = float(tempstring) #weird error when cs is lower than like 8?? but I guess that shouldn't happen
print("Cross Section Area: ", crossSection)

# create final dataframe with headers on line 4 (skip second header line to access load and strain arrays)
df = pd.read_csv('DAQ- Crosshead, … - (Timed).txt', header=[4], delimiter='\t')
x_strain = df['Strain 1 '][1:].astype(float)
y_load = df['Load '][1:].astype(float) #/crossSection (figure out math to just write F/A into stress array)
y_stress = y_load/crossSection

# create convolution kernel for calculating the smoothed second order derivative
smooth_width = 100

x1 = np.linspace(-3,3,smooth_width)
norm = np.sum(np.exp(-x1**2)) * (x1[1]-x1[0]) # ad hoc normalization
y1 = (4*x1**2 - 2) * np.exp(-x1**2) / smooth_width *2


# convolve stress data to form smoothed second order derivative 
y_conv = np.convolve(y_stress, y1, mode="same")

# find linear region from convolution
# bluetearstone_saveme
startflag = 0
endflag = 0
lincount = 0
for x in y_conv:
    lincount += 1
    #find start of linear region (must be true for 5 consecutive points)
    if (x < 0.01 and x > -0.01) and (startflag < 5):
        startflag += 1
        if startflag == 1:
            y_lin_start = lincount+5
    elif startflag == 5:
        pass
    elif (x > 0.01) or (x < -0.01):
        startflag = 0
            
    #find end of linear region (must be true for 5 consecutive points)
    if (x > 0.04 or x < -0.04) and (endflag < 5):
        endflag += 1
        if endflag == 1:
            y_lin_end = lincount-35
    elif endflag == 5:
        pass
    elif x > 0.01 or x < -0.01:
         endflag = 0
       
# make arrays for linear trendline
x_line = x_strain[y_lin_start:y_lin_end]
y_line = y_stress[y_lin_start:y_lin_end]
res = stats.linregress(x_line, y_line)
z = np.polyfit(x_line, y_line, 1)
p = np.poly1d(z)
print("Young's Modulus = ", "{:.4f}".format(res.slope))
print("R^2 = ", "{:.4f}".format(res.rvalue**2))

# plot data
# plt.plot(x_strain, y_stress,".k", markersize = 7)
plt.plot(x_strain, y_stress,"o", markersize = 2, label = "noisy data")
plt.plot(x_strain,y_conv, label = "second deriv")
# plt.plot(x_line, p(x_line), "0.3", label = "linear data")
plt.plot(x_strain, res.intercept + res.slope*x_strain, "0.3", label='linear fit')
# plt.hlines([0],-10, 20)
plt.axvspan(0,x_strain[y_lin_start], color="y", alpha=0.2)
plt.axvspan(x_strain[y_lin_start],x_strain[y_lin_end], color="b", alpha=0.2)
plt.axvspan(x_strain[y_lin_end],x_strain[lincount], color="y", alpha=0.2)
# plt.vlines([x_strain[y_lin_start], x_strain[y_lin_end], x_strain[lincount]],-0.01, 0.04)
plt.xlim(0, x_strain[lincount])
plt.ylim(-0.4, 1.5)
plt.legend(loc='upper left')
midgraph = floor((y_lin_end-y_lin_start)/2)
annotext = ("YM = ", "{:.4f}".format(res.slope), " MPa") #need to convert to str from tuple first
plt.annotate("YM", xy=(x_strain[midgraph-.1], y_stress[midgraph]+0.3))
#plt.show()