from numpy import pi, tan
import csv
import sys
import matplotlib.pyplot as plt
global reader, relP, relW, relP2W

# rel stands for revelant
# P ==> Power(W)
# W ==> Mass(kg)
# P2W ==> Power to Weight(W/kg)
relP = []
relW = []
relP2W = []

with open('powervsweight.csv', 'rb') as f:
    reader = csv.reader(f)
    for row in reader:
    	if row[7] != '' and float(row[7]) < 5000:
    		relP.append(float(row[7]))
    		relW.append(float(row[5]))
    		relP2W.append(float(row[7])/float(row[7]))

# the 5 index are the weights in kg
# the 7 index are the power output in W
        
plt.plot(relW,relP)
plt.ylabel('Power output (W)')
plt.xlabel('Engine weight (kg)')
plt.show()
            
        



        
        

