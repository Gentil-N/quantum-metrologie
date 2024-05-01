from qweather import QWeatherClient
import datetime
import time
import os
import numpy as np
import pyvisa
import random
import pandas as pd


cnx = QWeatherClient("tcp://10.209.64.56:5559")


osc = cnx.SR1LoggingOSC
print(cnx)
#osc = cnx.Sr1RTB2004OSC
print('Got oscilloscope')
#print(osc)
#print(osc.single_measurement)


RAMSEY = False
NOTHING= True
CPULLING= False
rm = pyvisa.ResourceManager()
print(rm.list_resources())
#
#rsosc2 = pyvisa.ResourceManager().open_resource('TCPIP::10.209.66.246::INSTR')
#print (rsosc2.query("*IDN?"))

#No detuning change
if NOTHING:
    for i in range(100):
        print(i)
#        rsosc2.single_measurement(4,save_filepath="Z:/Sr1 2022/10-07-2023/trctest.h5")
#        single_measurement(rsosc2,4,save_filepath="Z:/Sr1 2022/10-07-2023/trctest.h5")
        osc.single_measurement([1],save_filepath="Z:/Sr1 2024/2024-04-30/G1Measurement_{:d}.h5".format(i))
        #osc.single_measurement(1,save_filepath="Z:/Sr1 2022/10-07-2023/beat004.h5".format(i))
        print('Done')



#Ramsey Detuning change
if RAMSEY:
    
    print(rm.list_resources())
    rig = rm.open_resource('TCPIP0::10.209.96.227::INSTR')
    rig.write('*IDN?')
    time.sleep(0.01)
    rig.write('OUTPUT ON')
    time.sleep(0.01)
    rig.write(':SOUR1:FREQ 43.250e6')
    time.sleep(1.01)

    center_freq=43.185
    arr=np.linspace(0,0,1)
    detun=arr.tolist()
    random.shuffle(detun)

    detun=[round(elem, 3) for elem in detun]
    datasets_per_detuning = 5
    print(detun)
    #print(rig.read())
    w=0

    rig.write(':SOUR1:VOLT:UNIT DBM')
    rig.write(':SOUR1:VOLT 6.8')

    #fit RF Power to have even optical powers for different detunings
    #df = pd.read_csv('Powercalib0702.csv', usecols = ['Detuning','RF Power'])
    #z = np.polyfit(df['Detuning'],df['RF Power'],7)
    #poly = np.poly1d(z)
    
    ##
    for i in range(datasets_per_detuning):
        
        detun=arr.tolist()
        random.shuffle(detun)
        detun=[round(elem, 3) for elem in detun]
        for curr_detun in detun:
            new_freq = center_freq - curr_detun
            #new_amp = poly(curr_detun)
            rig.write(':SOUR1:FREQ '+str(new_freq)+'e6')
            #if new_amp>8.2:
            #    new_amp=8.0
            #rig.write(':SOUR1:VOLT '+str(new_amp))
            print(new_freq)
            #print(new_amp)
            print("changed center freq")
            print(w)
            
            osc.single_measurement([1],save_filepath="Z:/Sr1 2024/2024-04-30/Beattest2_detuning_"+str(curr_detun)+'_{:d}.h5'.format(i))
            #osc.single_measurement(3,save_filepath="Z:/Sr1 2022/03-08-2022/PulseLength1_detuning"+str(curr_detun)+'_{:d}.h5'.format(i))
            w=w+1
    rig.write(':SOUR1:FREQ '+str(center_freq)+'e6')
    rig.write(':SOUR1:VOLT 6.8')


if CPULLING:
    print(rm.list_resources())
    hp = rm.open_resource('ASRL3::INSTR')
    center_freq=345.083
    detun_start=-2.5
    detun_end=2.5
    detun_steps=50 #must be even
    detun_stepsize=abs(detun_end-detun_start)/detun_steps
    arr1=np.linspace(detun_start,detun_end-detun_stepsize,int(detun_steps/2))
    arr2=np.linspace(detun_end,detun_start+detun_stepsize,int(detun_steps/2))
    arr=np.hstack((arr1,arr2))
    detun=arr.tolist()
    #random.shuffle(detun)
    detun=[round(elem, 2) for elem in detun]
    datasets_per_detuning = 3
    print(detun)
    j=0
    for j in range(len(detun)):
        curr_detun=detun[j]

        new_freq = center_freq - curr_detun
        hp.write('CFRQ '+str(new_freq)+' MHz')
        print(new_freq)
        
        
        print("change center freq")
        time.sleep(.5)
        for i in range(datasets_per_detuning):
           osc.single_measurement([1],save_filepath="Z:/Sr1 2024/2024-04-30/NMScavDetuning_NoMolasses_"+str(2*curr_detun)+'_{:d}.h5'.format(i))
        
    hp.write('CFRQ '+str(center_freq)+' MHz')
    
   