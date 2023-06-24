#!/usr/bin/env python
# coding: utf-8

# In[1]:


##Sintela_ReadRaw

#Opens a raw time series recording (V3 or V4)
#Reads and Displays header
#Reads the data, extracts and plot a selection of data

#setup
import struct
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
import sys
from decimal import *

print("Python version %s" % sys.version)

#setup raw file to read
rawName = '/home/spri902/bc_cribs/'

#setup
datetimeFormat = "%Y-%m-%d %H:%M:%S"

maxTriggers = 8
triggerOffsets=[0]*maxTriggers


# In[2]:
## Read and check headers
fid = open(rawName, 'rb')

#get file size (used to calculate the number of samples in file)
fileContent = fid.read()
filesizeBytes = len(fileContent)
fid.seek(0)

dataType = -1

#read header
syncWord = struct.unpack("i", fid.read(4))[0]        #0
headerSize = struct.unpack("i", fid.read(4))[0]      #4
headerVersion = struct.unpack("i", fid.read(4))[0]   #8
    
if syncWord == 287454020 and headerSize == 64:
    print('Sintela Data Type 1  ** NO LONGER SUPPORTED **\n')
    exit()
elif syncWord == 287454020 and headerSize == 999:
    print("Sintela Data Type 2  ** NO LONGER SUPPORTED **\n")
    exit()
elif syncWord == 287454020 and headerSize == 96:
    print("Sintela Data Type 3 - Raw 32 Bit\n")
    dataType = 3;
elif syncWord == 287454020 and headerSize == 112:
    print("Sintela Data Type 4 - 16 Bit Differential\n")
    dataType = 4;
else:
    print("Unsupported File Format : Header = %X  DataType = %d\n" % (syncWord, headerSize))

    
if dataType == 3:
    #32 bit float
    sampleCount = struct.unpack("Q", fid.read(8))[0]     #12
    sampleTiming = struct.unpack("Q", fid.read(8))[0]    #20
    numChannels = struct.unpack("i", fid.read(4))[0]     #28
    numSamples = struct.unpack("i", fid.read(4))[0]      #32
    sampleRate = struct.unpack("f", fid.read(4))[0]      #36
    channelSpacing = struct.unpack("f", fid.read(4))[0]  #40
    gaugeLength = struct.unpack("f", fid.read(4))[0]     #44
    trigger = struct.unpack("i", fid.read(4))[0]         #48
    startChannel = struct.unpack("i", fid.read(4))[0]    #52
    channelStep = struct.unpack("i", fid.read(4))[0]     #56

    gpsPpsTime = struct.unpack("i", fid.read(4))[0]      #60
    gpsPpsOffset = struct.unpack("i", fid.read(4))[0]    #64
    gpsStatus = struct.unpack("i", fid.read(4))[0]       #68    

    triggerFlags = int.from_bytes(fid.read(1), byteorder='little')  #72
    flippedFlag = int.from_bytes(fid.read(1), byteorder='little')   #73
    spareFlags3 = int.from_bytes(fid.read(1), byteorder='little')   #74
    spareFlags4 = int.from_bytes(fid.read(1), byteorder='little')   #75
    spareFlags5 = int.from_bytes(fid.read(1), byteorder='little')   #76
    spareFlags6 = int.from_bytes(fid.read(1), byteorder='little')   #77
    spareFlags7 = int.from_bytes(fid.read(1), byteorder='little')   #78
    spareFlags8 = int.from_bytes(fid.read(1), byteorder='little')   #79

    for i in range(0, maxTriggers):                      #80..96
        triggerOffsets[i] = int.from_bytes(fid.read(2), byteorder='little')                

elif dataType == 4:
    #16 bit int differential
    codecType = struct.unpack("i", fid.read(4))[0]       #12
    sampleCount = struct.unpack("Q", fid.read(8))[0]     #16
    sampleTiming = struct.unpack("Q", fid.read(8))[0]    #24
    numChannels = struct.unpack("i", fid.read(4))[0]     #32
    numSamples = struct.unpack("i", fid.read(4))[0]      #36
    sampleRate = struct.unpack("f", fid.read(4))[0]      #40
    channelSpacing = struct.unpack("f", fid.read(4))[0]  #44
    gaugeLength = struct.unpack("f", fid.read(4))[0]     #48
    trigger = struct.unpack("i", fid.read(4))[0]         #52
    startChannel = struct.unpack("i", fid.read(4))[0]    #56
    channelStep = struct.unpack("i", fid.read(4))[0]     #60

    scaleFactor = struct.unpack("f", fid.read(4))[0]     #64
    align = struct.unpack("i", fid.read(4))[0]           #68
    
    gpsPpsTime = struct.unpack("i", fid.read(4))[0]      #72
    gpsPpsOffset = struct.unpack("i", fid.read(4))[0]    #76
    gpsStatus = struct.unpack("i", fid.read(4))[0]       #80    

    triggerFlags = int.from_bytes(fid.read(1), byteorder='little')  #84
    flippedFlag = int.from_bytes(fid.read(1), byteorder='little')   #85
    spareFlags3 = int.from_bytes(fid.read(1), byteorder='little')   #86
    spareFlags4 = int.from_bytes(fid.read(1), byteorder='little')   #87
    spareFlags5 = int.from_bytes(fid.read(1), byteorder='little')   #88
    spareFlags6 = int.from_bytes(fid.read(1), byteorder='little')   #89
    spareFlags7 = int.from_bytes(fid.read(1), byteorder='little')   #90
    spareFlags8 = int.from_bytes(fid.read(1), byteorder='little')   #91

    for i in range(0, maxTriggers):                      #92..104
        triggerOffsets[i] = int.from_bytes(fid.read(2), byteorder='little')     
        
    flippedFlag = int.from_bytes(fid.read(1), byteorder='little')   #108
    spareFlags22 = int.from_bytes(fid.read(1), byteorder='little')  #109
    spareFlags23 = int.from_bytes(fid.read(1), byteorder='little')  #110
    spareFlags24 = int.from_bytes(fid.read(1), byteorder='little')  #111
     
#rewind file
fid.seek(0)

# In[3]:
# Decode and Display 'first' header

#print('Sync Word = %s' % hex(syncWord))
print("Packet Chans x Samples : %d x %d  StartCh:%d  Fs:%8.3f  ChannelSpacing:%4.1f GaugeLength:%4.1f" %(numChannels, numSamples, startChannel, sampleRate, channelSpacing, gaugeLength)  )

startTime = datetime.fromtimestamp(sampleTiming/1e9).strftime(datetimeFormat)
startNsIntoSec = Decimal(sampleTiming) % Decimal(1e9)
print("StartTime: %s.%09d" % (startTime, startNsIntoSec))

prevFlippedFlag = flippedFlag; # initial value

# decode status flags
gpsAccuracy = gpsStatus & 0x00FFFFFF
gpsFlags = gpsStatus >> 24

accuracy =      gpsAccuracy;
accOverFlag =   (gpsFlags & 0x00000001) >> 0
valFlag =       (gpsFlags & 0x00000002) >> 1
syncFlag =      (gpsFlags & 0x00000004) >> 2
timestampFlag = (gpsFlags & 0x00000008) >> 3

if valFlag == 1:
    print("GPS:       Enabled (Accuracy=%dns  AccOver=%d  Valid=%d  Sync=%d  Timestamp=%d)" % (accuracy, accOverFlag, valFlag, syncFlag, timestampFlag) )
else:
    print("GPS:       Disabled") 
        
#convert trigger info
if triggerFlags > 0: 

    bit = 1
    for i in range(0, maxTriggers): 

        if (triggerFlags & bit) == bit:

            nsIntoPacket = triggerOffsets[i]* 1e9 / sampleRate
            triggerTime = datetime.fromtimestamp((sampleTiming + nsIntoPacket)/1e9).strftime(datetimeFormat)
            triggerNsIntoSec = (Decimal(sampleTiming) + Decimal(nsIntoPacket)) % Decimal(1e9)

            print("%s.%09d Trigger in #%d : Offset=%d  SampleIndex=%d" % (triggerTime, triggerNsIntoSec, i, triggerOffsets[i], sampleIndex + triggerOffsets[i]))

        bit = bit << 1

if dataType == 3:
    packetSizeBytes = (headerSize + (numChannels*numSamples*4));
    totalPackets = filesizeBytes / packetSizeBytes;
    totalSamples = totalPackets*numSamples;
elif dataType == 4:
    #TS02
    #packetSizeBytes = headerSize + (numChannels*numSamples*2); # int16
    #TS04
    packetSizeBytes = headerSize + (numChannels*4) + (numChannels*(numSamples-1)*2); # float32 + (n-1)int16
    
    totalPackets = filesizeBytes / packetSizeBytes;
    totalSamples = totalPackets*numSamples;

print("Total Samples Per Channel = %d   Total Packets = %d   NSamples Per Packet = %d\n" % (totalSamples, totalPackets, packetSizeBytes))


# In[4]:

##  Read and Extract Data Samples (Reading but ignoring headers blocks)

extractChannel = 0; # first channel in file
nChansToExtract = numChannels

if dataType == 3:
    #32 bit float
    chData = np.zeros((int(totalSamples), nChansToExtract))
    sampleIndex = 0
    
    for packet in range(0, int(totalPackets)):
        #read samples
        header = np.fromfile(fid, dtype=np.int8, count=headerSize)
        samplesIn = np.fromfile(fid, dtype=np.float32, count=numChannels*numSamples)
        
        samples = np.reshape(samplesIn, [numSamples, numChannels])
    
        for s in range(0, numSamples):
            for ch in range(0, nChansToExtract):
                chData[sampleIndex, ch] = samples[s, extractChannel+ch]
                    
            sampleIndex = sampleIndex+1
            
elif dataType == 4:
    #16 bit int differential
    integrationStore = np.zeros((1, nChansToExtract))
    chData = np.zeros((int(totalSamples), nChansToExtract))
    sampleIndex = 0
    
    for packet in range(0, int(totalPackets)):
        #read samples
        header = np.fromfile(fid, dtype=np.int8, count=headerSize)
        
        ##TS02
        #samplesIn = np.fromfile(fid, dtype=np.int16, count=numChannels*numSamples)
        #
        #samples = np.reshape(samplesIn, [numSamples, numChannels])
        # 
        #for s in range(0, numSamples):
        #    for ch in range(0, nChansToExtract):
        #        integrationStore[0,ch] = integrationStore[0, ch] + samples[s, extractChannel+ch]
        #        chData[sampleIndex, ch] = integrationStore[0,ch] / scaleFactor
        #            
        #    sampleIndex = sampleIndex+1
            
        #TS04   
        for ch in range(0, nChansToExtract):
            
            si = sampleIndex;
            
            integrationStore[0, ch] = np.fromfile(fid, dtype=np.float32, count=1)
            chData[si, ch] = integrationStore[0, ch]
            si = si + 1;
            
            data = np.fromfile(fid, dtype=np.int16, count=numSamples-1)

            for s in range(0, numSamples-1):
                integrationStore[0,ch] = integrationStore[0, ch] + (data[s] / scaleFactor)
                chData[si, ch] = integrationStore[0,ch]
                si = si + 1
        
        sampleIndex = si
            
fid.close()


# In[5]:


## process and plot data

#print(chData.shape)

plt.plot(chData)

fig, (ax) = plt.subplots(1)
ax.imshow(chData, aspect='auto' )
plt.show()