#! /bin/sh
# Generate Synthetic Seismograms

dt=0.00002                      # time sampling interval
lt=0.05                      # latest time modeled
fx=0                            # first x value
verbose=1                       # =1 chatty, =0 silent
snfile="snaps.su"                # output file for snapshots
hsz=0.0		                # z-position of horizontal line of geophones
vsx=180		                # x-position of vertical line of geophones
snaptime=0.002,0.004,0.006,0.008,0.01,0.012,0.014,0.016,0.018,0.02,0.022,0.024,0.026,0.028,0.03,	
# times of snapshots
bc=2,1,1,0                      # boundary conditions
qsw=0                           # =1 put in attenuation
asw=0                           # =1 anisotropy
sx=150                            # x-position of sources
sz=61                            # z-position of sources
favg=6.00                      # average frequency
ts=.05                         # source duration
wtype=sp                       # waveform type
stype=p                        # waveform type




# run suea2df to get first set of data
suea2df dt=$dt lt=$lt nz=$nz fx=$fx nx=$nx dx=$dx dz=$dz verbose=1 \
snfile=$snfile rhofile=$rhofile hsz=$hsz vsx=$vsx snaptime=$snaptime \
bc=$bc qsw=$qsw asw=$asw sx=$sx sz=$sz favg=$favg ts=$ts wtype=$wtype stype=$stype\
>out2

cp vsp.su vsp1.su

#run suea2df to get next set, with observation point moved 0.25 meters
vsx=180.25
suea2df dt=$dt lt=$lt nz=$nz fx=$fx nx=$nx dx=$dx dz=$dz verbose=1 \
snfile=$snfile rhofile=$rhofile hsz=$hsz vsx=$vsx snaptime=$snaptime \
bc=$bc qsw=$qsw asw=$asw sx=$sx sz=$sz favg=$favg ts=$ts wtype=$wtype stype=$stype\
>out2
