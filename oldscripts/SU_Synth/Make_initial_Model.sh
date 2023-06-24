#! /bin/sh
# make initial model for SU modeling

modelfile=model.unif            # input model file for unif2aniso
ninf=7                          # number of interfaces (surface counts)
x0=0,0,0,0,0,0,0,                   # x-position(s) for  vp00,vs00,rho00, etc.
z0=0,10,54,57,97,100,         	# z-position(s) for  vp00,vs00,rho00, etc.
nz=1200                         # size of z (depth) dimension  of model
nx=1200                         # size of x (horizontal) dimension of model
dz=0.25                          # increment in z direction
dx=0.25                          # increment in x direction
vp00=1800,4900,4900,4900,4900,4900,       # P-wavespeed(s) at (z0,x0)
vs00=1250,3188,3188,3188,3188,3188,    # S-wavespeed(s) at (z0,x0)
q00=30,500,30,500,500,500,500,
rho00=1800,2400,2000,2400,2400,2400       # density(s) at (z0,x0)
rhofile="rho_file"		# density file, if not specified rho=1000
method=linear			# linear, akima, spline, mono interpolation
				#   for boundaries in unif2aniso models


# build stiffness and density files
#unif2aniso < $modelfile ninf=$ninf x0=$x0 z0=$z0 nz=$nz nx=$nx \
#dx=$dx dz=$dz vp00=$vp00 vs00=$vs00 rho00=$rho00 q00=$q00 method=$method

# the files c11_file, c13_file c15_file c33_file c35_file c55_file rho_file
# are generated by unif2aniso
# transpose stiffness and density
#xbox=0
#ybox=100
for i in c11 c13 c15 c33 c35 c55 rho q
do


#	echo $xbox $ybox
#	ximage <  ${i}_file n1=$nz n2=$nx perc=99 xbox=$xbox ybox=$ybox \
#	 wbox=$nx hbox=$nz  legend=1 title=" ${i} parameter file  " &
#	sleep 2

#	xbox=`expr $xbox + 110 `
#	ybox=`expr $ybox + 5 `
	
        mv ${i}_file tmp.file
        transp n1=$nz < tmp.file > ${i}_file
done


rm tmp.file