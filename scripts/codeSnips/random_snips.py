### Switch to directory with Pengcheng's x,y,z data for the fractures OT-P, W-N, W-S, E-N, E-S, and W-L
# os.chdir('/home/spri902/EGS_Collab/4850/fractures/')
# vtkFile = 'fitted_fracs.vtk'
# reader = vtk.vtkUnstructuredGridReader()
# reader.SetFileName(vtkFile)
# reader.Update()
# data = reader.GetOutput()
# points = data.GetPoints()
# npts = points.GetNumberOfPoints()
# x = vtk_to_numpy(points.GetData())

### Plot stereonet for OT-P connector plane 
# fracs = pd.read_csv('/home/spri902/EGS_Collab/4850/fractures/fracs.csv',header=0)
# fig,ax = mpl.subplots(figsize=(9,6))
# ax.grid(color='k',alpha=0.2)

# types = {'One':'red','Four':'blue'}
# for i in range(len(fracs)):
#     strike, dip = fracs['Azimuth(deg)'][i]-90, fracs['Dip(deg)'][i]
#     ax.plane(strike,dip,markersize=4)
#     # ax.pole(strike,dip,markersize=4,alpha=0.5)
#     ax.grid()