'''
    Notes: Rotor 37 is already an airfoil designed using a different design philosopy than this code. To use rotor37 we need to import the points into airfoil 3D or airfoil 2D directly. 

    Why use airofil2D or 3D? 
    1. This enables the generation of intermediate profile.
    2. If you wanted to change the hub and shroud curves, this would fit the geometry in between the curves. 
    3. It also allows you to position the blade anywhere along the hub, so if you wanted to increase or decrease the gap between stator and rotor

    Limitations:
    1. When converting to Airfoil3D, you don't have access to stacking. Stacking uses the camberline 
'''
from microbladegen.aero import Airfoil2D,Airfoil3D,Passage2D
from microbladegen.helper import bezier, bezier3, csapi
import numpy as np
import copy
import matplotlib.pyplot as plt 
#%% R37 Blade by reading from hub and shroud profiles 
# hub_profile = np.loadtxt("hub.csv", delimiter=",",skiprows=1)
# tip_profile = np.loadtxt("tip.csv", delimiter=",",skiprows=1)

# # From the excel file 7in seems to be the radius of the bottom profile (could be mistaken)
# x = hub_profile[:,0]
# ss_y = hub_profile[:,1]
# ps_y = hub_profile[:,2]
# new_hub = np.vstack([np.vstack([x[:-1],ps_y[:-1],7+0*ss_y[:-1]]).transpose(), np.flipud(np.vstack([x,ss_y,7+0*ss_y]).transpose())])

# # plt.figure(1,clear=True)
# # plt.plot(x,ss_y,x,ps_y)
# # plt.axis('equal')
# # plt.show()

# # From the excel file 7.2in seems to be the radius of the top profile (could be mistaken)
# x = tip_profile[:,0]
# ss_y = tip_profile[:,1]
# ps_y = tip_profile[:,2]
# new_tip = np.vstack([np.vstack([x[:-1],ps_y[:-1],7.2+0*ss_y[:-1]]).transpose(), np.flipud(np.vstack([x,ss_y,7.2+0*ss_y]).transpose())])

# os.makedirs('temp',exist_ok=True)
# np.savetxt('temp/rotor37_0.txt',new_hub, delimiter=' ')
# np.savetxt('temp/rotor37_1.txt',new_tip, delimiter=' ')
# rotor37 = Airfoil3D.import_geometry(folder='temp',npoints=100,nspan=50,axial_chord=x.max(),ss_ps_split=len(ps_y)-1)
# rotor37.plot3D()

#%% R37 Blade by reading the 5 profiles from CSV 
# In this example the profiles are well defined - they have a radial component or z value.
files = ['R37_profile01.csv','R37_profile02.csv','R37_profile03.csv','R37_profile04.csv','R37_profile05.csv','R37_profile06.csv']
xss = np.zeros(shape=(len(files),151)) # 5 profiles, each one is 151 points
yss = np.zeros(shape=(len(files),151))
zss = np.zeros(shape=(len(files),151))

xps = np.zeros(shape=(len(files),151))
yps = np.zeros(shape=(len(files),151))
zps = np.zeros(shape=(len(files),151))
data = list()
cx = np.zeros(shape=(len(files),1)); cy = np.zeros(shape=(len(files),1)); cz = np.zeros(shape=(len(files),1))
airfoils = list()
for i,f in enumerate(files):
    data = np.loadtxt(f,skiprows=1,delimiter=',')
    xss[i,:] = data[:151,0]
    yss[i,:] = data[:151,1]
    zss[i,:] = data[:151,2]
    xps[i,:] = np.flipud(data[150:,0])
    yps[i,:] = np.flipud(data[150:,1])
    zps[i,:] = np.flipud(data[150:,2])
    
    # Airfoil 3D - Centroid 
    cx[i] = np.hstack([xps[i,:],xss[i,:]]).mean()
    cz[i] = np.hstack([zps[i,:],zss[i,:]]).mean()
    cy[i] = np.hstack([yps[i,:],yss[i,:]]).mean()

nspan = 100
R37 = Airfoil3D([],[],0)
R37.b3 = bezier3(cx,cy,cz)
R37.shft_ss = np.zeros(shape=(nspan,151,3))
R37.shft_ps = np.zeros(shape=(nspan,151,3))

for i in range(151):        # Lets get a higher resolution airfoil
    z = np.linspace(zss[0,i],zss[-1,i],nspan)
    R37.shft_ss[:,i,0] = csapi(zss[:,i],xss[:,i],z)
    R37.shft_ss[:,i,1] = csapi(zss[:,i],yss[:,i],z)
    R37.shft_ss[:,i,2] = z

    z = np.linspace(zps[0,i],zps[-1,i],nspan)
    R37.shft_ps[:,i,0] = csapi(zps[:,i],xps[:,i],z)
    R37.shft_ps[:,i,1] = csapi(zps[:,i],yps[:,i],z)
    R37.shft_ps[:,i,2] = z

# Build control_ss and control_ps as 3D arrays (nfiles, npts, 3)
R37.control_ss = np.zeros((len(files),151,3))
R37.control_ps = np.zeros((len(files),151,3))
R37.control_ss[:,:,0] = xss
R37.control_ss[:,:,1] = yss
R37.control_ss[:,:,2] = zss
R37.control_ps[:,:,0] = xps
R37.control_ps[:,:,1] = yps
R37.control_ps[:,:,2] = zps

# Build centroid-based stacking bezier
centroid = np.hstack([cx, cy, cz])
R37.stack_bezier_ctrl_pts = centroid
R37.stack_bezier = bezier3(centroid[:,0], centroid[:,1], centroid[:,2])

R37.ss = copy.deepcopy(R37.shft_ss)
R37.ps = copy.deepcopy(R37.shft_ps)
R37.zz = R37.shft_ss[:,0,2]
R37.te_center = None
R37.bImportedBlade = True
R37.stackType=2 # Centroid
R37.span = max(z)-min(z)
R37.spanwise_spline_fit()
R37.nspan = nspan



# Reading hub and shroud 
hub = np.loadtxt('hub_R37.dat')
shroud = np.loadtxt('shroud_R37.dat')
passage = Passage2D(airfoil_array=[R37],spacing_array=[0])
passage.add_endwalls(zhub=hub[:,0],rhub=hub[:,2],zshroud=shroud[:,0],rshroud=shroud[:,2])
passage.blade_fit(hub[:,0].min())
passage.plot2D()