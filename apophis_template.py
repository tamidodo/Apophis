# coding: utf-8

# In[1]:

import rebound as rb
import numpy as np
import sys
import reboundx
from reboundx import constants as rbxConstants

SIG=3.0 # strength of normal perturbations for exploring error elipse
npert=10000 # number of LOV samples

vpert = -0.01/(npert) # strength and direction of incremental perturbation for LOV in m/s
vpert0 =0.0 # centre of perturbation (used to step forward with finer resolution)

SEED=3205

# set our initial conditions
date="2022-01-21 00:00"
DTDIRECT=1

#solution last updated 2021-Jun-29 11:09:44
# small body pert ephem used: SB441-N16, planetary ephem used: DE441
# epoch 2459600.5 (2022-Jan-21.0) TDB
A2 = -2.901085583204654E-14	# Yarkovsky parameter
A2SIG=1.942E-16
betaAsteroid = 5E-13*(365.25/(2*np.pi))**2  # takes into account radiation effects other than Yarkovsy (A1)

twopi=2*np.pi

aum = rb.units.lengths_SI['au']
au=aum*1e2
aukm=aum/1e3 # au in km
msun=rb.units.masses_SI['msun']*1e3  # grams
code2sec = rb.units.times_SI['yr2pi']
sec2code = 1./code2sec
vcode2cmps=au*sec2code

A2PARAM = A2/( (1./365.25)*np.pi*2 )**2     # convert to code units
A2SIGPARAM=A2SIG/( (1./365.25)*np.pi*2 )**2

Rearth=6371. # km
Vesc = 11.186 # km/s

# used for Earth's J2 component
RE_eq = 6378.135/aukm
J2=1.0826157e-3
J4=-1.620e-6

dmin=4.326e-5 # Earth radius in au

tsimend=100 # simulation endtime in years
dtout=100 # end time in years
Noutputs = 100
tbfreeze=7.3*2*np.pi # Time of the 2029 flyby (years after the IC date * 2 pi)
# in between tbfreeze_start and tbfreeze_end, the b-plane location is frozen for the closest approach
# (so we don't keep recording b-plane coordinates on subsequent close approaches)
tbfreeze_start = 6*2*np.pi # We use 6 for the 2029 b-plane, and 13 for the 2036 b-plane
tbfreeze_end = 8*2*np.pi # We use 8 for the 2029 b-plane and 15 for the 2036 b-plane

TMAX=-2.3e9

vpert = vpert/vcode2cmps
vpert0 = vpert0/vcode2cmps

# add particles to sim

sim=rb.Simulation()
sim.integrator = "ias15" # IAS15 is the default integrator, so we actually don't need this line but just in case that changes...
sim.dt=0.01*DTDIRECT

print(rb.units.times_SI['yr2pi'])

# units of perturbations from JPL are in au and au/d so we need to convert the speed part. A2 in last
pertvals=np.array([2.52210990E-09,3.80951833E-09,2.24702572E-09,5.21392275E-11,2.66039622E-11,5.29135900E-11,A2SIGPARAM])

pertvals[3:6]*=code2sec/(24*3600)
pertvals=pertvals*SIG # scale 1 sigma perturbations to strength set in SIG variable

'''
sim.add("Sun",date=date)
sim.add("Mercury",date=date)
sim.add("Venus",date=date)
sim.add("Geocenter",date=date)
sim.add("Luna",date=date)
sim.add("Mars",date=date)
sim.add("Jupiter",date=date)
sim.add("Saturn",date=date)
sim.add("Uranus",date=date)
sim.add("Neptune",date=date)
sim.add("Apophis",date=date)
sim.add("Ceres",date=date)
sim.add("Vesta",date=date)
sim.add("Pallas",date=date)
sim.add("Hygiea",date=date)
sim.add("Davida",date=date)
sim.add("Interamnia",date=date)
sim.add("Eunomia",date=date)
sim.add("A804 RA",date=date)  # Juno
sim.add("Psyche",date=date)
sim.add("A858 CA",date=date)  # 52 Europa
sim.add("Thisbe",date=date)
sim.add("Sylvia",date=date)
sim.add("A847 NA",date=date) # 6 Hebe
sim.add("Cybele",date=date)
sim.add("A847 PA",date=date)  # 7 Iris
sim.add("Amphitrite",date=date)
sim.add("Herculina",date=date)
sim.add("Bamberga",date=date)
sim.add("Aletheia",date=date)
sim.add("Parthenope",date=date)
sim.add("Melete",date=date)
sim.add("Irene",date=date)
sim.add("Aurelia",date=date)
sim.add("Ausonia",date=date)
sim.add("Hertha",date=date)
'''

masteroid=[0,63.13,17.29,13.73,5.78,2.26,2.19,2.10,1.82,1.81,1.59,1.02,0.99,0.93,0.91,0.86,0.86,0.77,0.69,0.52,0.39,0.31,0.19,0.12,0.10,0.08]
masteroid=np.array(masteroid)
masteroid=masteroid/(aukm**3*sec2code**2)
print(masteroid)
names=["Sun","Mercury","Venus","Earth","Luna","Mars","Jupiter","Saturn","Uranus","Neptune","Apophis","Ceres","Vesta","Pallas","Hygiea","Davida","Interamnia","Eunomia","A804 RA","Psyche","A858 CA","Thisbe","Sylvia","Hebe","Cybele","A847 PA","Amphitrite","Herculina","Bamberga","Aletheia","Parthenope","Melete","Irene","Aurelia","Ausonia","Hertha"]

idEarth=3
idLuna=4
idAst=10


ps=sim.particles
for i in range(idAst,len(ps)):
   ps[i].m=masteroid[i-(idAst)]

'''
for iter,p in enumerate(ps):
   print("  sim.add(m={}, x={}, y={}, z={}, vx={}, vy={}, vz={}, hash='{}')".format(p.m, p.x, p.y, p.z, p.vx, p.vy, p.vz,names[iter]))

sys.exit(-1)
'''

np.random.seed(SEED)

for iloop in range(npert):

  sim=rb.Simulation()
  sim.integrator = "ias15" # IAS15 is the default integrator, so we actually don't need this line but just in case that changes...
  sim.dt=0.01*DTDIRECT

  # epoch 2022 Jan 21

  sim.add(m=1.0, x=-0.008645752226103389, y=0.0031777393722842036, z=0.0001759518809276825, vx=-0.0001826161440618218, vy=-0.0004932536835693154, vz=8.221380098083738e-06, hash='Sun')
  sim.add(m=1.6601141530543488e-07, x=-0.1092191120401082, y=0.30020659148499557, z=0.03367390382288566, vx=-1.8776538733671397, vy=-0.46454764356467837, vz=0.13430158634776393, hash='Mercury')
  sim.add(m=2.4478382877847715e-06, x=-0.4489559156593811, y=0.5699234947139606, z=0.03336213483046784, vx=-0.9327841828310308, vy=-0.7282494599620154, vz=0.04383336656148341, hash='Venus')
  sim.add(m=3.0034896149157645e-06, x=-0.5096379321842845, y=0.8501219972517323, z=0.00013411610325436964, vx=-0.8770808558990811, vy=-0.5131311801283556, vz=4.21861346502036e-05, hash='Earth')
  sim.add(m=3.694303310687701e-08, x=-0.5120048246663336, y=0.8512284348796718, z=0.0003643285922741746, vx=-0.8901995819351388, vy=-0.5441859907279303, vz=-0.00036268121378890114, hash='Luna')
  sim.add(m=3.2271560375549977e-07, x=-0.6208275045580232, y=-1.3777901085572537, z=-0.013749215822740257, vx=0.7740763312319248, vy=-0.2603460539307298, vz=-0.024430191817414834, hash='Mars')
  sim.add(m=0.0009547919152112404, x=4.699719663231086, y=-1.6424715721399088, z=-0.09833025824389338, vx=0.13947154994535904, vy=0.43469102695640327, vz=-0.004923764640054954, hash='Jupiter')
  sim.add(m=0.0002858856727222417, x=7.024350486954772, y=-6.984805543682808, z=-0.15821877383547872, vx=0.21052327889909733, vy=0.2293035798114121, vz=-0.01237017327084789, hash='Saturn')
  sim.add(m=4.36624373583127e-05, x=14.334571989635778, y=13.535698628409701, z=-0.13543496853663906, vx=-0.15865359681653518, vy=0.15558854767611008, vz=0.0026332475900192804, hash='Uranus')
  sim.add(m=5.151383772628674e-05, x=29.632786271397276, y=-4.024706676604673, z=-0.6000365686873408, vx=0.023350902936417825, vy=0.1819081175658396, vz=-0.0042841797886757895, hash='Neptune')
  sim.add(m=0.0, x=-1.0456881321299374, y=0.35094258894981506, z=-0.04293624008458587, vx=-0.23893409441254912, vy=-0.8280243776160886, vz=0.03847012036003141, hash='Apophis')
  sim.add(m=4.756901461539691e-10, x=0.6134651293198201, y=2.6326861851172394, z=-0.03135685721199767, vx=-0.5974636226513389, vy=0.09470968635305152, vz=0.11305505428040578, hash='Ceres')
  sim.add(m=1.302816826707132e-10, x=-0.2537961833901161, y=-2.1355824168206285, z=0.09394477256841491, vx=0.6993562098038364, vy=-0.09279166386497639, vz=-0.08235041357465825, hash='Vesta')
  sim.add(m=1.0345676709478846e-10, x=2.812401045506295, y=0.36637704796270953, z=-0.4944082885219968, vx=-0.24010383548267825, vy=0.43946002692246533, vz=-0.2841929607041755, hash='Pallas')
  sim.add(m=4.355281236765312e-11, x=-2.7300381318866234, y=-0.8343425932426698, z=-0.19004238743922028, vx=0.21888140147239332, vy=-0.5788909710135058, vz=0.005456775134758909, hash='Hygiea')
  sim.add(m=1.7029300337525267e-11, x=-1.5800461164183033, y=-3.268728168345612, z=0.7103432907032581, vx=0.39220783972942, vy=-0.252158458193233, vz=-0.08510938961799432, hash='Davida')
  sim.add(m=1.650184413238068e-11, x=0.6373623438433746, y=-2.929534673640189, z=0.03511695266785033, vx=0.5162115040293073, vy=0.20741068738203866, vz=0.16992767076506538, hash='Interamnia')
  sim.add(m=1.582368615433764e-11, x=-2.9772548882483507, y=-0.6267852318419523, z=-0.6196716972094697, vx=0.0736014757862834, vy=-0.5116841604257653, vz=-0.02728106698265056, hash='Eunomia')
  sim.add(m=1.3713861333759286e-11, x=1.1156844377669035, y=-2.6439267249610037, z=0.5556261239049757, vx=0.4462017237551057, vy=0.3168336232280958, vz=-0.09021188624168808, hash='A804 RA')
  sim.add(m=1.363851044731006e-11, x=-2.908804762541061, y=1.3100432325664344, z=0.017300720281454903, vx=-0.2681926163392616, vy=-0.46321912944428917, vz=0.028938957206562185, hash='Psyche')
  sim.add(m=1.1980790945427069e-11, x=-2.938162076018296, y=-0.08719409578776852, z=0.3081222191826396, vx=-0.036378911314455996, vy=-0.5933927222251281, vz=0.052274176956182215, hash='A858 CA')
  sim.add(m=7.685790417821139e-12, x=-2.751744351691551, y=-0.8866164507391334, z=-0.25761869754336764, vx=0.25850692329107355, vy=-0.5127823275890946, vz=0.01821827616355195, hash='Thisbe')
  sim.add(m=7.459737758473457e-12, x=-3.0629083161534103, y=-1.9355287680765512, z=0.4526321343496967, vx=0.2970393318575025, vy=-0.4089626453704056, vz=-0.07754206231866533, hash='Sylvia')
  sim.add(m=7.007632439778097e-12, x=1.9149189465907406, y=-0.06227664322461985, z=-0.321270381348591, vx=-0.04247766663231282, vy=0.7676620604298795, vz=-0.14431524444876684, hash='Hebe')
  sim.add(m=6.856930666879643e-12, x=2.6919594925709767, y=-1.7966120477366687, z=0.0327589338014858, vx=0.36162644826168855, vy=0.43869003677103396, vz=-0.03423694570104768, hash='Cybele')
  sim.add(m=6.480176234633508e-12, x=-0.9207818973068471, y=1.8950670934502356, z=-0.11969335020179461, vx=-0.7032716214880328, vy=-0.18303156277337448, vz=-0.06357303224217088, hash='A847 PA')
  sim.add(m=6.480176234633508e-12, x=-1.8678794400564265, y=-1.9923865576435447, z=-0.22472872657131027, vx=0.42053481113981384, vy=-0.4013705600558867, vz=-0.039753212906407265, hash='Amphitrite')
  sim.add(m=5.8020182565904676e-12, x=2.986807314578001, y=0.8747649858004135, z=-0.91164038905349, vx=-0.1630303877151722, vy=0.47690921464257724, vz=0.0035777525699662513, hash='Herculina')
  sim.add(m=5.199211164996652e-12, x=0.7745758800597393, y=-2.103050564443327, z=-0.26809171000888066, vx=0.5578847042858924, vy=0.4267494558567802, vz=0.12928793368771285, hash='Bamberga')
  sim.add(m=3.918246095359796e-12, x=2.5475442208760923, y=2.28336873297592, z=-0.46276967528073654, vx=-0.2989591130525284, vy=0.40712548971105067, vz=0.06132363365128795, hash='Aletheia')
  sim.add(m=2.938684571519847e-12, x=-1.9909658878201417, y=1.8311355477306643, z=0.044817455803290336, vx=-0.39324912711139887, vy=-0.4216749851723008, vz=0.045751964951978266, hash='Parthenope')
  sim.add(m=2.3358774799260323e-12, x=0.5370055543077894, y=2.939912931491588, z=-0.38886740287608285, vx=-0.4938054960669484, vy=0.18223714172047067, vz=-0.04098567423472462, hash='Melete')
  sim.add(m=1.43166684253531e-12, x=4.6396562248373225, y=-1.5409864493395842, z=-0.07620838377117704, vx=0.2257090517230397, vy=0.4738445610930855, vz=-0.032712103574278936, hash='Irene')
  sim.add(m=9.042106373907221e-13, x=2.2225072985634218, y=-0.4540649986495798, z=0.13642841975144457, vx=0.29031356453196816, vy=0.6375030806392412, vz=-0.01362021107405816, hash='Aurelia')
  sim.add(m=7.535088644922685e-13, x=1.9124048322172198, y=-1.0746214984791171, z=-0.027046583693988706, vx=0.4003684970505014, vy=0.5698794244013353, vz=0.06876374181030659, hash='Ausonia')
  sim.add(m=6.028070915938148e-13, x=-2.7550708535302744, y=-0.2230525198207123, z=-0.03981419623592202, vx=0.13397311346842286, vy=-0.544140198907304, vz=-0.019443065197370057, hash='Hertha')

  ps=sim.particles

  Nbod=sim.N

# Drawing perturbations evenly from the uncertainty ellipsoid
  multivar = []
  rad_u = (np.random.uniform(0.0,1.0))**(1/3)
  phi_u = np.random.uniform(0.0, 2*np.pi)
  theta_u = np.random.uniform(0.0,np.pi)
  multivar.append(rad_u*np.cos(phi_u)*np.sin(theta_u)*pertvals[0])
  multivar.append(rad_u*np.sin(phi_u)*np.sin(theta_u)*pertvals[1])
  multivar.append(rad_u*np.cos(theta_u)*pertvals[2])
  multivar.append(rad_u*np.cos(phi_u)*np.sin(theta_u)*pertvals[3])
  multivar.append(rad_u*np.sin(phi_u)*np.sin(theta_u)*pertvals[4])
  multivar.append(rad_u*np.cos(theta_u)*pertvals[5])
  multivar.append(np.random.uniform(-pertvals[6],pertvals[6]))
  print("perturbations",multivar)
  ps[idAst].x+=multivar[0]
  ps[idAst].y+=multivar[1]
  ps[idAst].z+=multivar[2]
  ps[idAst].vx+=multivar[3]
  ps[idAst].vy+=multivar[4]
  ps[idAst].vz+=multivar[5]
  A2PERT = A2PARAM+multivar[6]

  def rotateX(x,eps):
    xr=x[0]*1
    yr=x[1]*np.cos(eps)-x[2]*np.sin(eps)
    zr=x[1]*np.sin(eps)+x[2]*np.cos(eps)
    return np.array([xr,yr,zr])

  eps = 23.43651*twopi/360.
  for p in ps:
    x=[p.x,p.y,p.z]
    xr=rotateX(x,eps)
    p.x=xr[0]*1
    p.y=xr[1]*1
    p.z=xr[2]*1
    v=[p.vx,p.vy,p.vz]
    vr=rotateX(v,eps)
    p.vx=vr[0]*1
    p.vy=vr[1]*1
    p.vz=vr[2]*1

  p=ps[idAst]
  print("sim.add(m={}, x={}, y={}, z={}, vx={}, vy={}, vz={}, A2={} hash='{}')".format(p.m,p.x,p.y,p.z,p.vx,p.vy,p.vz,A2PERT,names[idAst]))

  rebx=reboundx.Extras(sim)
  gr=rebx.load_force("gr") # gr is the correction for just the Sun, gr_full is the correction for all bodies in the sim
  rebx.add_force(gr)
  gr.params["c"] = rbxConstants.C

  mig=rebx.load_force("modify_orbits_forces")
  rebx.add_force(mig)

  gh = rebx.load_force("gravitational_harmonics")
  rebx.add_force(gh)
  ps[idEarth].params["J2"] = J2
  ps[idEarth].params["J4"] = J4
  ps[idEarth].params["R_eq"] = RE_eq

  rf = rebx.load_force("radiation_forces")
  rebx.add_force(rf)
  rf.params["c"] = rbxConstants.C
  ps[idAst].params["beta"]=betaAsteroid

  mindis=rebx.load_operator("track_bplane")
  rebx.add_operator(mindis)

  sim.particles[idEarth].r=dmin # set size of Earth
  sim.particles[idLuna].r=dmin/4 # set size of Earth

  ps=sim.particles

  print(ps[idEarth].hash)

# for line of variations, we add an additional perturbation. This selects the direction of the asteroid and amplifies the vectors along track
  vxp=ps[idAst].vx
  vyp=ps[idAst].vy
  vzp=ps[idAst].vz

  vmag = np.sqrt(vxp*vxp+vyp*vyp+vzp*vzp)

  ps[idAst].vx = vxp + vxp/vmag*(vpert*(iloop)+vpert0)
  ps[idAst].vy = vyp + vyp/vmag*(vpert*(iloop)+vpert0)
  ps[idAst].vz = vzp + vzp/vmag*(vpert*(iloop)+vpert0)

  print("iloop {} Vpert {}cm/s vx {} vy {} vz {}".format(iloop,vcode2cmps*(vpert*(iloop)+vpert0),ps[idAst].vx,ps[idAst].vy,ps[idAst].vz))

  ps[idAst].params["tau_a"]=TMAX
  ps[idAst].params["min_distance"]=1e6
  ps[idAst].params["min_distance_from_id"] = 0
  ps[idAst].params["min_distance_from_id2"] = idEarth
  ps[idAst].params["min_distance_time"]=0.
  ps[idAst].params["trackx"]=0.
  ps[idAst].params["tracky"]=0.
  ps[idAst].params["trackz"]=0.
  ps[idAst].params["trackvx"]=0.
  ps[idAst].params["trackvy"]=0.
  ps[idAst].params["trackvz"]=0.
  ps[idAst].params["s1x"]=0.
  ps[idAst].params["s1y"]=0.
  ps[idAst].params["s1z"]=0.
  ps[idAst].params["s1vx"]=0.
  ps[idAst].params["s1vy"]=0.
  ps[idAst].params["s1vz"]=0.
  ps[idAst].params["s2x"]=0.
  ps[idAst].params["s2y"]=0.
  ps[idAst].params["s2z"]=0.
  ps[idAst].params["s2vx"]=0.
  ps[idAst].params["s2vy"]=0.
  ps[idAst].params["s2vz"]=0.
  ps[idAst].params["atractor"]=0.

  sim.move_to_com()
# We always move to the center of momentum frame before an integration

  year = 2.*np.pi # One year in units where G=1
  times = np.linspace(0.,year*tsimend, Noutputs)
  dtout = dtout*2*np.pi
  x = np.zeros((Nbod,Noutputs))
  y = np.zeros((Nbod,Noutputs))
  z = np.zeros((Nbod,Noutputs))
  vx = np.zeros((Nbod,Noutputs))
  vy = np.zeros((Nbod,Noutputs))
  vz = np.zeros((Nbod,Noutputs))

  sim.collision="direct"
  sim.collision_resolve="merge"
  sim.collision_resolve_keep_sorted = 1
  sim.track_energy_offset = 1

  ps = sim.particles       # ps is now an array of pointers and will change as the simulation runs

  tout=0
  iout=0

  dminafter=[]
  tdminafter=[]
  a_asteroid=[]
  dlast=1e6
  timeCloseApp=0.
  # flybypassed=0
  for i,time in enumerate(times):
      # Uncomment line 331 and the block below to modify Yarkovsky effect in 2029 by instantaneously scaling A2PERT
      # Note that we can't scale it to 0 because we can't divide by 0 but we can make it really small
      '''
      if flybypassed == 0:
          if time > tbfreeze:
             A2PERT = 1.1*A2PERT
             flybypassed=1
      '''
      dadt = 2.*A2PERT/(ps[idAst].n*ps[idAst].a**2*(1.-ps[idAst].e**2))
      TMAX = ps[idAst].a/dadt
      ps[idAst].params["tau_a"]=TMAX

      sim.integrate(time)
      if(time%1):sim.save("checkpoint.bin")
      if(time>=tout and PRINTOUT==1):
          output=sim.particles_ascii()
          fh=open("output.ascii."+repr(iout),"w")
          fh.write(output)
          fh.close()
          tout+=dtout
          iout+=1

      for p in range(len(ps)):
          x[p][i]=ps[p].x
          y[p][i]=ps[p].y
          z[p][i]=ps[p].z
          vx[p][i]=ps[p].vx
          vy[p][i]=ps[p].vy
          vz[p][i]=ps[p].vz


      try:
         tcol=ps[idEarth].lastcollision/(year)
         if tcol>0:
            dminafter.append(RE_eq)
            tdminafter.append(0.)
            break
      except: tcol=0

      if time > tbfreeze_end:
           dminafter.append(ps[idAst].params["min_distance"])
           tdminafter.append(ps[idAst].params["min_distance_time"])

      if time > tbfreeze_start and time < tbfreeze_end:
         if ps[idAst].params["min_distance"] < dlast:
           dlast = ps[idAst].params["min_distance"]
           timeCloseApp=ps[idAst].params["min_distance_time"]
           Xastr=np.array([ps[idAst].params["trackx"],ps[idAst].params["tracky"], ps[idAst].params["trackz"]])
           Vastr=np.array([ps[idAst].params["trackvx"],ps[idAst].params["trackvy"],ps[idAst].params["trackvz"]])

           Xearth=np.array([ps[idAst].params["s2x"],ps[idAst].params["s2y"],ps[idAst].params["s2z"]])
           Vearth=np.array([ps[idAst].params["s2vx"],ps[idAst].params["s2vy"],ps[idAst].params["s2vz"]])

           Xsun=np.array([ps[idAst].params["s1x"],ps[idAst].params["s1y"],ps[idAst].params["s1z"]])
           Vsun=np.array([ps[idAst].params["s1vx"],ps[idAst].params["s1vy"],ps[idAst].params["s1vz"]])

      ps[idAst].params["min_distance"]=1e6
      a_asteroid.append(ps[idAst].a)

  dminafter=np.array(dminafter)
  tdminafter=np.array(tdminafter)

  impact=-1
  try: tcol=ps[idEarth].lastcollision/(year)
  except: tcol=0
  if tcol > 0:
    print("Collision at time {} yr".format(tcol))
    impact=1.
    closeEncounterDist = Rearth
    closeEncounterTime = tcol
  else:
    closeEncounterDist = ps[idAst].params["min_distance"]*aukm
    closeEncounterTime = ps[idAst].params["min_distance_time"]/year
    print("Minimum Earth-asteroid Distance {} km at time {} yr".format(dlast*aukm, timeCloseApp/(2*np.pi )))


  distance=np.sqrt(np.square(x[idEarth]-x[idAst])+np.square(y[idEarth]-y[idAst]) + np.square(z[idEarth]-z[idAst]))*aukm
  idmin=np.argmin(distance)
  print("Close approach of {} from snap shots at time {}".format(distance[idmin], times[idmin]/year)) 

  #if distance[idmin] < 10*Rearth: idmin-=1
  print("Using distance for asymptote of {} from snap shots at time {}".format(distance[idmin],times[idmin]/year)) 

#Xsun  =np.array([x[0,idmin],y[0,idmin],z[0,idmin]])
#Vsun  =np.array([vx[0][idmin],vy[0][idmin],vz[0][idmin]])
#Xearth=np.array([x[idEarth,idmin],y[idEarth,idmin],z[idEarth,idmin]])
#Vearth=np.array([vx[idEarth][idmin],vy[idEarth][idmin],vz[idEarth][idmin]])
#Xastr =np.array([x[idAst,idmin],y[idAst,idmin],z[idAst,idmin]])
#Vastr =np.array([vx[idAst,idmin],vy[idAst,idmin],vz[idAst,idmin]])


  def dotpro(a,b):
    s=0.
    for i in range(len(a)): s+=a[i]*b[i]
    return s

# now we get the B-plane coordinate

  VREL = (Vastr-Vearth)
  Xdiff = (Xastr-Xearth)
  radial = np.sqrt(dotpro(Xdiff,Xdiff))

  h = np.cross(Xdiff,VREL)
  nhat = h/np.sqrt(dotpro(h,h))
  VREL2 = dotpro(VREL,VREL)

  evec = np.cross(VREL,h)/ps[idEarth].m-Xdiff/radial #(VREL2/ps[idEarth].m-1/radial)*Xdiff-dotpro(Xdiff,VREL)/ps[idEarth].m*VREL

  aenc = 1./(2/radial-VREL2/ps[idEarth].m)

  vinf = np.sqrt(ps[idEarth].m/(-aenc))

  ecc = np.sqrt(dotpro(evec,evec))

  benc = abs(aenc)*np.sqrt(ecc**2-1.)

  alpha = np.arccos(1/ecc)

  ehat = evec/ecc

  Qhat = np.cross(nhat,ehat)

  Shat = ehat/ecc + np.sqrt(ecc**2-1.)*Qhat/ecc

  lam = np.sqrt(1 + 2*ps[idEarth].m/(vinf**2*6378/aukm))

  B = np.cross(Shat,h)/vinf*aukm # lam
  minDist=np.sqrt(dotpro(B,B))

#V = np.array([0,0,-1])
  V=Vearth-Vsun

  Xi0 = np.cross(V,Shat)
  xihat = Xi0/np.sqrt(dotpro(Xi0,Xi0))

  zetahat = -np.cross(Shat,xihat)

  Xi = dotpro(B,xihat)
  Zeta = dotpro(B,zetahat)

  Vinf = np.sqrt(dotpro(VREL,VREL)-2*ps[idEarth].m/radial)*vcode2cmps/1e5
  Rimpact = Rearth*np.sqrt(1. + (Vesc/Vinf)**2)


  CorrectPrediction=1
  if impact>0 and minDist > Rimpact: CorrectPrediction=0
  elif impact<0 and Rimpact > minDist: CorrectPrediction=0

  id=np.argmin(dminafter)
  minDistAfter=dminafter[id]*aukm
  tminDistAfter=tdminafter[id]*code2sec/(365.25*24*3600)

  print(f"B-plane distance {minDist} Zeta {Zeta} Xi {Xi} Rimpact {Rimpact} minDistAfter {minDistAfter} tminDistAfter {tminDistAfter} Vinf {Vinf} [km or km/s] CorrectPrediction{CorrectPrediction}")

