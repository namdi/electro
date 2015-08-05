""" This program holds function used for loading data
"""

#-----------------------------------------------------------
# imports
#-----------------------------------------------------------
import numpy as np

#-----------------------------------------------------------
# Function Declarations
#-----------------------------------------------------------
def load_pos_new(fpath, pos_src, pos_targ, dosource, dotarget, dofmm, dodirect):
  
  """ load the data """
  for i in range(pos_src.shape[0]):
    fpath2 = fpath + "/t_" + str(i)
    
    if dosource:
      if dofmm:
	pos_src[i,:,:] 	= np.loadtxt( fpath2 + "/pos_source.txt").view('f8')
      elif dodirect:
	pos_src[i,:,:] 	= np.loadtxt( fpath2 + "/d_pos_source.txt").view('f8')
	
    if dotarget:
      if dofmm:
	pos_targ[i,:,:]	= np.loadtxt( fpath2 + "/pos_target.txt").view('f8')
      elif dodirect:
	pos_targ[i,:,:]	= np.loadtxt( fpath2 + "/d_pos_target.txt").view('f8')
  
  return
  
def load_vel_new(fpath, vel_src, vel_targ, dosource, dotarget, dofmm, dodirect):
  
  """ load the data """
  for i in range(vel_src.shape[0]):
    fpath2 = fpath + "/t_" + str(i)
    
    if dosource:
      if dofmm:
	vel_src[i,:,:] 	= np.loadtxt( fpath2 + "/vel_source.txt").view('f8')
      elif dodirect:
	vel_src[i,:,:] 	= np.loadtxt( fpath2 + "/d_vel_source.txt").view('f8')
	
    if dotarget:
      if dofmm:
	vel_targ[i,:,:]	= np.loadtxt( fpath2 + "/vel_target.txt").view('f8')
      elif dodirect:
	vel_targ[i,:,:]	= np.loadtxt( fpath2 + "/d_vel_target.txt").view('f8')
  
  return  

def load_pos_vel_accel_txt(fpath, pos_src, vel_src, accel_src, pos_targ, vel_targ, accel_targ):
  
  """ load the data """
  for i in range(pos_src.shape[0]):
    fpath2 = fpath + "/t_" + str(i)
    pos_src[i,:,:] 	= np.loadtxt( fpath2 + "/pos_source.txt").view('f8')
    vel_src[i,:,:] 	= np.loadtxt( fpath2 + "/vel_source.txt").view('f8')
    accel_src[i,:,:] 	= np.loadtxt( fpath2 + "/accel_source.txt").view('f8')
    
    pos_targ[i,:,:]	= np.loadtxt( fpath2 + "/pos_target.txt").view('f8')
    vel_targ[i,:,:]	= np.loadtxt( fpath2 + "/vel_target.txt").view('f8')
    accel_targ[i,:,:]	= np.loadtxt( fpath2 + "/accel_target.txt").view('f8')

def load_pot_field_txt(fpath, pot_src, field_src, pot_targ, field_targ):
  """ load the data
  It is good to load the data as floats even if the data is complex. 
  If a vector field is complex the first 3 columns correspond to the 
  real components of (x,y,z). The 2nd 3 columns correpsond to the imaginary 
  part of (x,y,z) components"""
  
  dtype = 'f8'
  
  for i in range(pot_src.shape[0]):
    fpath2 = fpath + "/t_" + str(i)
    
    temp 		= np.loadtxt( fpath2 + "/pot_source.txt").view(dtype)
    #temp		= temp.reshape(temp.shape[0],)
    pot_src[i,:]	= temp[:,0]
    
    temp		= np.loadtxt( fpath2 + "/field_source.txt").view(dtype)
    field_src[i,:,:] 	= temp[:,0:3]
    
    temp 		= np.loadtxt( fpath2 + "/pot_target.txt").view(dtype)
    #temp		= temp.reshape(temp.shape[0],)
    pot_targ[i,:]	= temp[:,0]
    
    temp		= np.loadtxt( fpath2 + "/field_target.txt").view(dtype)
    field_targ[i,:,:]	= temp[:,0:3]
    

    
def load_pot_field(fpath, pot_src, field_src, pot_targ, field_targ, dofmm, dodirect, dosource, dotarget, dopot, dofield):
  """ load the data
  This loads the data from a directory fpath. It correpsonds to only 1 timestep.
  It is good to load the data as floats even if the data is complex. 
  If a vector field is complex the first 3 columns correspond to the 
  real components of (x,y,z). The 2nd 3 columns correpsond to the imaginary 
  part of (x,y,z) components"""
  
  dtype = 'f8'
    
  
  if dofmm:
    pot_src_name = "/pot_source.txt"
    pot_targ_name = "/pot_target.txt"
    
    field_src_name = "/field_source.txt"
    field_targ_name = "/field_target.txt"
  elif dodirect:
    pot_src_name = "/d_pot_source.txt"
    pot_targ_name = "/d_pot_target.txt"
    
    field_src_name = "/d_field_source.txt"
    field_targ_name = "/d_field_target.txt"
    
  if dosource:
    if dopot:
      # this is done because we load the Nx1 complex  array as a N x 2 real array
      temp 	= np.loadtxt( fpath + pot_src_name).view(dtype)
      n 	= temp.shape[0]
      temp	= np.array(temp[:,0])
      temp 	= np.reshape(temp, (n,1))
      
      pot_src[:,0]	= temp[:,0]
    
    if dofield:
      # this is done because we load the Nx3 complex  array as a N x 6 real array
      # the first 3 columns correspond to the real values (x,y,z)
      temp		= np.loadtxt( fpath + field_src_name).view(dtype)
      field_src[:,:] 	= temp[:,0:3]
    
  if dotarget:
    if dopot:
      temp 		= np.loadtxt( fpath + pot_targ_name).view(dtype)
      n 	= temp.shape[0]
      temp	= np.array(temp[:,0])
      temp 	= np.reshape(temp, (n,1))
      
      pot_targ[:,0]	= temp[:,0]
    
    if dofield:
      temp		= np.loadtxt( fpath + field_targ_name).view(dtype)
      field_targ[:,:]	= temp[:,0:3]
      
  return
  
def load_pot_field_multi(fpath, pot_far_src, pot_mpole_src, pot_near_src, field_far_src, field_mpole_src, field_near_src,
	      pot_far_targ, pot_mpole_targ, pot_near_targ, field_far_targ, field_mpole_targ, field_near_targ, 
	      dosource, dotarget, dopot, dofield):
  """ load the data
  This loads the data from a directory fpath. It correpsonds to only 1 timestep.
  It is good to load the data as floats even if the data is complex. 
  If a vector field is complex the first 3 columns correspond to the 
  real components of (x,y,z). The 2nd 3 columns correpsond to the imaginary 
  part of (x,y,z) components"""
  
  dtype = 'f8'
        
  if dosource:
    if dopot:
      # this is done because we load the Nx1 complex  array as a N x 2 real array
      temp 		= np.loadtxt( fpath + "/pot_far_source.txt").view(dtype)
      n 		= temp.shape[0]
      temp		= np.array(temp[:,0])
      temp 		= np.reshape(temp, (n,1))
      
      pot_far_src[:,0]	= temp[:,0]
      
      temp 		= np.loadtxt( fpath + "/pot_mpole_source.txt").view(dtype)
      n 		= temp.shape[0]
      temp		= np.array(temp[:,0])
      temp 		= np.reshape(temp, (n,1))
      
      pot_mpole_src[:,0]= temp[:,0]
      
      temp 		= np.loadtxt( fpath + "/pot_near_source.txt").view(dtype)
      n 		= temp.shape[0]
      temp		= np.array(temp[:,0])
      temp 		= np.reshape(temp, (n,1))
      
      pot_near_src[:,0]	= temp[:,0]
      
    if dofield:
      # this is done because we load the Nx3 complex  array as a N x 6 real array
      # the first 3 columns correspond to the real values (x,y,z)
      temp			= np.loadtxt( fpath + "/field_far_source.txt").view(dtype)
      field_far_src[:,:] 	= temp[:,0:3]
    
      temp			= np.loadtxt( fpath + "/field_mpole_source.txt").view(dtype)
      field_mpole_src[:,:] 	= temp[:,0:3]
      
      temp			= np.loadtxt( fpath + "/field_near_source.txt").view(dtype)
      field_near_src[:,:] 	= temp[:,0:3]
      
  if dotarget:
    if dopot:
      temp 		= np.loadtxt( fpath + "/pot_far_target.txt").view(dtype)
      n 		= temp.shape[0]
      temp		= np.array(temp[:,0])
      temp 		= np.reshape(temp, (n,1))
      
      pot_far_targ[:,0]	= temp[:,0]
      
      temp 		= np.loadtxt( fpath + "/pot_mpole_target.txt").view(dtype)
      n 		= temp.shape[0]
      temp		= np.array(temp[:,0])
      temp 		= np.reshape(temp, (n,1))
      
      pot_mpole_targ[:.0]= temp[:,0]
      
      temp 		= np.loadtxt( fpath + "/pot_near_target.txt").view(dtype)
      n 		= temp.shape[0]
      temp		= np.array(temp[:,0])
      temp 		= np.reshape(temp, (n,1))
      
      pot_near_targ[:,0] =  temp[:,0]
      
    if dofield:
      temp			= np.loadtxt( fpath + "/field_far_target.txt").view(dtype)
      field_far_targ[:,:] 	= temp[:,0:3]
    
      temp			= np.loadtxt( fpath + "/field_mpole_target.txt").view(dtype)
      field_mpole_targ[:,:] 	= temp[:,0:3]
      
      temp			= np.loadtxt( fpath + "/field_near_target.txt").view(dtype)
      field_near_targ[:,:] 	= temp[:,0:3]
  return


""" this loads the direct data   """
def load_direct_pos_vel_accel_txt(fpath, pos_src, vel_src, accel_src, pos_targ, vel_targ, accel_targ):
  
  """ load the data """
  for i in range(pos_src.shape[0]):
    fpath2 = fpath + "/t_" + str(i)
    pos_src[i,:,:] 	= np.loadtxt( fpath2 + "/d_pos_source.txt").view('f8')
    vel_src[i,:,:] 	= np.loadtxt( fpath2 + "/d_vel_source.txt").view('f8')
    accel_src[i,:,:] 	= np.loadtxt( fpath2 + "/d_accel_source.txt").view('f8')
    
    pos_targ[i,:,:]	= np.loadtxt( fpath2 + "/d_pos_target.txt").view('f8')
    vel_targ[i,:,:]	= np.loadtxt( fpath2 + "/d_vel_target.txt").view('f8')
    accel_targ[i,:,:] 	= np.loadtxt( fpath2 + "/d_accel_target.txt").view('f8')
    
def load_direct_pot_field_txt(fpath, pot_src, field_src, pot_targ, field_targ):
  """ load the data """
  
  dtype = "f8"
  for i in range(pot_src.shape[0]):
    fpath2 = fpath + "/t_" + str(i)
    
    
    temp 		= np.loadtxt( fpath2 + "/d_pot_source.txt").view(dtype)
    #temp		= temp.reshape(temp.shape[0],)  
    pot_src[i,:]	= temp[:,0]
    
    temp		= np.loadtxt( fpath2 + "/d_field_source.txt").view(dtype)
    field_src[i,:,:] 	= temp[:,0:3]
    
    
    temp 		= np.loadtxt( fpath2 + "/d_pot_target.txt").view(dtype)
    #temp		= temp.reshape(temp.shape[0],)
    pot_targ[i,:]	= temp[:,0]
    
    temp		= np.loadtxt( fpath2 + "/d_field_target.txt").view(dtype)  
    field_targ[i,:,:]	= temp[:,0:3]
    
  return

""" This function loads a potential"""
def load_pot(fpath, fname, pot):
  """ load the data """
  
  for i in range(pot.shape[0]):
    
    fpath2 = fpath + "/t_" + str(i)
    fname2 = fpath2 + "/" + fname
    
    pot[i,:]	= np.loadtxt(fname2).view('f8')
  return
  
""" This function loads vector data from a files"""  
def load_vector(fpath, fname, data):
  
  nsteps = data.shape[0]
  npart  = data.shape[1]
  
  for i in range(nsteps):
    fpath2 = fpath + "/t_" + str(i)
    fname2 = fpath2 + "/" + fname
    
    data[i,:,:] = np.loadtxt(fname2).view('f8')    
    
  return
  
def load_far_near_pot_field_targ_txt(fpath, pot_far_targ, field_far_targ, 
				      pot_mpole_targ, field_mpole_targ, pot_near_targ, field_near_targ):
  """ load the data """
  
  dtype = "f8"
  for i in range(pot_far_targ.shape[0]):
    fpath2 = fpath + "/t_" + str(i)
    
    temp 			= np.loadtxt( fpath2 + "/pot_far_target.txt").view(dtype)
    #temp			= temp.reshape(temp.shape[0],)
    pot_far_targ[i,:]		= temp[:,0]
    
    temp			= np.loadtxt( fpath2 + "/field_far_target.txt").view(dtype)
    field_far_targ[i,:,:]	= temp[:,0:3]
    
    temp 			= np.loadtxt( fpath2 + "/pot_mpole_target.txt").view(dtype)
    #temp			= temp.reshape(temp.shape[0],)
    pot_mpole_targ[i,:]		= temp[:,0]
    
    temp			= np.loadtxt( fpath2 + "/field_mpole_target.txt").view(dtype)
    field_mpole_targ[i,:,:]	= temp[:,0:3]
    
    temp 			= np.loadtxt( fpath2 + "/pot_near_target.txt").view(dtype)
    #temp			= temp.reshape(temp.shape[0],)
    pot_near_targ[i,:]		= temp[:,0]
    
    temp			= np.loadtxt( fpath2 + "/field_near_target.txt").view(dtype)
    field_near_targ[i,:,:]	= temp[:,0:3]
    
  return

def load_far_near_pot_field_src_txt(fpath, pot_far_src, field_far_src, pot_mpole_src, field_mpole_src,
				     pot_near_src, field_near_src):
  """ load the data """
  
  dtype = "f8"
  for i in range(pot_far_src.shape[0]):
    fpath2 = fpath + "/t_" + str(i)
    
    temp 			= np.loadtxt( fpath2 + "/pot_far_source.txt").view(dtype)
    #temp			= temp.reshape(temp.shape[0],)
    pot_far_src[i,:]		= temp[:,0]
    
    temp			= np.loadtxt( fpath2 + "/field_far_source.txt").view(dtype)
    field_far_src[i,:,:]	= temp[:,0:3]
    
    temp 			= np.loadtxt( fpath2 + "/pot_mpole_source.txt").view(dtype)
    #temp			= temp.reshape(temp.shape[0],)
    pot_mpole_src[i,:]		= temp[:,0]
    
    temp			= np.loadtxt( fpath2 + "/field_mpole_source.txt").view(dtype)
    field_mpole_src[i,:,:]	= temp[:,0:3]
    
    temp 			= np.loadtxt( fpath2 + "/pot_near_source.txt").view(dtype)
    #temp			= temp.reshape(temp.shape[0],)
    pot_near_src[i,:]		= temp[:,0]
    
    temp			= np.loadtxt( fpath2 + "/field_near_source.txt").view(dtype)
    field_near_src[i,:,:]	= temp[:,0:3]
    
  return
  
  
def load_cplex_potential(data_path, fname, npart, nsteps):
  
  """ This function loads the data for a potential over all timesteps"""
  
  potential	= np.zeros((npart, nsteps+1), dtype=np.dtype('c16'))
  
  for i in range(nsteps+1):
    ftemp = data_path + "/t_" + str(i) + "/" + fname
    
    # have to use temp and reshape it because of array dimension inssues  
    temp 		= np.loadtxt(ftemp).view(potential.dtype)
    
    potential[:,i]	= np.reshape(temp, npart)
    
  return potential