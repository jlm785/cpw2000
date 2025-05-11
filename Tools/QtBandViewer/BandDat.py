import numpy as np

class DosDat_t:
  def ReadFile(self,fname):
    
    HARTREE = 27.21138386
    
    print ("fname is", fname)
    f = open(fname)

    rec      = np.fromfile(f, dtype=np.float32, count=1)
    self.nband,self.norb,self.n = np.fromfile(f, dtype=np.int32, count=3)
    rec      = np.fromfile(f, dtype=np.float32, count=1)

    print(self.nband)
    print(self.norb)
    print(self.n)
  
#    quit()
    
    rec    = np.fromfile(f, dtype=np.float32, count=1)
    self.energy  = np.fromfile(f, dtype=np.float64, count=self.n)
    rec    = np.fromfile(f, dtype=np.float32, count=1)
    
    rec    = np.fromfile(f, dtype=np.float32, count=1)
    self.total_dos  = np.fromfile(f, dtype=np.float64, count=self.n)
    rec    = np.fromfile(f, dtype=np.float32, count=1)

    for i in range(self.n):
      self.energy[i] = self.energy[i]*HARTREE
      self.total_dos[i] = self.total_dos[i]/HARTREE
  
    self.projected_dos = np.empty([self.norb, self.n],dtype=np.float64)
    
    for iorb in range(self.norb):
      rec    = np.fromfile(f, dtype=np.float32, count=1)
      self.projected_dos[iorb,:]  = np.fromfile(f, dtype=np.float64, count=self.n)
      rec    = np.fromfile(f, dtype=np.float32, count=1)

      for i in range(self.n):
        self.projected_dos[iorb,i] = self.projected_dos[iorb,i]/HARTREE

  
class BandDat_t:

  def ReadFile(self,fname):
  
    # Reads the binary file and gets all the needed info
    # BEGIN FILE READ eigenvalues and unfolding
    # open binary file for reading
    f = open(fname)

    rec   = np.fromfile(f, dtype=np.float32, count=1)
    self.title = np.fromfile(f,np.uint8,50).tobytes()
    rec   = np.fromfile(f, dtype=np.float32, count=1)

    rec   = np.fromfile(f, dtype=np.float32, count=1)
    self.subtitle = np.fromfile(f,np.uint8,140).tobytes()  #this needs to be fixed
    rec   = np.fromfile(f, dtype=np.float32, count=1)

    rec   = np.fromfile(f, dtype=np.float32, count=1)
    self.eref  = np.fromfile(f, dtype=np.float64, count=1)[0]
    rec   = np.fromfile(f, dtype=np.float32, count=1)

    rec   = np.fromfile(f, dtype=np.float32, count=1)
    self.xmin,self.xmax  = np.fromfile(f, dtype=np.float64, count=2)
    rec   = np.fromfile(f, dtype=np.float32, count=1)

    rec   = np.fromfile(f, dtype=np.float32, count=1)
    self.ymin,self.ymax  = np.fromfile(f, dtype=np.float64, count=2)
    rec   = np.fromfile(f, dtype=np.float32, count=1)

    rec   = np.fromfile(f, dtype=np.float32, count=1)
    self.nticks = np.fromfile(f, dtype=np.int32, count=1)[0]
    rec   = np.fromfile(f, dtype=np.float32, count=1)

    rec    = np.fromfile(f, dtype=np.float32, count=1)
    self.xklab  = np.fromfile(f, dtype=np.float64, count=self.nticks)
    rec    = np.fromfile(f, dtype=np.float32, count=1)


    self.label    = np.empty(self.nticks, dtype="S6")

    rec     = np.fromfile(f, dtype=np.float32, count=1)
    for i in range(self.nticks):
      self.label[i]   = np.fromfile(f,np.uint8,6).tobytes()
    rec     = np.fromfile(f, dtype=np.float32, count=1)

    self.labels = []

    for i in range(self.nticks):
      self.labels.append(self.label[i].decode("UTF-8").rstrip().lstrip())


    rec      = np.fromfile(f, dtype=np.float32, count=1)
    self.nrk,self.neig = np.fromfile(f, dtype=np.int32, count=2)
    rec      = np.fromfile(f, dtype=np.float32, count=1)


### Insert code here

    rec    = np.fromfile(f, dtype=np.float32, count=1)
    self.xk     = np.fromfile(f, dtype=np.float64, count=self.nrk)
    rec    = np.fromfile(f, dtype=np.float32, count=1)

    rec    = np.fromfile(f, dtype=np.float32, count=1)
    self.bands    = np.fromfile(f, dtype=np.float64,count=self.neig*self.nrk).reshape(self.nrk,self.neig)
    rec    = np.fromfile(f, dtype=np.float32, count=1)

    rec    = np.fromfile(f, dtype=np.float32, count=1)
    self.pkn    = np.fromfile(f, dtype=np.float64,count=self.neig*self.nrk).reshape(self.neig,self.nrk)
    rec    = np.fromfile(f, dtype=np.float32, count=1)

    print(self.title)
    print(self.subtitle)
    print(self.eref)
    print(self.xmin,self.xmax)
    print(self.ymin,self.ymax)
    print(self.neig,self.nrk)
    print(self.xk)

    # BEGIN FILE READ orbital information

    rec   = np.fromfile(f, dtype=np.float32, count=1)
    self.nbaslcao = np.fromfile(f, dtype=np.int32, count=1)[0]
    rec   = np.fromfile(f, dtype=np.float32, count=1)

    print (self.nbaslcao)

    rec   = np.fromfile(f, dtype=np.float32, count=1)
    self.infolcao = np.fromfile(f, dtype=np.int32, count=5*self.nbaslcao).reshape(self.nbaslcao,5)
    rec   = np.fromfile(f, dtype=np.float32, count=1)

    rec   = np.fromfile(f, dtype=np.float32,  count=1)
    self.basxpsi = np.fromfile(f, dtype=np.float64,count=self.nbaslcao*self.neig*self.nrk).reshape(self.nrk,self.neig,self.nbaslcao)
    rec   = np.fromfile(f, dtype=np.float32,  count=1)

    rec   = np.fromfile(f, dtype=np.float32,  count=1)
    self.ntype = np.fromfile(f, dtype=np.int32, count=1)[0]
    rec   = np.fromfile(f, dtype=np.float32,  count=1)
    
    print ("ntype=", self.ntype)
    
    self.nameats    = np.empty(self.ntype, dtype="S2")

    rec     = np.fromfile(f, dtype=np.float32, count=1)
    for i in range(self.ntype):
      self.nameats[i]   = np.fromfile(f,np.uint8,2).tobytes()
    rec     = np.fromfile(f, dtype=np.float32, count=1)

    self.nameat = []

    for i in range(self.ntype):
      self.nameat.append(self.nameats[i].decode("UTF-8").rstrip().lstrip())
    
    print ("nameat=")
    print (self.nameat)
    

    for i in range(self.label.size):
      if (self.labels[i] == "Gamma" or self.labels[i] == "GAMMA" or self.labels[i] == "gamma" ):
        self.labels[i] = r"$\Gamma$"
      if (self.labels[i] == "Lambda" or self.labels[i] == "LAMBDA" or self.labels[i] == "Lambda" ):
        self.labels[i] = r"$\Lambda$"
      if (self.labels[i] == "Sigma" or self.labels[i] == "SIGMA" or self.labels[i] == "sigma" ):
        self.labels[i] = r"$\Sigma$"
      if (self.labels[i] == "Delta" or self.labels[i] == "DELTA" or self.labels[i] == "delta" ):
        self.labels[i] = r"$\Delta$"
      if (self.labels[i] == "--" or self.labels[i] == "--" or self.labels[i] == "--" ):
        self.labels[i] = " "

  def ReadFileFl32(self,fname):
  
    # Reads the binary file and gets all the needed info
    # BEGIN FILE READ eigenvalues and unfolding
    # open binary file for reading
    f = open(fname)

    rec   = np.fromfile(f, dtype=np.float32, count=1)
    self.title = np.fromfile(f,np.uint8,50).tobytes()
    rec   = np.fromfile(f, dtype=np.float32, count=1)

    rec   = np.fromfile(f, dtype=np.float32, count=1)
    self.subtitle = np.fromfile(f,np.uint8,140).tobytes()  #this needs to be fixed
    rec   = np.fromfile(f, dtype=np.float32, count=1)

    rec   = np.fromfile(f, dtype=np.float32, count=1)
    self.eref  = np.fromfile(f, dtype=np.float64, count=1)[0]
    rec   = np.fromfile(f, dtype=np.float32, count=1)

    rec   = np.fromfile(f, dtype=np.float32, count=1)
    self.xmin,self.xmax  = np.fromfile(f, dtype=np.float64, count=2)
    rec   = np.fromfile(f, dtype=np.float32, count=1)

    rec   = np.fromfile(f, dtype=np.float32, count=1)
    self.ymin,self.ymax  = np.fromfile(f, dtype=np.float64, count=2)
    rec   = np.fromfile(f, dtype=np.float32, count=1)

    rec   = np.fromfile(f, dtype=np.float32, count=1)
    self.nticks = np.fromfile(f, dtype=np.int32, count=1)[0]
    rec   = np.fromfile(f, dtype=np.float32, count=1)

    rec    = np.fromfile(f, dtype=np.float32, count=1)
    self.xklab  = np.fromfile(f, dtype=np.float64, count=self.nticks)
    rec    = np.fromfile(f, dtype=np.float32, count=1)


    self.label    = np.empty(self.nticks, dtype="S6")

    rec     = np.fromfile(f, dtype=np.float32, count=1)
    for i in range(self.nticks):
      self.label[i]   = np.fromfile(f,np.uint8,6).tobytes()
    rec     = np.fromfile(f, dtype=np.float32, count=1)

    self.labels = []

    for i in range(self.nticks):
      self.labels.append(self.label[i].decode("UTF-8").rstrip().lstrip())


    rec      = np.fromfile(f, dtype=np.float32, count=1)
    self.nlines = np.fromfile(f, dtype=np.int32, count=1)[0]
    rec      = np.fromfile(f, dtype=np.float32, count=1)

#    print("self.nlines=", self.nlines)

    rec    = np.fromfile(f, dtype=np.float32, count=1)
    self.xk_start  = np.fromfile(f, dtype=np.float64, count=self.nlines)
    rec    = np.fromfile(f, dtype=np.float32, count=1)

    rec    = np.fromfile(f, dtype=np.float32, count=1)
    self.xk_end  = np.fromfile(f, dtype=np.float64, count=self.nlines)
    rec    = np.fromfile(f, dtype=np.float32, count=1)




    rec      = np.fromfile(f, dtype=np.float32, count=1)
    self.nrk,self.neig = np.fromfile(f, dtype=np.int32, count=2)
    rec      = np.fromfile(f, dtype=np.float32, count=1)




### Insert code here

    rec    = np.fromfile(f, dtype=np.float32, count=1)
    self.xk     = np.fromfile(f, dtype=np.float32, count=self.nrk)
    rec    = np.fromfile(f, dtype=np.float32, count=1)

    rec    = np.fromfile(f, dtype=np.float32, count=1)
    self.bands    = np.fromfile(f, dtype=np.float32,count=self.neig*self.nrk).reshape(self.nrk,self.neig)
    rec    = np.fromfile(f, dtype=np.float32, count=1)

    rec    = np.fromfile(f, dtype=np.float32, count=1)
    self.pkn    = np.fromfile(f, dtype=np.float32,count=self.neig*self.nrk).reshape(self.neig,self.nrk)
    rec    = np.fromfile(f, dtype=np.float32, count=1)

    print(self.title)
    print(self.subtitle)
    print(self.eref)
    print(self.xmin,self.xmax)
    print(self.ymin,self.ymax)
    print(self.neig,self.nrk)
    print(self.xk)

    # BEGIN FILE READ orbital information

    rec   = np.fromfile(f, dtype=np.float32, count=1)
    self.nbaslcao = np.fromfile(f, dtype=np.int32, count=1)[0]
    rec   = np.fromfile(f, dtype=np.float32, count=1)

    print (self.nbaslcao)

    rec   = np.fromfile(f, dtype=np.float32, count=1)
    self.infolcao = np.fromfile(f, dtype=np.int32, count=5*self.nbaslcao).reshape(self.nbaslcao,5)
    rec   = np.fromfile(f, dtype=np.float32, count=1)

    rec   = np.fromfile(f, dtype=np.float32,  count=1)
    self.basxpsi = np.fromfile(f, dtype=np.float32,count=self.nbaslcao*self.neig*self.nrk).reshape(self.nrk,self.neig,self.nbaslcao)
    rec   = np.fromfile(f, dtype=np.float32,  count=1)

    rec   = np.fromfile(f, dtype=np.float32,  count=1)
    self.ntype = np.fromfile(f, dtype=np.int32, count=1)[0]
    rec   = np.fromfile(f, dtype=np.float32,  count=1)
    
    print ("ntype=", self.ntype)
    
    self.nameats    = np.empty(self.ntype, dtype="S2")

    rec     = np.fromfile(f, dtype=np.float32, count=1)
    for i in range(self.ntype):
      self.nameats[i]   = np.fromfile(f,np.uint8,2).tobytes()
    rec     = np.fromfile(f, dtype=np.float32, count=1)

    self.nameat = []

    for i in range(self.ntype):
      self.nameat.append(self.nameats[i].decode("UTF-8").rstrip().lstrip())
    
    print ("nameat=")
    print (self.nameat)
    

    for i in range(self.label.size):
      if (self.labels[i] == "Gamma" or self.labels[i] == "GAMMA" or self.labels[i] == "gamma" ):
        self.labels[i] = r"$\Gamma$"
      if (self.labels[i] == "Lambda" or self.labels[i] == "LAMBDA" or self.labels[i] == "Lambda" ):
        self.labels[i] = r"$\Lambda$"
      if (self.labels[i] == "Sigma" or self.labels[i] == "SIGMA" or self.labels[i] == "sigma" ):
        self.labels[i] = r"$\Sigma$"
      if (self.labels[i] == "Delta" or self.labels[i] == "DELTA" or self.labels[i] == "delta" ):
        self.labels[i] = r"$\Delta$"
      if (self.labels[i] == "--" or self.labels[i] == "--" or self.labels[i] == "--" ):
        self.labels[i] = " "

#   additional information

    rec    = np.fromfile(f, dtype=np.float32, count=1)
    self.rkpt     = np.fromfile(f, dtype=np.float64, count=3*self.nrk).reshape(self.nrk,3)
    rec    = np.fromfile(f, dtype=np.float32, count=1)

    rec    = np.fromfile(f, dtype=np.float32, count=1)
    self.rkpt_fld     = np.fromfile(f, dtype=np.float64, count=3*self.nrk).reshape(self.nrk,3)
    rec    = np.fromfile(f, dtype=np.float32, count=1)

    for irk in range(self.nrk):
      print("%5d %8.5f unf %8.5lf %8.5lf %8.5lf %8.5lf %8.5lf %8.5lf " %( 
      irk,self.xk[irk], 
      self.rkpt[irk,0], 
      self.rkpt[irk,1], 
      self.rkpt[irk,2], 
      self.rkpt_fld[irk,0],
      self.rkpt_fld[irk,1],
      self.rkpt_fld[irk,2]
      )
      )


