import cp2k_interface_low,os,copy

cp2kDefaultBase=os.getenv("HOME",os.getcwd())+"/cp2k"
defaultRunBaseDir="/tmp"

def calc_default_input(cp2kBase=cp2kDefaultBase):
    default_input={
        'CP2K':{'PROGRAM':'Quickstep', 'PROJECT':'Default', 'IOLEVEL':'10',
                'FFTLIB':'FFTSG', 'RUN_TYPE':'GEO_OPT',
                'BASIS_SET_FILE_NAME':cp2kBase+"/tests/QS/GTH_BASIS_SET",
                'POTENTIAL_FILE_NAME':cp2kBase+"/tests/QS/POTENTIAL" },
        'DFT':{'FUNCTIONAL':'Pade', 'XC_SMOOTH_RHO':'NN50','XC_DERIV':'NN50_SMOOTH'},
        'QS':{'CUTOFF':'200'},
        'SCF':{'GUESS':'ATOMIC'},
        'PRINT MEDIUM':{}
        'force_eval':{
         'dft':{
          'BASIS_SET_FILE_NAME':cp2kBase+"/tests/QS/GTH_BASIS_SET",
          'POTENTIAL_FILE_NAME':cp2kBase+"/tests/QS/POTENTIAL",
         }
        }
        }
    return default_input

def completeInput(input,defaults):
    """Completes missing values from input from the dictionary defaults"""
    for k,v in defaults.items():
        if not input.has_key(k):
            input[k]=copy.deepcopy(v)
        elif isinstance(input[k],dict):
            completeInput(input[k],v)

def writeInput(input,outfile,intent=0):
    """writes out the dictionary input in the file outFile using the
    cp2k input format"""
    outfile.write(indent*" ")
    for k,v in input.items():
        if isinstance(v,dict):
            outfile.write("&"+k+"\n")
            writeInput(v,outFile,indent=indent+1)
            outfile.write((indent*" ")+"&END "+k+"\n")
        elif isinstance(v,list):
            outfile.write(k)
            for i in range(len(v)):
                el=v[i]
                if isinstance(el,list):
                    for el2 in el:
                        outfile.write(" "+str(el2))
                    if i!=len(v)-1:
                        outfile.write("\n")
                else:
                    outfile.write(" "+str(el))
            outfile.write("\n")
        else:
            outfile.write(k+" "+str(v))

def initialSetup():
    cp2k_interface_low.cp_init_cp2k(1)

def finalCleanup():
    cp2k_interface.cp_finalize_cp2k(1)

def createPath(p):
    "creates the path if it does not exist"
    if not os.path.exits(p):
        s=os.path.split(p)
        if s[0]:
            createPath(s[0])
        if s[1]:
            os.mkdir(p)

def readInput(inFile):
    "reads the input in cp2k format from the file inFile"
    root={'__name__':'__root__'}
    path=[root]
    while 1:
        line=inFile.readline()
        l=line.split()
        if l:
            if l[0].lower()=='&end':
                
            if l[0][0]=='&':
                sect={'__name__':l[0][1:],
                      '__args__':l[1:]}
                path.append(sect)
            

def setInitialPos(inputDict,atoms,pos):
    "set the initial position in the given input"
    
    inputDict['force_eval']
    
class cp2k_env:
    def __init__(self,input,runDir=None,from_scratch=1):
        self.input=input
        self.runDir=runDir
        self.outFilePath
        self.from_scratch=from_scratch
        self.force_env_id=None
        self.atom_labels=None
        self.v=None
        self.f=None
        if self.force_env_id!=None:
            throw Error("start called twice on the same object")
        if not self.runDir:
            self.runDir=defaultRunBaseDir+"/cp2kRun"
        if self.from_scratch:
            rDir=self.runDir
            i=0
            while os.path.exits(rDir):
                i+=1
                rDir=self.runDir+str(i)
            self.runDir=rDir
        self.runDir=os.path.normpath(os.path.abspath(self.runDir))
        createPath(self.runDir)
        self.inputPath=os.path.join(self.runDir,"input.inp")
        iFile=file(self.inputPath,"w")
        writeInput(iFile)
        iFile.close()
        self.outputPath=os.path.join(self.runDir,"output.out")
        self.force_env_id=self.cp2k_interface_low.cp_create_f_env(self.inputPath,self.outputPath)
    def set_pos(newPos):
        self.pos=newPos
        if 
        
        
initialSetup()
