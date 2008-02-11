import subprocess,os,os.path,re
from numpy import *

class Cp2kEnv:
    def __init__(self,cpshell,env_id):
        self.cpshell=cpshell
        self.env_id=env_id
    def eval_e(self):
        'evaluates the energy in background'
        nerr=self.cpshell.nErrors
        self.cpshell.exec_cmd('eval_e %d\n'%(self.env_id),async=1)
        self.cpshell.raiseLastError(nerr)
    def eval_ef(self):
        'evaluates the energy and forces in background'
        nerr=self.cpshell.nErrors
        self.cpshell.exec_cmd('eval_ef %d\n'%(self.env_id),async=1)
        self.cpshell.raiseLastError(nerr)
    def get_pos(self):
        'returns the actual position of the atoms'
        nerr=self.cpshell.nErrors
        lines=self.cpshell.exec_cmd('get_pos %d\n'%(self.env_id))
        self.cpshell.raiseLastError(nerr)
        ndim=int(lines[0])
        pos=array(map(float,("".join(lines[1:])).split()[:ndim]))
        return reshape(pos,(ndim/3,3))
    def get_f(self):
        'returns the force of the atoms'
        nerr=self.cpshell.nErrors
        lines=self.cpshell.exec_cmd('get_f %d\n'%(self.env_id))
        self.cpshell.raiseLastError(nerr)
        ndim=int(lines[0])
        f=array(map(float,("".join(lines[1:])).split()[:ndim]))
        return reshape(f,(ndim/3,3))
    def get_natom(self):
        '''returns the number of atoms in the environment'''
        nerr=self.cpshell.nErrors
        lines=self.cpshell.exec_cmd('natom %d\n'%(self.env_id))
        self.cpshell.raiseLastError(nerr)
        return int(lines[0])
    def get_e(self):
        '''returns the energy, of the last configuaration evaluated,
         0.0 if none was evaluated'''
        nerr=self.cpshell.nErrors
        lines=self.cpshell.exec_cmd('get_e %d\n'%(self.env_id))
        self.cpshell.raiseLastError(nerr)
        return float(lines[0])
    def set_pos(self,pos):
        """changes the position of the environement
        returns the maximum change in coordinates"""
        nerr=self.cpshell.nErrors
        fin=self.cpshell.st.stdin
        fin.write('set_pos %d\n'%(self.env_id))
        p=ravel(pos)
        fin.write("%d\n"%(p.size))
        for pp in p:
            fin.write(' %18.13f'%(pp))
        fin.write("\n* END\n")
        fin.flush()
        lines=self.cpshell.getReady()
        self.cpshell.raiseLastError(nerr)
        return float(lines[0])
    def destroy(self):
        "destroys this environment"
        nerr=self.cpshell.nErrors
        lines=self.cpshell.exec_cmd('destroy %d\n'%(self.env_id))
        self.cpshell.envs[self.env_id]
        self.env_id=-1
        self.cpshell.raiseLastError(nerr)

class Cp2kShell:
    def __init__(self,task,taskDir):
    	"starts a cp2k task"
    	myDir=os.getcwd()
    	os.chdir(taskDir)
    	self.st=subprocess.Popen(task,shell=True,
    		stdin=subprocess.PIPE,
    		stdout=subprocess.PIPE,
    		stderr=subprocess.STDOUT,
    		close_fds=1, env=os.environ)
    	self.isReady=0
    	self.envs={}
    	self.nErrors=0
    	self.lastErrors=[]
    	os.chdir(myDir)
    
    def raiseLastError(self,nErrors=0):
        if nErrors!=self.nErrors:
            raise Exception("".join(self.lastErrors))
    
    def addErrs(self,errs):
        'adds the given errors to the task errors'
        if errs:
            self.nErrors+=1
            self.lastErrors=errs
    
    def clearErrs(self):
        self.nErrors=0
        self.lastErrors=[]
    
    def getReady(self):
        "reads the readiness message of the subtask"
        readyRe=re.compile(r"^ *\* *(?:ready).*",re.IGNORECASE)
        errorRe=re.compile(r"^ *\* *error.*",re.IGNORECASE)
        lines=[]
        errors=[]
        while 1:
            line=self.st.stdout.readline()
            if not line: break
            if errorRe.match(line):
                errors.append(line)
            if readyRe.match(line):
                self.isReady=1
                self.addErrs(errors)
                return lines
            lines.append(line)
        self.addErrs(errors)        
        raise Exception('cp2k_shell did not get ready')
    
    def exec_cmd(self,cmdStr,async=0):
        "executes a command and (if async) returns the output"
        if not self.isReady:
            lines=self.getReady()
        self.st.stdin.write(cmdStr)
        if cmdStr[-1]!='\n':
            self.st.stdin.write('\n')            
        self.st.stdin.flush()
        if async:
            self.isReady=0
        else:
            return self.getReady()
    
    def newEnv(self,inputFile,runDir=""):
        """creates a new environment
        changes the current working directory to it"""
        nerr=self.nErrors
        if runDir and runDir!=".":
            self.exec_cmd('cd '+runDir)
        if self.nErrors==nerr:
            lines=self.exec_cmd('load '+inputFile)
            if self.nErrors==nerr:
                env_id=int(lines[0])
                self.envs[env_id]=Cp2kEnv(self,env_id)
                return self.envs[env_id]
        return None
    
    def newBgEnv(self,inputFile,runDir=""):
        "creates a new environment in background, use lastEnv to get it"
        nerr=self.nErrors
        if runDir and runDir!=".":
            self.exec_cmd('cd '+runDir)
        if self.nErrors==nerr:
            self.exec_cmd('bg_load '+inputFile,async=1)
    
    def lastEnv(self):
        "returns the last environment created"
        lines=self.exec_cmd('last_env_id')
        env_id=int(lines[0])
        if env_id<=0:
            return None
        if not self.envs.has_key(env_id):
            self.envs[env_id]=Cp2kEnv(self,env_id)
        return self.envs[env_id]
    
    def chdir(self,newDir):
        "changes the working directory"
        self.exec_cmd('cd '+newDir)
    
    def getcwd(self):
        "returns the working directory"
        return self.exec_cmd('pwd')[0][:-1]
        
        