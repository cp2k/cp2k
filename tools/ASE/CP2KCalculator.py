"""Simple version of CP2K calculator."""

import weakref
import Numeric as num

from ASE.Calculators.Calculator import Calculator
import cp2k_interface
import numarray as numa


class CP2KCalculator(Calculator):
    #To determine whether or not CP2K was initialized
    cp2kWasInit = 0
    
    #kindParams is a dictionary of the settings to be changed
    kindParams = {}
    
    #Default dictionary of parameters to build input file for CP2K    
    params = {'CP2K':{'PROGRAM':'Quickstep', 'PROJECT':'Default', 'IOLEVEL':'10',
                      'FFTLIB':'FFTSG', 'RUN_TYPE':'GEO_OPT'},
              'DFT':{'FUNCTIONAL':'Pade', 'XC_SMOOTH_RHO':'NN50','XC_DERIV':'NN50_SMOOTH'},
              'QS':{'CUTOFF':'200'},
              'SCF':{'GUESS':'ATOMIC'},
              'PRINT MEDIUM':{}}
    
    #Default atom paramters
    atomParams = []
    
    #Forces
    forces = []

    #Positions
    positions = []

    #ID for qs environment
    env_id = 0
    
    def __init__(self, inputParams):
        self.ready = False
        
        #Assign user defined params to class variable
        self.kindParams.update(inputParams)
                
                
    def _SetListOfAtoms(self, atoms):
        Calculator._SetListOfAtoms(self, atoms)

        #Integrate user defined params, and build other params
        self.UpdateParams()
        
        
        
    def GetCartesianForces(self):
        #Control from ASE API?

        #Make array of positions
        #Get number of atoms
        numAtoms = len(self.GetListOfAtoms())

        #Debugging
        print 'there are %i atoms' %(numAtoms)

        #Make a numarray to hold the atoms
        positions = numa.zeros((numAtoms, 3,), numa.Float64)

        #Get the info from each atom and copy it to the numarray
        count = 0
        tempList = self.GetListOfAtoms()
        for k in tempList:
            positions[count] = k.GetCartesianPosition()
            count = count+1

        #Debugging
        print positions    

        #Run a set positions
        cp2k_interface.c_set_pos(positions, self.env_id) 
      
        #Run a get to check
        positions = numa.zeros((numAtoms, 3,), numa.Float64)
        cp2k_interface.c_get_pos(positions, self.env_id)

        #Debugging
        print positions
        
        #Calculate forces (update_forces)
        cp2k_interface.c_update_forces(self.env_id)       

        #Now run a get force


        #Use parse tuple to get pyobject

        #Junk return
        self.forces = num.zeros((numAtoms, 3), num.Float)
        return self.forces
    

        
    #Looks through ListOfAtoms to generate new parameters
    def UpdateParams(self):
        list = self.GetListOfAtoms()
        cell = list.GetUnitCell()

        #Add cell info to dictionary
        cellParams = {'CELL':{'UNIT':'ANGSTROM', 'ABC':'%s %s %s'
                              %(cell[0][0], cell[1][1], cell[2][2])}}
        
        self.params.update(cellParams)

        #Add user specified kinds to the dictionary
        self.params.update(self.kindParams)
        
        #Add coordinate info to list
        for k in list:
            self.atomParams.append([k.GetChemicalSymbol(), k.GetCartesianPosition()[0],
                                    k.GetCartesianPosition()[1], k.GetCartesianPosition()[2]])

        #Do some error checking on the kinds
        for k in self.atomParams:
            flag = 0
            if self.params.has_key('KIND '+self.atomParams[self.atomParams.index(k)][0]):
                flag = 1
            if flag == 0:
                print 'ERROR: no kind specified for',self.atomParams[self.atomParams.index(k)][0]
                break
            
           
        #Print params to screen    
        self.PrintParams()

        #Print params to file
        self.PrintParamsToFile()

        #Init CP2K if it's the first run
        #Create a QS Environment if it hasn't been done
        if self.cp2kWasInit == 0:
            cp2k_interface.c_init_cp2k()
            print 'cp2k initialized'
            self.env_id = cp2k_interface.c_create_qs_env()
            print 'envid %i' %(self.env_id)
            self.cp2kWasInit = 1
            


    #Prints all of the parameters out
    def PrintParams(self):
        #First print params list
        for k in self.params.keys():
            print "&"+k
            for j in self.params[k].keys():
                print ' %s %s' %(j, self.params[k][j])
            print "&END\n"


        print "&COORD"
        #Then print atomParams list    
        for k in self.atomParams:
            print '%s %s %s %s' %(self.atomParams[self.atomParams.index(k)][0],
                                  self.atomParams[self.atomParams.index(k)][1],
                                  self.atomParams[self.atomParams.index(k)][2],
                                  self.atomParams[self.atomParams.index(k)][3])
        print "&END"        
                
                
    #Prints all of the parameters out
    def PrintParamsToFile(self):
        f=file('input.inp', 'w')
        
        #First print params list
        for k in self.params.keys():
            f.write("&"+k+"\n")
            for j in self.params[k].keys():
                f.write(' %s %s\n' %(j, self.params[k][j]))
            f.write("&END\n\n")


        f.write("&COORD\n")
        #Then print atomParams list    
        for k in self.atomParams:
            f.write('%s %s %s %s\n' %(self.atomParams[self.atomParams.index(k)][0],
                                  self.atomParams[self.atomParams.index(k)][1],
                                  self.atomParams[self.atomParams.index(k)][2],
                                  self.atomParams[self.atomParams.index(k)][3]))
        f.write("&END\n")   
