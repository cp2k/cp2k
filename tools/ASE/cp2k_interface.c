/*-----------------------------------------------------------------------------!
 *   CP2K: A general program to perform molecular dynamics simulations         !
 *   Copyright (C) 2002,2003,2004  CP2K developers group                       !
 *-----------------------------------------------------------------------------!
 *
 * python bridge of the f77 interface of cp2k
 *
 * 8.2004 created [Johanna & fawzi]
 * 11.2004 adapted to the new f77 interface, better error checking
 * 
 */


#include <Python.h>
#include "numarray.h"
#include "libnumeric.h"
#include "cp2k_c_bridge.h"

#define tFReal tFloat64

static PyObject *cp2kError;
/* python interface functions */

static PyObject* cp2k_interface_cp_init_cp2k(PyObject *self, PyObject *args){ 
  int ierr,init_mpi;
  if (!PyArg_ParseTuple(args, "i",&init_mpi))
    return NULL;
      
  ierr=ccp_init_cp2k(init_mpi);
  if (ierr) {
    PyErr_Format(cp2kError,"error %d initializing cp2k",ierr);
    return NULL;
  }
  return Py_BuildValue("");
}

static PyObject* cp2k_interface_cp_create_f_env(PyObject *self, PyObject *args){
  int ierr,new_id;
  const char *in_path, *out_path;

  if (!PyArg_ParseTuple(args, "ss",&in_path,&out_path))
    return NULL;
  ierr=ccp_create_fenv(&new_id,in_path,out_path);
  if (ierr) {
    PyErr_Format(cp2kError,"error %d creating a new force env",ierr);
    return NULL;
  }
  return Py_BuildValue("i", new_id);
}

static PyObject* cp2k_interface_cp_set_pos(PyObject *self, PyObject *args){
  PyObject* positions;
  int env_id,dim,ierr;
  PyArrayObject *c_positions;

  if (!PyArg_ParseTuple(args, "iO", &env_id, &positions))
    return NULL;
  
  c_positions = NA_InputArray(positions,tFReal,NUM_C_ARRAY);
  if (!c_positions) {
    PyErr_SetString(cp2kError,"error converting the input array");
    return NULL;
  }
  if (c_positions->nd!=2) return NULL;
    
  dim = c_positions->dimensions[0]*c_positions->dimensions[1];

  ierr=ccp_set_pos(env_id,NA_OFFSETDATA(c_positions), dim);

  Py_XDECREF(c_positions);

  if (ierr) {
    PyErr_SetString(cp2kError,"error setting the atomic positions");
    return NULL;
  }
  return Py_BuildValue("");
}

static PyObject* cp2k_interface_cp_get_pos(PyObject *self, PyObject *args){
  PyObject* positions;
  int env_id,dim,ierr;
  PyArrayObject *c_positions;
  f_real *array;

  if (!PyArg_ParseTuple(args, "i|O", &env_id, &positions))
    return NULL;

  /* Just to get dimension */
  if (!positions){
    ierr=ccp_get_natom(env_id,&dim);
    dim=dim*3;
    array=malloc(sizeof(f_real)*dim);
    if (!array) {
      PyErr_NoMemory();
      return NULL;
    }
  }
  else{
    c_positions = NA_OutputArray(positions,tFReal,NUM_C_ARRAY);
    dim = c_positions->dimensions[0]*c_positions->dimensions[1];
    array=NA_OFFSETDATA(c_positions);
  }

  /* Get the positions and put them in c_positions */
  ierr=ccp_get_pos(env_id, array, dim);
  if (ierr) {
    PyErr_SetString(cp2kError,"error getting the atomic positions");
    return NULL;
  }

  if (!positions){
    c_positions=NA_NewArray(array,tFReal,2,dim/3,3); /* swap dimensions? */
    free(array);
    return PyArray_Return(c_positions);
  }
  else {
    Py_XDECREF(c_positions);
    return positions;
  }
}

static PyObject* cp2k_interface_cp_get_force(PyObject *self, PyObject *args){
  PyObject* positions;
  int env_id,dim,ierr;
  PyArrayObject *c_forces;

/*  if (!PyArg_ParseTuple(args, "iO", &env_id, &forces))
    return NULL;

  / * Just to get dimension * /
  c_forces = NA_InputArray(forces,tFReal,NUM_C_ARRAY);
  dim = c_forces->dimensions[0]*c_forces->dimensions[1];

  / * Get the force and put them in c_force * /
  ierr=ccp_get_force(env_id, NA_OFFSETDATA(c_forces), dim);
  if (ierr) {
    PyErr_SetString(cp2kError,"error getting the forces");
    return NULL;
  }

  / * Make a suitable python array * /
  return NA_OutputArray(c_forces, tFReal, NUM_C_ARRAY);*/
  return NULL;
}

static PyObject* cp2k_interface_cp_get_energy(PyObject *self, PyObject *args){
  int env_id,ierr;
  f_real e_pot;

  if (!PyArg_ParseTuple(args, "i", &env_id))
    return NULL;

  ierr=ccp_get_energy(env_id,&e_pot);
  if (ierr) {
    PyErr_SetString(cp2kError,"error getting the potential energy");
    return NULL;
  }
  if (ierr) return NULL;
  return Py_BuildValue("d",e_pot);
}

static PyObject* cp2k_interface_cp_update_forces(PyObject *self, PyObject *args){
  int env_id,ierr,calcForce;

  if (!PyArg_ParseTuple(args, "ii", &env_id,&calcForce))
    return NULL;

  ierr=ccp_calc_energy_force(env_id,calcForce);
  if (ierr) {
    PyErr_SetString(cp2kError,"error updating the forces");
    return NULL;
  }
  return Py_BuildValue("");
}


/* python extension stuff */
static PyMethodDef cp2k_interfaceMethods[] = {
  {"cp_init_cp2k", cp2k_interface_cp_init_cp2k, METH_VARARGS, "Initialize CP2K."},
  {"cp_create_f_env", cp2k_interface_cp_create_f_env, METH_VARARGS, "Create a Force Environment"},
  {"cp_set_pos", cp2k_interface_cp_set_pos, METH_VARARGS, "Set Positions"},
  {"cp_get_pos", cp2k_interface_cp_get_pos, METH_VARARGS, "Get Positions"},
  {"cp_get_force", cp2k_interface_cp_get_force, METH_VARARGS, "Get Forces"},
  {"cp_get_energy", cp2k_interface_cp_get_energy, METH_VARARGS, "Get Positions"},
  {"cp_update_forces", cp2k_interface_cp_update_forces, METH_VARARGS, "Update Forces"},
  {NULL, NULL, 0, NULL} /* Sentinel */
};


PyMODINIT_FUNC initcp2k_interface(void){
  PyObject *m;

  m=Py_InitModule("cp2k_interface", cp2k_interfaceMethods);
  import_libnumarray();

  cp2kError = PyErr_NewException("cp2k_interface_low.error", NULL, NULL);
  Py_INCREF(cp2kError);
  PyModule_AddObject(m, "error", cp2kError);
  if (tFReal==tFloat64 && sizeof(Float64)!=sizeof(f_real)) {
    PyErr_Format(cp2kError,"*ERROR* fortran real (%d) is not c double (%d)",
		 (int)sizeof(f_real),(int)sizeof(Float64));
  }
}

