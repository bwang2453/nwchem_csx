/*
 $Id: nwchem_wrap.c 26319 2014-10-14 22:54:32Z edo $
*/
#include <Python.h>
#if defined(DECOSF)
#include <alpha/varargs.h>
#endif
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#include <abstract.h>

#include "rtdb.h"
#include "macdecls.h"
#include "typesf2c.h"

#if PY_MAJOR_VERSION >= 3
#define PyString_FromString PyBytes_FromString
#define PyString_Check      PyBytes_Check
#define PyInt_FromLong      PyLong_FromLong
#define PyInt_Check         PyLong_Check
#define PyInt_AS_LONG       PyLong_AS_LONG
#endif

static PyObject *NwchemError;

static Integer rtdb_handle;            /* handle to the rtdb */
extern void task_(Integer *);


static PyObject *nwwrap_integers(int n, Integer a[])
{
  PyObject *sObj;
  int i;

  if (n == 1)
    return PyInt_FromLong(a[0]);

  if (!(sObj=PyList_New(n))) return NULL;
  for(i=0; i<n; i++) {
    PyObject *oObj = PyInt_FromLong(a[i]);
    if (!oObj) {
      Py_DECREF(sObj);
      return NULL;
    }
    if (PyList_SetItem(sObj,i,oObj)) {
      Py_DECREF(sObj);
      Py_DECREF(oObj);
      return NULL;
    }
  }
  return sObj;
}

static PyObject *nwwrap_doubles(int n, double a[])
{
  PyObject *sObj;
  int i;

  if (n == 1)
    return PyFloat_FromDouble(a[0]);

  if (!(sObj=PyList_New(n))) return NULL;
  for(i=0; i<n; i++) {
    PyObject *oObj = PyFloat_FromDouble(a[i]);
    if (!oObj) {
      Py_DECREF(sObj);
      return NULL;
    }
    if (PyList_SetItem(sObj,i,oObj)) {
      Py_DECREF(sObj);
      Py_DECREF(oObj);
      return NULL;
    }
  }
  return sObj;
}

static PyObject *nwwrap_strings(int n, char *a[])
{
  PyObject *sObj;
  int i;

  if (n == 1)
    return PyString_FromString(a[0]);

  if (!(sObj=PyList_New(n))) return NULL;
  for(i=0; i<n; i++) {
    PyObject *oObj = PyString_FromString(a[i]);
    if (!oObj) {
      Py_DECREF(sObj);
      return NULL;
    }
    if (PyList_SetItem(sObj,i,oObj)) {
      Py_DECREF(sObj);
      Py_DECREF(oObj);
      return NULL;
    }
  }
  return sObj;
}


static PyObject *wrap_rtdb_open(PyObject *self, PyObject *args)
{
   const char *filename, *mode;
   int inthandle;

   if (PyArg_ParseTuple(args, "ss", &filename, &mode)) {
       if (!rtdb_seq_open(filename, mode, &inthandle)) {
           PyErr_SetString(NwchemError, "rtdb_seq_open failed");
           return NULL;
       }
       rtdb_handle = inthandle;
   }
   else {
      PyErr_SetString(PyExc_TypeError, "Usage: rtdb_open(filename, mode)");
      return NULL;
   }
  (void) MA_init(MT_CHAR, -1, -1);
  (void) MA_set_auto_verify(1);
   Py_INCREF(Py_None);
   return Py_None;
}

static PyObject *wrap_rtdb_close(PyObject *self, PyObject *args)
{
   const char *mode;
   int  result;

   if (PyArg_ParseTuple(args, "s", &mode)) {
       if (!(result = rtdb_seq_close(rtdb_handle, mode))) {
           PyErr_SetString(NwchemError, "rtdb_close failed");
           return NULL;
       }
   }
   else {
       PyErr_SetString(PyExc_TypeError, "Usage: rtdb_close(mode)");
       return NULL;
   }
   Py_INCREF(Py_None);
   return Py_None;
}


static PyObject *wrap_pass_handle(PyObject *self, PyObject *args)
{
  int inthandle;
   if (!(PyArg_ParseTuple(args, "i", &inthandle))) {
      PyErr_SetString(PyExc_TypeError, "Usage: pass_handle(rtdb_handle)");
      return NULL;
   }
   rtdb_handle = inthandle;
   Py_INCREF(Py_None);
   return Py_None;
}


static PyObject *wrap_rtdb_print(PyObject *self, PyObject *args)
{
   int flag;

//   if (!rtdb_seq_print(rtdb_handle, 0)) 
//        PyErr_SetString(NwchemError, "rtdb_print failed");
   if (PyArg_ParseTuple(args, "i", &flag)) {
      if (!rtdb_seq_print(rtdb_handle, flag)) 
           PyErr_SetString(NwchemError, "rtdb_print failed");
   }
   else {
      PyErr_SetString(PyExc_TypeError, "Usage: rtdb_print(flag)");
      return NULL;
   }
   Py_INCREF(Py_None);
   return Py_None;
//   Py_RETURN_NONE
}


static PyObject *wrap_rtdb_put(PyObject *self, PyObject *args)
{
    int i, list, list_len;
    int ma_type = -1;
    char *name;
    Integer* int_array;
    double *dbl_array;
    char *char_array;
    char cbuf[8192], *ptr;
    void *array = 0;
    PyObject *obj;

     if ((PyTuple_Size(args) == 2) ) {
      PyArg_ParseTuple(args, "sO", &name,&obj);
      

      if (PyList_Check(obj)) 
        list = 1; 
      else 
        list = 0;
      
      if (list) {
        list_len = PyList_Size(obj);
        if (   PyInt_Check(PyList_GetItem(obj, 0)))  
          ma_type = MT_F_INT;
        else if ( PyFloat_Check(PyList_GetItem(obj, 0)))  
          ma_type = MT_F_DBL;
        else if (PyString_Check(PyList_GetItem(obj, 0))) 
          ma_type = MT_CHAR;
        else {
          printf("ERROR A\n");
          ma_type = -1;
        }
      } else {
        list_len = 1;
        if (   PyInt_Check(obj))  
          ma_type = MT_F_INT;
        else if ( PyFloat_Check(obj))  
          ma_type = MT_F_DBL;
        else if (PyString_Check(obj))  
          ma_type = MT_CHAR; 
        else {
          printf("ERROR B\n");
          ma_type = -1;
        }
      } 

      if (ma_type == -1) {
          PyErr_SetString(PyExc_TypeError, 
                          "Usage: rtdb_put - ma_type is confused");
          return NULL;
      }
      
     
      if (ma_type != MT_CHAR) {
        if (!(array = malloc(MA_sizeof(ma_type, list_len, MT_CHAR)))) {
          PyErr_SetString(PyExc_MemoryError,
                          "rtdb_put failed allocating work array");
          return NULL;
        }
      }
      
      switch (ma_type) {
      case MT_INT:
      case MT_F_INT:  
      case MT_BASE + 11:        /* Logical */
        int_array = array;
        for (i = 0; i < list_len; i++) {
          if (list) 
            int_array[i] = PyInt_AS_LONG(PyList_GetItem(obj, i));
          else 
            int_array[i] = PyInt_AS_LONG(obj);
        }
        break;
        
      case MT_DBL:  
      case MT_F_DBL:
        dbl_array = array;
        for (i = 0; i < list_len; i++) {
          if (list) 
            PyArg_ParseTuple(PyList_GetItem(obj, i), "d", dbl_array+i);
          else 
            PyArg_ParseTuple(obj, "d", dbl_array+i);
        }
        break;
        
      case MT_CHAR: 
        ptr = cbuf;
        *ptr = 0;
        for (i = 0; i < list_len; i++) {
          if (list) 
            PyArg_ParseTuple(PyList_GetItem(obj, i), "s", &char_array); 
          else 
            PyArg_ParseTuple(obj, "s", &char_array); 
          /*printf("PROCESSED '%s'\n", char_array);*/
          if ((ptr+strlen(char_array)) >= (cbuf+sizeof(cbuf))) {
             PyErr_SetString(PyExc_MemoryError,"rtdb_put too many strings");
             return NULL;
           }
          strcpy(ptr,char_array);
          ptr = ptr+strlen(char_array);
          strcpy(ptr,"\n");
          ptr = ptr + 1;
        }                
        list_len = strlen(cbuf) + 1;
        array = cbuf;
        break;
        
      default:
        PyErr_SetString(NwchemError, "rtdb_put: ma_type is incorrect");
        if (array) free(array);
        return NULL;
        break;
      }                
      
      if (!(rtdb_seq_put(rtdb_handle, name, ma_type, list_len, array))) {
        PyErr_SetString(NwchemError, "rtdb_seq_put failed");
        if ((ma_type != MT_CHAR) && array) free(array);
        return NULL;
      }
      
    } else {
      PyErr_SetString(PyExc_TypeError, 
                      "Usage: rtdb_put(value or values,[optional type])");
      if ((ma_type != MT_CHAR) && array) free(array);
      return NULL;
    }
    Py_INCREF(Py_None);
    if ((ma_type != MT_CHAR) && array) free(array);
    return Py_None;
   Py_RETURN_NONE;
}

PyObject *test_parse(PyObject *self, PyObject *args)
{
  char *name;
  if (PyArg_ParseTuple(args, "s", &name)) {
      printf("it worked \n");
      printf("%s + \n",name);
      }
  else {
      printf("it did not \n");
  }
   Py_RETURN_NONE;

}
PyObject *wrap_rtdb_get(PyObject *self, PyObject *args)
{
  int nelem, ma_type;
  char *name;
#define MAXPTRS 2048
  char *ptrs[MAXPTRS];
  PyObject *returnObj = NULL;
  char format_char;
  void *array=NULL;
  int ma_handle;
  
  if (PyArg_ParseTuple(args, "s", &name)) {
    if (!rtdb_seq_ma_get(rtdb_handle, name, &ma_type, &nelem, &ma_handle)) {
      PyErr_SetString(NwchemError, "rtdb_ma_get failed");
      return NULL;
    }
    if (!MA_get_pointer(ma_handle, &array)) {
      PyErr_SetString(NwchemError, "rtdb_ma_get failed");
      return NULL;
    }
//    printf("name=%s ma_type=%d nelem=%d ptr=%x\n",name, ma_type, 
//      nelem, array);
    
    switch (ma_type) {
    case MT_F_INT:
    case MT_INT  : 
    case MT_BASE + 11  : 
      format_char = 'i'; break;
    case MT_F_DBL: 
    case MT_DBL  : 
      format_char = 'd'; break;
      break;
    case MT_CHAR : 
      format_char = 's'; break;
    default:
      PyErr_SetString(NwchemError, "rtdb_get: ma type incorrect");
      (void) MA_free_heap(ma_handle);
      return NULL;
      break;
    }
    
    /* For character string need to build an array of pointers */
    
    if (ma_type == MT_CHAR) {
      char *ptr, *next;
      nelem = 0;
      next = ptr = array;
      while (1) {
        int eos = (*ptr == 0);
        if ((*ptr == '\n') || (*ptr == 0)) {
          *ptr = 0;
          if (nelem >= MAXPTRS) {
            PyErr_SetString(PyExc_MemoryError,"rtdb_get too many strings");
            (void) MA_free_heap(ma_handle);
            return NULL;
          }
          if (strlen(next) > 0) {
            ptrs[nelem] = next;
            nelem++;
          }
          next = ptr+1;
          if (eos) break;
        }
        ptr++;
      }
    }
           
    switch (format_char) {
    case 'i':
      returnObj = nwwrap_integers(nelem, array); break;
    case 'd':
      returnObj = nwwrap_doubles(nelem, array); break;
    case 's':
      returnObj = nwwrap_strings(nelem, ptrs); break;
    }

    (void) MA_free_heap(ma_handle);

    if (!returnObj) {
      PyErr_SetString(PyExc_TypeError, "rtdb_get: conversion to python object failed.");
    }
  }
  else {
    PyErr_SetString(PyExc_TypeError, "Usage: value = rtdb_get(name)");
  }

  return returnObj;
}

PyObject *wrap_rtdb_delete(PyObject *self, PyObject *args)
{
   char *name;
   PyObject *returnObj = NULL;

   if (PyArg_Parse(args, "s", &name)) {
       if (rtdb_seq_delete(rtdb_handle, name)) {
         returnObj = Py_None;
         Py_INCREF(Py_None);
       }
       else {
           PyErr_SetString(NwchemError, "rtdb_delete failed");
       }
   }
   else {
       PyErr_SetString(PyExc_TypeError, "Usage: value = rtdb_delete(name)");
   }
   return returnObj;
}

PyObject *wrap_rtdb_get_info(PyObject *self, PyObject *args)
{
   int nelem, ma_type;
   char *name;
   PyObject *returnObj = 0;
   char date[26];

   if (PyArg_Parse(args, "s", &name)) {
       if (!rtdb_seq_get_info(rtdb_handle, name, &ma_type, &nelem, date)) {
           PyErr_SetString(NwchemError, "rtdb_seq_get_info failed");
           return NULL;
       }
       if (!(returnObj = PyTuple_New(3))) {
           PyErr_SetString(NwchemError, "rtdb_seq_get_info failed with pyobj");
           return NULL;
       }
       PyTuple_SET_ITEM(returnObj, 0, PyInt_FromLong((long) ma_type)); 
       PyTuple_SET_ITEM(returnObj, 1, PyInt_FromLong((long) nelem)); 
       PyTuple_SET_ITEM(returnObj, 2, PyString_FromString(date)); 
   }
   else {
       PyErr_SetString(PyExc_TypeError, "Usage: value = rtdb_get_info(name)");
       return NULL;
   }
   return returnObj;
}


PyObject *wrap_rtdb_first(PyObject *self, PyObject *args)
{
   char name[256];
   PyObject *returnObj = NULL;

   if (rtdb_seq_first(rtdb_handle, sizeof(name), name)) {
     returnObj = PyString_FromString(name); /*Py_BuildValue("s#", name, 1); */
   }
   else {
       PyErr_SetString(NwchemError, "rtdb_first: failed");
       return NULL;
   }
   return returnObj;
}

PyObject *wrap_rtdb_next(PyObject *self, PyObject *args)
{
   char name[256];
   PyObject *returnObj = NULL;

   if (rtdb_seq_next(rtdb_handle, sizeof(name), name)) {
     returnObj = PyString_FromString(name); /*Py_BuildValue("s#", name, 1); */
   }
   else {
       PyErr_SetString(NwchemError, "rtdb_next: failed");
       return NULL;
   }
   return returnObj;
}




/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
static char second_doc[] = 
"This module is just a simple example.  It provides one function: func().";

static PyObject*
second_func(PyObject *self, PyObject *args)
{
	PyObject *a, *b;
	
	if (!PyArg_UnpackTuple(args, "func", 2, 2, &a, &b)) {
		return NULL;
	}

	return PyNumber_Add(a, b);
}

static char second_func_doc[] = 
"func(a, b)\n\
\n\
Return the sum of a and b.";

static PyMethodDef nwchem_methods[] = {
	{"func",          second_func,        METH_VARARGS, second_func_doc},
        {"rtdb_open",     wrap_rtdb_open,     METH_VARARGS, second_func_doc}, 
        {"rtdb_close",    wrap_rtdb_close,    METH_VARARGS, second_func_doc}, 
        {"rtdb_print",    wrap_rtdb_print,    METH_VARARGS, second_func_doc}, 
        {"pass_handle",   wrap_pass_handle,   METH_VARARGS, 0}, 
        {"rtdb_put",      wrap_rtdb_put,      METH_VARARGS, 0}, 
        {"rtdb_get",      wrap_rtdb_get,      METH_VARARGS, 0}, 
        {"rtdb_delete",   wrap_rtdb_delete,   METH_VARARGS, 0}, 
        {"rtdb_get_info", wrap_rtdb_get_info, METH_VARARGS, 0}, 
        {"rtdb_first",    wrap_rtdb_first,    METH_VARARGS, 0}, 
        {"rtdb_next",     wrap_rtdb_next,     METH_VARARGS, 0}, 
        {"test_parse",    test_parse,         METH_VARARGS, 0}, 
	{NULL, NULL}
};

#if PY_MAJOR_VERSION >= 3
static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,
        "nwchem",
        "module to access the NWChem Runtime DataBase (RTDB)",
        -1,
        nwchem_methods,
        NULL,
        NULL,
        NULL,
        NULL
};
#endif

#if PY_MAJOR_VERSION >= 3
PyMODINIT_FUNC
PyInit_nwchem(void)
{
        PyObject *m, *d;

        m =  PyModule_Create(&moduledef);

        d = PyModule_GetDict(m);

        NwchemError = PyErr_NewException("nwchem.NwchemError",NULL,d);
        Py_INCREF(NwchemError);
        PyModule_AddObject(m,"NwchemError",NwchemError);

        return m;
}
#else
PyMODINIT_FUNC
initnwchem(void)
{
        PyObject *m, *d;

	m = Py_InitModule3("nwchem", nwchem_methods, second_doc);

        d = PyModule_GetDict(m);

        NwchemError = PyErr_NewException("nwchem.NwchemError",NULL,d);
        Py_INCREF(NwchemError);
        PyModule_AddObject(m,"NwchemError",NwchemError);

}
#endif


