/*
Implementation of dl functions using WIN API.
(c) J. Saroun, 2017, saroun@ujf.cas.cz

IMPORTANT: This is not full implementation of dl, the API is simplified 
just to fit the needs of the SIMRES project. A lot of error handling stuff is missing ...

On Windows, use windows DLL API: LoadLibrary, GetProcAddress, ...
Otherwise assume dlfnc and just wrap the dl functions.

*/

#if defined(_WIN32) || defined(_WIN64) || defined(__CYGWIN__)

#include <stdio.h>
#include <string.h>
#include <windows.h>
#include "dlapi.h"


/* Simple replacement for dlerror() 
Note the different delcaration and return value!
*/
int mydlerror() {
  int ie;
  ie=GetLastError();
  if (ie !=0 ) {
     SetLastError(0); // clear error info
     return ie;
  } else {
     SetLastError(0); // clear error info
     return 0;
  };
}

/* Simple replacement for dlclose() 
Note the different delcaration and return value!
*/
int mydlclose(void *handle) {
	HMODULE hModule = (HMODULE) handle;
	BOOL res;
	res = FreeLibrary(handle);
	return (! res);
}

/* Load symbol from dll 
Unlike loadsym, returns error number as a parameter ires.
if ires<0 or handle==NULL, does nothing.
Returned ires vales:
	ires=1 if successful
	ires=<0 if error 
*/
void *mydlsym(void *handle, const char *name, int *ires) {
  int error;
  int *sh=NULL;
  mydlerror();    /* Clear any existing error */
  if (*ires >= 0 && handle != NULL ) {
    *(void **) (&sh) = GetProcAddress((HMODULE) handle, name);
    if ((error = mydlerror()) != 0)  {
      fprintf (stderr, "Error: cannot read symbol %s \n", name);
      *ires=-4;
    } else {
    *ires=1;
    };
  };
  return sh;
};

/* Simple replacement for dlopen()  */
void *mydlopen( const char *libname) {
	HMODULE hModule;
	hModule =  LoadLibrary(libname);
	return (void *) hModule;
}

#else
#include <cstddef>
#include <stdio.h>
#include <dlfcn.h>	
#include "dlapi.h"

int mydlerror() {
  char *error;
  error = dlerror();
  if (error==NULL) {
	  return 0;
  } else {
	  return 1;
  }
}


int mydlclose(void *handle) {
	return dlclose(handle);
}

void *mydlsym(void *handle, const char *name, int *ires) {
  void * res;	
  res = dlsym(handle, name);
  if (res==NULL) {
	  *ires=1;
  } else {
	   *ires=-4;
  }
  return res;
};


void *mydlopen( const char *libname) {
	return dlopen(libname, RTLD_LAZY);
}

#endif