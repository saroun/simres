/*
$Id: mcplhandle.c,v 1.8 2019/08/15 15:02:08 saroun Exp $
C-wrappers for dynamically linked procedures from mcplio.dll
(c) J. Saroun, 2017, saroun@ujf.cas.cz
*/


//define ABSOFT in order to use windows API: LoadLibrary, GetProcAddress, ...

/*
 Windows, Absoft  - loaded at runtime using LoadLibrary ...
*/
#ifdef ABSOFT

#include <stdio.h>
#include <string.h>
#include <windows.h>

static char *lastname=NULL;
static void *mcplhandle=NULL;
void (*mcplopenwrite_)();
void (*mcplclosewrite_)();
void (*mcpladd_)();
void (*mcplopenread_)();
void (*mcplcloseread_)();
void (*mcplget_)();
void (*mcplnmax_)();

#define RTLD_NOW	0
#define RTLD_LAZY	1
#define RTLD_GLOBAL	2
#define RTLD_LOCAL	4


/* Simple replacement for dlerr() */
int mcplerr() {
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

void mcplclose() {
	FreeLibrary(mcplhandle);
}

/* Load symbol from dll */
void *loadsym(const char *name, int *ires) {
  int error;
  int *sh=NULL;
  mcplerr();    /* Clear any existing error */
  if (*ires >= 0 && mcplhandle != NULL ) {
    *(void **) (&sh) = GetProcAddress(mcplhandle, name);
    if ((error = mcplerr()) != 0)  {
      fprintf (stderr, "Error: cannot read symbol %s from %s \n", name,lastname);
      RELEASEMCPL();
      *ires=-4;
    } else {
    *ires=1;
    };
  };
  return sh;
};

/*  equivalent of the gcc dlopen 
void * dlopen( const char *libname, int code) {
	return LoadLibrary(libname);
}
*/

void * dlopen( const char *libname) {
	return LoadLibrary(libname);
}

/* Copy string src to lastname, takes care about memory allocation
 returns number of copied characters */
int cpstring(const char *src) {
  char *cpt;
  if (lastname != NULL) {
    free(lastname);
    lastname=NULL;
  };
  if (src == NULL) {
    lastname=NULL;
    return 0;
  };

  lastname = (char *) malloc(strlen(src) + 1);
  if (lastname == NULL) return -1;
  cpt=strcpy(lastname,src);
  if (cpt == NULL) {
        printf("Error in strcpy: <%s> <%s>\n",lastname,src);
        return -2;
  };
  return strlen(lastname);
}

/* Release EXCI library */
void RELEASEMCPL() {
    int error;
    if (mcplhandle != NULL) {
      mcplerr();    /* Clear any existing error */
      mcplclose();
      if ((error = mcplerr()) != 0)  {
         fprintf (stderr, "Error in RELEASEMCPL :%d\n", error);
         return;
      };
      printf("MCPL released: %s\n",lastname);
      if (lastname != NULL) {
        free(lastname);
        lastname=NULL;
      };
      mcplhandle=NULL;
      mcplopenwrite_=NULL;
	    mcplclosewrite_=NULL;
      mcpladd_=NULL;
      mcplopenread_=NULL;
	    mcplcloseread_=NULL;
      mcplget_=NULL;
    };
}


/* Load the EXCI library dynamically
Return:
<0 if error occured (EXCI released)
=0 if nothing changed
=1 if new EXCI was loaded
*/
int LOADMCPL(char *libname) {
    int ires=0;
    int icom;
    int ilen;

/* libname=="" => free previously loaded libraray and exit*/
    if (libname=="") {
      RELEASEMCPL();
      return -1;
    };

 /* compare new library name with the current one */
   if (lastname != NULL) {
      icom=strcmp(lastname,libname);
    } else { icom=-1; };

/* Load a new library named *libname if
 (a) there is no library loaded, or
 (b) the names of the current and requested libraries are different, or
 (c) mcplhandle == NULL
*/
    if ((icom!=0) || (!mcplhandle)) {

/* release previously loaded libraray if necessary */
      RELEASEMCPL();
      mcplerr(0);
      // mcplhandle = dlopen(libname, RTLD_LAZY);
	  mcplhandle = dlopen(libname);
      ires=mcplerr();
      if ((ires != 0) || !mcplhandle) {
        fprintf (stderr, "Cannot load library: %s, error=%d\n", libname,ires);
        return -abs(ires);
      };
      ilen=cpstring(libname);
      ires=0;
      /* load pointers to subroutines */
      *(void **) (&mcplopenwrite_)=loadsym("mcplopenwrite",&ires);
      *(void **) (&mcplclosewrite_)=loadsym("mcplclosewrite",&ires);
      *(void **) (&mcpladd_)=loadsym("mcpladd",&ires);
      *(void **) (&mcplopenread_)=loadsym("mcplopenread",&ires);
      *(void **) (&mcplcloseread_)=loadsym("mcplcloseread",&ires);
      *(void **) (&mcplget_)=loadsym("mcplget",&ires);
      *(void **) (&mcplnmax_)=loadsym("mcplnmax",&ires);
      return ires;
      fprintf (stdout, "Loaded library %s \n", libname);
    };
    return 0;
};


/* C wrapper for dynamically loaded fortran subroutines
*/

void MCPLADD(double *r[],double *k[],double *t, double *p)
{
  if (mcpladd_ !=NULL ) (*mcpladd_)(r,k,*t,*p);
};

void MCPLGET(double *r[],double *k[], double *s[], double *t, double *p)
{
  if (mcpladd_ !=NULL ) (*mcplget_)(r,k,s,t,p);
};

void MCPLOPENWRITE(char *fnam)
{
 if (mcplopenwrite_ !=NULL ) (*mcplopenwrite_)(fnam,strlen(fnam));
};

void MCPLOPENREAD(char *fnam)
{
 if (mcplopenread_ !=NULL ) (*mcplopenread_)(fnam,strlen(fnam));
};

void MCPLCLOSEWRITE()
{
 if (mcplclosewrite_ !=NULL ) (*mcplclosewrite_)();
};

void MCPLCLOSEREAD()
{
 if (mcplcloseread_ !=NULL ) (*mcplcloseread_)();
};

void MCPLNMAX(double *n)
{
  (*mcplnmax_)(n);
};

#else
/* ---------------------------------------------
 Linux, GCC - loaded at startup
---------------------------------------------*/
#include "mcplio.h"
#include <stdio.h>
int loadmcpl_(char *libname) {
	return 0;
};

void releasemcpl_() {
}

void mcpladd_(double *r, double *k, double *t, double *p)
{
  double tt = *t;
  double pp =*p;
   mcpladd(r, k, tt, pp);
};

void mcplget_(double *r, double *k, double *s, double *t, double *p)
{
   mcplget(r, k, s, t, p);
};

void mcplopenwrite_(char *fnam)
{
  mcplopenwrite(fnam);
};

void mcplclosewrite_()
{
  mcplclosewrite();
};

void mcplopenread_(char *fnam)
{
  mcplopenread(fnam);
};

void mcplcloseread_()
{
  mcplcloseread();
};

int mcplnmax_(double *n)
{
  mcplnmax(n);
};

#endif
