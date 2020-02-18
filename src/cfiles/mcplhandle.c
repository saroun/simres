/*
$Id$
Wrappers callable by fortran for dynamically linked procedures from mcplio.dll
(c) J. Saroun, 2017, saroun@ujf.cas.cz

Define LINK_DYNAMIC to link dynamically (using dlopen, loadsym, ...). Otherwise,
the library is loaded at startup.

*/

/*
 Dynamic loading  ...
*/
#ifdef LINK_DYNAMIC
#include "dlapi.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


/* Define library name*/
#if defined(_WIN32) || defined(_WIN64) || defined(__CYGWIN__)
static const char LIBMCPLIO[] = "libmcplio.dll";
#else
static const char LIBMCPLIO[] = "libmcplio.so";
#endif


static char *lastname=NULL;
static void *mcplhandle=NULL;
void (*dll_mcplopenwrite)();
void (*dll_mcplclosewrite)();
void (*dll_mcpladd)();
void (*dll_mcplopenread)();
void (*dll_mcplcloseread)();
void (*dll_mcplget)();
void (*dll_mcplnmax)();

/* Copy string src to lastname, takes care about memory allocation
 returns number of copied characters */
int cpname(const char *src) {
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

/* Release library */
void releasemcpl_() {
    int error;
    if (mcplhandle != NULL) {
      mydlerror();    /* Clear any existing error */
      mydlclose(mcplhandle);
      if ((error = mydlerror()) != 0)  {
         fprintf (stderr, "Error in releasemcpl_ :%d\n", error);
         return;
      };
      printf("MCPL released: %s\n",lastname);
      if (lastname != NULL) {
        free(lastname);
        lastname=NULL;
      };
      mcplhandle=NULL;
      dll_mcplopenwrite=NULL;
	  dll_mcplclosewrite=NULL;
      dll_mcpladd=NULL;
      dll_mcplopenread=NULL;
	  dll_mcplcloseread=NULL;
      dll_mcplget=NULL;
	  dll_mcplnmax=NULL;
    };
}


/* Load the library dynamically
Return:
<0  library was released, either on request or due to an error
=0  already loaded, nothing changed
=1  newly loaded
on error, return negative value of the error code
*/
int loadmcpl_() {
    int ires=0;
    int ilen;
/* Load the library if not yeat loaded */
    if (!mcplhandle) {
      mydlerror();
	  mcplhandle = mydlopen(LIBMCPLIO);
      ires=mydlerror();
      if ((ires != 0) || !mcplhandle) {
        fprintf (stderr, "Cannot load library: %s, error=%d\n", LIBMCPLIO,ires);
        return -abs(ires);
      };
      ilen=cpname(LIBMCPLIO);
      ires=0;
      /* load pointers to subroutines */
      *(void **) (&dll_mcplopenwrite)=mydlsym(mcplhandle, "mcplopenwrite",&ires);
      *(void **) (&dll_mcplclosewrite)=mydlsym(mcplhandle, "mcplclosewrite",&ires);
      *(void **) (&dll_mcpladd)=mydlsym(mcplhandle, "mcpladd",&ires);
      *(void **) (&dll_mcplopenread)=mydlsym(mcplhandle, "mcplopenread",&ires);
      *(void **) (&dll_mcplcloseread)=mydlsym(mcplhandle, "mcplcloseread",&ires);
      *(void **) (&dll_mcplget)=mydlsym(mcplhandle, "mcplget",&ires);
      *(void **) (&dll_mcplnmax)=mydlsym(mcplhandle, "mcplnmax",&ires);
	  if (ires > 0) {
         fprintf (stdout, "Loaded library %s \n", LIBMCPLIO);
	  };
	  return ires;
    };
    return 0;
};


/* C wrapper for dynamically loaded fortran subroutines
*/

void mcpladd_(double *r[],double *k[],double *t, double *p)
{
  if (dll_mcpladd !=NULL ) (*dll_mcpladd)(r,k,*t,*p);
};

void mcplget_(double *r[],double *k[], double *s[], double *t, double *p)
{
  if (dll_mcplget !=NULL ) (*dll_mcplget)(r,k,s,t,p);
};

void mcplopenwrite_(char *fnam)
{
 if (dll_mcplopenwrite !=NULL ) (*dll_mcplopenwrite)(fnam,strlen(fnam));
};

void mcplopenread_(char *fnam)
{
 if (dll_mcplopenread !=NULL ) (*dll_mcplopenread)(fnam,strlen(fnam));
};

void mcplclosewrite_()
{
 if (dll_mcplclosewrite !=NULL ) (*dll_mcplclosewrite)();
};

void mcplcloseread_()
{
 if (dll_mcplcloseread !=NULL ) (*dll_mcplcloseread)();
};

void mcplnmax_(double *n)
{
  if (dll_mcplnmax !=NULL ) (*dll_mcplnmax)(n);
};

#else
/* ---------------------------------------------
 loaded at startup
---------------------------------------------*/
#include "mcplio.h"
#include <stdio.h>
int loadmcpl_() {
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
