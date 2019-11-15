#ifndef MCPLIO_H
#define MCPLIO_H

#ifdef __cplusplus
extern "C" {
#endif

#if defined _WIN32 || defined __CYGWIN__
  #ifdef COMPILE_DLL
    #define MY_DLL  __declspec(dllexport)
  #else
    #define MY_DLL  __declspec(dllimport)
  #endif
#else
  #ifdef COMPILE_DLL
    #define MY_DLL
  #else
    #define MY_DLL
  #endif
#endif

MY_DLL void mcplopenwrite(char *fname);
MY_DLL void mcpladd(double r[], double k[], double t, double p);
MY_DLL void mcplclosewrite();
MY_DLL void mcplopenread(char *fname);
MY_DLL void mcplget(double r[], double k[], double s[], double *t, double *p);
MY_DLL void mcplcloseread();
MY_DLL void mcplnmax(double *n);


#undef MY_DLL
#ifdef __cplusplus
}
#endif

#endif

