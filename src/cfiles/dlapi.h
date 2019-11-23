#ifndef DLAPI_H
#define DLAPI_H
#ifdef __cplusplus
extern "C" {
#endif

int mydlerror();
int mydlclose(void *handle); 
void *mydlsym(void *handle, const char *name, int *ires);
void *mydlopen( const char *libname);

#ifdef __cplusplus
}
#endif
#endif /* DLAPI_H */

