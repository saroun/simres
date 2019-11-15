BIN?=.
LIB?=.
MCDIR?=.
MCIMP?=libmcplio.a
MCLIB?=libmcplio.so
MCOPT?= -std=c99 -fPIC -I$(MCDIR)/mcpl -I$(MCDIR)
default: $(LIB)/$(MCLIB)
$(BIN)/mcpl.o: $(MCDIR)/mcpl/mcpl.c $(MCDIR)/mcpl/mcpl.h 
	gcc $(MCOPT) -c $(MCDIR)/mcpl/mcpl.c -o $(BIN)/mcpl.o 
$(BIN)/mcplio.o: $(MCDIR)/mcplio.c $(MCDIR)/mcplio.h 
	gcc $(MCOPT) -DCOMPILE_DLL -c $(MCDIR)/mcplio.c -o $(BIN)/mcplio.o 
$(LIB)/$(MCLIB): $(BIN)/mcpl.o $(BIN)/mcplio.o 
#	gcc -shared -static-libgcc  -static-libgcc_s_dw2-1 $(MCOPT) $(BIN)/mcpl.o $(BIN)/mcplio.o -Wl,--out-implib,$(LIB)/$(MCIMP) -o $(LIB)/$(MCLIB)
	gcc -shared $(MCOPT) $(BIN)/mcpl.o $(BIN)/mcplio.o -o $(LIB)/$(MCLIB)
 
clean:
	rm -r *.so *.lib *.o
