#! /usr/bin/perl
###########################################################
#  $Id: configure.pl,v 1.112 2019/09/06 19:32:11 saroun Exp $
#  Creates  'makefile' for SIMRES
#  (c) J. Saroun, 2006-2018
#  Uses: Restrax.pm
#  Usage: perl configure.pl <sys> [options] 
#  Takes system-dependent information from the 'config.<sys>' files.
#  Options:
#      -noparse    ... do not process *.in templates
#      -parseonly  ... only process *.in templates, no configuration of makefile etc.\n";
#      -nomake     ... do not create makefile
#      -m64        ... build 64 bit version
###########################################################
BEGIN {
	push @INC, '.';
}
use Restrax;
use strict;
use Cwd;

use vars qw/*dbg *SYSNAME *HOST/;
*dbg      = *Restrax::dbg;
*SYSNAME  = *Restrax::SYSNAME;
*HOST     = *Restrax::HOST;
printf("System: %s\n",$SYSNAME);

my $CD=cwd();       # current directory
$dbg=0;          # set 1 for debug: no system commands will be executed
my $ThisScript=$0;  # name of this script
my $PARSE="yes";    # parse *.in files by default (change via cmd options)
my $DOMAKE="yes";   # set to "no" if you don't want to create makefile
#------------  DEFINE VERSION HERE  ------------
my $PGMNAME="simres";      # program name
my $VERSION="6.4.0";       # version
my $LPGPLOT="";            # PGPLOT link option (if empty,  it will be set by this script)
my $LMCPL="";              # MCPL link option (if empty,  it will be set by this script)
#--------------------------------------------

# global declarations
my @cmdopt=();        # command line options

my $CFG;           # filename with system dependent options (./config/config.<sys>)
my $OBJFILES;      # list of common objects
my @MODFILES=();   # list of module objects (full path names with *.mod extension)
my %VARS;          # list of variables passed to the makefile
my %MODULES;       # list of modules (keys are module names, values are source files)
my %FSOURCES;      # Fortran source files (all)
my %CSOURCES;      # common C source files
my @INCFILES;      # list of *.inc files
my @HFILES=();      # list of *.h files
my @SRCDIR=();     # source directories
my @INFILES=();    # array with *.in files

# system dependent variables
# Linux:
my $rm="rm -f ";      # remove comman
my $cp="cp -f -p ";   # copy command
my $ext="";           # executable file extension
my $SHEXT="so";       # shared library extension
my $PGLIB="libpgplot.so"; # PGPLOT library name
my $PGSERV="pgxwin_server"; # PGPLOT server executable
my $MCPLLIB="libmcplio.so"; # MCPL shared library name
my $PGIMP="";
my $MCPLIMP="";

# Windows:
if ($SYSNAME eq 'win32') {
  $rm="ERASE ";       #  remove command
  $cp="COPY /Y ";     # copy command  
  $ext=".exe"; # we are using restrax as DLL in win32
  $SHEXT="dll"; # shared library extension
  $PGLIB="libpgplot.dll"; # PGPLOT library name
  $PGIMP="libpgplot.lib"; # PGPLOT import library name (to be linked width)
  $PGSERV="jsdriv_server.exe"; # PGPLOT server executable   
  $MCPLLIB="libmcplio.dll"; # MCPL shared library name
  $MCPLIMP="libmcplio.lib"; # MCPL import library name (to be linked width)
};


#---------------------------------- PROJECT DEFINITIONS -------------------------------------
# format version ID
my $BUILDNUM=$VERSION;
my $PACKNUM=$VERSION;
$PACKNUM =~ s/[.]//g;
$BUILDNUM =~ s/(\d*)[.](\d*)[.](\d*)/$1.$2$3/;

# get VARS hash table
# $VARS{'CD'}=UX2DOS($CD);   # current directory
$VARS{'BITS'}=32;
$VARS{'HOMEPAGE'}="http://neutron.ujf.cas.cz/restrax";  # RESTRAX hommepage
$VARS{'PGMNAME'}=$PGMNAME;        # program name
$VARS{'VERSION'}=$VERSION;        # version number
$VARS{'BUILDNUM'}=$BUILDNUM;      # build number (derived from version number)
$VARS{'PACKNUM'}=$PACKNUM;        # package number (derived from version number)
$VARS{'CFGDIR'}="./config";       # Path to the 'config.<sys> files
$VARS{'BDATE'}=localtime();       # build time stamp
$VARS{'EXEC'}="$PGMNAME$PACKNUM$ext"; # Name of the executable
$VARS{'SRC'}="src";               # Path to the source files
$VARS{'BIN'}="bin";               # Path to executable files
$VARS{'LIB'}="lib";               # Path to library files
$VARS{'MKFILE'}="makefile";            # makefile with full path
$VARS{'INSTDIR'}=".";             # default installation directory
$VARS{'SHEXT'}="$SHEXT";          # shared library extension, definned above
$VARS{'J3DCLS'}="j3d-jre/lib/ext"; # location of J3D jar classes (relative to ./GUI)
$VARS{'PGSRC'}="3rdparty/pgplot";        # location of PGPPLOT source files 
$VARS{'MCPLIO'}="mcplio";          # location of MCPLIO source files - binding of MCPL to SIMRES
$VARS{'MCPLSRC'}="3rdparty/mcpl/src/mcpl";        # location of MCPL source files 
$VARS{'SYSNAME'}="$SYSNAME";       # OS name

# Add some system-dependent items
if ($HOST eq 'x86_64') {           # location of J3D shared libraries (relative to ./GUI)
  $VARS{'BITS'}=64; 
} else {
  $VARS{'BITS'}=32;
}
# Windows:
if ($SYSNAME eq 'win32') {
  $VARS{'PGTGT'}="pgplot/windows";  # location of PGPLOT binding files
  $VARS{'PGPLOT_DEV'}="/jsdriv";     # default pgplot device
  $VARS{'J3DLIB'}="j3d-jre/bin";     # location of J3D shared libraries (relative to ./GUI)	
} else {
# Linux:
  $VARS{'PGTGT'}="pgplot/linux";  # location of PGPLOT binding files
  $VARS{'PGPLOT_DEV'}="/xserve";         # default pgplot device
  if ($HOST eq 'x86_64') {         # location of J3D shared libraries (relative to ./GUI)
	$VARS{'J3DLIB'}="j3d-jre/lib/amd64"; 
  } else {
	$VARS{'J3DLIB'}="j3d-jre/lib/i386"; 
  } 	
 };

# Executable files (will be given 0755 permissions)
my @XFILES=();
push @XFILES,glob("lib/*.$VARS{'SHEXT'}");
push @XFILES,glob("bin/*.dll");
push @XFILES,glob("*.pl");
push @XFILES,glob("*.bat");

# 
#-------------------------- FORTRAN MODULES ----------------------
# base system
$MODULES{"CONSTANTS"}="src/constants.f";
$MODULES{"SIMATH"}="src/simath.f";
$MODULES{"FILETOOLS"}="src/filetools.f";
$MODULES{"IO"}="src/io.f";
$MODULES{"MESSAGES"}="src/messages.f";
$MODULES{"LATTICE"}="src/lattice.f";
$MODULES{"mtmod"}="src/mtmod.f";
$MODULES{"RNDGEN"}="src/rndgen.f";
# primitives
$MODULES{"ENUMDEF"}="src/primitives/enumdef.f";
$MODULES{"FIELDDEF"}="src/primitives/fielddef.f";
$MODULES{"CLASSDEF"}="src/primitives/classdef.f";
$MODULES{"FIELDDATA"}="src/primitives/fielddata.f";
$MODULES{"FIELDDEF"}="src/primitives/fielddef.f";
$MODULES{"CLASSDATA"}="src/primitives/classdata.f";
$MODULES{"TRACELIB"}="src/primitives/tracelib.f";
# components
$MODULES{"FRAMES"}="src/components/frames.f";
$MODULES{"FRAMES_TRACE"}="src/components/frames_trace.f";
# $MODULES{"MIRRORS"}="src/components/mirrors.f";
$MODULES{"SOURCES"}="src/components/sources.f";
$MODULES{"SOURCES_PULSE"}="src/components/sources_pulse.f";
$MODULES{"SOURCES_TABLE"}="src/components/sources_table.f";
$MODULES{"SOURCES_TRACE"}="src/components/sources_trace.f";
$MODULES{"SOURCES_ISIS"}="src/components/sources_isis.f";
$MODULES{"DETECTORS"}="src/components/detectors.f";
$MODULES{"DETECTORS_TRACE"}="src/components/detectors_trace.f";
$MODULES{"SAMPLES"}="src/components/samples.f";
$MODULES{"SAMPLES_TRACE"}="src/components/samples_trace.f";
$MODULES{"SAMPLE_SINGLE"}="src/components/sample_single.f";
$MODULES{"SAMPLE_POLY_TABLE"}="src/components/sample_poly_table.f";
$MODULES{"SAMPLE_POLY"}="src/components/sample_poly.f";
$MODULES{"XTALS"}="src/components/xtals.f";
$MODULES{"XTALS_REF"}="src/components/xtals_ref.f";
$MODULES{"XTALS_TRACE"}="src/components/xtals_trace.f";
$MODULES{"GUIDES"}="src/components/guides.f";
$MODULES{"GUIDES_TRACE"}="src/components/guides_trace.f";
$MODULES{"CRYSTALS"}="src/components/crystals.f";
$MODULES{"CRYSTALS_TRACE"}="src/components/crystals_trace.f";
$MODULES{"DCHOPPERS"}="src/components/tof/dchoppers.f";
$MODULES{"DCHOPPERS_TRACE"}="src/components/tof/dchoppers_trace.f";
$MODULES{"MONITORS"}="src/components/monitors.f";
$MODULES{"MONITORS_TRACE"}="src/components/monitors_trace.f";
$MODULES{"IGROUPS"}="src/components/igroups.f";
$MODULES{"FOCARRAYS"}="src/components/focarrays.f";
$MODULES{"SGUIDES"}="src/components/sguides.f";
$MODULES{"SGUIDES_TRACE"}="src/components/sguides_trace.f";
# tables
$MODULES{"TABLES"}="src/tables/tables.f";
$MODULES{"TABLE_CRYSTALS"}="src/tables/table_crystals.f";
$MODULES{"TABLE_ATOMS"}="src/tables/table_atoms.f";
$MODULES{"MIRROR_TABLE"}="src/tables/mirror_table.f";
$MODULES{"STRAIN_TABLE"}="src/tables/strain_table.f";
$MODULES{"NITI"}="src/tables/niti.f";
# instruments
$MODULES{"SPECTROMETER"}="src/instruments/spectrometer.f";
# virtual classes
$MODULES{"VCTABLE_SAMPLES"}="src/vctable_samples.f";
$MODULES{"VCTABLE_COMPONENTS"}="src/vctable_components.f";
$MODULES{"VCTABLE_INSTRUMENTS"}="src/vctable_instruments.f";
$MODULES{"VCTABLE_OPTIONS"}="src/vctable_options.f";
$MODULES{"VCTABLE_COMMANDS"}="src/vctable_commands.f";
$MODULES{"VCTABLE"}="src/vctable.f";
# NESS
$MODULES{"TRACINGOPT"}="src/ness/tracingopt.f";
$MODULES{"REPORTSOPT"}="src/ness/reportsopt.f";
$MODULES{"TRACINGDATA"}="src/ness/tracingdata.f";
$MODULES{"CLASSES"}="src/ness/classes.f";
$MODULES{"NSTORE"}="src/ness/nstore.f";
$MODULES{"COMPONENTS"}="src/ness/components.f";
$MODULES{"COMPONENTS_IO"}="src/ness/components_io.f";
$MODULES{"EVENTLOGS"}="src/ness/eventlogs.f";
$MODULES{"REPORTS"}="src/ness/reports.f";
$MODULES{"GENERATOR"}="src/ness/generator.f";
$MODULES{"RESMAT"}="src/ness/resmat.f";
$MODULES{"TRACING"}="src/ness/tracing.f";
$MODULES{"INSTCONTROL"}="src/ness/instcontrol.f";
$MODULES{"DATASETS"}="src/ness/datasets.f";
$MODULES{"ARRAY3D"}="src/ness/array3d.f";
$MODULES{"EVENTMONITOR"}="src/ness/eventmonitor.f";
# XML
$MODULES{"XMLPARSE"}="src/xml/xmlparse.f";
$MODULES{"XMLSIMRES"}="src/xml/xmlsimres.f";
$MODULES{"XMLHANDLER"}="src/xml/xmlhandler.f";
$MODULES{"XMLWRITER"}="src/xml/xmlwriter.f";
$MODULES{"XMLINFO"}="src/xml/xmlinfo.f";
$MODULES{"XMLUTILS"}="src/xml/xmlutils.f";
# Console interface
$MODULES{"DIALOGS"}="src/dialogs.f";
$MODULES{"COMMANDS"}="src/commands.f";
# Graphics
$MODULES{"GRFDATA"}="src/graphics/grfdata.f";
$MODULES{"GRFPLOT"}="src/graphics/grfplot.f";
$MODULES{"GRFEXEC"}="src/graphics/grfexec.f";

# SIMRES specific (high-level)
$MODULES{"RESULTS"}="src/results.f";
$MODULES{"BEAM1D"}="src/beam1d.f";
$MODULES{"RESOL1D"}="src/resol1d.f";
$MODULES{"BEAM2D"}="src/beam2d.f";
$MODULES{"RESOL2D"}="src/resol2d.f";
$MODULES{"SIMEXEC"}="src/simexec.f";
$MODULES{"CMDHANDLER"}="src/cmdhandler.f";
# Optimization 
$MODULES{"OPTIMIZATION"}="src/optimization.f";
$MODULES{"MIRRLOG"}="src/mirrlog.f";
# Matrix optimization 
# excluded $MODULES{"DYNARRAY"}="src/matrix/dynarray.f";
# excluded $MODULES{"MATRIX"}="src/matrix/matrix.f";

# ---------------------------  directories with source files ------------------------
@SRCDIR=("$VARS{'SRC'}");
push @SRCDIR,("$VARS{'SRC'}/ness");
push @SRCDIR,("$VARS{'SRC'}/components");
push @SRCDIR,("$VARS{'SRC'}/components/tof");
push @SRCDIR,("$VARS{'SRC'}/instruments");
push @SRCDIR,("$VARS{'SRC'}/tables");
push @SRCDIR,("$VARS{'SRC'}/xml");
push @SRCDIR,("$VARS{'SRC'}/graphics");
push @SRCDIR,("$VARS{'SRC'}/primitives");
# excluded push @SRCDIR,("$VARS{'SRC'}/matrix");


# Target directories - will be created if they do not exist
my @TARDIR=("$VARS{'BIN'}/", "$VARS{'LIB'}/", "$VARS{'LIB'}/pgplot","doc/");

# Directories with *.inc and *.h files
my @INCDIR=("$VARS{'SRC'}", "$VARS{'SRC'}/ness", "$VARS{'MCPLIO'}");

#
# template *.in files to be processed at configiration time
my @TEMPLATES = ("templates");

# in windows, place *.dll in the same directory as binaries
# if ($SYSNAME eq 'win32') {$VARS{'LIB'}=$VARS{'BIN'};};

#---------------------------------- SUBROUTINES  --------------------------------------------

#-------------------------------
# handle command line parameters
#-------------------------------
sub CmdParam {
  # at least the argument with system configuration is required
  if ( $#ARGV >= 0) {
	  if ($ARGV[0] ne "-parseonly") {
			$CFG="$VARS{'CFGDIR'}/config.$ARGV[0]";
			if (! -e $CFG) { die "Can't open $CFG:\n $!\n";};
		};
  };
  if ( $#ARGV < 0) {
    print "\n";
    print "Usage: perl $0 <sys> [options]\n";
    print "Arguments:\n";
    print "    <sys> ... defines the configuration file $VARS{'CFGDIR'}/config.<sys>\n";
    print "Options:\n";
    print "    -noparse    ... do not process *.in templates\n";
    print "    -parseonly  ... only process *.in templates, no configuration of makefile etc.\n";
    print "    -nomake     ... do not create makefile\n";
    print "\n";
    die;
  };
  foreach my $a (@ARGV) {
    push @cmdopt,$a;
    if ($a =~ m/[-]*noparse\z/) {$PARSE="no";};
    if ($a =~ m/[-]*parseonly\z/) {$PARSE="only";};
    if ($a =~ m/[-]*nomake\z/) {$DOMAKE="no";};
    if ($a =~ m/[-]*m(32|64)\z/) {
	$VARS{'BITS'}=$1;
        printf("set bits=%d\n",$VARS{'BITS'});
    };
  };
};

#------------------------------
# find linked libraries
#------------------------------
sub FindLinks {
 # PGPLOT
 if ($LPGPLOT eq "") {
   if ($SYSNAME eq 'win32') {
   # windows: target is import library, install to lib
   # dll will go to bin
	$LPGPLOT=UX2DOS("\$(BIN)/$PGIMP");
   } else {
   # linux: target is so, install to lib
	$LPGPLOT="\$(LIB)/$PGLIB";
   };
 } else {
	if ( ! -e $LPGPLOT) {
      die "PGPLOT library '$LPGPLOT' does not exist\n $!";
    };
 };
 printf("PGPLOT link is %s\n",$LPGPLOT);
 # MCPL
 if ($LMCPL eq "") {
   if ($SYSNAME eq 'win32') {
   # windows: target is import library, install to lib
   # dll will go to bin
	$LMCPL=UX2DOS("\$(BIN)/$MCPLIMP");
   } else {
   # linux: target is so, install to lib
	$LMCPL="\$(LIB)/$MCPLLIB";
   };
 } else {
	if ( ! -e $LMCPL) {
      die "MCPL library '$LMCPL' does not exist\n $!";
    };
 };
 printf("MCPL link is %s\n",$LMCPL);  
};

#------------------------------------
# Check source and target directories
#------------------------------------
sub CheckDirectories {
# Check source directories
  my @dirlist=("$VARS{'CFGDIR'}");
  foreach my $d (@SRCDIR) {
    push @dirlist,$d;
  };
  foreach my $F (@dirlist) {
    print "testing source directory $F\n";
    if ( ! -e $F) {
      die "Inconsistent distribution, directory '$F' does not exist\n $!";
    };
  };
# check target directories, create when needed
  foreach my $F (@TARDIR) {
    print "testing target directory $F\n";
    if (! -d $F) {
      MkSubDirCmd($F);
    };
  };
};

#------------------------------
# Process configuration file
#------------------------------
sub ProcessConfig {
my $lino; my $LINE;
# set defaults:
  $VARS{'MAKE'}="make";             # make executable
  $VARS{'FC'}="gfortran";
  $VARS{'DBGFLAGS'}= "-Og -fcheck=all -Waliasing -Wampersand -Wsurprising -Wc-binding-type -Wintrinsics-std -Wintrinsic-shadow -Wline-truncation -Wtarget-lifetime -Winteger-division -Wreal-q-constant -Wundefined-do-loop -Wmaybe-uninitialized -Werror";
  $VARS{'FCFLAGS'}="-O -funderscoring -ffixed-form -ffixed-line-length-132 -Jbin";
  $VARS{'CC'}="gcc";
  $VARS{'CCFLAGS'}="-O";
  $VARS{'USER_LIBS'}="";
  $VARS{'LPATHS'}="";               # Path for include files required by some make versions
  $VARS{'restrax_LDADD'}="-lX11";
  $VARS{'restrax_LD'}="\$(FC)";
  $VARS{'restrax_LDFLAGS'}="";

# replace defaults by whatever is found in the config file
  open(INFILE,  "<$CFG") or die "Can't open $CFG:\n $!\n";
  while (<INFILE>) {
    $lino = $lino + 1;
    $LINE=$_;
    for my $key (keys %VARS) {
      if ($LINE =~ m/^$key=[ ]*(.*)/) {
#        printf("replace %s(%s) with %s\n",$key,"$VARS{$key}","$1");
        $VARS{$key}="$1";
      };
    };
  };
  close(INFILE);
};

#------------------------------------------------------------------------------------------------
# Collect source files
# look for all *.f and *.c files in directories listed in @SRCDIR
# fill %FSOURCES table with *.f files (key=filename without extension, value - full pathname)
# fill %CSOURCES table with *.c files (key=filename without extension, value - full pathname)
# fill @OBJFILES array with the list *.o files
# fill @INCFILES array with the list *.inc files
# fill @HFILES array with the list *.h files
#------------------------------------------------------------------------------------------------
sub CollectSources {
  my @LST;
  my $D;
  my $F;
  my $ff;
  my $file;
  my @splname;

# collect modules first
   $OBJFILES="";
   for my $key (keys %MODULES) {
      if ($dbg != 0) {printf("Modules: %s\n","$key")};
      @splname=SplitFileName("$MODULES{$key}");
      $OBJFILES="$OBJFILES $VARS{'BIN'}/$splname[2].o";
      #$F="$VARS{'BIN'}/$key.mod";
      $F="$VARS{'BIN'}/$splname[2].o";      
      push @MODFILES,$F;
   };
# $MODOBJFILES="";
# foreach $F (@MODFILES) {
#   $MODOBJFILES="$MODOBJFILES $F";
# };

# collect Fortran and C sources in %FSOURCES and %CSOURCES, respectively
# key=filename without extension, without path
# value=full filename with extension
  foreach my $D (@SRCDIR) {
    printf("Collecting sources from %s\n",$D);
    @LST=glob("$D/*.f");
    foreach $F (@LST) {
      if ($F =~ m/(.*[\w]+)([.]f)/i) {
        if ($dbg != 0) {printf("Fortran source: %s\n","$1$2")};
        $file=StripPath("$1");
        $FSOURCES{$file}="$1$2";
        $ff="$VARS{'BIN'}/$file.o";
        if ($OBJFILES !~ m/$ff/i) { # don't add objects for modules again
          $OBJFILES="$OBJFILES $ff";
        };
      };
    };
    @LST=glob("$D/*.c");
    foreach $F (@LST) {
      if ($F =~ m/(.*[\w]+)([.]c)/i) {
        if ($dbg != 0) {printf("C source: %s\n","$1$2")};
        $file=StripPath("$1");
        $CSOURCES{$file}="$1$2";
        $OBJFILES="$OBJFILES $VARS{'BIN'}/$file.o";
      };
    };
  };
# collect include and header files
  @INCFILES=();
  foreach my $D (@INCDIR) {
    printf("Collecting header files from %s\n",$D);
    @LST=glob("$D/*.inc");
    foreach $F (@LST) {
      if ($F =~ m/.*[.]inc/) {
        push @INCFILES,$F;
        if ($dbg != 0) {printf("Include file: %s\n",$F)};
      };
    };
  };
  @HFILES=();
  foreach my $D (@INCDIR) {
    @LST=glob("$D/*.h");
    foreach $F (@LST) {
      if ($F =~ m/.*[.]h/) {
        push @INCFILES,$F;
	push @HFILES,$F;
        if ($dbg != 0) {printf("Header file: %s\n",$F)};
      };
    };
  };
};

#-------------------------------------------
# Create config.inc
# Fortran include file with release + system info
#-------------------------------------------
sub CreateConfigHdr {
  open(OUTFILE,">$VARS{'SRC'}/config.inc") or die "Cannot create $VARS{'SRC'}/config.inc:\n $!\n";
  printf(OUTFILE "C! generated by configure:\n");
  printf(OUTFILE "      CHARACTER*32 PACKAGE,PACKAGE_VERSION,PACKAGE_DATE,SYSNAME\n");
#  printf(OUTFILE "      REAL*4 BUILD_NUMBER\n");
  printf(OUTFILE "      CHARACTER*3 SHEXT\n");
  printf(OUTFILE "      PARAMETER (\n");
  printf(OUTFILE "     &  PACKAGE=\"%s\",\n",$VARS{'PGMNAME'});
  printf(OUTFILE "     &  PACKAGE_VERSION=\"%s\",\n",$VARS{'VERSION'});
  printf(OUTFILE "     &  PACKAGE_DATE=\"%s\",\n",$VARS{'BDATE'});
#  printf(OUTFILE "     &  BUILD_NUMBER=\"%s\",\n",$VARS{'BUILDNUM'});
  printf(OUTFILE "     &  SYSNAME=\"%s\",\n",$SYSNAME);
  printf(OUTFILE "     &  SHEXT=\"%s\"\n",$VARS{'SHEXT'});
  printf(OUTFILE "     &  )\n");
  close(OUTFILE);
};


# ===================================================================
# Create complete makefile
# ===================================================================
sub CreateMakefile {
  my $f;
  my $key;
  my $depends;
  my $moddepends;
  my $incldir;
  my $object;
  my $mobject;

  open(OUTFILE,">$VARS{'MKFILE'}") or die "Cannot create $VARS{'MKFILE'}:\n $!\n";
  printf(OUTFILE "#  Automatically created makefile for %s %s\n",$VARS{'PGMNAME'}, $VARS{'VERSION'});
  printf(OUTFILE "#  %s\n",$VARS{'BDATE'});
  printf(OUTFILE "#  Configuration: %s\n",$CFG);
  printf(OUTFILE "#------------------------------------------------------\n#\n");
# Pass VARS to makefile
  for $key (sort (keys(%VARS))) {
    printf(OUTFILE "%s=%s\n",$key,$VARS{$key});
    # printf("%s=%s\n",$key,$VARS{$key});
  };
  printf(OUTFILE "\n\n");

# Default target
  printf(OUTFILE "# ----------------------------- Targets -------------------------------------\n");
  printf(OUTFILE "default: \$(BIN)/\$(EXEC) \n\n"); 
  printf(OUTFILE "\n",);

# Executable target
  # format options for linking
  my $libraries="$LPGPLOT $LMCPL \$(USER_LIBS) \$(restrax_LDADD)";
  if ($SYSNAME eq 'win32') {
	# Windows: don't link with MCPLIO, it will be loaded dynamically at runtime.
	$libraries="$LPGPLOT \$(USER_LIBS) \$(restrax_LDADD)";
  };
  my $options="\$(LPATHS) \$(restrax_LDFLAGS)";
  my $output="\$(BIN)/\$(EXEC)";
  printf(OUTFILE "# Executable target  \n");
  printf(OUTFILE "%s: libs objfiles \$(USER_LIBS) \n",$output);
  printf(OUTFILE "\t\@echo 'Linking \$(PGMNAME) \$(VERSION)' \n");
  printf(OUTFILE "\t\$(restrax_LD) \$(OBJFILES) %s  %s  -o %s \n",$libraries,$options,$output);
  if (! $SYSNAME eq "win32") {printf(OUTFILE "\tchmod 0755 $output \n");};
  printf(OUTFILE "\t\@echo ----------------------------------------------------\n",$f);
  printf(OUTFILE "\t\@echo Before running the program, execute Install.pl:\n",$f);
  printf(OUTFILE "\t\@echo perl Install.pl [installation directory]\n",$f);
  printf(OUTFILE "\t\@echo ----------------------------------------------------\n",$f);
  printf(OUTFILE "\n\n");

# clean
  $f="$VARS{'BIN'}/*.o  $VARS{'BIN'}/*.mod $VARS{'BIN'}/*.def";
  if ($SYSNAME eq "win32") {
    $f =~ s/\//\\/g;
  };
  printf(OUTFILE "clean: \n");
  printf(OUTFILE "\t$rm %s \n",$f);
  printf(OUTFILE "\n");

# cleandist
  $f="$VARS{'BIN'}/*.o $VARS{'BIN'}/*.dll $VARS{'BIN'}/*.lib $VARS{'BIN'}/*.mod ./*.mod \$(BIN)/\$(EXEC) ZipBin.pl ZipSrc.pl ZipDoc.pl ";
  $f=$f." $VARS{'LIB'}/pgplot/grfont.dat";
  $f=$f." $VARS{'LIB'}/*.so";
  if ($SYSNAME eq "win32") {
	$f=$f." $VARS{'BIN'}/$MCPLIMP $VARS{'BIN'}/$MCPLLIB";
	$f=$f." $VARS{'BIN'}/$PGIMP $VARS{'BIN'}/$PGLIB";
    $f =~ s/\//\\/g;
  } else {
    $f=$f." $VARS{'LIB'}/$MCPLLIB";
    $f=$f." $VARS{'LIB'}/$PGLIB";
	$f=$f." $VARS{'LIB'}/pgplot/pgxwin_server";
  };
  printf(OUTFILE "cleandist: \n");
  printf(OUTFILE "\t$rm %s \n",$f);
  printf(OUTFILE "\t\$(MAKE) -C \$VARS{'PGTGT'} -f makefile_gfortran  erase \n");
  printf(OUTFILE "\n");

# install
  printf(OUTFILE "install: \n");
  printf(OUTFILE "\t\@echo To install, you have to execute script: \n",$f);
  printf(OUTFILE "\t\@echo perl Install.pl [installation directory]\n",$f);
  printf(OUTFILE "\t\@echo \n",$f);
  printf(OUTFILE "\n");

# distbin
  printf(OUTFILE "distbin: ZipBin.pl\n");
  printf(OUTFILE "\tperl ZipBin.pl -Q\n",$f);
  printf(OUTFILE "\n");

# distsrc
  printf(OUTFILE "distsrc: ZipSrc.pl\n");
  printf(OUTFILE "\tperl ZipSrc.pl -Q\n",$f);
  printf(OUTFILE "\n");

# distutil
  printf(OUTFILE "distutil: ZipSrc.pl\n");
  printf(OUTFILE "\tperl ZipSrc.pl -util -Q\n",$f);
  printf(OUTFILE "\n");

# distdoc
  printf(OUTFILE "distdoc: ZipDoc.pl\n");
  printf(OUTFILE "\tperl ZipDoc.pl -Q\n",$f);
  printf(OUTFILE "\n");

# objfiles
  printf(OUTFILE "objfiles: \$(OBJFILES)\n");
  printf(OUTFILE "\n");

# libs
  printf(OUTFILE "libs: $LPGPLOT $LMCPL\n");
  printf(OUTFILE "\n");

#--------------------
# Build libmcplio
#--------------------
  printf(OUTFILE "# Compile %s \n",$MCPLLIB);
  printf(OUTFILE "$LMCPL: \n");
  if ($SYSNAME eq 'win32') {
   # windows
    $options="BIN=$CD/\$(BIN) MCPLSRC=$CD/\$(MCPLSRC)";
	$options =~ s/\//\\/g;
    printf(OUTFILE "\t\$(MAKE) -C \$(MCPLIO) -f makefile_windows  %s\n", $options);
    printf(OUTFILE "\t\$(MAKE) -C \$(MCPLIO) -f makefile_windows  clean\n");	
  } else { 
   # linux
    $options="BIN=$CD/\$(LIB) MCPLSRC=$CD/\$(MCPLSRC)";
    printf(OUTFILE "\t\$(MAKE) -C \$(MCPLIO) -f makefile_linux  %s\n", $options);
    printf(OUTFILE "\t\$(MAKE) -C \$(MCPLIO) -f makefile_linux  clean\n");		
  };    
  printf(OUTFILE "\n\n");
 
#----------------  
# Build libpgplot
#----------------
  printf(OUTFILE "# Compile %s\n",$PGLIB);
  printf(OUTFILE "$LPGPLOT:   \n");
  my $opt1="SRC=$CD/\$(PGSRC)";
  my $opt2="BIN=$CD/\$(LIB)  PGD=$CD/\$(LIB)/pgplot";	
# move result to the appropriate targets
  if ($SYSNAME eq 'win32') {
    # windows: install to ./bin, swap slashes
	$opt2="BIN=$CD/\$(BIN)  PGD=$CD/\$(LIB)/pgplot";
	$opt1 =~ s/\//\\/g;	
	$opt2 =~ s/\//\\/g;	
  };
  printf(OUTFILE "\t\$(MAKE) -C \$(PGTGT) -f makefile_gfortran  all %s\n",  $opt1);
  printf(OUTFILE "\t\$(MAKE) -C \$(PGTGT) -f makefile_gfortran  clean \n"); 	
  printf(OUTFILE "\t\$(MAKE) -C \$(PGTGT) -f makefile_gfortran  install %s\n", $opt2); 
  printf(OUTFILE "\n\n");

#------------------------------------------- 
# format options for compiling with Fortran
#-------------------------------------------
  $incldir="";
  foreach $f (@INCDIR) {$incldir="$incldir -I$f";};
  $options="\$(FCFLAGS)";

# Dependences for common Fortran sources
  printf(OUTFILE "# Dependences for Fortran sources \n");
  for $key (keys %FSOURCES) {
#    printf("key='%s'  value='%s' \n",$key,$FSOURCES{$key});
    $depends=FindDependences("$FSOURCES{$key}","f",@INCFILES);
    $moddepends=FindDependences("$FSOURCES{$key}","mod",@MODFILES);
    $object="$VARS{'BIN'}/$key.o";
    printf(OUTFILE "%s: %s %s %s\n",$object,$FSOURCES{$key}, $moddepends, $depends);
    printf(OUTFILE "\t\$(FC)  -c %s  %s  %s  -o %s \n",$FSOURCES{$key},$options,$incldir,$object);
  };
  printf(OUTFILE "\n\n");

# format options for compiling with C
  $options="\$(CCFLAGS)";

# Dependences for common C sources
  printf(OUTFILE "# Dependences for C sources \n");
  for $key (keys %CSOURCES) {
    $depends=FindDependences("$CSOURCES{$key}","c","h",@HFILES);
    $object="$VARS{'BIN'}/$key.o";
    printf(OUTFILE "%s: %s %s \n",$object,$CSOURCES{$key},$depends);
    printf(OUTFILE "\t\$(CC)  -c %s %s  %s  -o %s \n",$CSOURCES{$key},$options,$incldir,$object);
  };
  printf(OUTFILE "\n\n");

# close makefile
  close(OUTFILE);
};

# ===================================================================


# Scan *.in files in directories defined in @_ and substitute for
# variables marked as , using the %VARS hash table
sub ProcessInFiles {
  if ($PARSE ne "no") {
    printf("Processing templates: \n");
    foreach my $F (@INFILES) {
      SubstituteInFile($F,$VARS{'INSTDIR'},\%VARS);
    };
  };
};

# message at the end
sub EndTasks {
  my $mode=0755;
# define exec permission for executables
  foreach my $F (@XFILES) {
    my $LF=UX2DOS("$F");
    if ( -f $LF ) {chmod $mode,"$LF";};
  };
# print final messages
  printf("----------------------------------------------------------------\n ");
  printf("Use following command to compile RESTRAX:\n");
  printf("%s -f %s \n",$VARS{'MAKE'},$VARS{'MKFILE'});
  printf("----------------------------------------------------------------\n ");
};

#--------------------------------------------------------------
# MAIN
#--------------------------------------------------------------

# ---------------  Handle command line parameters  ---------------
CmdParam;

# Find PGPLOT
FindLinks;

# collect and process templates
@INFILES=CollectResources("in","",@TEMPLATES);
ProcessInFiles;
if ($PARSE eq "only") { exit;};

# ---------------  Check source and target directories  ---------------
CheckDirectories;

# ---------------  Process configuration file  ---------------
ProcessConfig;

# ---------------  create config.inc ----------------------
CreateConfigHdr;

# ---------------  Collect source files ----------------------
CollectSources;
$VARS{'OBJFILES'}=$OBJFILES;

#  ---------------  create makefile for RESTRAX  --------------------------
if ($DOMAKE eq "yes") {CreateMakefile;};

# finalization
EndTasks;
exit;
