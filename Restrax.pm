package Restrax;
###########################################################
#  $Id: Restrax.pm,v 1.8 2019/08/15 15:22:19 saroun Exp $
#  Perl package used in the RESTRAX script system.
#  Used for configuration, installation and packaging.
#  (c) J. Saroun, 2006
#
###########################################################
use base 'Exporter';
use Cwd;
use Config;

our @EXPORT = ('UX2DOS','SplitFileName','MkSubDirCmd','FileCopyCmd','getFileCopyCmd',
    'dosystem','CollectResources','RmFileCmd','RmDirCmd','StripPath','getTUGZip',
    'ZipDirCmd','FindDependences','SubstituteInFile');

our $dbg=0;
our $SYSNAME= "$Config{'osname'}";
our $ARCHNAME="$Config{'archname'}";
our $PROMPT="yes";   # set to no if you dont want any user confirmation on file delete
our $PATHDEL="\\";
our $HOST="i586";
use strict;
use File::Find ();
# for the convenience of &wanted calls, including -eval statements:
use vars qw/*name *dir *prune/;
*name   = *File::Find::name;
*dir    = *File::Find::dir;
*prune  = *File::Find::prune;

# set the pathname of the TUGZip tzscript.exe utility
# if not installed in the default location in "ProgramFiles\tugzip"
my $TUGZIP="";

my @INFILES=();     # array with *.in files
my @SRCFILES=();    # array with files to install
my @SRCDIRS=();     # array with directories to install


# Get unified OS name
$SYSNAME=$^O;   # system name
if ($SYSNAME =~ m/MSWin32/i) {
  $SYSNAME='win32';
  if ($ARCHNAME =~ m/.*x64.*/ or $ARCHNAME =~ m/.*x86_64.*/) {
    $HOST='x86_64';
  } else {    
    $HOST='x86_32';
  }
} else {
  if ($ARCHNAME =~ m/.*x86_64.*/) {
    $HOST='x86_64';
  } else {
    $HOST=`uname -m`;
	$HOST =~ s/\n//;
  }
};

# get TUGZip executable for win32
sub getTUGZip {
  if ( ! -e $TUGZIP) {
    $TUGZIP="$ENV{'ProgramFiles'}\\tugzip\\tzscript.exe";
    if ( ! -e $TUGZIP) {
      die "Packaging in Win32 requires TUGZip. \nSee http://www.tugzip.com/ and path setting in Restrax.pm.";
    };
  };
  return $TUGZIP;
};

# UX2DOS - replace / delimiters with \
sub UX2DOS {
 my $f=$_[0];
 if ($SYSNAME eq 'win32') {
  # skip PGPLOT_DEV which uses slash in device names ...
	if ( $f !~ m/^PGPLOT_DEV=\//) {
		$f =~ s/\/\//#b#/g;
		$f =~ s/\//\\/g;
		$f =~ s/#b#/\//g;
	};
 };
 return $f;
};

# if $dbg=0, executes system command
# if $dbg!=0, prints the command, but don't execute it
sub dosystem {
  if ($dbg ==0) {
    if (system(@_)==0) {return 0;} else {return 1;};}
  else {
    printf("  ");
    foreach my $f (@_) {printf("%s ",$f);};
    printf("\n");
    return 0;
  };
};

# Strip path from filename
sub StripPath {
  my $f=$_[0]; # full filename
  my $fstrip=$f;
  if ($f =~ m/^(.*)\/([\w.]*)/) {$fstrip="$2";};
  return $fstrip;
};

# Split file name to drive and path
# Except for drive, no path delimiters are included
# ./ prefix is removed from path
sub SplitDirName {
  my $f=UX2DOS("$_[0]");# full filename
  my $drv="";
  my $path="";

# get pathname without initial ./
  if ($f =~ m/^[.][\/\\](.*)/ ) {$f="$1"};
# get drive/root
  if ($f =~ m/^([A-Z][:][\/\\]+)(.*)[\/\\]*\z/i) { $drv="$1"; $path="$2";}
  elsif ($f =~ m/^([\/\\]+)(.*)[\/\\]*\z/i) { $drv="$1"; $path="$2";}
  else { $drv=""; $path="$f";};
  my @comp = ("$drv","$path");

#  printf("drive=     '%s'\n",$drv);
#  printf("path=      '%s'\n",$path);
  return @comp;
};

# Split file name to drive, path, name and extension
# Except for drive, no path delimiters are included
# ./ prefix is removed from path
sub SplitFileName {
  my $f=UX2DOS("$_[0]");# full filename
  my $drv="";
  my $path="";
  my $ext="";
  my $name="";

#  printf("splitting  '%s'\n",$f);

  my $filechars=qr/[\w\s.\-\&{}]+/; # reg exp. for allowed filename characters
# get pathname without initial ./
  if ($f =~ m/^[.][\/\\](.*)/ ) {$f="$1"};
# get drive/root
  if ($f =~ m/^([A-Z][:][\/\\]+)(.*)\z/i) { $drv="$1"; $f="$2";}
  elsif ($f =~ m/^([\/\\]+)(.*)\z/i) { $drv="$1"; $f="$2";};
# get path
  if ($f =~ m/^(.*)[\/\\](${filechars})\z/) { $path="$1"; $f="$2";}
  elsif ($f =~ m/^(.*)[\/\\]\z/) { $path="$1"; $f="";};
# get name and extension
  if ($f =~ m/^(${filechars})([.][\w]*)\z/) { $name="$1"; $ext="$2";}
  # no name
  elsif ($f =~ m/^([.][\w]*)\z/) { $ext="$1";}
  # no extension
  else { $name="$f";}
  my @comp = ("$drv","$path","$name","$ext");

#  printf("drive=     '%s'\n",$drv);
#  printf("path=      '%s'\n",$path);
#  printf("name=      '%s'\n",$name);
#  printf("extension= '%s'\n",$ext);

  return @comp;
};

# Copy subdirectory $_[0] to another parent directory, $_[1]
sub DirCopyCmd {
  my $SD=UX2DOS("$_[0]");
  my $PD=UX2DOS("$_[1]");
  my @cmd=();
  if ( -d "$PD") {
    if ($SYSNAME eq 'win32') {
      @cmd=("XCOPY","/E","/C","/Y","/I","$SD","$PD\\$SD");
    } else {
      @cmd=("cp","-rfp","$SD","$PD");
    };
    dosystem(@cmd);
  };
};

# create specified directory
# set privileges to 0755
# the specified directory can also include a filename
sub MkSubDirCmd {
  my $LD=UX2DOS("$_[0]"); # directory name
  my $mode=0755;
  my @path=SplitFileName("$LD"); # split pathname to drive, path, filename and extension
  my $PD="$path[0]";      # drive name for subsequent use
# split path to individual subdirectories
#                                        printf("MkSubDirCmd:  %s \n","$LD");
  my @dlist=split(/[\/\\]/,"$path[1]");
# create the subdirctory branch if necessary
  foreach my $D (@dlist) {
  if ($D ne "") {
# append next directory
    if ($PD eq "") {$PD = "$D";} else {$PD = "$PD/$D";};
# get required target directory, convert delimiters
    my $PTD=UX2DOS("$PD");
# create when necessary, set privileges to $mode
    if ( ! -d "$PTD" ) {
      printf("Creates missing directory:  %s \n","$PTD");
      if ($dbg == 0) {
        mkdir("$PTD");
        chmod $mode,"$PTD";
      };
    };
  };
  };
};

# return command for file to file copy, system dependent
sub getFileCopyCmd {
  my $source=UX2DOS("$_[0]"); # source file
  my $target=UX2DOS("$_[1]"); # target directory
  my @cmd=();
  if ($SYSNAME eq 'win32') {
    @cmd=("COPY","/Y","$source","$target");
  } else {
    @cmd=("cp","-fp","$source","$target");
  };
  return @cmd;
};

# Copy file $_[0] to another directory, $_[1]
# Set access permissions to $_[2]
# $_[0] may start with a path
# $_[1] may end with a filename (would be ignored)
sub FileCopyCmd {
  my $source="$_[0]"; # source file
  my $target="$_[1]"; # target directory
#                           printf("FileCopyCmd:  %s --> %s\n",$source,$target);
  my $mode=oct("$_[2]"); # permissions
# split source file to drive, path, filename and extension
  my @input=SplitFileName("$source");
#          printf("input  :");foreach my $f (@input) {printf("%s:",$f)};printf("\n");
# split target path to drive, path, filename and extension
  my @output=SplitDirName("$target");
#          printf("output :");foreach my $f (@output) {printf("%s:",$f)};printf("\n");
# get subdirectory to be appended to the target
  my $path="";
  if ($input[1] ne "") {$path="$input[1]/"};
# define target filename (with full path)
  my $prefix="$output[0]$output[1]/";
  if ($prefix eq "") {$prefix="./"};
  $target="$prefix$path$input[2]$input[3]";
#          printf("target: %s\n",$target);
# copy the file using a command appropriate for given system
  my $result="";
  printf("Restrax.pm: [%s] [%s]\n",$source, $target);
  if ( -f "$source") {
    MkSubDirCmd("$target"); # create directories if necessary
    my @cmd=getFileCopyCmd("$source","$target");
    dosystem(@cmd);
    if ($dbg == 0) {chmod $mode,"$target";};
  } else {$result="$source"};
  return $result;
};


# remove specified directory (CAUTION!!)
sub RmDirCmd {
  my $SD=UX2DOS("$_[0]");
  my @cmd=();
  my $answ = "yes";
  if ( -d "$SD") {
    if ($PROMPT eq "yes") {
      printf("WARNING ! Existing directory %s will be removed. Type 'yes' to continue: ",$SD);
      $answ=readline(*STDIN);
    };
    if ($answ =~ m/^yes[\n]*\z/) {
      if ($SYSNAME eq 'win32') {
        @cmd=("RMDIR","/Q","/S","$SD");
      } else {
        @cmd=("rm","-rf","$SD");
      };
      dosystem(@cmd);
    };
  };
};

# remove specified file (CAUTION!!)
sub RmFileCmd {
  my $SD=UX2DOS("$_[0]");
  my @cmd=();
  my $answ = "yes";
  if ( -f "$SD") {
    if ($PROMPT eq "yes") {
      printf("WARNING ! Existing file %s will be removed. Type 'yes' to continue: ",$SD);
      $answ=readline(*STDIN);
    };
    if ($answ =~ m/^yes[\n]*\z/) {
      if ($SYSNAME eq 'win32') {
        @cmd=("ERASE","/Q","/S","/F","$SD");
      } else {
        @cmd=("rm","-f","$SD");
      };
      dosystem(@cmd);
    };
  };
};

# Collects all files from the given list.
# The list may include both files and directories.
# Apply filters on allowed and denied file extensions.
# Return :
# @SRCFILES  ... all source files (incl. those in subdirectories)
# Both files and subdirectories remain listed in the global
# arrays @SRCFILES and @SRCDIRS, respectively.

sub CollectResources {
  my $incl= shift @_; # include extensions
  my $excl= shift @_; # exclude extensions
  my @SRCLIST= @_; # list of sources (files or directories)
  my @lst=();
  while ($#SRCFILES > -1) {shift @SRCFILES;};
  while ($#SRCDIRS > -1) {shift @SRCDIRS;};
  my $CD=cwd();       # current directory
#  printf("Current directory: %s\n",$CD);
  if ($dbg != 2) {printf("Collecting resources (allow:%s, deny:%s)\n",$incl,$excl)};

# scan input list and separate files from directories
  foreach my $D (@SRCLIST) {
# list directory + subdirectories
    my $DD=UX2DOS("$D");
#  printf("%s    \n",$DD);
    if ( -d $DD ) {
#  printf("directory: %s    \n",$DD);
      File::Find::find(\&wantedSrcDirs, "$D");}
# push files through filter
    elsif ( -f $DD ) {
#  printf("file: %s    \n",$DD);
      addSrcFile($DD,$incl,$excl);}
# glob * entries and push files through filter
    elsif ($DD =~ m/\*/) {
#  printf("glob: %s    \n",$DD);
      @lst=glob("$DD");
      foreach my $f (@lst) {
        addSrcFile($f,$incl,$excl);
      };
    };
  };
# scan directories and push files through filter
  foreach my $D (@SRCDIRS) {
    if ($D ne "") {
      my @list=glob("$D/*");
      foreach my $f (@list) {
        my $ff=UX2DOS("$f");
        if ( -f $ff ) {addSrcFile($ff,$incl,$excl);};
      };
    };
  };
  return @SRCFILES;
};

# Collect directories to be installed
sub wantedSrcDirs {
    my ($dev,$ino,$mode,$nlink,$uid,$gid);
    (($dev,$ino,$mode,$nlink,$uid,$gid) = lstat($_)) &&
    -d _ &&
    addSrcDir($dir,$_);
};

sub addSrcDir {
  my $D="$_[0]";
  my $F="$_[1]";
  my $DD="";
  if ($F ne ".") {$DD=UX2DOS("$D/$F")} else { $DD=UX2DOS("$D")};
# exclude CVS from distributions
  if ( $DD !~ m/.*CVS\z/)  {
    push @SRCDIRS, ("$DD");
  };
  return 1;
};

# Put given file on the list, apply specified filter on file extensions
sub addSrcFile {
  my $F="$_[0]";     # file name to be added
  my $incl="$_[1]";  # include extensions
  my $excl="$_[2]";  # exclude extensions
  my $regin=qr/.+/;  # default: allow all (except of empty strings)
  my $regout=qr/\|/; # default: deny none
  if ($incl ne "") {$regin=qr/^.+\.($incl)\z/};
  if ($excl ne "") {$regout=qr/^.+\.($excl)\z/};
#  printf("%s    \n",$F);
  if (( $F =~ m/$regin/i) && ( $F !~ m/$regout/i)) {
    push @SRCFILES, ("$F");
    if ($dbg != 2) {printf("   added %s\n",$F)};
  };
  return 1;
};

# create archive from a subdirectory
# returns archive name, when created, or ""
sub ZipDirCmd {
  my $SD="$_[0]";
  my @cmd=();
# script name for TUGZIP
  my $scrname="$SD.tzs";
# archive filename
  my $ARCH="$SD.tar.gz";
  if ($SYSNAME eq 'win32') {$ARCH="$SD.tgz";};

  if ( -d "$SD") {
# remove existing archive
    RmFileCmd("$ARCH");
# in win32, TUGZip is required
    if ($SYSNAME eq 'win32') {
# create and execute script for TUGZip
      my $TZP=getTUGZip;
      printf("Creating TZS script: %s \n",$scrname);
      my $SCR="";
      $SCR = $SCR . "function main() {\n";
      $SCR = $SCR . "var Comp = new Compress();\n";
      $SCR = $SCR . "Comp.Archive = \"$ARCH\";\n";
      $SCR = $SCR . "Comp.Type = \"TGZ\";\n";
      $SCR = $SCR . "Comp.Compression = 3;\n";
      $SCR = $SCR . "Comp.WorkingDir = \".\\\\\";\n";
      $SCR = $SCR . "Comp.Data = \"$SD\";\n";
      $SCR = $SCR . "Comp.Password = \"\";\n";
      $SCR = $SCR . "Comp.DateExtension = 0;\n";
      $SCR = $SCR . "Comp.TimeExtension = 0;\n";
      $SCR = $SCR . "Comp.Overwrite = 1;\n";
      $SCR = $SCR . "Comp.Recurse = 1;\n";
      $SCR = $SCR . "Comp.StoreFolderNames = 1;\n";
      $SCR = $SCR . "Comp.IncludeHiddenFiles = 0;\n";
      $SCR = $SCR . "Comp.DoCompress(); }\n";
      if ($dbg == 0) {
        open(SCRFILE, ">$scrname") or die "Cannot open TZS file: $scrname\n $! \n";
        printf(SCRFILE "%s",$SCR);
        close(SCRFILE);
      } else {
        printf("%s \n",$SCR);
      };
      @cmd=("$TZP","$scrname");
      printf("Executing %s \%s\n",$TZP,$scrname);
      dosystem(@cmd);

# in other systems, use GNU tar, gzip
    } else {
      @cmd=("tar","-cf","$SD.tar","$SD");
      dosystem(@cmd);
      @cmd=("gzip","$SD.tar");
      dosystem(@cmd);
    };
  } else {
    printf("Can't create archive from %s: directory does not exist \n",$SD);
  };
  if ( ! -e $ARCH) {$ARCH=""};
  return $ARCH;
};

# Find dependence in a file
sub FindDepInFile {
  my $f=$_[0];     # include file without path to search for
  my $LINE=$_[1];  # line to search
  my $exs=$_[2];   # source extension
  my $dep="";
  my @splname;
  if ($LINE =~ m/(.*)include(.*)/i) {
# Fortran *.inc
    if ($exs =~ m/f/i) {
      if ($LINE =~ m/^([\t ]*)include[ ]*[']($f)[']/i ) {
        $dep= $f;
      };
# C *.h
    } elsif ($exs =~ m/c/i) {
      if ($LINE =~ m/^([#])include[ ]*["]($f)["]/i ) {
        $dep= $f;
      };
    };
  } elsif ($LINE =~ m/^[\t ]*use[ ]+/i) {
# Fortran *.mod
    if ($exs =~ m/mod/i) {
      @splname=SplitFileName($f);
      if ($LINE =~ m/^([\t ]*)use[ ]*($splname[2])([\W]+|$)/i ) {
        $dep= $f;
      };
    };
  };
  return $dep;
};

# Find dependences for given source file
# Handles *.f files for Fortran and *.c files for C
sub FindDependences {
# input parameters
  my $FN= shift @_;    # source filename
  my $exs= shift @_;   # source extension
  my @LOOKUP=@_;       # list of all available include files
# local declarations
  my $inf="";
  my $LINE="";
  my $dep="";
  my $res="";

#  if ($FN =~ m/^(.*)[.]$exs/) {
    printf("$FN\n");
    open(SRCFILE, "<$FN") or die "$FN: $! \n";
    while (<SRCFILE>) {
      $LINE=$_;
      foreach my $f (@LOOKUP) {
        $inf=StripPath($f);
        if ($dep !~ m/^$inf$/i) { # add each include file only once
          if (($res=FindDepInFile($inf,$LINE,$exs)) ne "") {$dep="$dep $f";};
        };
      };
    };
    close(SRCFILE);
#  };
  return $dep;
};

# Substitutes for @var@ in a *.in file, using given hash table
# Saves the file with .in stripped off to a given target directory
# call as SubstituteInFile($inp,$TD,\%VAR)
# Update 11/2019: handle templates directory: ignore "templates" parent on targets
sub SubstituteInFile {
  my $inp= shift  @_;  # input file
  my $TD=  shift  @_;  # target directory
  my $addr= shift  @_;  # address of the %VAR hash
  my %VAR= %$addr;
  my $lin="";
  my @path=SplitFileName("$inp"); # Split file name to drive, path, name and extension
  my $file="$path[2]$path[3]";
  my $subpath="$path[1]"; # no drive nor leading /
  my $doUX2DOS= "no";
  if ($file =~ m/([\w.]*)([.]in)\z/) {
    my $fname="$1";
# in *.ini files, set correct path delimiters for filename entries
    if ($fname =~ m/^.*\.ini\z/) {$doUX2DOS= "yes"};
# Can not process itself
    if ($0 !~ m/.*[\/\\]$fname\z/) {
# construct source and target names
      my $source = "$file";
      if ($subpath ne "") {$source = UX2DOS("$subpath/$file")};
	  my $tsubpath = $subpath;
	  # in templates subdirectory: exclude "templaes/" from target file name 
	  if ( $subpath =~ m/^templates[\/\\](.*)/) {
		 $tsubpath = "$1";
	  } elsif ( $subpath =~ m/^templates/) {
		 $tsubpath = "";
	  };
	  #printf("subpath=%s, tsubpath=%s\n",$subpath,$tsubpath);
      my $target=UX2DOS("$TD/$fname");
      if ($tsubpath ne "") {$target=UX2DOS("$TD/$tsubpath/$fname")};
# do parsing
      if ($dbg != 2) {printf("    %s\n",$target)};
      if ($dbg == 0) {
        open(INFILE,"<$source") or die "Cannot open input file $source:\n $!\n";
        open(OUTFILE,">$target") or die "Cannot create output file $target:\n $!\n";
        while (<INFILE>) {
          if ($doUX2DOS eq "no") { $lin=$_;} else { $lin=UX2DOS($_);};
		  # printf("line: [%s]\n",$lin);
          for my $key (%VAR) {
		  if ( $lin =~ m/^.*[@]$key[@].*/ ) {
		     #printf("substitute: [%s] for [%s]\n",$key,$VAR{$key});
			 $lin =~ s/[@]$key[@]/$VAR{$key}/g;
		  };
		  };
          # $lin =~ s/\@[\w]*\@//g; # undefined variables replace with ""
          printf(OUTFILE "%s",$lin);
        };
        close(OUTFILE);
        close(INFILE);
      };
    };
  };
};
