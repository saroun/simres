#! /usr/bin/perl
# Attempts to finds Java VM for given system.
# This is a test for the findJRE function used by Install.pl script.
#
BEGIN {
  push @INC, '.';
}
use strict;
my %VARS;
my $SYSNAME="win32";
my $BITS=64; # 32 or 64
my $TD=".";

sub findJRE {
	my $abovever="1.7"; # version newer than 1.7 
	my @javas = ();
	my %goodJRE;
	my $JRE="";
	my $JRE_VER="0.0";
	my $jexe="jre/bin/java"; 
	my $pgm="/usr/lib/jvm";
	my $scrname="GUI/startGUI";
	if ($SYSNAME eq "win32") {
	   $jexe="bin/java.exe";  	
		if ($BITS eq 32) {
			$pgm=sprintf("%s/Java",$ENV{'PROGRAMFILES(X86)'});
		} else {
			$pgm=sprintf("%s/Java",$ENV{'PROGRAMFILES'});
		};
		$scrname="GUI/startGUI_win32.bat";
	};
	#printf("pgm=%s\n",$pgm);
	my @DIRS = glob("'$pgm/*'"); 	
	printf("Searching for Java VM, %d-bit:\n",$BITS);
	foreach my $d (@DIRS) {
	  #printf("dir=%s\n",$d);
	  if ( -d $d ) {
		my $jre=$d."/".$jexe;
		if ( -f $jre) {
		  push @javas ,$jre;
		};
	  };
	};

	foreach my $f (@javas) {
	  #printf("jre=%s\n",$f);
	  my $out=`"$f" -version 2>&1`;
	  my $b=32;
	  #printf("%s\n\n",$out);
	  if ( $out =~ m/version[\s\t]*"([0-9_.]+)"/) {
		 my $vstr=$1;
		 if ( $out =~ m/64[-]bit/i) {
		   $b=64;
		 };
		 if (($vstr > $abovever) and ($b==$BITS)) {
		   $goodJRE{$vstr}=$f;
		 };
	  };  
	};

	#printf("found Java VM: \n");
	for my $key (sort keys %goodJRE) {
	  printf("version=%-11s: %s\n",$key,$goodJRE{$key});
	  $JRE=$goodJRE{$key};
	  $JRE_VER=$key;
	};
	printf("\n");
  
  if ($JRE_VER > $abovever) {
    printf("Found Java VM version=%s\n%s\n\n",$JRE_VER,$JRE);
    printf("You can change this choice in %s/%s.\n",$TD,$scrname);
    printf("Press any key to continue ...");
    <STDIN>;
	$VARS{'JRE'}=$JRE;
  } else {
    printf("No suitable Java VM found for %s\n",$SYSNAME);
    printf("You will need to define JRE maunally in %s/%s.\n",$TD,$scrname);
    printf("Press any key to continue ...");
    <STDIN>;
  };
};

findJRE;
