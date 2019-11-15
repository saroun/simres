#include "mcpl.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mcplio.h" 

static mcpl_outfile_t mcplout;
static mcpl_file_t mcplin;
static mcpl_particle_t * particle;
static long long nparticles;
const double hsqov2m = 2.072124655;

void mcplopenwrite(char *fname) {
	mcplout = mcpl_create_outfile(fname);
	mcpl_hdr_set_srcname(mcplout,"SIMRES");
        mcpl_enable_universal_pdgcode(mcplout,2112);//all particles are neutrons
	// mcpl_hdr_add_comment(mcplout,"Test output.");
	particle = (mcpl_particle_t*)calloc(sizeof(mcpl_particle_t),1);
	mcpl_enable_doubleprec(mcplout);
}

void mcplopenread(char *fname) {
	mcplin = mcpl_open_file(fname);
	printf("Opened MCPL file produced with %s\n",mcpl_hdr_srcname(mcplin));
    for (unsigned i = 0; i < mcpl_hdr_ncomments(mcplin); ++i)
		printf("file had comment: '%s'\n",mcpl_hdr_comment(mcplin,i));
	nparticles = mcpl_hdr_nparticles(mcplin);
	printf("File containts %lu particles\n",nparticles);
	particle = (mcpl_particle_t*)calloc(sizeof(mcpl_particle_t),1);
}

void mcplnmax(double *n) {
	*n = (double)(1.0*nparticles);
}

void mcpladd(double r[], double k[], double t, double p)
{
    double nrm;
    /*positions are in cm*/
    particle->position[0]=r[0]*0.1;
    particle->position[1]=r[1]*0.1;
    particle->position[2]=r[2]*0.1;
    nrm =sqrt(k[0]*k[0] + k[1]*k[1] + k[2]*k[2]);
    /*ekin is in MeV*/
    particle->ekin = hsqov2m*nrm*nrm/1e9;
    particle->direction[0] = k[0]/nrm;
    particle->direction[1] = k[1]/nrm;
    particle->direction[2] = k[2]/nrm;
    /*time in ms:*/
    particle->time = t*1e-3;
    /*weight in unspecified units:*/
    particle->weight = p;
    mcpl_add_particle(mcplout,particle);
}

void mcplget(double r[], double k[], double s[], double *t, double *p)
{
	const mcpl_particle_t *particle; 
	double nrm;
    particle = mcpl_read(mcplin);
	if (!particle) {
		mcpl_rewind(mcplin);
		particle = mcpl_read(mcplin);
	}
    /*positions are in cm*/
    r[0]=particle->position[0]*10;
    r[1]=particle->position[1]*10;
	r[2]=particle->position[2]*10;
    /*ekin is in MeV*/
	nrm = sqrt(particle->ekin/hsqov2m*1e9);
	/* k-vector */ 
    k[0]=particle->direction[0]*nrm;
    k[1]=particle->direction[1]*nrm;
	k[2]=particle->direction[2]*nrm;
    /*time in ms: convert to us*/
    *t = (double)(particle->time*1e3);
    /*weight in unspecified units:*/
    *p = (double)(particle->weight);
	s[0]=particle->polarisation[0];
    s[1]=particle->polarisation[1];
	s[2]=particle->polarisation[2];
	if (s[0]+s[1]+s[2] == 0) {
		s[0]=1.0;
	}
}

void mcplclosewrite() {	
	mcpl_close_outfile(mcplout);
	free(particle);
}


void mcplcloseread() {	
	mcpl_close_file(mcplin);
	free(particle);
}	

