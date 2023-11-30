// copyright ############################### //
// This file is part of the Xtrack Package.  //
// Copyright (c) CERN, 2021.                 //
// ######################################### //
#define POW2(X) ((X)*(X))
#define POW3(X) ((X)*(X)*(X))
#define POW4(X) ((X)*(X)*(X)*(X))
#ifndef XTRACK_IONLASERIP_H
#define XTRACK_IONLASERIP_H

/*gpufun*/
void CWLaser_track_local_particle(CWLaserData el, LocalParticle* part0){

    //The algorithm is partially from https://anaconda.org/petrenko/psi_beam_vs_laser

    double nx  = CWLaserData_get_laser_direction_nx(el);
    double ny  = CWLaserData_get_laser_direction_ny(el);
    double nz  = CWLaserData_get_laser_direction_nz(el);
    
    double laser_x = CWLaserData_get_laser_x(el);
    double laser_y = CWLaserData_get_laser_y(el);
    double laser_z = CWLaserData_get_laser_z(el);
    double w0 = CWLaserData_get_laser_waist_radius(el);

    double laser_intensity = CWLaserData_get_laser_intensity(el);   // W/m^2
    double laser_wavelength = CWLaserData_get_laser_wavelength(el); // Hz
    
    double ion_excited_lifetime = CWLaserData_get_ion_excited_lifetime(el); // sec
    double ion_excitation_energy = CWLaserData_get_ion_excitation_energy(el); // eV

    double cooling_section_length = CWLaserData_get_cooling_section_length(el); // m
    
    double p0c = LocalParticle_get_p0c(part0); // eV
    double m0  = LocalParticle_get_mass0(part0); // eV/c^2
    double hbar = 1.054571817e-34; // J*sec
        
    double gamma0 = sqrt(1.0 + p0c*p0c/(m0*m0));
    double beta0  = sqrt(1.0 - 1.0/(gamma0*gamma0));
    double OmegaTransition = ion_excitation_energy*QELEM/hbar; // rad/sec
    // printf("OmegaTransition = %e\n", OmegaTransition);
    
    double number_of_excitations = cooling_section_length/(beta0*gamma0*C_LIGHT*ion_excited_lifetime);

    // double kick_energy = number_of_excitations*-2.0*ion_excitation_energy*2.0*gamma0; //eV
    // double kick_momentum = kick_energy/C_LIGHT; //eV/c
    // double kick_delta = kick_momentum/p0c; //eV/c
    //printf("n_scattered = %e\n", n_scattered);
    //printf("kick_strength = %e\n", 1e10*kick_strength/(C_LIGHT*p0c));
    //printf("kick_delta = %e\n", kick_delta);


    //start_per_particle_block (part0->part)
    
        double state = LocalParticle_get_state(part);
        double delta = LocalParticle_get_delta(part);
        double z     = LocalParticle_get_zeta(part);
        double x     = LocalParticle_get_x(part);
        double y     = LocalParticle_get_y(part);
        double px    = LocalParticle_get_px(part);
        double py    = LocalParticle_get_py(part);
    
        double pc = p0c*(1.0+delta); // eV
        double gamma = sqrt(1.0 + pc*pc/(m0*m0));
        double beta  = sqrt(1.0 - 1.0/(gamma*gamma));
        double beta_x  = px*p0c/m0/gamma;
        double beta_y  = py*p0c/m0/gamma;
        double beta_z  = sqrt(beta*beta - beta_x*beta_x -beta_y*beta_y);

        double vx  = C_LIGHT*beta_x; // m/sec
        double vy  = C_LIGHT*beta_y; // m/sec
        double vz  = C_LIGHT*beta_z; // m/sec
           
        // Collision of ion with the laser pulse:
        // The position of the laser beam center is rl=rl0+ct*n. We can find the moment
        // when a particle with a position r=r0+vt collides with the laser as the moment
        // when r−rl is perpendicular to n. Then (r−rl,n)=0, which yields the equation
        // (r0,n)+(v,n)t−(rl0,n)−ct(n,n)=0. Hence
        // tcol=(r0−rl0,n)/[c−(v,n)]

        double tcol = ( (x-laser_x)*nx + (y-laser_y)*ny + (z-laser_z)*nz ) / (C_LIGHT - (vx*nx+vy*ny+vz*nz)); // sec

	    double xcol = x + vx*tcol; // m
	    double ycol = y + vy*tcol; // m
	    double zcol = z + vz*tcol; // m

        // r^2 to the laser center = |r-rl| at the moment tcol:
        double r2 = (\
                POW2(xcol - (laser_x+C_LIGHT*nx*tcol)) + \
                POW2(ycol - (laser_y+C_LIGHT*ny*tcol)) + \
                POW2(zcol - (laser_z+C_LIGHT*nz*tcol)) \
             ); // m

        double I = 4.0*gamma*gamma * laser_intensity;
        
        double OmegaRabi = C_LIGHT*sqrt(6*PI)*sqrt(I)/
        sqrt(hbar*ion_excited_lifetime*POW3(OmegaTransition));
                    
        // Detuning from the ion transition resonance in the ion rest frame:        
        double cos_theta = -(nx*vx + ny*vy + nz*vz)/(beta*C_LIGHT);
        double laser_omega_ion_frame = (2.0*PI*C_LIGHT/laser_wavelength)*(1.0+beta*cos_theta)*gamma;

        
        double DeltaDetuning = (OmegaTransition - laser_omega_ion_frame);
        double DeltaDetuningTau = DeltaDetuning/(2.0*gamma);

        double gamma_decay=1/ion_excited_lifetime;
        double k1=OmegaRabi*OmegaRabi/(gamma_decay*gamma_decay);
        double ratio_detuning_gamma = DeltaDetuning/gamma_decay;
        //printf("k1=%f\n\n",k1);
        //printf("ratio_detuning_gamma=%f\n\n",ratio_detuning_gamma);
        
        //printf("state: %lf\n", state);
        if (r2 < POW2(w0))
            {
            if (state > 0)
            {
                    double excitation_probability = 0.5*k1 / (4*ratio_detuning_gamma * ratio_detuning_gamma + k1 + 1);
                    //printf("excitation_probability=%f\n\n",excitation_probability);
                                    
                    double rnd = (float)rand()/(float)(RAND_MAX);
                    if ( rnd < excitation_probability )
                        {
                        LocalParticle_set_state(part, 2); // Excited particle
                        // photon recoil (from emitted photon!):
                        double rnd = (float)rand()/(float)(RAND_MAX);
                        
                        //LocalParticle_add_to_energy(part,-2.0*ion_excitation_energy*rnd*2.0*gamma, 0); // eV
                        LocalParticle_add_to_energy(part,number_of_excitations*-2.0*ion_excitation_energy*rnd*2.0*gamma, 0); // eV
                        }	
                    else
                        {
                        LocalParticle_set_state(part, 1); // Still particle
                        }
                    
                
            }
            }    
	//end_per_particle_block
    
}

#endif

