// copyright ############################### //
// This file is part of the Xtrack Package.  //
// Copyright (c) CERN, 2021.                 //
// ######################################### //
#define POW2(X) ((X)*(X))
#define POW3(X) ((X)*(X)*(X))
#define POW4(X) ((X)*(X)*(X)*(X))
#ifndef XTRACK_ChirpedLaser_H
#define XTRACK_ChirpedLaser_H

/*gpufun*/
void ChirpedLaser_track_local_particle(ChirpedLaserData el, LocalParticle* part0){

    //The algorithm is partially from https://anaconda.org/petrenko/psi_beam_vs_laser

    double nx  = ChirpedLaserData_get_laser_direction_nx(el);
    double ny  = ChirpedLaserData_get_laser_direction_ny(el);
    double nz  = ChirpedLaserData_get_laser_direction_nz(el);
    
    double laser_x = ChirpedLaserData_get_laser_x(el);
    double laser_y = ChirpedLaserData_get_laser_y(el);
    double laser_z = ChirpedLaserData_get_laser_z(el);
    double w0 = ChirpedLaserData_get_laser_waist_radius(el);

    double laser_intensity = ChirpedLaserData_get_laser_intensity(el);   // W/m^2
    double laser_wavelength = ChirpedLaserData_get_laser_wavelength(el); // Hz
    double chirp_rate = ChirpedLaserData_get_chirp_rate(el); // Hz^2
    
    double ion_excited_lifetime = ChirpedLaserData_get_ion_excited_lifetime(el); // sec
    double ion_excitation_energy = ChirpedLaserData_get_ion_excitation_energy(el); // eV

    //double cooling_section_length = ChirpedLaserData_get_cooling_section_length(el); // m
    
    double p0c = LocalParticle_get_p0c(part0); // eV
    double m0  = LocalParticle_get_mass0(part0); // eV/c^2
    double hbar = 1.054571817e-34; // J*sec
        
    //double gamma0 = sqrt(1.0 + p0c*p0c/(m0*m0));
    //double beta0  = sqrt(1.0 - 1.0/(gamma0*gamma0));
    double OmegaTransition = ion_excitation_energy*QELEM/hbar; // rad/sec

    //number of excitations that will occur over the entire cooling section:        
    //double number_of_excitations = cooling_section_length/(beta0*gamma0*C_LIGHT*ion_excited_lifetime);

    // constants for the Map_of_Excitation vs OmegaRabi and Detuning:
    int64_t N_Omega_star_values = ChirpedLaserData_get_N_Omega_star_values(el);
    int64_t N_Detuning_star_values = ChirpedLaserData_get_N_Detuning_star_values(el);
    double  Omega_star_max = ChirpedLaserData_get_Omega_star_max(el);
    double  Detuning_star_max = ChirpedLaserData_get_Detuning_star_max(el);
    double  dOmega_star = Omega_star_max/(N_Omega_star_values-1.0);
    double  dDetuning_star = Detuning_star_max/(N_Detuning_star_values-1.0);


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
        // when r-rl is perpendicular to n. Then (r-rl,n)=0, which yields the equation
        // (r0,n)+(v,n)t-(rl0,n)-ct(n,n)=0. Hence
        // tcol=(r0-rl0,n)/[c-(v,n)]

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

        double OmegaRabi_star = OmegaRabi/chirp_rate;
                    
        // Detuning from the ion transition resonance in the ion rest frame:        
        double cos_theta = -(nx*vx + ny*vy + nz*vz)/(beta*C_LIGHT);
        double laser_omega_ion_frame = (2.0*PI*C_LIGHT/laser_wavelength)*(1.0+beta*cos_theta)*gamma;

        
        double DeltaDetuning = fabs(OmegaTransition - laser_omega_ion_frame);
        double Detuning_star = DeltaDetuning/chirp_rate;
        printf("OmegaRabi_star = %e\n", OmegaRabi_star);
        //printf("Detuning_star = %e\n", Detuning_star);

        //double gamma_decay=1/ion_excited_lifetime;
        //compute saturation parameter and normalized detuning:
       
        //1. only apply laser cooling to particles that are within the laser radius
        //2. only apply laser cooling to particles that have not been lost. In Xsuite, this means positive state
        if (r2 < POW2(w0))
            {
                    //only apply laser cooling to particles that have not been lost. In Xsuite, this means positive state
        if (state > 0)
            {    
 
        if (Detuning_star < Detuning_star_max && OmegaRabi_star > Omega_star_max)
            {
                //printf("high field case: DeltaDetuningTau = %f\n", DeltaDetuningTau);
                // printf("high field case: OmegaRabiTau = %e\n", OmegaRabiTau);
                // In case of a very high laser field:
                LocalParticle_set_state(part, 2); // Excited particle
                double rnd = RandomUniform_generate(part); //between 0 and 1
                LocalParticle_add_to_energy(part,-ion_excitation_energy*2.0*rnd*2.0*gamma, 0); // eV
            }
    
        else if (Detuning_star < Detuning_star_max &&
            OmegaRabi_star > dOmega_star/10.0)
            {
            // N_OmegaRabiTau_values  N_DeltaDetuningTau_values
            //   OmegaRabiTau_max       DeltaDetuningTau_max
            int64_t row = (int)floor(OmegaRabi_star/dOmega_star);
            int64_t col = (int)floor(Detuning_star/dDetuning_star);
            int64_t idx = row*N_Detuning_star_values + col;
            
            double excitation_probability = ChirpedLaserData_get_Map_of_Excitation(el, idx);

            //printf("excitation_probability = %e\n", excitation_probability);
                                        
            double rnd = RandomUniform_generate(part); //between 0 and 1
            if ( rnd < excitation_probability )
                {
                LocalParticle_set_state(part, 2); // Excited particle
                // photon recoil (from emitted photon!):
                double rnd = RandomUniform_generate(part); //between 0 and 1

                // If particle is excited, reduce its energy by, on average, the excitation energy with Lorentz boost
                // 2.0*rnd ensures that the average energy lost is the excitation energy
                // 2.0*gamma is the Lorentz boost
                LocalParticle_add_to_energy(part,-ion_excitation_energy*2.0*rnd*2.0*gamma, 0); // eV

                // //Also do transverse emission. This doesnt have a Lorentz boost
                // double rnd_x = ((double)rand() / (RAND_MAX)) * 2.0 - 1.0;//between -1 and +1
                // LocalParticle_add_to_px(part,rnd_x*ion_excitation_energy/p0c);
                // double rnd_y = ((double)rand() / (RAND_MAX)) * 2.0 - 1.0;//between -1 and +1
                // LocalParticle_add_to_py(part,rnd_y*ion_excitation_energy/p0c);
                }	
                else
                {
                LocalParticle_set_state(part, 1); // Still particle
                }

            }
            }
            
	//end_per_particle_block
    
}
}
#endif

