/*
###############################################################################
# If you use PhysiCell in your project, please cite PhysiCell and the version #
# number, such as below:                                                      #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1].    #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# See VERSION.txt or call get_PhysiCell_version() to get the current version  #
#     x.y.z. Call display_citations() to get detailed information on all cite-#
#     able software used in your PhysiCell application.                       #
#                                                                             #
# Because PhysiCell extensively uses BioFVM, we suggest you also cite BioFVM  #
#     as below:                                                               #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1],    #
# with BioFVM [2] to solve the transport equations.                           #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# [2] A Ghaffarizadeh, SH Friedman, and P Macklin, BioFVM: an efficient para- #
#     llelized diffusive transport solver for 3-D biological simulations,     #
#     Bioinformatics 32(8): 1256-8, 2016. DOI: 10.1093/bioinformatics/btv730  #
#                                                                             #
###############################################################################
#                                                                             #
# BSD 3-Clause License (see https://opensource.org/licenses/BSD-3-Clause)     #
#                                                                             #
# Copyright (c) 2015-2021, Paul Macklin and the PhysiCell Project             #
# All rights reserved.                                                        #
#                                                                             #
# Redistribution and use in source and binary forms, with or without          #
# modification, are permitted provided that the following conditions are met: #
#                                                                             #
# 1. Redistributions of source code must retain the above copyright notice,   #
# this list of conditions and the following disclaimer.                       #
#                                                                             #
# 2. Redistributions in binary form must reproduce the above copyright        #
# notice, this list of conditions and the following disclaimer in the         #
# documentation and/or other materials provided with the distribution.        #
#                                                                             #
# 3. Neither the name of the copyright holder nor the names of its            #
# contributors may be used to endorse or promote products derived from this   #
# software without specific prior written permission.                         #
#                                                                             #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" #
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE   #
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE  #
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE   #
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR         #
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF        #
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS    #
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN     #
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)     #
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE  #
# POSSIBILITY OF SUCH DAMAGE.                                                 #
#                                                                             #
###############################################################################
*/
#include <iostream>
#include <fstream>
#include "./custom.h"

void create_cell_types( void )
{
    // set the random seed
    SeedRandom( parameters.ints("random_seed") );
    
    /*
       Put any modifications to default cell definition here if you
       want to have "inherited" by other cell types.
       
       This is a good place to set default functions.
    */
    
    initialize_default_cell_definition();
    cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment );
    
    cell_defaults.functions.volume_update_function = standard_volume_update_function;
    cell_defaults.functions.update_velocity = standard_update_cell_velocity;

    cell_defaults.functions.update_migration_bias = NULL;
    cell_defaults.functions.update_phenotype = NULL; // update_cell_and_death_parameters_O2_based;
    cell_defaults.functions.custom_cell_rule = NULL;
    cell_defaults.functions.contact_function = NULL;
    
    cell_defaults.functions.add_cell_basement_membrane_interactions = NULL;
    cell_defaults.functions.calculate_distance_to_membrane = NULL;
    
    /*
       This parses the cell definitions in the XML config file.
    */
    
    initialize_cell_definitions_from_pugixml();
    
    /*
       Put any modifications to individual cell definitions here.
       
       This is a good place to set custom functions.
    */
    
    cell_defaults.functions.update_phenotype = phenotype_function;
    cell_defaults.functions.custom_cell_rule = custom_function;
    cell_defaults.functions.contact_function = contact_function;
    
    Cell_Definition* pCD = find_cell_definition( "tumor" );
    pCD->functions.update_phenotype = tumor_phenotype;
    
    /*
       This builds the map of cell definitions and summarizes the setup.
    */
        
    build_cell_definitions_maps();
    display_cell_definitions( std::cout );
    
    return;
}

void setup_microenvironment( void )
{
    // set domain parameters
    
    // put any custom code to set non-homogeneous initial conditions or
    // extra Dirichlet nodes here.
    
    // initialize BioFVM
    
    initialize_microenvironment();
    
    return;
}

void setup_tissue( void )
{
    double Xmin = microenvironment.mesh.bounding_box[0];
    double Ymin = microenvironment.mesh.bounding_box[1];
    double Zmin = microenvironment.mesh.bounding_box[2];

    double Xmax = microenvironment.mesh.bounding_box[3];
    double Ymax = microenvironment.mesh.bounding_box[4];
    double Zmax = microenvironment.mesh.bounding_box[5];
    
    if( default_microenvironment_options.simulate_2D == true )
    {
        Zmin = 0.0;
        Zmax = 0.0;
    }
    
    double Xrange = Xmax - Xmin;
    double Yrange = Ymax - Ymin;
    double Zrange = Zmax - Zmin;
    
    // create some of each type of cell
    
    Cell* pC;
    
    for( int k=0; k < cell_definitions_by_index.size() ; k++ )
    {
        Cell_Definition* pCD = cell_definitions_by_index[k];
        std::cout << "Placing cells of type " << pCD->name << " ... " << std::endl;
        for( int n = 0 ; n < parameters.ints("number_of_cells") ; n++ )
        {
            std::vector<double> position = {0,0,0};
            position[0] = Xmin + UniformRandom()*Xrange;
            position[1] = Ymin + UniformRandom()*Yrange;
            position[2] = Zmin + UniformRandom()*Zrange;
            
            pC = create_cell( *pCD );
            pC->assign_position( position );
        }
    }    
    // load cells from your CSV file (if enabled)
    load_cells_from_pugixml();
    
    return;
}

std::vector<std::string> my_coloring_function( Cell* pCell )
{ return damage_coloring(pCell); }

void phenotype_function( Cell* pCell, Phenotype& phenotype, double dt )
{ return; }

void custom_function( Cell* pCell, Phenotype& phenotype , double dt )
{ return; }

void contact_function( Cell* pMe, Phenotype& phenoMe , Cell* pOther, Phenotype& phenoOther , double dt )
{ return; }

//void tumor_phenotype( Cell* pC, Phenotype& p, double dt)
//{
//    // find my cell definition
//    Cell_Definition* pCD = find_cell_definition( pC->type );
//    // find index of O2 in the microenvironment
//    static int nAO = microenvironment.find_density_index( "antioxygen" );
//    // find index of necrosis death model
//    static int nApop = p.death.find_death_model_index( "apoptosis" );
//    // sample O2
//    double pAO = pC->nearest_density_vector()[nAO];
//    // set birth rate
//    // get base rate from cell definition
//    double base_rate = pCD->phenotype.cycle.data.transition_rate(0,0);
//    // set multiplier to 1.0
//    double multiplier = 1.0;
//
//    multiplier = pC->custom_data["pAO_min_prolif_factor"] + ( 1 - pC->custom_data["pAO_min_prolif_factor"] )
//    /( 1 + pow(pAO / pC->custom_data["pAO_half_change"], pC->custom_data["pAO_hill_power"]) );
//
//    // transition rate = base_rate * multiplier
//    p.cycle.data.transition_rate(0,0) = base_rate * multiplier;
//
//    if( p.death.dead == true )
//    {
//        p.secretion.set_all_secretion_to_zero();
//        p.secretion.set_all_uptake_to_zero();
//        pC->functions.update_phenotype = NULL;
//    }
//    return;
//}

void tumor_phenotype( Cell* pC, Phenotype& p, double dt)
{
    // find my cell definition
    Cell_Definition* pCD = find_cell_definition( pC->type );
    // find index of aO2 in the microenvironment
    static int nAO = microenvironment.find_density_index( "antioxygen" );
    // find index of necrosis death model
    static int nApop = p.death.find_death_model_index( "apoptosis" );
    // find index of damage variable
    static int nD = pC->custom_data.find_variable_index( "damage" );
    //find index of damage accumulation rate
    static int nDA = pC->custom_data.find_variable_index( "damage_accumulation_rate" );
    //find index of antioxygen metabolism within cells
    static int nM = pC->custom_data.find_variable_index( "AO_metabolism_rate" );
    // find index of repair parameter
    static int nR = pC->custom_data.find_variable_index( "repair_rate" );
    
    // use pressure to arrest proliferation
    if( pC->state.simple_pressure < pC->custom_data["pressure_threshold"] )
    {
        p.cycle.data.transition_rate(0,0) = pCD->phenotype.cycle.data.transition_rate(0,0);
        pC->custom_data["arrested"] = 0;
    }
    else
    {
        p.cycle.data.transition_rate(0,0) = 0;
        pC->custom_data["arrested"] = 1;
    }
    
    // internalized aO2 causes damage
    double pAO = p.molecular.internalized_total_substrates[nAO];
    pAO -= pC->custom_data[nM] * pAO * dt; // create degradation within cell to clear antioxygen
    if(pAO<0)
    {
        pAO=0;
    }
    p.molecular.internalized_total_substrates[nAO] = pAO; // set antioxygen based on this
    if(pAO>0)
    {
        pC->custom_data[nD] += pC->custom_data[nDA] * pAO * dt; // later can add proportionality constant, but this is just an abstract quantity, so why? could add to PD model linking to cell fate decisions
    }
    
    // set apoptosis rate
    // get base rate from cell definition
    double base_rate = pCD->phenotype.death.rates[nApop];
    
    // cell repairs damage
    pC->custom_data[nD] -= pC->custom_data[nR] * dt;
    if(pC->custom_data[nD]<=0)
    {
        pC->custom_data[nD]=0; // handle negative damage becuase of "too much" repair
        p.death.rates[nApop] = base_rate;
    }
    else // update apoptosis rate if there is damage
    {
        p.death.rates[nApop] = base_rate * ( 1 + pC->custom_data[nD] );
    }
    
    

    if( p.death.dead == true )
    {
        p.secretion.set_all_secretion_to_zero();
        p.secretion.set_all_uptake_to_zero();
        pC->functions.update_phenotype = NULL;
    }
    return;
}


static double tolerance = 0.01 * diffusion_dt; // using this in PK_model and write_cell_data_for_plots
static int dose_count = 0;

void PK_model( double current_time ) // update the Dirichlet boundary conditions as systemic circulation decays and/or new doses given
{
    static double next_dose_time = 0;
    static double systemic_circulation_concentration = 0.0;
    static double periphery_concentration = 0.0;

    // update systemic circulation and Dirichlet boundary conditions
    if( current_time > next_dose_time - tolerance && dose_count < parameters.ints("max_number_doses") )
    {
        systemic_circulation_concentration += parameters.doubles("systemic_circulation_increase_on_dose");
        for( int n=0; n < microenvironment.number_of_voxels(); n++ )
        {
            if( microenvironment.is_dirichlet_node( n ) )
            {
                microenvironment.update_dirichlet_node( n, 0, systemic_circulation_concentration * parameters.doubles("biot_number") );
            }
        }
        
        next_dose_time += parameters.doubles("dose_interval");
        dose_count++;
    }
    else
    {
        double systemic_circulation_change_rate = - parameters.doubles("systemic_circulation_elimination_rate");
        double periphery_change_rate = 0;
        double k = parameters.doubles("drug_flux_across_capillaries");
        double R = parameters.doubles("systemic_circulation_to_periphery_volume_ratio");
        
        systemic_circulation_change_rate += k * ( periphery_concentration/R -systemic_circulation_concentration );
        periphery_change_rate += k * ( systemic_circulation_concentration*R - periphery_change_rate );
        
        systemic_circulation_concentration +=  systemic_circulation_change_rate * diffusion_dt;
        periphery_concentration +=  periphery_change_rate * diffusion_dt;
        
        for( int n=0; n < microenvironment.number_of_voxels(); n++ )
        {
            if( microenvironment.is_dirichlet_node( n ) )
            {
                microenvironment.update_dirichlet_node( n, 0, systemic_circulation_concentration * parameters.doubles("biot_number") );
            }
        }
        
        
        
    }
    
    return;
 
}

void write_cell_data_for_plots( double current_time, char delim = ',') {
	// Write cell number data to a CSV file format time,tumor_cell_count
	// Can add different classes of tumor cells - apoptotic, necrotic, hypoxic, etc to this
	
	// NEED TO FIX - get the time in terms of minutes from start
	static double next_write_time = 0;
	if( current_time > next_write_time - tolerance ) {
		char dataFilename [256];
		sprintf(dataFilename, "%s/cell_counts.csv", PhysiCell_settings.folder.c_str());

		int tumorCount = 0;
		Cell* pC = NULL;
		
		for( int i=0; i < (*all_cells).size(); i++ ) {
			pC = (*all_cells)[i];
			if ( pC->type == 0 ) {
				tumorCount += 1;
			}
		}
		char dataToAppend [1024];
		sprintf(dataToAppend, "%d%c%d", next_write_time, delim, tumorCount);

		// append to file
		std::ofstream file_out;

		file_out.open(dataFilename, std::ios_base::app);
		file_out << dataToAppend << std::endl;
		next_write_time += parameters.doubles("csv_data_interval");
	}
	return;

}

std::vector<std::string> damage_coloring( Cell* pCell )
{
	// Update color of the cell based on damage in the cell
	// Initial color: grey, damaged cell: gradient of red, dead cell: black
	std::vector< std::string > output( 4 , "black" );
	// nucleus and both outlines are already black.

	// cytoplasm color
	// first, get the damage and convert to a value between 0 to 1
	
	// determine damage index
	int d_index = pCell->custom_data.find_variable_index( "damage" );
	double damage_value = pCell->custom_data[d_index];
	double norm_damage_value = damage_value/(1 + damage_value);
	
	if (pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::apoptotic ) { // apoptotic - black
		return output;
	}

	char colorTempString [128];
	if( pCell->phenotype.death.dead == false ) { // live cells
		if ( norm_damage_value <= 0 ) {
			sprintf(colorTempString, "rgb(128, 128, 128)");
		} else if ( norm_damage_value >= 1 ) {
			sprintf(colorTempString, "rgb(0, 0, 0)");
		} else {
			// Red gradient goes from (255, 200, 200) to (51, 0, 0) 
			int color = (int) round(norm_damage_value*100)*2; // converting to percentage to make the color gradient easier to set
			if ( color < 256) {
				sprintf(colorTempString, "rgb(%u, %u, %u)", 255-color, 200-color, 200-color); 
			} else {
			sprintf(colorTempString, "rgb(0, 0, 0)"); }
		}

		output[0].assign( colorTempString );
		output[2].assign( colorTempString );
		output[3].assign( colorTempString );
	}		
	return output;

}
