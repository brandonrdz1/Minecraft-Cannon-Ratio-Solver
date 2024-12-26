#include <iostream>
#include <iomanip>
#include <cmath>

#include "optimization.h"
namespace hammer {
	double explosion_height = 0.061250001192092896;
	double gravity = -0.04;
	double drag = 0.98;

	void compute_booster_to_power_x_exposure(double power_x0, double* booster_x_values_G, double* booster_to_power_x_exposures) {
		std::cout << "Booster to Hammer's Power Exposure: " << std::endl;
		for (int i = 0; i < 12; i++) {
			double distance = power_x0 - booster_x_values_G[i];
			double f = (1.0 - distance / 8.0) * 1.0 / distance; // Just to keep correct calculations TODO: simplify
			booster_to_power_x_exposures[i] = distance * f;			// Output

			std::cout << "  booster[" << i << "] (blocks/tick): " << booster_to_power_x_exposures[i] << std::endl;
		}
		return;
	}

	double position_after_17gt_L(double* booster_to_power_x_exposures, unsigned short int* booster_amounts, double power_x0_L, unsigned int power_amount) {
		// X Position of Hammer's Power
		double power_x1_L = power_x0_L;
		// Affected by boosters
		for (int i = 0; i < 12; i++) {
			power_x1_L += booster_amounts[i] * booster_to_power_x_exposures[i];
		}

		double hammer_x0_L = 0.0;							// Local initial position of the hammer
		double distance = hammer_x0_L - power_x1_L;		// Distance between Hammer's power and Hammer
		double f = (1.0 - distance / 8.0) * 1.0 / distance; // Force

		// Freefall of the hammer for 17 gt
		double hammer_u0 = f * power_amount * distance;	// Initial velocity of the Hammer
		double hammer_x17_L = hammer_x0_L;
		for (int i = 0; i < 17; i++) {
			hammer_x17_L += hammer_u0;
			hammer_u0 *= 0.98;
		}

		// Return's Hammer's local x position
		return hammer_x17_L;
	}

	void compute_booster_x_Hrange_virtual_exposure(unsigned short int* booster_amounts, double* booster_toHP_x_exposures, double power_H_x1i_L, unsigned int power_H_amount, double* booster_x_Hrange_virtual_exposure) {
		// Compute the reference position
		std::cout << "Booster to Hammer's Final Position Exposure: " << std::endl;
		for (int i = 0; i < 12; i++) { booster_amounts[i] = 6; }
		double hammer_x_ref_final_L = position_after_17gt_L(booster_toHP_x_exposures, booster_amounts, power_H_x1i_L, power_H_amount);

		// Find the effects of pertubations
		for (int i = 0; i < 12; i++) {
		// Decrement 1
		booster_amounts[i] = 5;
		double hammer_x_neg_perturbation_final_L = position_after_17gt_L(booster_toHP_x_exposures, booster_amounts, power_H_x1i_L, power_H_amount);
		double decrement_rel_change = hammer_x_ref_final_L - hammer_x_neg_perturbation_final_L;

		// Increment 1
		booster_amounts[i] = 7;
		double hammer_x_pos_perturbation_final_L = position_after_17gt_L(booster_toHP_x_exposures, booster_amounts, power_H_x1i_L, power_H_amount);
		double increment_rel_change = hammer_x_pos_perturbation_final_L - hammer_x_ref_final_L;

		// Compute Inderect relation of incrementing / decrementing the booster[i]
		booster_x_Hrange_virtual_exposure[i] = (decrement_rel_change + increment_rel_change) / 2.0;
		// Print out virtual exposure
		std::cout << "  booster[" << i << "] (blocks): " << booster_x_Hrange_virtual_exposure[i] << std::endl;

		// reset
		booster_amounts[i] = 6;
		}

	}

	double compute_virtual_initial_position(double* booster_x_range_virtual_exposure, double x_ref_final_L) {
		double ref_booster_dist = 0.0;
		for (int i = 0; i < 12; i++) {
			ref_booster_dist += 6 * booster_x_range_virtual_exposure[i];
		}
		double x0virtual_L = x_ref_final_L - ref_booster_dist;
		std::cout << "Initial Virtual Hammer's Position (local): " << x0virtual_L << std::endl;
		return x0virtual_L;
	}

	double compute_objective(
		double* booster_toHP_x_exposures,			// To check if the Hammer's power is in range
		unsigned short int* booster_amounts,		// To check if the Hammer's power is in range
		double power_H_x1i_L,						// To check if the Hammer's power is in range
		double* booster_x_Hrange_virtual_exposure,	// To compute Hammer's final position
		double hammer_x0virtual_L,					// Hammer's virtual initial position
		double target_x_L							// To minimize
	) {
		// Check if this will place the hammer's power within range
		// Position of hammer's power with current boosters:
		double power_H_x1_L = power_H_x1i_L;
		for (int i = 0; i < 12; i++) {
			power_H_x1_L += booster_amounts[i] * booster_toHP_x_exposures[i];
		}
		// Check if the Hammer's power is within explosion range
		if (power_H_x1_L > -0.5 || power_H_x1_L < -7.5) { // TODO: Play with bounds
			return 10000000.0; // should be infinity
		}
		// Find the final position of the Hammer
		double hammer_x17_L = hammer_x0virtual_L;
		for (int i = 0; i < 12; i++) {
			hammer_x17_L += booster_amounts[i] * booster_x_Hrange_virtual_exposure[i];
		}
		// Compute the Error
		double error = abs(hammer_x17_L - target_x_L);
		return error;
	}

	double debug_objective(
		double* booster_toHP_x_exposures,				// To check if the Hammer's power is in range
		std::vector<unsigned short>& booster_amounts,	// To check if the Hammer's power is in range
		double power_H_x1i_L,							// To check if the Hammer's power is in range
		double* booster_x_Hrange_virtual_exposure,		// To compute Hammer's final position
		double hammer_x0virtual_L,						// Hammer's virtual initial position
		double target_x_L								// To minimize
	) {
		// Check if this will place the hammer's power within range
		// Position of hammer's power with current boosters:
		double power_H_x1_L = power_H_x1i_L;
		for (int i = 0; i < 12; i++) {
			power_H_x1_L += booster_amounts[i] * booster_toHP_x_exposures[i];
		}
		// Check if the Hammer's power is within explosion range
		if (power_H_x1_L > -0.5 || power_H_x1_L < -7.5) { // TODO: Play with bounds
			// return invalid
			// std::cout << "  Hammer's power is out of range/ahead of it" << std::endl; // Penalty funciton
			return 10000000.0; // should be infinity
		}
		// else
		std::cout << "  Hammer's Power x position (local): " << power_H_x1_L << std::endl;
		std::cout << "  Hammer's Power x position (global): " << power_H_x1_L + -86214.5 << std::endl;
		// Find the final position of the Hammer
		double hammer_x17_L = hammer_x0virtual_L;
		for (int i = 0; i < 12; i++) {
			hammer_x17_L += booster_amounts[i] * booster_x_Hrange_virtual_exposure[i];
		}
		std::cout << "  Hammer at 17gt x position (local): " << hammer_x17_L << std::endl;
		std::cout << "  Hammer at 17gt x position (global): " << hammer_x17_L + -86214.5 << std::endl;
		// Compute the Error
		double error = abs(hammer_x17_L - target_x_L);
		return error;
	}
}

namespace arrow {

	double arrow_eye = 0.12999999523162842;
	double gravity = -0.05f;
	double drag = 0.99f;

	void compute_booster_to_power_x_exposure(double power_x0, double* booster_x_values_G, double* booster_to_power_x_exposures) {
		std::cout << "Booster to Arrows's Power Exposure: " << std::endl;
		for (int i = 0; i < 12; i++) {
			double distance = power_x0 - booster_x_values_G[i];
			double f = (1.0 - distance / 8.0) * 1.0 / distance; // Just to keep correct calculations TODO: simplify
			booster_to_power_x_exposures[i] = distance * f;			// Output

			std::cout << "  booster[" << i << "] (blocks/tick): " << booster_to_power_x_exposures[i] << std::endl;
		}
		return;
	}

	void compute_arrow_x17gt_L_map(double* arrow_x17_map, double* arrow_u17_map, double power_A_amount) {
		// double x_min = -7.5, x_max = -0.5, step = 0.05;
		// double power_x_arr[141];
		// for (int i = 0; i <= 140; ++i) { power_x_arr[i] = x_min + i * step; }
		// for (int i = 0; i <= 140; ++i) { std::cout << "power_loc_arr[" << i << "] = " << power_x_arr[i] << std::endl; }
		
		double dy_squared = arrow_eye * arrow_eye;
		for (int i = 0; i <= 140; i++) {
			double dx = -7.5 + i * 0.05;
			double distance = std::sqrt(dx * dx + dy_squared);
			double f = (1.0 - distance / 8.0);
			double f0 = f;

			// elaborate way to compute initial velocity, similar to game
			for (int j = 1; j < power_A_amount; j++) { f0 += f; }
			
			// freefall
			double x = 0.0;		// initial local position of the arrow
			double v = f0;		// initial velocity of the arrow
			for (int j = 0; j < 17; j++) {
				x += v;
				v *= drag;
			}
			arrow_x17_map[i] = x;
			arrow_u17_map[i] = v;

			//std::cout << dx << ", " << arrow_u17_map[i] << std::endl;
		}
		return;
	}

	double compute_objective(
		double* booster_toAP_x_exposures,			// To check if the Hammer's power is in range
		unsigned short int* booster_amounts,		// To check if the Hammer's power is in range
		double power_A_x1i_L,						// To check if the Hammer's power is in range
		double* arrow_x17_map,	
		double target_x_L							// To minimize
	) {
		// Check if this will place the hammer's power within range
		// Position of hammer's power with current boosters:
		double power_A_x1_L = power_A_x1i_L;
		for (int i = 0; i < 12; i++) {
			power_A_x1_L += booster_amounts[i] * booster_toAP_x_exposures[i];
		}
		// Check if the Hammer's power is within explosion range
		if (power_A_x1_L > -1 || power_A_x1_L < -7.5) { // TODO: Play with bounds
			return 10000000.0; // should be infinity
		}
		
		// Map Lookup
		// Find the final position of the Hammer
		int lower_index = static_cast<int>((power_A_x1_L - -7.5) / 0.05);
		int upper_index = lower_index + 1;

		// Fetch corresponding values from the map
		double lower_value = arrow_x17_map[lower_index];
		double upper_value = arrow_x17_map[upper_index];

		// Linear interpolation
		double fraction = (power_A_x1_L - (-7.5 + lower_index * 0.05)) / 0.05;
		double arrow_x17_L = lower_value + fraction * (upper_value - lower_value);

		// Compute the error
		double error = std::abs(arrow_x17_L - target_x_L);
		return error;
	}

	double debug_objective(
		double* booster_toAP_x_exposures,				// To check if the Hammer's power is in range
		std::vector<unsigned short>& booster_amounts,	// To check if the Hammer's power is in range
		double power_A_x1i_L,							// To check if the Hammer's power is in range
		double* arrow_x17_map,
		double target_x_L								// To minimize
	) {
		// Check if this will place the hammer's power within range
		// Position of hammer's power with current boosters:
		double power_A_x1_L = power_A_x1i_L;
		for (int i = 0; i < 12; i++) {
			power_A_x1_L += booster_amounts[i] * booster_toAP_x_exposures[i];
		}
		// Check if the Hammer's power is within explosion range
		if (power_A_x1_L > -1 || power_A_x1_L < -7.5) { // TODO: Play with bounds
			return 10000000.0; // should be infinity
		}

		// Arrow's Power initial position
		std::cout << "  Arrow's Power x position (local): " << power_A_x1_L << std::endl;
		std::cout << "  Arrow's Power x position (global): " << power_A_x1_L + -86214.5 << std::endl;

		// Map Lookup
		// Find the final position of the Hammer
		int lower_index = static_cast<int>((power_A_x1_L - -7.5) / 0.05);
		int upper_index = lower_index + 1;

		// Fetch corresponding values from the map
		double lower_value = arrow_x17_map[lower_index];
		double upper_value = arrow_x17_map[upper_index];

		// Linear interpolation
		double fraction = (power_A_x1_L - (-7.5 + lower_index * 0.05)) / 0.05;
		double arrow_x17_L = lower_value + fraction * (upper_value - lower_value);
		std::cout << "  Arrow at 17gt x position (local): " << arrow_x17_L << std::endl;
		std::cout << "  Arrow at 17gt x position (global): " << arrow_x17_L + -86214.5 << std::endl;

		// Compute the error
		double error = std::abs(arrow_x17_L - target_x_L);
		return error;
	}
}

int main() {
	/************* Output Formatting *************/
	std::cout << std::setprecision(17);

	/* In-game Conditions (Global) */
	//// Positions/Velocities
	double center_x = -86214.5;
	unsigned short int reference_booster_amounts[12] = { 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6 }; // Reference values, Limited to Between 0-12

	/** Hammer Input **/
	double hammer_booster_x_values_G[12] = { 
		-86234.78586161736, -86235.03129429705, -86235.29123018472,
		-86235.56509943453, -86235.85226417152, -86236.15201904518,
		-86233.91205295271, -86234.25911295973, -86234.62067852802,
		-86234.99617983929, -86235.38497905637, -86235.78637088025
	};
	// State before/after afected by ratio and constant booster and amount of hammer's power amount
	double power_H_x0_G = -86228.49000000954;	// -86228.49000000954 278.06125000119226 87930.5
	double power_H_u0_G = 0.0;					// 0.0 13.581624504779418 -2.4523938435550006E-15
	double power_H_x1i_G = -86232.08075066289;	// -86232.08075066289 316.3339737459176 87930.5
	double power_H_u1i_G = -3.5189356402947953;	// -3.5189356402947953 37.507269269830815 -1.158739770801276E-14
	unsigned int power_H_amount = 36;			// minimum power
	// Computed Quantities (local)
	double power_H_x0_L = power_H_x0_G - center_x;		// Local Starting location of the hammer
	double power_H_x1i_L = power_H_x1i_G - center_x;	// Local Starting location of the hammer
	double booster_toHP_x_exposures[12];				// The effect of the hammer's power's boosters to the hammer's power
	hammer::compute_booster_to_power_x_exposure(power_H_x0_G, hammer_booster_x_values_G, booster_toHP_x_exposures);
	double booster_x_Hrange_virtual_exposure[12];	// The inderect affect of the hammer's power's boosters onto the final position of the hammer
	hammer::compute_booster_x_Hrange_virtual_exposure(reference_booster_amounts, booster_toHP_x_exposures, power_H_x1i_L, power_H_amount, booster_x_Hrange_virtual_exposure);
	double hammer_x_ref_final_L = hammer::position_after_17gt_L(booster_toHP_x_exposures, reference_booster_amounts, power_H_x1i_L, power_H_amount); // Uses the reference values of {6 ... 6}
	// Finding the hammer's virtual starting position
	double hammer_x0virtual_L = hammer::compute_virtual_initial_position(booster_x_Hrange_virtual_exposure, hammer_x_ref_final_L);

	/** Arrow Input **/
	double arrow_booster_x_values_G[12] = {
		-86234.66023655233, -86234.92907985931, -86235.21242460351,
		-86235.50970078005, -86235.82027033666, -86236.14342772469,
		-86234.53621773126, -86234.8832512168, -86235.2447882072,
		-86235.62025871217, -86236.00902470213, -86236.41038066245
	};
	// State before/after afected by ratio and constant booster and amount of arrow's power amount
	double power_A_x0_G = -86228.49000000954;	// -86228.49000000954 278.06125000119226 87930.5
	double power_A_u0_G = 0.0;					// 0.0 13.581624504779418 -2.4523938435550006E-15
	double power_A_x1i_G = -86230.15683194347;	// -86230.15683194347 298.0187499494105 87930.5
	double power_A_u1i_G = -1.6334952952582336;	// -1.6334952952582336 19.55834994925387 -5.9841021027295936E-15
	unsigned int power_A_amount = 36;			// minimum power
	// Computed Quantities (local)
	double power_A_x0_L = power_A_x0_G - center_x;		// Local Starting location of the arrow
	double power_A_x1i_L = power_A_x1i_G - center_x;	// Local Starting location of the arrow
	double booster_toAP_x_exposures[12];				// The effect of the arrow's power's boosters to the arrow's power
	arrow::compute_booster_to_power_x_exposure(power_A_x0_G, arrow_booster_x_values_G, booster_toAP_x_exposures);
	std::cout << "Starting loc of arrow power w/out booster " << power_A_x1i_L << std::endl;
	double arrow_x17_map[141]; // based on 36 tnt
	double arrow_u17_map[141];
	arrow::compute_arrow_x17gt_L_map(arrow_x17_map, arrow_u17_map, power_A_amount);


	/************* Optimization problem *************/
	double target_x_L = 490.0; // we want to find hammer/arrow that is as close to x blocks away.
	double error;

	// Define the hammer objective function lambda
	auto hammer_objective_function = [&](const std::vector<unsigned short>& booster_amounts) {
		unsigned short booster_amounts_array[12];
		std::copy(booster_amounts.begin(), booster_amounts.end(), booster_amounts_array);
		return hammer::compute_objective(
			booster_toHP_x_exposures,
			booster_amounts_array,
			power_H_x1i_L,
			booster_x_Hrange_virtual_exposure,
			hammer_x0virtual_L,
			target_x_L
		);
	};
	
	// Define the arrow objective function lambda
	auto arrow_objective_function = [&](const std::vector<unsigned short>& booster_amounts) {
		unsigned short booster_amounts_array[12];
		std::copy(booster_amounts.begin(), booster_amounts.end(), booster_amounts_array);
		return arrow::compute_objective(
			booster_toAP_x_exposures,
			booster_amounts_array,
			power_A_x1i_L,
			arrow_x17_map,
			target_x_L
		);
	};

	/* Initialize booster amounts */
	std::vector<unsigned short> hammer_booster_amounts = { 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6 };
	std::vector<unsigned short> arrow_booster_amounts = { 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6 };

	/* Perform Simulated Annealing optimization */
	double best_error_SA = optimization::optimize_boosters_simulated_annealing(
		hammer_objective_function,
		hammer_booster_amounts,	// Just to get the size of the vector
		0,						// Lower bound
		12,						// Upper bound
		1000000,				// Max iterations
		1000.0,					// Initial temperature
		0.999999				// Cooling rate
	);
	// Display results 
	std::cout << "\nOptimal Simulated Annealing Booster Configuration: { ";
	for (size_t i = 0; i < hammer_booster_amounts.size(); ++i) { std::cout << hammer_booster_amounts[i] << (i < hammer_booster_amounts.size() - 1 ? ", " : " "); }
	std::cout << "}\n";
	std::cout << "  Optimal Error: " << best_error_SA << "\n";
	std::cout << "Debug Hammer: " << std::endl;
	error = hammer::debug_objective(booster_toHP_x_exposures, hammer_booster_amounts, power_H_x1i_L, booster_x_Hrange_virtual_exposure, hammer_x0virtual_L, target_x_L);
	std::cout << "  Error: " << error << std::endl;


	//////////////////////////////////////////////

	/* Perform Simulated Annealing optimization */
	double best_error_SA_arrow = optimization::optimize_boosters_simulated_annealing(
		arrow_objective_function,
		arrow_booster_amounts,   // Vector size
		0,                       // Lower bound
		12,                      // Upper bound
		1000000,                 // Max iterations
		1000.0,                  // Initial temperature
		0.999999                 // Cooling rate
	);

	// Display results
	std::cout << "\nOptimal Arrow Simulated Annealing Booster Configuration: { ";
	for (size_t i = 0; i < arrow_booster_amounts.size(); ++i) {
		std::cout << arrow_booster_amounts[i] << (i < arrow_booster_amounts.size() - 1 ? ", " : " ");
	}
	std::cout << "}\n";
	std::cout << "  Optimal Arrow Error: " << best_error_SA_arrow << "\n";
	error = arrow::debug_objective(booster_toAP_x_exposures, arrow_booster_amounts, power_A_x1i_L, arrow_x17_map, target_x_L);
	std::cout << "  Error: " << error << std::endl;
	return 0;
}