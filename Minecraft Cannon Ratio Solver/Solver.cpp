#include <iostream>
#include <iomanip>
#include <cmath>


double objective_function(
	double* booster_toHP_x_exposures,			// To check if the Hammer's power is in range
	unsigned short int* booster_amounts,		// To check if the Hammer's power is in range
	double power_H_x1i_L,						// To check if the Hammer's power is in range
	double* booster_x_Hrange_virtual_exposure,	// To compute Hammer's final position
	double power_H_x0virtual_L,					// Hammer's virtual initial position
	double target_x_L							// To minimize
) {
	// Check if this will place the hammer's power within range
	// Position of hammer's power with current boosters:
	double power_H_x1_L = power_H_x1i_L;
	for (int i = 0; i < 12; i++) {
		power_H_x1_L += booster_amounts[i] * booster_toHP_x_exposures[i];
	}
	// Check if the Hammer's power is within explosion range
	if (power_H_x1_L > -1 || power_H_x1_L < -7.5) { // TODO: Play with bounds
		// return invalid
		std::cout << "  Hammer's power is out of range/ahead of it" << std::endl;
	}
	// else
	std::cout << "  Hammer's Power x position (local): " << power_H_x1_L << std::endl;
	std::cout << "  Hammer's Power x position (global): " << power_H_x1_L + -86214.5 << std::endl;
	// Find the final position of the Hammer
	double hammer_x17_L = power_H_x0virtual_L;
	for (int i = 0; i < 12; i++) {
		hammer_x17_L += booster_amounts[i] * booster_x_Hrange_virtual_exposure[i];
	}
	std::cout << "  Hammer at 17gt x position (local): " << hammer_x17_L << std::endl;
	std::cout << "  Hammer at 17gt x position (global): " << hammer_x17_L + -86214.5 << std::endl;
	// Compute the Error
	double error = hammer_x17_L - target_x_L;
	return error;
}

double hammer_projectile_after_17gt_L(double* booster_toHP_x_exposures, unsigned short int* booster_amounts, double power_H_x0_L, unsigned int power_H_amount) {
	// X Position of Hammer's Power
	double power_H_x1_L = power_H_x0_L;
	// Affected by boosters
	for (int i = 0; i < 12; i++) {
		power_H_x1_L += booster_amounts[i] * booster_toHP_x_exposures[i];
	}

	// std::cout << "  Hammer's power position (Global): " << power_H_x1_L + -86214.5 << std::endl;

	double hammer_x0_L = 0.0;							// Local initial position of the hammer
	double distance = hammer_x0_L - power_H_x1_L;		// Distance between Hammer's power and Hammer
	double f = (1.0 - distance / 8.0) * 1.0 / distance; // Force

	// Freefall of the hammer for 17 gt
	double hammer_u0 = f * power_H_amount * distance;	// Initial velocity of the Hammer
	double hammer_x17_L = hammer_x0_L;
	for (int i = 0; i < 17; i++) {
		hammer_x17_L += hammer_u0;
		hammer_u0 *= 0.98;
	}

	// Return's Hammer's local x position
	return hammer_x17_L;
}

int main() {
	/* Output Formatting */
	std::cout << std::setprecision(17);

	/* In-game Conditions (Global) */
	//// Positions/Velocities
	double center_x = -86214.5;
	double booster_x_values_G[12] = { 
		-86234.78586161736,
		-86235.03129429705,
		-86235.29123018472,
		-86235.56509943453,
		-86235.85226417152,
		-86236.15201904518,
		-86233.91205295271,
		-86234.25911295973,
		-86234.62067852802,
		-86234.99617983929,
		-86235.38497905637,
		-86235.78637088025
	};
	// State before afected by ratio and constant booster
	double power_H_x0_G = -86228.49000000954;	// -86228.49000000954 278.06125000119226 87930.5
	double power_H_u0_G = 0.0;					// 0.0 13.581624504779418 -2.4523938435550006E-15
	// State after constant booster only
	double power_H_x1i_G = -86232.08075066289;	// -86232.08075066289 316.3339737459176 87930.5
	double power_H_u1i_G = -3.5189356402947953;	// -3.5189356402947953 37.507269269830815 -1.158739770801276E-14
	// Amounts 
	unsigned int power_H_amount = 36; // minimum power
	
	/* Computated Quantities (Local) */
	// Positions
	double power_H_x0_L = power_H_x0_G - center_x; // Local Starting location of the hammer
	std::cout << "power_H_x0_L: " << power_H_x0_L << std::endl;
	double power_H_x1i_L = power_H_x1i_G - center_x; // Local Starting location of the hammer
	// Exposure from Booster to Hammer's power 
	double booster_toHP_x_exposures[12]; // Is also the distance it pushes the hammer along x
	std::cout << "Booster to Hammer's Power Exposure: " << std::endl;
	for (int i = 0; i < 12; i++) {
		double distance = power_H_x0_G - booster_x_values_G[i];
		double f = (1.0 - distance / 8.0) * 1.0 / distance; // Just to keep correct calculations TODO: simplify
		booster_toHP_x_exposures[i] = distance * f;
		std::cout << "  booster[" << i << "] (blocks/tick): " << booster_toHP_x_exposures[i] << std::endl;
	}
	// We can figure out how each booster affects the Hammer's position after 17 gameticks by the following:
	unsigned short int booster_ref_amounts[12] = { 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6 }; // Limited to Between 0-12
	double hammer_x_ref_final_L;

	// Reference Case
	std::cout << "Reference Booster Configuration { 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6 }: " << std::endl;
	hammer_x_ref_final_L = hammer_projectile_after_17gt_L(booster_toHP_x_exposures, booster_ref_amounts, power_H_x1i_L, power_H_amount);
	std::cout << "  Position of reference hammer x (Global): " << hammer_x_ref_final_L + center_x << std::endl;

	// We will increment/decrement booster[i] by 1 tnt and take the average affect one tnt change is
	std::cout << "Booster's indirect(virtual) affect to Hammer's Range after 17gt: " << std::endl;
	// Store the inderect relation of incrementing/decrementing the booster[i]
	double booster_x_Hrange_virtual_exposure[12]; // the effect of the hammer's range by adding one tnt to this booster
	for (int i = 0; i < 12; i++) {
		// Decrement 1
		booster_ref_amounts[i] = 5;
		double hammer_x_neg_perturbation_final_L = hammer_projectile_after_17gt_L(booster_toHP_x_exposures, booster_ref_amounts, power_H_x1i_L, power_H_amount);
		double decrement_rel_change = hammer_x_ref_final_L - hammer_x_neg_perturbation_final_L;

		// Increment 1
		booster_ref_amounts[i] = 7;
		double hammer_x_pos_perturbation_final_L = hammer_projectile_after_17gt_L(booster_toHP_x_exposures, booster_ref_amounts, power_H_x1i_L, power_H_amount);
		double increment_rel_change = hammer_x_pos_perturbation_final_L - hammer_x_ref_final_L;

		// Compute Inderect relation of incrementing / decrementing the booster[i]
		booster_x_Hrange_virtual_exposure[i] = (decrement_rel_change + increment_rel_change) / 2.0;
		// Print out virtual exposure
		std::cout << "  booster[" << i << "] (blocks): " << booster_x_Hrange_virtual_exposure[i] << std::endl;

		// reset
		booster_ref_amounts[i] = 6;
	}

	// We find the virtual starting location of the hammer range if no booster are applied
	double ref_booster_dist = 0.0;
	for (int i = 0; i < 12; i++) {
		ref_booster_dist += 6 * booster_x_Hrange_virtual_exposure[i];
	}
	double power_H_x0virtual_L = hammer_x_ref_final_L - ref_booster_dist;
	std::cout << "Virtual initial Hammer Position (local): " << power_H_x0virtual_L << std::endl;
	std::cout << "Virtual initial Hammer Position (global): " << power_H_x0virtual_L + center_x << std::endl;

	// Function Target
	double target_x_L = 100.0; // we want to find Hammer that is as close to 100.0 blocks away.
	double error;

	// What is the error of the reference case?
	std::cout << "\nReference Booster Configuration { 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6 }: " << std::endl;
	error = objective_function(booster_toHP_x_exposures, booster_ref_amounts, power_H_x1i_L, booster_x_Hrange_virtual_exposure, power_H_x0virtual_L, target_x_L);
	std::cout << "  Error: " << error << std::endl;

	// Test different ratios for validity
	unsigned short int booster_amounts[12] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 };
	std::cout << "\nTest Booster Configuration { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 }: " << std::endl;
	error = objective_function(booster_toHP_x_exposures, booster_amounts, power_H_x1i_L, booster_x_Hrange_virtual_exposure, power_H_x0virtual_L, target_x_L);
	std::cout << "  Error: " << error << std::endl;


	return 0;
}