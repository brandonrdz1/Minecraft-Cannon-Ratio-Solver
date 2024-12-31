#include <iostream>
#include <iomanip>
#include <cmath>
#include <sstream>

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

class Tnt {
public:
	double x = 0.0, y = 0.0, z = 0.0;
	double u = 0.0, v = 0.0, w = 0.0;
	unsigned int amount = 1;
	
	static const double explosion_height, gravity, drag;

	Tnt(const std::string& attributes, unsigned int explosivePower = 1.0) : amount(explosivePower) {
		std::istringstream stream(attributes);
		stream >> x >> y >> z >> u >> v >> w;
	};

	void explosion(Tnt source, double exposure = 1.0) {
		double dx = x - source.x;
		double dy = y - (source.y + explosion_height);
		double dz = z - source.z;
		double distance = std::sqrt(dx * dx + dy * dy + dz * dz);
		double f = (1.0 - distance / 8.0) * 1.0 / distance * source.amount * exposure;
		
		u += f * dx ;
		v += f * dy ;
		w += f * dz ;
	}

	void freefall(int ticks) {
		// print(" 0 ");
		for (int i = 0; i < ticks; i++) {
			v += Tnt::gravity;

			x += u;
			y += v;
			z += w;
			// print(" " + std::to_string(i + 1) + " ");
			u *= Tnt::drag;
			v *= Tnt::drag;
			w *= Tnt::drag;
		}
	}

	std::vector<double> getVel() {
		return { u, v, w };
	}

	void print(std::string str = "") {
		std::cout << str << "x: " << x << "\ty: " << y << "\tz: " << z << "\t| u: " << u << "\tv: " << v << "\tw: " << w << std::endl;
	}
};
// Define static const members outside the class
const double Tnt::explosion_height = 0.061250001192092896;
const double Tnt::gravity = -0.04;
const double Tnt::drag = 0.98;

bool permutate(int* current, const int* start, const int* end, int size) {
	int index = size - 1;

	// Increment the current array to the next permutation
	while (index >= 0) {
		if (current[index] < end[index]) {
			current[index]++;
			return true; // Next permutation exists
		}
		else {
			current[index] = start[index]; // Reset this position and carry over
			index--;
		}
	}

	// If we've carried out of bounds, there are no more permutations
	return false;
}

// Function to check if two arrays are equal
bool isEqual(const int* arr1, const int* arr2, int size) {
	for (int i = 0; i < size; ++i) {
		if (arr1[i] != arr2[i]) {
			return false;
		}
	}
	return true;
}

int main() {
	/************* Output Formatting *************/
	std::cout << std::setprecision(17);

	Tnt arrowPower0("-132125.49000000954 250.0199999809265 1008.5");
	const int size = 8;
	Tnt* aPBoosters[size];

	aPBoosters[0] = new Tnt("-132125.50999999046 250.0 1000.5099999904633",3 * 2);
	aPBoosters[1] = new Tnt("-132125.50999999046 249.0 1002.4275000095367",2 * 2);
	aPBoosters[2] = new Tnt("-132125.50999999046 248.0 1004.4275000095367",10 * 2);
	aPBoosters[3] = new Tnt("-132125.50999999046 247.0 1006.4275000095367",5 * 2);
	aPBoosters[4] = new Tnt("-132125.49000000954 250.0 1000.5099999904633",1 * 2);
	aPBoosters[5] = new Tnt("-132128.50999999046 249.0 1002.4900000095367",9 * 2);
	aPBoosters[6] = new Tnt("-132128.50999999046 248.0 1004.4900000095367",3 * 2);
	aPBoosters[7] = new Tnt("-132128.50999999046 247.0 1006.4900000095367",9 * 2);

	for (int i = 0; i < size; i++) {
		arrowPower0.explosion(*aPBoosters[i]);
		delete aPBoosters[i]; // clear old boosters
	}

	aPBoosters[0] = new Tnt("-132145.49000000954 250.0 989.6775000095367");
	aPBoosters[1] = new Tnt("-132145.49000000954 250.0 991.2400000095367");
	aPBoosters[2] = new Tnt("-132145.49000000954 250.0 993.4900000095367");
	aPBoosters[3] = new Tnt("-132145.49000000954 250.0 995.4900000095367");
	aPBoosters[4] = new Tnt("-132147.50999999046 250.0 997.5099999904633");
	aPBoosters[5] = new Tnt("-132149.50999999046 250.0 997.5099999904633");
	aPBoosters[6] = new Tnt("-132151.75999999046 250.0 997.5099999904633");
	aPBoosters[7] = new Tnt("-132153.32249999046 250.0 997.5099999904633");

	Tnt preciseBooster0("-132145.50999999046 250.0 997.4900000095367");
	Tnt preciseBooster1("69420.0 250.0 1008.5"); 
	double boosterInfluence[size];
	std::cout << "Precise Exposure: " << std::endl;
	for (int i = 0; i < size; i++) {
		Tnt temp = preciseBooster0;
		temp.explosion(*aPBoosters[i], 0.9f); // snow slow-down
		boosterInfluence[i] = temp.u;
		std::cout << "  u[" << i << "]: " << boosterInfluence[i] << std::endl;
		delete aPBoosters[i];
	}

	int boosterAmountStart[size] = { 0, 0, 0, 0, 0, 0, 0, 0};
	int boosterAmount[size];
	int boosterAmountFinal[size] = { 34, 23, 23, 23, 23, 23, 23, 34};
	int targetArr[size] = { 1, 1, 1, 1, 1, 1, 1, 1 };
	// Initialize boosterAmount to the start values
	std::copy(boosterAmountStart, boosterAmountStart + size, boosterAmount);
	double diff;
	double closestValue = std::numeric_limits<double>::max();
	double startX = preciseBooster0.x;
	double target = 278 + Tnt::explosion_height;
	std::vector<int> bestBoosterAmount(size, 0);
	do {
		double finalX = startX;
		for (int i = 0; i < size; ++i) {
			finalX += boosterAmount[i] * boosterInfluence[i];
		}

		if (finalX > arrowPower0.x) { // check if its too far forward; continue the script and skip calculations

			permutate(boosterAmount, boosterAmountStart, boosterAmountFinal, size);
			if (!isEqual(boosterAmount, boosterAmountFinal, size)) {
				continue;
			}
		}
		
		preciseBooster1.x = finalX;
		Tnt arrowPower1 = arrowPower0;
		arrowPower1.explosion(preciseBooster1);
		arrowPower1.freefall(2);
		diff = std::abs(target - arrowPower1.y);
		if ((diff <= closestValue && arrowPower1.y >= target)){
			closestValue = diff;
			bestBoosterAmount = std::vector<int>(boosterAmount, boosterAmount + size);

			std::cout << "Current Combination: ";
			for (int i = 0; i < 8; ++i) {
				std::cout << boosterAmount[i] << " ";
			}
			std::cout << std::endl;
			std::cout << " x: " << finalX << std::endl;
			std::cout << " y: " << arrowPower1.y << std::endl;
			std::cout << " Error: " << diff << std::endl;
			
		}

		// Generate the next permutation
		permutate(boosterAmount, boosterAmountStart, boosterAmountFinal, size);
	} while (!isEqual(boosterAmount, boosterAmountFinal, size));

	std::cout << "End" << std::endl;
	while (true) {};
	return 0;
}
