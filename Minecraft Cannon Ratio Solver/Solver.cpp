#include <iostream>
#include <iomanip>
#include <cmath>
#include <sstream>

#include "optimization.h"

#define center -132111.5 // Where everything starts from
#define center_x -132111.5 // Where everything starts from

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
		for (int i = 0; i < 17; i++) { // 17
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
		std::cout << "  Hammer's Power x position (global): " << power_H_x1_L + center << std::endl;
		// Find the final position of the Hammer
		double hammer_x17_L = hammer_x0virtual_L;
		for (int i = 0; i < 12; i++) {
			hammer_x17_L += booster_amounts[i] * booster_x_Hrange_virtual_exposure[i];
		}
		std::cout << "  Hammer at 17gt x position (local): " << hammer_x17_L << std::endl;
		std::cout << "  Hammer at 17gt x position (global): " << hammer_x17_L + center << std::endl;
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
		std::cout << "  Arrow's Power x position (global): " << power_A_x1_L + center << std::endl;

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
		std::cout << "  Arrow at 17gt x position (global): " << arrow_x17_L + center << std::endl;

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

		u += f * dx;
		v += f * dy;
		w += f * dz;
	}

	void swing(Tnt source) {
		Tnt origin(source);
		std::vector<std::vector<double>> locations;
		double dx, dy, dz, distance;
		float f;
		for (int i = 0; i < source.amount; i++) {
			Tnt next = origin;
			next.freefall(1);
			if (next.y > 254.0199999809265) { next.y = 254.0199999809265; }
			std::vector<double> loc(3);
			loc[0] = next.x;
			loc[1] = next.y;
			loc[2] = next.z;
			locations.push_back(loc);
			dx = origin.x - locations[i][0], dy = origin.y - locations[i][1], dz = origin.z - locations[i][2];
			distance = sqrt(dx * dx + dy * dy + dz * dz);
			if (distance > 8.0 || distance == 0.0) {
				continue;
			}
			f = (1.0 - distance / 8.0) * 1.0 / distance;
			origin.u += dx * f;
			origin.v += dy * f;
			origin.w += dz * f;
		}

		for (int i = 0; i < locations.size(); i++) {
			std::cout << "x: " << locations[i][0] << " y: " << locations[i][1] << " z: " << locations[i][2] << std::endl;
		}
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

int main() {
	/************* Output Formatting *************/
	std::cout << std::setprecision(17);

	/* In-game Conditions (Global) */
	//// Positions/Velocities
	
	unsigned short int reference_booster_amounts[12] = { 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6}; // Reference values, Limited to Between 0-12

	/** Hammer Input **/
	double hammer_booster_x_values_G[12] = { 
		-132131.78586161736, -132132.03129429705, -132132.29123018472,
		-132132.56509943452, -132132.85226417152, -132133.1520190452,
		-132130.91192610186, -132131.2589431057, -132131.62046530127,
		-132131.99592300472, -132132.3846785224, -132132.78602670677
	};
	// State before/after afected by ratio and constant booster and amount of hammer's power amount
	double power_H_x0_G = -132125.49000000954;	// -132125.49000000954 278.06125000119226 1008.5
	double power_H_u0_G = 0.0;					// 0.0 13.581624504779414 9.596323735650002E-16
	double power_H_x1i_G = -132129.0807506629;	// -132129.0807506629 316.33397374591755 1008.5
	double power_H_u1i_G = -3.5189356402947953;	// -3.5189356402947953 37.5072692698308 -7.235323451482145E-15
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
		-132131.66023655233, -132131.9290798593, -132132.21242460352,
		-132132.50970078004, -132132.82027033667, -132133.1434277247,
		-132131.53607811127, -132131.8830685851, -132132.24456219366,
		-132132.61998908114, -132133.00871136136, -132133.41002367114
	};
	// State before/after afected by ratio and constant booster and amount of arrow's power amount
	double power_A_x0_G = -132125.49000000954;	// -132125.49000000954 278.06125000119226 1008.5
	double power_A_u0_G = 0.0;					// 0.0 13.581624504779414 9.596323735650002E-16
	double power_A_x1i_G = -132127.1568319435;	// -132127.1568319435 298.0187499494105 1008.5
	double power_A_u1i_G = -1.6334952952582336;	// -1.6334952952582336 19.558349949253866 -2.5024426975051027E-15
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


	//////////////////////////////////////////////
	/************* Optimization problem *************/
	double target_x_L = 200.0; // we want to find hammer/arrow that is as close to x blocks away.
	double error;

	// Define the hammer and arrow objective function lambdas
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

	//////////////////////////////////////////////
	const double error_threshold = 1e-5; // Desired error threshold
	int max_retries = 100; // Maximum retries to prevent infinite loops
	int retry_count = 0;

	/* Perform Simulated Annealing optimization */
	double best_error_SA = std::numeric_limits<double>::infinity();
	while (best_error_SA > error_threshold && retry_count < max_retries) {
		best_error_SA = optimization::optimize_boosters_simulated_annealing(
			hammer_objective_function,
			hammer_booster_amounts,
			0,                        // Lower bound
			12,                       // Upper bound
			1000000,                  // Max iterations
			100000.0,                 // Initial temperature
			0.99999                   // Cooling rate
		);
		retry_count++;
	}
	if (best_error_SA <= error_threshold) {
		std::cout << "\nOptimal Hammer Booster Configuration Found: { ";
		for (size_t i = 0; i < hammer_booster_amounts.size(); ++i) {
			std::cout << hammer_booster_amounts[i] << (i < hammer_booster_amounts.size() - 1 ? ", " : " ");
		}
		std::cout << "}\n";
		std::cout << "  Optimal Hammer Error: " << best_error_SA << "\n";

		std::cout << "Debug Hammer: " << std::endl;
		error = hammer::debug_objective(booster_toHP_x_exposures, hammer_booster_amounts, power_H_x1i_L, booster_x_Hrange_virtual_exposure, hammer_x0virtual_L, target_x_L);
		std::cout << "  Error: " << error << std::endl;
	}
	else {
		std::cout << "Failed to find a configuration with error < " << error_threshold << " after " << retry_count << " retries.\n";
	}

	// Reset retry counter for the arrow optimization
	retry_count = 0;
	/* Perform Simulated Annealing optimization */
	double best_error_SA_arrow = std::numeric_limits<double>::infinity();
	while (best_error_SA_arrow > error_threshold && retry_count < max_retries) {
		best_error_SA_arrow = optimization::optimize_boosters_simulated_annealing(
			arrow_objective_function,
			arrow_booster_amounts,
			0,                        // Lower bound
			12,                       // Upper bound
			1000000,                  // Max iterations
			100000.0,                 // Initial temperature
			0.99999                   // Cooling rate
		);
		retry_count++;
	}
	if (best_error_SA_arrow <= error_threshold) {
		std::cout << "\nOptimal Arrow Booster Configuration Found: { ";
		for (size_t i = 0; i < arrow_booster_amounts.size(); ++i) {
			std::cout << arrow_booster_amounts[i] << (i < arrow_booster_amounts.size() - 1 ? ", " : " ");
		}
		std::cout << "}\n";
		std::cout << "  Optimal Arrow Error: " << best_error_SA_arrow << "\n";

		std::cout << "Debug Arrow: " << std::endl;
		error = arrow::debug_objective(booster_toAP_x_exposures, arrow_booster_amounts, power_A_x1i_L, arrow_x17_map, target_x_L);
		std::cout << "  Error: " << error << std::endl;
	}
	else {
		std::cout << "Failed to find a configuration with error < " << error_threshold << " after " << retry_count << " retries.\n";
	}
	///////////////////////////////
	// Brute Force Optimization for Hammer
	std::vector<unsigned short> best_hammer_configuration = hammer_booster_amounts;
	double best_hammer_error = optimization::optimize_boosters_brute_force(
		hammer_objective_function,
		best_hammer_configuration,
		0,   // Lower bound
		12   // Upper bound
	);

	std::cout << "\nBrute Force Hammer Booster Configuration Found: { ";
	for (size_t i = 0; i < best_hammer_configuration.size(); ++i) {
		std::cout << best_hammer_configuration[i] << (i < best_hammer_configuration.size() - 1 ? ", " : " ");
	}
	std::cout << "}\n";
	std::cout << "  Optimal Hammer Error: " << best_hammer_error << "\n";

    std::cout << "Final Best Value: " << global_best_value << std::endl;
    std::cout << "End of computation." << std::endl;
    while (true) {

	std::cout << "\nBrute Force Arrow Booster Configuration Found: { ";
	for (size_t i = 0; i < best_arrow_configuration.size(); ++i) {
		std::cout << best_arrow_configuration[i] << (i < best_arrow_configuration.size() - 1 ? ", " : " ");
	}
	std::cout << "}\n";
	std::cout << "  Optimal Arrow Error: " << best_arrow_error << "\n";
	return 0;
}