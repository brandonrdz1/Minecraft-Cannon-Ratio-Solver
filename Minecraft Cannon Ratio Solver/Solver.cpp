#include <iostream>
#include <iomanip>
#include <cmath>

double objective_function(double* booster_exposure_x, unsigned short int* booster_amount, double hammer_x0) {
	// not yet goning to return minimum or max
	double result = hammer_x0;
	for (int i = 0; i < 12; i++) {
		result += booster_amount[i] * booster_exposure_x[i];
	}
	return result;
}

int main() {

	std::cout << std::setprecision(17);
	//// In-game Conditions
	double booster_x_values[12] = { 
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
	double power_H_x0_value = -86228.49000000954;	// -86228.49000000954 278.06125000119226 87930.5
	double power_H_u0_value = 0.0;					// 0.0 13.581624504779418 -2.4523938435550006E-15

	// State after constant booster only
	double power_H_x1_value = -86232.08075066289;	// -86232.08075066289 316.3339737459176 87930.5
	double power_H_u1_value = -3.5189356402947953;		// -3.5189356402947953 37.507269269830815 -1.158739770801276E-14
	
	// Function parameters
	double booster_x_exposure[12]; // Is also the distance it pushes the hammer along x
	std::cout << "Booster Exposure: " << std::endl;
	for (int i = 0; i < 12; i++) {
		double distance = power_H_x0_value - booster_x_values[i];
		double f = (1.0 - distance / 8.0) * 1.0 / distance;
		booster_x_exposure[i] = distance * f;
		std::cout << "  booster[" << i << "]: " << booster_x_exposure[i] << std::endl;
	}

	double power_H_x0 = power_H_x1_value; // Virtual Starting location of the hammer
	
	// Result?
	unsigned short int booster_amount[12] = { 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6 }; // Limited to Between 0-12
	double result = objective_function(booster_x_exposure, booster_amount, power_H_x0);
	
	std::cout << "Hammer's Power location: " << result << std::endl;
	std::cout << "End!" << std::endl;
	return 0;
}


// all 6
// -86219.89860780213 316.33397374591794 87930.5