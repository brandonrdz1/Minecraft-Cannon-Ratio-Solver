#include <iostream>
#include <iomanip>
#include <cmath>
#include <sstream>
#include <vector>
#include <thread>
#include <mutex>

const int ratio_size = 8;
class Tnt {
public:
    double x = 0.0, y = 0.0, z = 0.0;
    double u = 0.0, v = 0.0, w = 0.0;
    unsigned int amount = 1;

    static const double explosion_height, gravity, drag;
    static const double distance(const Tnt& a, const Tnt& b) {
        double dx = a.x - b.x;
        double dy = a.y - b.y;
        double dz = a.z - b.z;
        return std::sqrt(dx * dx + dy * dy + dz * dz);
    }

    Tnt(const std::string& attributes, unsigned int explosivePower = 1.0) : amount(explosivePower) {
        std::istringstream stream(attributes);
        stream >> x >> y >> z >> u >> v >> w;
    };

    Tnt(double x, double y, double z, double u = 0.0, double v = 0.0, double w = 0.0, unsigned int explosivePower = 1.0) : x(x), y(y), z(z), u(u), v(v), w(w), amount(explosivePower) {};

    void explosion(Tnt source, double exposure = 1.0) {
        double dx = x - source.x;
        double dy = y - source.y; // double dy = y - (source.y + explosion_height);
        double dz = z - source.z;
        float distance = std::sqrt(dx * dx + dy * dy + dz * dz); // double
        double f = (-1.0 / 8.0 + 1.0 / distance) * source.amount * exposure;
        u += f * dx;
        v += f * dy;
        w += f * dz;
    }

    void swing(Tnt source) {
        Tnt origin(source);
        std::vector<std::vector<double>> locations(source.amount);
        double dx, dy, dz;
        float distance;
        double f;
        for (int i = 0; i < source.amount; i++) {
            Tnt next = origin;
            next.freefall(1);
            if (next.y > 254.0199999809265) { next.y = 254.0199999809265; }
            std::vector<double> loc(3);
            loc[0] = next.x;
            loc[1] = next.y;
            loc[2] = next.z;
            locations[i]=loc;
            dx = origin.x - locations[i][0], dy = origin.y - locations[i][1], dz = origin.z - locations[i][2];
            distance = std::sqrt(dx * dx + dy * dy + dz * dz);
            if (distance > 8.0 || distance == 0.0) {
                continue;
            }
            f = -1.0 / 8.0 + 1.0 / distance;
            origin.u += dx * f;
            origin.v += dy * f;
            origin.w += dz * f;
        }

        for (int i = 0; i < locations.size(); i++) {
            Tnt source(locations[i][0], locations[i][1], locations[i][2]);
            //source.print("Source[" + std::to_string(i) + "]: ");
            this->explosion(source);
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
    void freefall_print(int ticks, std::string str = "") {
        print(str+ " 0 ");
        for (int i = 0; i < ticks; i++) {
            v += Tnt::gravity;

            x += u;
            y += v;
            z += w;
            print(str + " " + std::to_string(i + 1) + " ");
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

    void printPos(std::string str = "") {
        std::cout << str << "x: " << x << "\ty: " << y << "\tz: " << z << std::endl;
    }
};
// Define static const members outside the class
const double Tnt::explosion_height = 0.061250001192092896;
const double Tnt::gravity = -0.04f;
const double Tnt::drag = 0.98f;

// Hash function for the vector
struct VectorHash {
    std::size_t operator()(const std::vector<int>& vec) const {
        std::size_t seed = 0;
        for (int num : vec) {
            seed ^= std::hash<int>()(num) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        }
        return seed;
    }
};

bool generateNextPermutation(short int currentArray[], const short int startArray[], const short int finalArray[]) {
    size_t i = ratio_size;

    while (i > 0) {
        --i;
        if (currentArray[i] < finalArray[i]) {
            ++currentArray[i];
            for (size_t j = i + 1; j < ratio_size; ++j) {
                currentArray[j] = startArray[j]; // Reset subsequent values to their start
            }
            return true;
        }
    }
    return false;
}

double fastff(double vel) {
    double start = 253.97999998182058;
    vel += Tnt::gravity;
    start += vel;
    vel *= Tnt::drag;
    vel += Tnt::gravity;
    start += vel;
    vel *= Tnt::drag;
    vel += Tnt::gravity;
    start += vel;
    vel *= Tnt::drag;
    vel += Tnt::gravity;
    start += vel;
    vel *= Tnt::drag;
    vel += Tnt::gravity;
    start += vel;
    vel *= Tnt::drag;
    vel += Tnt::gravity;
    start += vel;
    vel *= Tnt::drag;
    vel += Tnt::gravity;
    start += vel;
    return start;
}

std::mutex output_mutex;
double global_best_error = std::numeric_limits<double>::max(); // Start with the worst possible error

void process_range(int threadNum, const std::vector<int>& assigned_numbers) {
    std::cout << "Thread[" << threadNum << "] assigned numbers: ";
    unsigned long long int outputModNumber = 0;
    for (int num : assigned_numbers) {
        std::cout << num << " ";
    }
    std::cout << std::endl;

    for (int param : assigned_numbers) {
		// Reset/initialize has map for each thread
        std::cout << "Processing parameter: " << param << std::endl;
		std::this_thread::sleep_for(std::chrono::seconds(3)); // Add delay to stagger threads
        
        Tnt swing0("-200001.4037412894 229.0199999809265 -333579.5 -0.4335096811347792 0.0 0.0");


        double decimal = 0.0;
        int forward = 0;
        short int ratio_start[ratio_size] = { param, 0, 0, 0, 0, 0, 0, 0 };
        short int ratio_curr[ratio_size] = { param, 0, 0, 0, 0, 0, 0, 0 };
        short int ratio_final[ratio_size] = { param, 12, 12, 12, 12, 12, 12, 12 };
		unsigned long long int counter = 0;
        do {
            counter++;
            if (counter % 1000000 == 0) {
                std::cout << "Thread[" << threadNum << "] | is at: ";
                for (int i = 0; i < ratio_size; i++) {
                    std::cout << ratio_curr[i] << (i < ratio_size - 1 ? ", " : " ");
                }
                std::cout << std::endl;
            }
            Tnt swing1 = swing0;
            swing1.freefall(ratio_curr[0]);
			swing1.amount = ratio_curr[1];
            swing1.swing(swing1);
            if (swing1.x > -200000) { continue; }

			swing1.freefall(ratio_curr[2]);
			swing1.amount = ratio_curr[3];
            swing1.swing(swing1);
            if (swing1.x > -200000) { continue; }

			swing1.freefall(ratio_curr[4]);
			swing1.amount = ratio_curr[5];
            swing1.swing(swing1);
            if (swing1.x > -200000) { continue; }

			swing1.freefall(ratio_curr[6]);
			swing1.amount = ratio_curr[7];
            swing1.swing(swing1);
            if (swing1.x > -200000) { continue; }

            decimal = swing1.x - floor(swing1.x);
            forward = (swing1.x * -1) - 200000;
            if (forward == 3 && decimal < 0.51f && decimal > 0.49f && swing1.y > 229.0 && swing1.y < 230.0 && abs(swing1.u) < 0.01) {
                std::cout << "Thread[" << threadNum << "] | ";
                for (int i = 0; i < ratio_size; i++) {
                    std::cout << ratio_curr[i] << (i < ratio_size - 1 ? ", " : " ");
                }
                swing1.print();

                Tnt swing2 = swing0;

                std::cout << "  FF0: " << ratio_curr[0] << std::endl;
                swing2.freefall_print(ratio_curr[0]);
                std::cout << "  Amt1: "  << ratio_curr[1] << std::endl;
                swing2.amount = ratio_curr[1];
                swing2.swing(swing2);
                
                std::cout << "  FF2: " << ratio_curr[2] << std::endl;
                swing2.freefall_print(ratio_curr[2]);
                std::cout << "  Amt3: " << ratio_curr[3] << std::endl;
                swing2.amount = ratio_curr[3];
                swing2.swing(swing2);

                std::cout << "  FF4: " << ratio_curr[4] << std::endl;
                swing2.freefall_print(ratio_curr[4]);
                std::cout << "  Amt5: " << ratio_curr[5] << std::endl;
                swing2.amount = ratio_curr[5];
                swing2.swing(swing2);

                std::cout << "  FF6: " << ratio_curr[6] << std::endl;
                swing2.freefall_print(ratio_curr[6]);
                std::cout << "  Amt7: " << ratio_curr[7] << std::endl;
                swing2.amount = ratio_curr[7];
                swing2.swing(swing2);

                swing2.print(" Final State Vector: ");
                std::cout << std::endl;
            }

        } while (generateNextPermutation(ratio_curr, ratio_start, ratio_final));
    }

    std::cout << "Thread[" << threadNum << "] finished processing." << std::endl;
}

int main() {
    std::cout << std::setprecision(17);
    unsigned int thread_count = std::thread::hardware_concurrency();
    std::cout << "Available threads: " << thread_count << std::endl;

    const int start = 0, end = 12;
    std::vector<std::thread> threads;
    std::vector<std::vector<int>> thread_assignments(thread_count);

    int current_thread = 0;
    for (int i = start; i <= end; i ++) {
        thread_assignments[current_thread].push_back(i);
        current_thread = (current_thread + 1) % thread_count;
    }

    for (unsigned int i = 0; i < thread_count; i++) {
        threads.emplace_back([i, &thread_assignments]() {
            std::this_thread::sleep_for(std::chrono::seconds(1) * i * 0.05); // Add delay to stagger threads
            process_range(i, thread_assignments[i]);
            });
    }

    for (auto& thread : threads) {
        thread.join();
    }

    std::cout << "Final Best Error: " << global_best_error << std::endl;
    std::cout << "End of computation." << std::endl;
    while (true) {

    }
    return 0;
}