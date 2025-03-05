#include <iostream>
#include <iomanip>
#include <cmath>
#include <sstream>
#include <vector>
#include <thread>
#include <mutex>

const int ratio_size = 13;
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
        if (distance > 8.0 || distance == 0.0) {
            return;
        }
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
double global_best_value = std::numeric_limits<double>::max(); // Start with the worst possible error

double min_height, max_height;
unsigned int final_hammer_amount;

void process_range(int threadNum, const std::vector<int>& assigned_numbers) {
    std::cout << "Thread[" << threadNum << "] assigned numbers: ";
    for (int num : assigned_numbers) {
        std::cout << num << " ";
    }
    std::cout << std::endl;

    for (int param : assigned_numbers) {
        std::cout << "Processing parameter: " << param << std::endl;

        double hx, hy, hz, hu, hv, hw;
        double sx, sy, sz, su, sv, sw;
        int ff_max, tnt_max;
        std::cout << "Enter hammer x, y, z, u, v, w: ";
        std::cin >> hx >> hy >> hz >> hu >> hv >> hw;
        std::cout << "Enter slabbust x, y, z, u, v, w: ";
        std::cin >> sx >> sy >> sz >> su >> sv >> sw;
        std::cout << "Enter final TNT amount: ";
        std::cin >> final_hammer_amount;
        std::cout << "Enter min and max height constraints: ";
        std::cin >> min_height >> max_height;
        std::cout << "Enter max freefall constraints: ";
        std::cin >> ff_max;
        std::cout << "Enter max tnt constraints: ";
		std::cin >> tnt_max;

        Tnt hammer0(hx, hy, hz, hu, hv, hw);
        Tnt slabbust0(sx, sy, sz, su, sv, sw);

        double y_vel = DBL_MAX, y_vel_record = DBL_MAX;
        short int ratio_start[ratio_size] = { 0, 0, 0, 0, 0, 0, 0 , 0, 0, 0, 0, 0, 1 };
        short int ratio[ratio_size] = { 0, 0, 0, 0, 0, 0, 0 , 0, 0, 0, 0, 0, 1 };
        short int ratio_final[ratio_size] = { 0, 1, 20, 0, 1, 20, 0, 1, 20, 0, 1, 20, 1};
        unsigned long long int counter = 0;
        double ff;
        double ff_sum, tnt_sum;
        do {
            Tnt hammer1 = hammer0;
            Tnt slabbust1 = slabbust0;

            ff = ratio[0];
            ff_sum = ff;
            tnt_sum = 0;
            hammer1.freefall(ff);
            slabbust1.freefall(ff);

            for (int i = 1; i + 2 < ratio_size; i += 3) {
                if (ratio[i] == 0) {
                    slabbust1.amount = ratio[i + 1];
                    hammer1.swing(slabbust1);
                    slabbust1.swing(slabbust1);
                }
                else {
                    hammer1.amount = ratio[i + 1];
                    slabbust1.swing(hammer1);
                    hammer1.swing(hammer1);
                }
				tnt_sum += ratio[i + 1];
                ff = ratio[i + 2];
                ff_sum += ff;
                hammer1.freefall(ff);
                slabbust1.freefall(ff);
            }

            if (hammer1.y < min_height || hammer1.y > max_height || slabbust1.y < min_height || slabbust1.y > max_height || ff_sum > ff_max || tnt_sum > tnt_max) {
                continue;
            }

            hammer1.amount = final_hammer_amount;
            slabbust1.swing(hammer1);
            y_vel = slabbust1.getVel()[1];
            if (y_vel < y_vel_record) {
                y_vel_record = y_vel;
                std::cout << "Thread[" << threadNum << "] New best vel: " << y_vel_record << " with parameters: ";
                for (int i = 0; i < ratio_size; i++) {
                    std::cout << ratio[i] << " ";
                }
                std::cout << std::endl;

                // Print sequence of actions
                std::cout << "freefall: " << ratio[0] << std::endl;
                for (int i = 1; i + 2 < ratio_size; i += 3) {
                    if (ratio[i] == 0) {
                        std::cout << "slabbust swing: " << ratio[i + 1] << std::endl;
                    }
                    else {
                        std::cout << "hammer swing: " << ratio[i + 1] << std::endl;
                    }
                    std::cout << "freefall: " << ratio[i+2] << std::endl;
                }

                // Print history
                hammer1 = hammer0;
                slabbust1 = slabbust0;

                ff = ratio[0];
                std::cout << "Initial freefall: " << ff << std::endl;
                hammer1.freefall_print(ff, " H:");
                slabbust1.freefall_print(ff, " S:");

                for (int i = 1; i + 2 < ratio_size; i += 3) {
                    if (ratio[i] == 0) {
                        std::cout << " slabbust swing: " << ratio[i + 1] << std::endl;
                        slabbust1.amount = ratio[i + 1];
                        hammer1.swing(slabbust1);
                        slabbust1.swing(slabbust1);
                    }
                    else {
                        std::cout << " hammer swing: " << ratio[i + 1] << std::endl;
                        hammer1.amount = ratio[i + 1];
                        slabbust1.swing(hammer1);
                        hammer1.swing(hammer1);
                    }
                    std::cout << " freefall: " << ratio[i + 2] << std::endl;
                    ff = ratio[i + 2];
                    hammer1.freefall_print(ff, " H:");
                    slabbust1.freefall_print(ff, " S:");
                }
                hammer1.amount = final_hammer_amount;
                slabbust1.swing(hammer1);

                hammer1.print("Hammer Position: ");
                slabbust1.print("Slabbust Position: ");
                std::cout << std::endl;
            }
            counter++;
        } while (generateNextPermutation(ratio, ratio_start, ratio_final));
    }
    std::cout << "Thread[" << threadNum << "] finished processing." << std::endl;
}
int main() {
    std::cout << std::setprecision(17);
    unsigned int thread_count = 1; // std::thread::hardware_concurrency();
    std::cout << "Available threads: " << thread_count << std::endl;

    const int start = 0, end = 0;
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

    std::cout << "Final Best Value: " << global_best_value << std::endl;
    std::cout << "End of computation." << std::endl;
    while (true) {

    }
    return 0;
}