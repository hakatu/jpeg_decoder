#include <iostream>
#include <cstring>  // For std::memcpy
#include "idct_fast.h"
#include "idct_aan.h"

int main() {
    // Initialize an 8x8 block with only the DC coefficient set
    int data_in[64] = {0};
    data_in[0] = 1024;  // DC coefficient

    // Output arrays
    int data_out_fast[64];
    int data_out_aan[64];

    // Temporary input buffers
    int temp_in_fast[64];
    int temp_in_aan[64];

    // Copy original data for idct_fast
    std::memcpy(temp_in_fast, data_in, sizeof(data_in));
    idct_fast idct_fast;
    idct_fast.process(temp_in_fast, data_out_fast);

    // Copy original data for idct_aan
    std::memcpy(temp_in_aan, data_in, sizeof(data_in));
    idct_aan idct_aan;
    idct_aan.process(temp_in_aan, data_out_aan);

    // Print outputs
    std::cout << "idct_fast output:" << std::endl;
    for (int i = 0; i < 64; ++i) {
        std::cout << data_out_fast[i] << " ";
        if ((i + 1) % 8 == 0) std::cout << std::endl;
    }

    std::cout << "\nidct_aan output:" << std::endl;
    for (int i = 0; i < 64; ++i) {
        std::cout << data_out_aan[i] << " ";
        if ((i + 1) % 8 == 0) std::cout << std::endl;
    }

    // Verify uniformity
    bool is_uniform_fast = true;
    int first_value_fast = data_out_fast[0];
    for (int i = 1; i < 64; ++i) {
        if (data_out_fast[i] != first_value_fast) {
            is_uniform_fast = false;
            break;
        }
    }

    bool is_uniform_aan = true;
    int first_value_aan = data_out_aan[0];
    for (int i = 1; i < 64; ++i) {
        if (data_out_aan[i] != first_value_aan) {
            is_uniform_aan = false;
            break;
        }
    }

    // Report results
    std::cout << "\nVerification:" << std::endl;
    if (is_uniform_fast) {
        std::cout << "idct_fast output is uniform with value: " << first_value_fast << std::endl;
    } else {
        std::cout << "idct_fast output is NOT uniform!" << std::endl;
    }

    if (is_uniform_aan) {
        std::cout << "idct_aan output is uniform with value: " << first_value_aan << std::endl;
    } else {
        std::cout << "idct_aan output is NOT uniform!" << std::endl;
    }

    return 0;
}