#include <iostream>
#include <iomanip>

void progressbar(Long64_t current, Long64_t total, int barWidth = 100) {
    double progress = (double)current / total;
    int pos = barWidth * progress;

    std::cout << "[";
    for (int i = 0; i < barWidth; ++i) {
        if (i < pos) std::cout << "█";
        else std::cout << "░";
    }
    std::cout << "] " << int(progress * 100.0) << " %\r";
    std::cout.flush();
}

