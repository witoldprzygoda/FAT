/**
 * @file progressbar.h
 * @brief Progress bar with time estimation and graceful termination
 *
 * @author Witold Przygoda (witold.przygoda@uj.edu.pl)
 * @date 2025
 */

#ifndef PROGRESSBAR_H
#define PROGRESSBAR_H

#include <iostream>
#include <iomanip>
#include <chrono>
#include <string>
#include <sstream>
#include <csignal>
#include <atomic>
#include <Rtypes.h>

// ============================================================================
// Global signal handling for graceful termination
// ============================================================================

namespace SignalHandler {
    // Atomic flag for thread-safe signal handling
    inline std::atomic<bool> g_interrupted{false};
    
    /**
     * @brief Check if termination was requested (via Ctrl+C)
     */
    inline bool wasInterrupted() {
        return g_interrupted.load();
    }
    
    /**
     * @brief Reset the interrupted flag
     */
    inline void reset() {
        g_interrupted.store(false);
    }
    
    /**
     * @brief Internal signal handler
     */
    inline void signalHandler(int signum) {
        if (signum == SIGINT) {
            g_interrupted.store(true);
            std::cout << "\n\n*** Ctrl+C detected - finishing current event and saving results... ***\n";
        }
    }
    
    /**
     * @brief Install signal handler for graceful termination
     * Call this at the beginning of main()
     */
    inline void install() {
        std::signal(SIGINT, signalHandler);
        reset();
    }
}

// ============================================================================
// ProgressBar class with time estimation
// ============================================================================

/**
 * @brief Progress bar class with time estimation
 *
 * Usage:
 *   ProgressBar progress(total_events);
 *   for (Long64_t i = 0; i < total_events; ++i) {
 *       if (SignalHandler::wasInterrupted()) break;  // Graceful exit
 *       progress.update(i + 1);
 *       // ... processing ...
 *   }
 *   progress.finish();
 */
class ProgressBar {
private:
    Long64_t total_;
    int bar_width_;
    std::chrono::steady_clock::time_point start_time_;
    double last_update_percent_;
    
    std::string formatTime(double seconds) const {
        if (seconds < 0) return "--:--";
        
        int mins = static_cast<int>(seconds) / 60;
        int secs = static_cast<int>(seconds) % 60;
        
        std::ostringstream oss;
        if (mins >= 60) {
            int hours = mins / 60;
            mins = mins % 60;
            oss << hours << ":" << std::setfill('0') << std::setw(2) << mins 
                << ":" << std::setw(2) << secs;
        } else {
            oss << std::setfill('0') << std::setw(2) << mins 
                << ":" << std::setw(2) << secs;
        }
        return oss.str();
    }
    
public:
    /**
     * @brief Construct progress bar
     * @param total Total number of items to process
     * @param bar_width Width of the progress bar in characters (default: 50)
     */
    ProgressBar(Long64_t total, int bar_width = 50) 
        : total_(total), bar_width_(bar_width), last_update_percent_(-1.0) {
        start_time_ = std::chrono::steady_clock::now();
    }
    
    /**
     * @brief Update progress bar
     * @param current Current item number (1-based)
     * 
     * Updates display:
     * - Every 0.1% for first 1%
     * - Every 1% after that
     */
    void update(Long64_t current) {
        if (total_ <= 0) return;
        
        double progress = static_cast<double>(current) / total_;
        double percent = progress * 100.0;
        
        // Determine update threshold based on progress
        double threshold;
        if (percent < 1.0) {
            threshold = 0.1;  // Update every 0.1% for first 1%
        } else {
            threshold = 1.0;  // Update every 1% after that
        }
        
        // Check if we should update (crossed a threshold)
        double last_bucket = std::floor(last_update_percent_ / threshold);
        double current_bucket = std::floor(percent / threshold);
        
        if (current_bucket == last_bucket && last_update_percent_ >= 0) return;
        last_update_percent_ = percent;
        
        int pos = static_cast<int>(bar_width_ * progress);
        
        // Calculate elapsed time
        auto now = std::chrono::steady_clock::now();
        double elapsed = std::chrono::duration<double>(now - start_time_).count();
        
        // Calculate ETA (start from 0.1% progress)
        std::string eta_str;
        if (percent >= 0.1 && progress > 0.0001) {
            double total_estimated = elapsed / progress;
            double remaining = total_estimated - elapsed;
            eta_str = formatTime(remaining);
        } else {
            eta_str = "--:--";
        }
        
        // Build progress bar
        std::cout << "\r[";
        for (int i = 0; i < bar_width_; ++i) {
            if (i < pos) std::cout << "█";
            else std::cout << "░";
        }
        
        // Show percent with decimal for early progress
        if (percent < 1.0) {
            std::cout << "] " << std::fixed << std::setprecision(1) << std::setw(4) << percent 
                     << "% | ETA: " << eta_str << "    ";
        } else {
            std::cout << "] " << std::setw(3) << static_cast<int>(percent) 
                     << "% | ETA: " << eta_str << "    ";
        }
        std::cout.flush();
    }
    
    /**
     * @brief Finish progress bar (prints newline and elapsed time)
     * @param interrupted If true, indicates early termination
     */
    void finish(bool interrupted = false) {
        auto now = std::chrono::steady_clock::now();
        double elapsed = std::chrono::duration<double>(now - start_time_).count();
        
        // Complete the bar (or show partial if interrupted)
        std::cout << "\r[";
        if (interrupted) {
            double progress = last_update_percent_ / 100.0;
            int pos = static_cast<int>(bar_width_ * progress);
            for (int i = 0; i < bar_width_; ++i) {
                if (i < pos) std::cout << "█";
                else std::cout << "░";
            }
            std::cout << "] " << static_cast<int>(last_update_percent_) 
                     << "% | Interrupted after: " << formatTime(elapsed) << "    \n";
        } else {
            for (int i = 0; i < bar_width_; ++i) {
                std::cout << "█";
            }
            std::cout << "] 100% | Done in: " << formatTime(elapsed) << "    \n";
        }
    }
    
    /**
     * @brief Reset progress bar for reuse
     * @param new_total New total count (optional, keeps old if 0)
     */
    void reset(Long64_t new_total = 0) {
        if (new_total > 0) total_ = new_total;
        start_time_ = std::chrono::steady_clock::now();
        last_update_percent_ = -1.0;
    }
    
    /**
     * @brief Get elapsed time in seconds
     */
    double elapsed() const {
        auto now = std::chrono::steady_clock::now();
        return std::chrono::duration<double>(now - start_time_).count();
    }
};

// ============================================================================
// Legacy function interface (for backwards compatibility)
// ============================================================================

/**
 * @brief Simple progress bar function (no time estimation)
 * @deprecated Use ProgressBar class instead for time estimation
 */
inline void progressbar(Long64_t current, Long64_t total, int barWidth = 50) {
    if (total <= 0) return;
    
    double progress = static_cast<double>(current) / total;
    int pos = static_cast<int>(barWidth * progress);

    std::cout << "\r[";
    for (int i = 0; i < barWidth; ++i) {
        if (i < pos) std::cout << "█";
        else std::cout << "░";
    }
    std::cout << "] " << std::setw(3) << static_cast<int>(progress * 100) << "%";
    std::cout.flush();
}

#endif // PROGRESSBAR_H
