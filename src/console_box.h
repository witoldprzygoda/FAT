/**
 * @file console_box.h
 * @brief Console output formatting utilities - boxes, frames, separators
 *
 * Provides clean, reusable functions for printing formatted console output.
 * Automatically handles padding and alignment for box frames.
 *
 * Usage:
 * @code
 *   ConsoleBox::printHeader("FAT Framework", "pp → npπ+ Analysis");
 *   ConsoleBox::printBox({"Processing complete!", "Events: 1000"});
 *   ConsoleBox::printInfo("Beam", "4500 MeV");
 * @endcode
 *
 * @author Witold Przygoda (witold.przygoda@uj.edu.pl)
 * @date 2025
 */

#ifndef CONSOLE_BOX_H
#define CONSOLE_BOX_H

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <iomanip>

/**
 * @class ConsoleBox
 * @brief Static utility class for formatted console output
 */
class ConsoleBox {
public:
    // Box styles
    enum class Style {
        DOUBLE,     // ╔═══╗ style (default for headers)
        SINGLE,     // ┌───┐ style (for info boxes)
        SIMPLE      // +---+ style (ASCII only)
    };

    // Default box width
    static constexpr int DEFAULT_WIDTH = 68;
    static constexpr int MIN_WIDTH = 20;

    // ========================================================================
    // Main Header - Large title box
    // ========================================================================
    
    /**
     * @brief Print a main program header with title and subtitle
     * @param title Main title text
     * @param subtitle Optional subtitle text
     * @param width Box width (default: 68)
     */
    static void printHeader(const std::string& title,
                           const std::string& subtitle = "",
                           int width = DEFAULT_WIDTH) {
        std::vector<std::string> lines;
        lines.push_back("");  // Empty line for spacing
        lines.push_back(title);
        if (!subtitle.empty()) {
            lines.push_back(subtitle);
        }
        lines.push_back("");  // Empty line for spacing
        
        printBox(lines, Style::DOUBLE, width);
    }

    /**
     * @brief Print a completion/status header
     * @param message Status message
     * @param width Box width (default: 68)
     */
    static void printStatus(const std::string& message, int width = DEFAULT_WIDTH) {
        printBox({message}, Style::DOUBLE, width);
    }

    // ========================================================================
    // Info Box - Single-line frame
    // ========================================================================
    
    /**
     * @brief Print an info box with single-line border
     * @param message Message to display
     * @param width Box width (default: 64)
     */
    static void printInfoBox(const std::string& message, int width = DEFAULT_WIDTH - 4) {
        printBox({message}, Style::SINGLE, width);
    }

    /**
     * @brief Print multiple lines in an info box
     * @param lines Vector of lines to display
     * @param width Box width (default: 64)
     */
    static void printInfoBox(const std::vector<std::string>& lines, int width = DEFAULT_WIDTH - 4) {
        printBox(lines, Style::SINGLE, width);
    }

    // ========================================================================
    // Generic Box Printing
    // ========================================================================
    
    /**
     * @brief Print a box with given lines and style
     * @param lines Vector of strings to display (each is centered)
     * @param style Box style (DOUBLE, SINGLE, SIMPLE)
     * @param width Total box width including borders
     */
    static void printBox(const std::vector<std::string>& lines,
                        Style style = Style::DOUBLE,
                        int width = DEFAULT_WIDTH) {
        // Get box characters based on style
        BoxChars chars = getBoxChars(style);
        
        // Ensure minimum width
        width = std::max(width, MIN_WIDTH);
        
        // Find max line length to auto-adjust width if needed
        int max_len = 0;
        for (const auto& line : lines) {
            max_len = std::max(max_len, static_cast<int>(displayWidth(line)));
        }
        
        // Auto-expand width if content is too wide (with padding)
        if (max_len + 4 > width) {
            width = max_len + 4;
        }
        
        int content_width = width - 2;  // Subtract border characters
        
        // Top border
        std::cout << chars.top_left;
        for (int i = 0; i < content_width; ++i) {
            std::cout << chars.horizontal;
        }
        std::cout << chars.top_right << "\n";
        
        // Content lines
        for (const auto& line : lines) {
            std::cout << chars.vertical;
            printCentered(line, content_width);
            std::cout << chars.vertical << "\n";
        }
        
        // Bottom border
        std::cout << chars.bottom_left;
        for (int i = 0; i < content_width; ++i) {
            std::cout << chars.horizontal;
        }
        std::cout << chars.bottom_right << "\n";
    }

    // ========================================================================
    // Separators and Simple Lines
    // ========================================================================
    
    /**
     * @brief Print a horizontal separator line
     * @param style Style for the line
     * @param width Line width (default: 68)
     */
    static void printSeparator(Style style = Style::DOUBLE, int width = DEFAULT_WIDTH) {
        BoxChars chars = getBoxChars(style);
        for (int i = 0; i < width; ++i) {
            std::cout << chars.horizontal;
        }
        std::cout << "\n";
    }

    // ========================================================================
    // Utility Functions
    // ========================================================================
    
    /**
     * @brief Print a blank line
     */
    static void newLine() {
        std::cout << "\n";
    }

    /**
     * @brief Calculate display width of a string (handles UTF-8 box-drawing chars)
     * 
     * UTF-8 box-drawing characters (like ═, ║) are 3 bytes but display as 1 char.
     * Greek letters (like π) are 2 bytes but display as 1 char.
     */
    static size_t displayWidth(const std::string& str) {
        size_t width = 0;
        size_t i = 0;
        
        while (i < str.size()) {
            unsigned char c = str[i];
            
            if ((c & 0x80) == 0) {
                // ASCII character (1 byte)
                width++;
                i++;
            } else if ((c & 0xE0) == 0xC0) {
                // 2-byte UTF-8 (e.g., Greek letters π, θ)
                width++;
                i += 2;
            } else if ((c & 0xF0) == 0xE0) {
                // 3-byte UTF-8 (e.g., box-drawing characters ═, ║)
                width++;
                i += 3;
            } else if ((c & 0xF8) == 0xF0) {
                // 4-byte UTF-8 (e.g., emojis)
                width++;
                i += 4;
            } else {
                // Invalid or continuation byte, skip
                i++;
            }
        }
        
        return width;
    }

private:
    // Box drawing character set
    struct BoxChars {
        std::string top_left;
        std::string top_right;
        std::string bottom_left;
        std::string bottom_right;
        std::string horizontal;
        std::string vertical;
    };

    static BoxChars getBoxChars(Style style) {
        switch (style) {
            case Style::DOUBLE:
                return {"╔", "╗", "╚", "╝", "═", "║"};
            case Style::SINGLE:
                return {"┌", "┐", "└", "┘", "─", "│"};
            case Style::SIMPLE:
            default:
                return {"+", "+", "+", "+", "-", "|"};
        }
    }

    /**
     * @brief Print text centered within given width
     */
    static void printCentered(const std::string& text, int width) {
        int text_width = static_cast<int>(displayWidth(text));
        int padding = width - text_width;
        
        if (padding < 0) {
            // Text too long, truncate (shouldn't happen with auto-expand)
            std::cout << text.substr(0, width);
            return;
        }
        
        int left_pad = padding / 2;
        int right_pad = padding - left_pad;
        
        // Print left padding
        for (int i = 0; i < left_pad; ++i) std::cout << ' ';
        
        // Print text
        std::cout << text;
        
        // Print right padding
        for (int i = 0; i < right_pad; ++i) std::cout << ' ';
    }
};

#endif // CONSOLE_BOX_H
