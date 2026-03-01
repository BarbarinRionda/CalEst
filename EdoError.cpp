// EdoError.cpp

#include <iostream>
#include <stdexcept>

class EdoError : public std::runtime_error {
public:
    explicit EdoError(const std::string &message)
        : std::runtime_error("Edo Error: " + message) {}
};

void handleError(const std::string &errorMessage) {
    throw EdoError(errorMessage);
}

int main() {
    try {
        handleError("An error occurred during processing.");
    } catch (const EdoError &e) {
        std::cerr << e.what() << std::endl;
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}