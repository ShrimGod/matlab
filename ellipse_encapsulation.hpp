#ifndef ELLIPSE_ENCAPSULATION_HPP
#define ELLIPSE_ENCAPSULATION_HPP

#include "Options.hpp"
#include "pipelinemonitor.hpp"
#include "PeriodicClusterContainer.hpp"
#include "FrequencyContainer.hpp"

#include <vector>
#include <string>
#include <utility>
#include <tuple>
#include <cmath>
#include <sstream>

class EllipseEncapsulation {
public:
    EllipseEncapsulation(std::shared_ptr<EBSOptions> o, std::shared_ptr<PeriodicClusterContainer> pcc);
    ~EllipseEncapsulation() {
        std::cout << "EllipseEncapsulation destructor called." << std::endl;
    }

    void start(StageConnector connector);

private:
    std::shared_ptr<EBSOptions> options;
    std::shared_ptr<PeriodicClusterContainer> periodicClusterStream;
    
    void computeMVEE(std::vector<std::pair<float, float>>& points, float& center_x, float& center_y, float& major_radius, float& minor_radius,
                     float& rot1, float& rot2, float& rot3, float& rot4,float tolerance);
    
    // Process clusters and compute ellipse encapsulation for each
    void encapsulateEllipse(const std::vector<std::reference_wrapper<const PeriodicCluster>>& clusterVector);
};

#endif