#include "ellipse_encapsulation.hpp"
#include <vector>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <limits>

// Matrix operations for 2D and 3D matrices needed for MVEE
struct Matrix2D {
    float data[2][2];
    
    Matrix2D() {
        data[0][0] = data[0][1] = data[1][0] = data[1][1] = 0.0f;
    }
    
    Matrix2D(float a, float b, float c, float d) {
        data[0][0] = a; data[0][1] = b;
        data[1][0] = c; data[1][1] = d;
    }
    
    float& operator()(int i, int j) { return data[i][j]; }
    const float& operator()(int i, int j) const { return data[i][j]; }
    
    // Matrix multiplication
    Matrix2D operator*(const Matrix2D& other) const {
        Matrix2D result;
        for(int i = 0; i < 2; i++) {
            for(int j = 0; j < 2; j++) {
                result(i,j) = 0;
                for(int k = 0; k < 2; k++) {
                    result(i,j) += data[i][k] * other(k,j);
                }
            }
        }
        return result;
    }
    
    // Matrix inversion for 2x2
    Matrix2D inverse() const {
        float det = data[0][0] * data[1][1] - data[0][1] * data[1][0];
        if(std::abs(det) < 1e-10) {
            // Singular matrix - return identity as fallback
            return Matrix2D(1, 0, 0, 1);
        }
        return Matrix2D(data[1][1]/det, -data[0][1]/det, -data[1][0]/det, data[0][0]/det);
    }
    
    // SVD for 2x2 symmetric positive definite matrices
    void svd(Matrix2D& U, float eigenvals[2]) const {
        // For symmetric 2x2 matrix, compute eigenvalues and eigenvectors
        float trace = data[0][0] + data[1][1];
        float det = data[0][0] * data[1][1] - data[0][1] * data[1][0];
        
        float lambda1 = (trace + std::sqrt(trace*trace - 4*det)) / 2.0f;
        float lambda2 = (trace - std::sqrt(trace*trace - 4*det)) / 2.0f;
        
        eigenvals[0] = std::max(lambda1, lambda2);
        eigenvals[1] = std::min(lambda1, lambda2);
        
        // Compute eigenvectors
        if(std::abs(data[0][1]) > 1e-10) {
            float v1_x = data[0][1];
            float v1_y = eigenvals[0] - data[0][0];
            float norm1 = std::sqrt(v1_x*v1_x + v1_y*v1_y);
            
            float v2_x = data[0][1]; 
            float v2_y = eigenvals[1] - data[0][0];
            float norm2 = std::sqrt(v2_x*v2_x + v2_y*v2_y);
            
            U(0,0) = v1_x / norm1; U(0,1) = v2_x / norm2;
            U(1,0) = v1_y / norm1; U(1,1) = v2_y / norm2;
        } else {
            // Diagonal matrix case
            U = Matrix2D(1, 0, 0, 1);
        }
    }
};

struct Matrix3D {
    float data[3][3];
    
    Matrix3D() {
        for(int i = 0; i < 3; i++) {
            for(int j = 0; j < 3; j++) {
                data[i][j] = 0.0f;
            }
        }
    }
    
    float& operator()(int i, int j) { return data[i][j]; }
    const float& operator()(int i, int j) const { return data[i][j]; }
    
    // Matrix inversion for 3x3 using Cramer's rule
    Matrix3D inverse() const {
        Matrix3D result;
        float det = determinant();
        
        if(std::abs(det) < 1e-10) {
            // Singular matrix - return identity as fallback
            result(0,0) = result(1,1) = result(2,2) = 1.0f;
            return result;
        }
        
        // Calculate adjugate matrix elements and divide by determinant
        result(0,0) = (data[1][1]*data[2][2] - data[1][2]*data[2][1]) / det;
        result(0,1) = (data[0][2]*data[2][1] - data[0][1]*data[2][2]) / det;
        result(0,2) = (data[0][1]*data[1][2] - data[0][2]*data[1][1]) / det;
        
        result(1,0) = (data[1][2]*data[2][0] - data[1][0]*data[2][2]) / det;
        result(1,1) = (data[0][0]*data[2][2] - data[0][2]*data[2][0]) / det;
        result(1,2) = (data[0][2]*data[1][0] - data[0][0]*data[1][2]) / det;
        
        result(2,0) = (data[1][0]*data[2][1] - data[1][1]*data[2][0]) / det;
        result(2,1) = (data[0][1]*data[2][0] - data[0][0]*data[2][1]) / det;
        result(2,2) = (data[0][0]*data[1][1] - data[0][1]*data[1][0]) / det;
        
        return result;
    }
    
private:
    float determinant() const {
        return data[0][0] * (data[1][1]*data[2][2] - data[1][2]*data[2][1]) -
               data[0][1] * (data[1][0]*data[2][2] - data[1][2]*data[2][0]) +
               data[0][2] * (data[1][0]*data[2][1] - data[1][1]*data[2][0]);
    }
};

EllipseEncapsulation::EllipseEncapsulation(std::shared_ptr<EBSOptions> o, std::shared_ptr<PeriodicClusterContainer> pcc) :
    options(o),
    periodicClusterStream(pcc) {
}

void EllipseEncapsulation::computeMVEE(std::vector<std::pair<float, float>>& points, 
                                       float& center_x, float& center_y,
                                       float& major_radius, float& minor_radius,
                                       float& rot1, float& rot2, float& rot3, float& rot4,
                                       float tolerance) {
    
    const int N = points.size();
    const int d = 2; // 2D points
    
    // Handle degenerate cases
    if(N < 3) {
        std::cerr << "Warning: Cannot compute MVEE with fewer than 3 points. Using default ellipse." << std::endl;
        center_x = center_y = 0.0f;
        major_radius = minor_radius = 1.0f;
        rot1 = rot4 = 1.0f; rot2 = rot3 = 0.0f; // Identity rotation
        return;
    }
    
    // TODO: Add convex hull optimization here for performance
    // For now, using all points directly
    
    // Create augmented matrix Q: each column is [x_i, y_i, 1]^T
    std::vector<std::vector<float>> Q(d + 1, std::vector<float>(N));
    for(int i = 0; i < N; i++) {
        Q[0][i] = points[i].first;   // x coordinate
        Q[1][i] = points[i].second;  // y coordinate  
        Q[2][i] = 1.0f;              // homogeneous coordinate
    }
    
    // Initialize uniform weights (Khachiyan algorithm step 1)
    std::vector<float> u(N, 1.0f / N);
    
    // Khachiyan iterative algorithm parameters
    const int max_iter = 100; // TODO: Could be moved to configuration
    int iter = 0;
    float err = 1.0f;
    
    std::cout << "Starting MVEE computation with " << N << " points, tolerance = " << tolerance << std::endl;
    
    // Main Khachiyan iteration loop
    while(err > tolerance && iter < max_iter) {
        // Compute weighted covariance matrix X = Q * diag(u) * Q^T
        Matrix3D X;
        for(int i = 0; i < d + 1; i++) {
            for(int j = 0; j < d + 1; j++) {
                X(i,j) = 0.0f;
                for(int k = 0; k < N; k++) {
                    X(i,j) += Q[i][k] * u[k] * Q[j][k];
                }
            }
        }
        
        // Compute M = diagonal of Q^T * X^(-1) * Q  
        Matrix3D X_inv = X.inverse();
        std::vector<float> M(N);
        for(int k = 0; k < N; k++) {
            M[k] = 0.0f;
            for(int i = 0; i < d + 1; i++) {
                for(int j = 0; j < d + 1; j++) {
                    M[k] += Q[i][k] * X_inv(i,j) * Q[j][k];
                }
            }
        }
        
        // Find the point that violates the ellipse constraint most
        int j = 0;
        float maximum = M[0];
        for(int k = 1; k < N; k++) {
            if(M[k] > maximum) {
                maximum = M[k];
                j = k;
            }
        }
        
        // Compute optimal step size for weight update
        float step_size = (maximum - d - 1) / ((d + 1) * (maximum - 1));
        
        // Update weights: reduce all weights and increase weight of maximally violating point
        std::vector<float> new_u(N);
        for(int k = 0; k < N; k++) {
            new_u[k] = (1 - step_size) * u[k];
        }
        new_u[j] += step_size;
        
        // Check convergence based on weight change
        err = 0.0f;
        for(int k = 0; k < N; k++) {
            float diff = new_u[k] - u[k];
            err += diff * diff;
        }
        err = std::sqrt(err);
        
        u = new_u;
        iter++;
    }
    
    if(iter >= max_iter) {
        std::cout << "Warning: MVEE algorithm reached maximum iterations (" << max_iter << ") without convergence" << std::endl;
    } else {
        std::cout << "MVEE converged after " << iter << " iterations with error " << err << std::endl;
    }
    
    // Extract ellipse parameters from final solution
    
    // Compute weighted mean (ellipse center)
    center_x = center_y = 0.0f;
    for(int k = 0; k < N; k++) {
        center_x += u[k] * points[k].first;
        center_y += u[k] * points[k].second;
    }
    
    // Compute covariance matrix P * U * P^T - (P * u) * (P * u)^T
    Matrix2D P_U_Pt, center_outer;
    
    // P * U * P^T term  
    for(int i = 0; i < 2; i++) {
        for(int j = 0; j < 2; j++) {
            P_U_Pt(i,j) = 0.0f;
            for(int k = 0; k < N; k++) {
                P_U_Pt(i,j) += Q[i][k] * u[k] * Q[j][k];
            }
        }
    }
    
    // (P * u) * (P * u)^T term (outer product of center)
    center_outer(0,0) = center_x * center_x;
    center_outer(0,1) = center_x * center_y;
    center_outer(1,0) = center_y * center_x;
    center_outer(1,1) = center_y * center_y;
    
    // A = (1/d) * inverse(P * U * P^T - center_outer)
    Matrix2D cov_matrix;
    cov_matrix(0,0) = P_U_Pt(0,0) - center_outer(0,0);
    cov_matrix(0,1) = P_U_Pt(0,1) - center_outer(0,1);
    cov_matrix(1,0) = P_U_Pt(1,0) - center_outer(1,0);
    cov_matrix(1,1) = P_U_Pt(1,1) - center_outer(1,1);
    
    Matrix2D A = cov_matrix.inverse();
    // Scale by 1/d as in MATLAB code
    A(0,0) /= d; A(0,1) /= d; A(1,0) /= d; A(1,1) /= d;
    
    // Perform SVD to get rotation matrix and radii
    Matrix2D rotation;
    float eigenvals[2];
    A.svd(rotation, eigenvals);
    
    // Extract ellipse parameters
    major_radius = 1.0f / std::sqrt(eigenvals[0]);
    minor_radius = 1.0f / std::sqrt(eigenvals[1]);
    
    // Store rotation matrix elements
    rot1 = rotation(0,0); rot2 = rotation(0,1);
    rot3 = rotation(1,0); rot4 = rotation(1,1);
    
    std::cout << "MVEE result: center=(" << center_x << "," << center_y 
              << "), radii=(" << major_radius << "," << minor_radius << ")" << std::endl;
}

void EllipseEncapsulation::encapsulateEllipse(const std::vector<std::reference_wrapper<const PeriodicCluster>>& clusterVector) {
    std::cout << "Starting ellipse encapsulation for " << clusterVector.size() << " clusters" << std::endl;
    
    for (const auto& clusterRef : clusterVector) {
        try {
            // Get mutable reference to cluster for writing results
            auto& cluster = const_cast<PeriodicCluster&>(clusterRef.get());
            
            // Extract pixel coordinates from cluster's frequency array
            std::vector<std::pair<float, float>> points;
            for (const auto& freqPtr : cluster.freqArray) {
                if (freqPtr) {
                    points.emplace_back(static_cast<float>(freqPtr->x), static_cast<float>(freqPtr->y));
                } else {
                    std::cerr << "Warning: Null frequency pointer in cluster " << cluster.clusterID << std::endl;
                }
            }
            
            if (points.empty()) {
                std::cerr << "Warning: No valid points found in cluster " << cluster.clusterID << std::endl;
                continue;
            }
            
            std::cout << "Processing cluster " << cluster.clusterID << " with " << points.size() << " points" << std::endl;
            
            // Compute MVEE using appropriate tolerance
            float tolerance;
            if (cluster.classification == "blinking" || cluster.classification == "rotating" || cluster.classification == "vibrating") {
                // This is a periodic cluster
                tolerance = options->periodic_ellipse_enc.tol;
            } else {
                // This is likely a non-periodic cluster, need to refine, especially for other case
                tolerance = options->nonperiodic_ellipse_enc.tol;
            }
            
            float center_x, center_y, major_rad, minor_rad;
            float rot1, rot2, rot3, rot4;
            
            computeMVEE(points, center_x, center_y, major_rad, minor_rad, 
                       rot1, rot2, rot3, rot4, tolerance);
            
            // Store results in cluster structure
            cluster.x_c = static_cast<int>(std::round(center_x));
            cluster.y_c = static_cast<int>(std::round(center_y));
            cluster.majorRadius = static_cast<int>(std::round(major_rad));
            cluster.minorRadius = static_cast<int>(std::round(minor_rad));
            cluster.rotMat1 = rot1;
            cluster.rotMat2 = rot2;
            cluster.rotMat3 = rot3;
            cluster.rotMat4 = rot4;
            
            std::cout << "Cluster " << cluster.clusterID << " ellipse: center=(" 
                      << cluster.x_c << "," << cluster.y_c << "), radii=(" 
                      << cluster.majorRadius << "," << cluster.minorRadius << ")" << std::endl;
                      
        } catch (std::exception& e) {
            std::ostringstream msg;
            msg << "Error during ellipse encapsulation: " << e.what();
            throw std::runtime_error(msg.str());
        }
    }
}

void EllipseEncapsulation::start(StageConnector connector) {
    std::cout << "Starting Ellipse Encapsulation" << std::endl;
    
    // Wait for initial data from previous stage
    while(connector.leftLoc->getValue() == 0) {
        // Wait for data
    }

    uint64_t currentTimeBin = options->freq.time_bin_duration;

    // Process data until previous stage signals completion
    while(connector.leftLoc->getValue() > connector.rightLoc->getValue() || !connector.leftDone->getValue()) {
        
        if(connector.leftLoc->getValue() >= currentTimeBin || connector.leftDone->getValue()) {
            std::cout << "Processing ellipse encapsulation batch at time: " << currentTimeBin << std::endl;
            
            // Get clusters for current time bin
            auto periodicClusterVect = periodicClusterStream->find(currentTimeBin);
            
            if(!periodicClusterVect.empty()) {
                encapsulateEllipse(periodicClusterVect);
            } else {
                std::cout << "No clusters found for time bin " << currentTimeBin << std::endl;
            }
            
            // Signal completion of this batch to next stage
            connector.rightLoc->setValue(currentTimeBin);
            currentTimeBin += options->freq.time_bin_duration;
        }
    }
    
    // Signal completion to next stage
    connector.rightDone->setValue(true);
    std::cout << "Finished Ellipse Encapsulation" << std::endl;
}