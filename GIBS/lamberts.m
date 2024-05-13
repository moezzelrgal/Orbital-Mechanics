%***********************************************************************
%* This function is provided under the license: CC BY-NC-SA 4.0        * 
%* License details: https://creativecommons.org/licenses/by-nc-sa/4.0/ *
%* Author: Mohammed Ezzelrgal                                          *
%***********************************************************************
% This function computes the velocity vectors v1 and v2 for a given pair of
% position vectors r1 and r2, separated by a specified time interval delta_t,
% using Lambert's problem solver.
%
% Inputs:
%   - r1: Initial position vector (3x1) [km]
%   - r2: Final position vector (3x1) [km]
%   - delta_t: Time interval between r1 and r2 [s]
%
% Outputs:
%   - v1: Velocity vector at r1 (3x1) [km/s]
%   - v2: Velocity vector at r2 (3x1) [km/s]
function [v1, v2] = lamberts(r1, r2, delta_t)
    % Constants
    mu = 3.986e5; % Gravitational constant for Earth (km^3/s^2)
    
    % Step 1: Calculate r1 and r2
    r1_mag = norm(r1);
    r2_mag = norm(r2);
    
    % Step 2: Calculate Δθ
    delta_theta = acos(dot(r1, r2) / (r1_mag * r2_mag));
    if dot(cross(r1, r2), [0; 0; 1]) < 0
        delta_theta = 2*pi - delta_theta;
    end
    
    % Step 3: Calculate A
    A = sqrt(r1_mag * r2_mag) * sin(delta_theta) / (1 - cos(delta_theta));
    
    % Step 4: Iteration to solve for z
    z = 0;
    while true
        F = (1 - (r1_mag / sqrt(A)) * (1 - cos(sqrt(A) * delta_t)) - (r2_mag / sqrt(A)) * (1 - cos(sqrt(A) * delta_t)) + (r1_mag * r2_mag / A) * (sqrt(A) * delta_t - sin(sqrt(A) * delta_t))) / sqrt(A);
        F_prime = sqrt(2 / A) * (delta_t - (1 - (r1_mag / sqrt(A))) * sqrt(A) * (1 - cos(sqrt(A) * delta_t)) - (r2_mag / sqrt(A)) * (1 - cos(sqrt(A) * delta_t)));
        z_new = z - F / F_prime;
        if abs(z_new - z) < 1e-6
            z = z_new;
            break;
        end
        z = z_new;
    end
    
    % Step 5: Calculate y
    y = r1 + r2 + ((A / sqrt(z)) * ((z * sin(sqrt(z) * delta_t)) / sqrt(1 - cos(sqrt(z) * delta_t))));
    
    % Step 6: Calculate f, g, and g_ functions
    f = 1 - (y / r1_mag);
    g = A * sqrt(y / mu);
    g_prime = 1 - (y / r2_mag);
    
    % Step 7: Calculate v1 and v2
    v1 = (1 / g) * (r2 - f * r1);
    v2 = (1 / g) * (g_prime * r2 - r1);
end
