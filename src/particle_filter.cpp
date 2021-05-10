/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include "particle_filter.h"

#include <math.h>

#include <algorithm>
#include <iostream>
#include <iterator>
#include <map>
#include <numeric>
#include <random>
#include <set>
#include <string>
#include <vector>

#include "helper_functions.h"

using std::string;
using std::vector;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1.
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method
   *   (and others in this file).
   */
  num_particles = 0;  // TODO: Set the number of particles
  std::default_random_engine gen;
  std::normal_distribution<double> dist_x(x, std[0]);
  std::normal_distribution<double> dist_y(y, std[1]);
  std::normal_distribution<double> dist_theta(theta, std[2]);
  for (int i = 0; i < num_particles; ++i) {
    Particle particle;
    particle.x = dist_x(gen);
    particle.y = dist_y(gen);
    particle.theta = dist_theta(gen);
    particle.weight = 1.0;
    particles.push_back(particle);
  }
  is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[],
                                double velocity, double yaw_rate) {
  /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */
  for (int i = 0; i < num_particles; ++i) {
    particles[i].x += velocity / yaw_rate *
                      (sin(particles[i].theta + yaw_rate * delta_t) -
                       sin(particles[i].theta));
    particles[i].y += velocity / yaw_rate *
                      (cos(particles[i].theta) -
                       cos(particles[i].theta + yaw_rate * delta_t));
    particles[i].theta += yaw_rate * delta_t;
  }
}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted,
                                     vector<LandmarkObs>& observations) {
  /**
   * TODO: Find the predicted measurement that is closest to each
   *   observed measurement and assign the observed measurement to this
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will
   *   probably find it useful to implement this method and use it as a helper
   *   during the updateWeights phase.
   */
  std::set<int> paired_landmarks;  // assume landmark only corresponds with one
                                   // observation
  std::map<float, int> distance_map;

  for (unsigned int i = 0; i < observations.size(); ++i) {
    for (unsigned int j = 0; j < predicted.size(); ++j) {
      // id not in paired landmarks
      if (!paired_landmarks.count(predicted[j].id)) {
        double distance = dist(observations[i].x, observations[i].y,
                               predicted[j].x, predicted[j].y);
        distance_map.insert(std::pair<float, int>(distance, j));
      }
    }
    // insert the closest landmark
    int closest_index = distance_map.begin()->second;
    paired_landmarks.insert(predicted[closest_index].id);
    // exchange the closest landmark with the ith element
    int id_temp = predicted[i].id;
    double x_temp = predicted[i].x;
    double y_temp = predicted[i].y;
    predicted[i].id = predicted[closest_index].id;
    predicted[i].x = predicted[closest_index].x;
    predicted[i].y = predicted[closest_index].y;
    predicted[closest_index].id = id_temp;
    predicted[closest_index].x = x_temp;
    predicted[closest_index].y = y_temp;

    distance_map.clear();
  }
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[],
                                   const vector<LandmarkObs>& observations,
                                   const Map& map_landmarks) {
  /**
   * TODO: Update the weights of each particle using a mult-variate Gaussian
   *   distribution. You can read more about this distribution here:
   *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
   * NOTE: The observations are given in the VEHICLE'S coordinate system.
   *   Your particles are located according to the MAP'S coordinate system.
   *   You will need to transform between the two systems. Keep in mind that
   *   this transformation requires both rotation AND translation (but no
   * scaling). The following is a good resource for the theory:
   *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
   *   and the following is a good resource for the actual equation to implement
   *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
   */
  bool closest_landmark_association = true;
  bool print_particle_association = true;

  vector<double>().swap(weights);
  double weight_sum = 0.0;
  for (int i = 0; i < num_particles; ++i) {
    // transform each observation from local frame into map frame
    vector<LandmarkObs> trans_observations;
    for (unsigned int j = 0; j < observations.size(); ++j) {
      trans_observations.push_back(homogenousTransform(
          particles[i].x, particles[i].y, particles[i].theta, observations[j].x,
          observations[j].y, observations[j].id));
    }

    // filter out the map_landmarks with current guessed position with
    // sensor_range
    vector<LandmarkObs> predicted;
    for (unsigned int k = 0; k < map_landmarks.landmark_list.size(); ++k) {
      if (dist(particles[i].x, particles[i].y,
               map_landmarks.landmark_list[k].x_f,
               map_landmarks.landmark_list[k].y_f) < sensor_range) {
        LandmarkObs landmark;
        landmark.id = map_landmarks.landmark_list[k].id_i;
        landmark.x = map_landmarks.landmark_list[k].x_f;
        landmark.y = map_landmarks.landmark_list[k].y_f;
        predicted.push_back(landmark);
      }
    }

    // for each observation find a closest landmark
    if (closest_landmark_association) {
      dataAssociation(predicted, trans_observations);

      // optional print msg
      if (print_particle_association) {
        vector<int> associated_id;
        vector<double> associated_x;
        vector<double> associated_y;

        for (unsigned int m = 0; m < predicted.size(); ++m) {
          associated_id.push_back(predicted[m].id);
          associated_x.push_back(predicted[m].x);
          associated_y.push_back(predicted[m].y);
        }
        SetAssociations(particles[i], associated_id, associated_x,
                        associated_y);
      }
    }

    // update weight
    double weight = 1.0;
    for (unsigned int j = 0; j < trans_observations.size(); ++j) {
      // calc multi gauss
      double gauss_norm = 1 / (2 * M_PI * std_landmark[0] * std_landmark[1]);
      double exponent = (pow(trans_observations[j].x - predicted[j].x, 2) /
                         (2 * pow(std_landmark[0], 2))) +
                        (pow(trans_observations[j].y - predicted[j].y, 2) /
                         (2 * pow(std_landmark[1], 2)));
      weight *= gauss_norm * exp(-exponent);
    }
    weights.push_back(weight);
    weight_sum += weight;
  }
  // normalize probability
  for (int k = 0; k < num_particles; ++k) {
    weights[k] /= weight_sum;
  }
}

LandmarkObs ParticleFilter::homogenousTransform(double origin_x,
                                                double origin_y, double theta,
                                                double obj_x, double obj_y,
                                                int id) {
  LandmarkObs parent_frame;
  parent_frame.x = origin_x + cos(theta) * obj_x - sin(theta) * obj_y;
  parent_frame.y = origin_y + sin(theta) * obj_x - cos(theta) * obj_y;
  parent_frame.id = id;
  return parent_frame;
}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional
   *   to their weight.
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */
  std::default_random_engine gen;
  std::uniform_int_distribution<int> distribution(0, num_particles);
  int index = distribution(gen);
  double beta = 0.0;
  vector<Particle> resampled;
  for (int i = 0; i < num_particles; ++i) {
    beta += distribution(gen) / num_particles * 2.0 *
            (*std::max_element(weights.begin(), weights.end()));
    while (weights[index] < beta) {
      beta -= weights[index];
      index = (index + 1) % num_particles;
    }
    resampled.push_back(particles[index]);
  }
  particles = resampled;
}

void ParticleFilter::SetAssociations(Particle& particle,
                                     const vector<int>& associations,
                                     const vector<double>& sense_x,
                                     const vector<double>& sense_y) {
  // particle: the particle to which assign each listed association,
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  particle.associations = associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length() - 1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord) {
  vector<double> v;

  if (coord == "X") {
    v = best.sense_x;
  } else {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length() - 1);  // get rid of the trailing space
  return s;
}