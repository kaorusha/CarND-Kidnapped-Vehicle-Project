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
  for (unsigned int i = 0; i < particles.size(); ++i) {
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
                                     vector<LandmarkObs>& observations,
                                     const Map& map_landmarks) {
  /**
   * TODO: Find the predicted measurement that is closest to each
   *   observed measurement and assign the observed measurement to this
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will
   *   probably find it useful to implement this method and use it as a helper
   *   during the updateWeights phase.
   */
  std::set<int> paired_landmarks;
  std::map<float, int> distance_map;

  for (unsigned int i = 0; i < observations.size(); ++i) {
    for (unsigned int j = 0; j < map_landmarks.landmark_list.size(); ++j) {
      // id not in paired landmarks
      if (!paired_landmarks.count(map_landmarks.landmark_list[j].id_i)) {
        float distance = (float)calcDistance(observations[i],
                                             map_landmarks.landmark_list[j]);
        distance_map.insert(std::pair<float, int>(
            distance, map_landmarks.landmark_list[j].id_i));
      }
    }
    // insert the closest landmark
    paired_landmarks.insert(distance_map.begin()->second);
    distance_map.clear();
    // map id starts with 1
    int index = distance_map.begin()->second - 1;
    LandmarkObs landmark;
    landmark.id = map_landmarks.landmark_list[index].id_i;
    landmark.x = map_landmarks.landmark_list[index].x_f;
    landmark.y = map_landmarks.landmark_list[index].y_f;
    predicted.push_back(landmark);
  }
}

double calcDistance(LandmarkObs obs, Map::single_landmark_s landmark) {
  return sqrt(pow(obs.x - landmark.x_f, 2) + pow(obs.y - landmark.y_f, 2));
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

  for (unsigned int i = 0; i < particles.size(); ++i) {
    vector<LandmarkObs> trans_observations;
    vector<LandmarkObs> predicted;
    for (unsigned int j = 0; j < observations.size(); ++j) {
      trans_observations.push_back(homogenousTransform(
          particles[i].x, particles[i].y, particles[i].theta, observations[j].x,
          observations[j].y, observations[j].id));
    }
    // for each observation find a closest landmark
    vector<LandmarkObs> predicted;
    dataAssociation(predicted, trans_observations, map_landmarks);
    for (unsigned int i = 0; i < trans_observations.size(); ++i) {
    }
    // copy landmark vector, if associated with observation, erase it from the
    // landmark list
  }
}

LandmarkObs homogenousTransform(double origin_x, double origin_y, double theta,
                                double obj_x, double obj_y, int id) {
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