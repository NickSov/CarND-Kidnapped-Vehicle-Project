//
// particle_filter.cpp
//
// Created on: Dec 12, 2016
// Last Editted: Jan 27, 2019
// Authors: Tiffany Huang & Nicholas Atanasov
//

#include "particle_filter.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <vector>
#include "Eigen/Dense"

using  Eigen::MatrixXd;
using  Eigen::VectorXd;

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
  num_particles = 100;  // TODO: Set the number of particles

  // Create a random number engine instance

  std::default_random_engine gen;

  // Create a normal Gaussian distribution for x, y and theta
  std::normal_distribution<double> distX(x, std[0]);
  std::normal_distribution<double> distY(y, std[1]); // y
  std::normal_distribution<double> distTheta(theta, std[2]); // theta

  // Initialize particles based on the Gaussian distributon
  for(int i = 0; i < num_particles; i++)
  {
    // Sample a Gaussian distribution, taking into account sensor uncertainty/noise
    Particle newParticle;
    newParticle.x = distX(gen);
    newParticle.y = distY(gen);
    newParticle.theta = distTheta(gen);
    newParticle.weight = 1; // initial weight is set to 1

    particles.push_back(newParticle); // add the new particle to the particle vector

    //Question: Do I need to add more noise?
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

   for(int i=0; i<particles.size(); i++)
   {
     double x0 = particles[i].x; // initial x position of particle
     double *xf = &particles[i].x; // x position of particle after update
     double y0 = particles[i].y; // initial y position of particle
     double *yf = &particles[i].y; // y position of particle after update
     double theta0 = particles[i].theta; // initial heading (theta) of particle
     double *thetaf = &particles[i].theta; // heading (theta) of particle after update
     double v = velocity; // velocity of the vehicle
     double thetaDot = yaw_rate; // yaw rate of the vehicle - the rate of change of theta
     double dt = delta_t; // change in time - time increment

     // Create a random number engine instance

     std::default_random_engine genPred;

     // Add measurements to each particle

     if(yaw_rate > 0.0001){
       // update x position, y position, and heading of a particle if the yaw rate is NOT close to zero
       *xf = x0 + (v/thetaDot)*(sin(theta0+thetaDot*dt)-sin(theta0));
       *yf = y0 + (v/thetaDot)*(cos(theta0)-cos(theta0+thetaDot*dt));
       *thetaf = theta0 + thetaDot*dt;

       // Create random noise values
       std::normal_distribution<double> distX(*xf, std_pos[0]); // x
       std::normal_distribution<double> distY(*yf, std_pos[1]); // y
       std::normal_distribution<double> distTheta(*thetaf, std_pos[2]); // theta

       // add noise

       *xf = distX(genPred);
       *yf = distY(genPred);
       *thetaf = distTheta(genPred);

     } else {
       // update x position, y position, and heading of a particle if the yaw rate is close to zero
       *xf = x0 + v*dt*cos(theta0);
       *yf = y0 + v*dt*sin(theta0);
       *thetaf = theta0;

       // Create random noise values
       std::normal_distribution<double> distX(x0, std_pos[0]); //x
       std::normal_distribution<double> distY(y0, std_pos[1]); // y
       std::normal_distribution<double> distTheta(theta0, std_pos[2]); // theta

       // add noise

       *xf = distX(genPred);
       *yf = distY(genPred);
       *thetaf = distTheta(genPred);
     }
   }

   // NOTE: Need to still add gaussian sensor/positioning noise

}

VectorXd ParticleFilter::transformMarker(LandmarkObs obs, int indexParticle)
{
  // TRANSFORMATION FUNCTION
  // Transforms the particles observations (x,y) from vehicle to map coordinates for a single observation

  double x_obs = obs.x; // observation distance, x component
  double y_obs = obs.y; // observation distance, y component

  double x = particles[indexParticle].x; // current x position of particle
  double y = particles[indexParticle].y; // current y position of particle
  double theta = particles[indexParticle].theta; // initial heading (theta) of particle

  MatrixXd transformMatrix(3,3);
  transformMatrix << cos(theta), -sin(theta), x,
                     sin(theta), cos(theta), y,
                     0, 0, 1;

  VectorXd markerVector(3);
  markerVector << x_obs, y_obs, 1;
  VectorXd mapVector(3);
  mapVector = transformMatrix * markerVector;

  return mapVector;
}

void ParticleFilter::dataAssociation(const Map &map_landmarks,
                                     vector<LandmarkObs> &observations) {
  /**
   * TODO: Find the predicted measurement that is closest to each
   *   observed measurement and assign the observed measurement to this
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will
   *   probably find it useful to implement this method and use it as a helper
   *   during the updateWeights phase.
   */

   for (int i = 0; i < observations.size(); ++i)
   {
    int closest_landmark = 0;
    int min_dist = 999999;
    int curr_dist;
    // Iterate through all landmarks to check which is closest
     for (int j = 0; j < map_landmarks.landmark_list.size(); ++j) {
       // Calculate Euclidean distance
       curr_dist = sqrt(pow(observations[i].x - map_landmarks.landmark_list[j].x_f, 2)+ pow(observations[i].y - map_landmarks.landmark_list[j].y_f, 2));
       //std::cout << " Landmark Id: " << map_landmarks.landmark_list[j].id_i << " | Landmark Val x_f: " << map_landmarks.landmark_list[j].x_f << " | curr_dist: " << curr_dist  << std::endl;
       // Compare to min_dist and update if closest
       if (curr_dist < min_dist) {
         min_dist = curr_dist;
         closest_landmark = map_landmarks.landmark_list[j].id_i;
         observations[i].id = closest_landmark;
       }
     }
    // std::cout << " Observation: " << i+1 << " | Landmark: " << observations[i].id << std::endl;
   }
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[],
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {//double sensor_range, double std_landmark[],
                                   //const vector<LandmarkObs> &observations,
                                   //const Map &map_landmarks) {
  /**
   * TODO: Update the weights of each particle using a mult-variate Gaussian
   *   distribution. You can read more about this distribution here:
   *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
   * NOTE: The observations are given in the VEHICLE'S coordinate system.
   *   Your particles are located according to the MAP'S coordinate system.
   *   You will need to transform between the two systems. Keep in mind that
   *   this transformation requires both rotation AND translation (but no scaling).
   *   The following is a good resource for the theory:
   *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
   *   and the following is a good resource for the actual equation to implement
   *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
   */


   VectorXd obsTransPos(2); // empty transform observation vector
   LandmarkObs obs;
   LandmarkObs transformObs;
   vector<LandmarkObs> transformObsVect;
   double stdDevX = std_landmark[0];
   double stdDevY = std_landmark[1];

   for(int i=0; i<particles.size(); i++)
   {
     particles[i].weight = 1.0;




     transformObsVect.clear(); //clear the transformed observation vector so it's ready for the next particle
     //std::cout << "observations for Particle " << i+1 << std::endl;
     for(int j=0; j<observations.size(); j++)
     {
       obsTransPos.fill(0);
       obs = observations[j];
       obsTransPos = transformMarker(obs,i);
       transformObs.id = 1; // temp id
       transformObs.x = obsTransPos[0];
       transformObs.y = obsTransPos[1];
       transformObsVect.push_back(transformObs);
       //std::cout << "x: " << transformObs.x << "| y: " << transformObs.y << std::endl;
     }

     /*
     // vector<LandmarkObs> landmarks;
     // for(int i=0; i<map_landmarks.landmark_list.size();i++)
     // {
     //  //vector<double> landmarkObjects = {map_landmarks.landmark_list[i].id_i,map_landmarks.landmark_list[i].x_f, map_landmarks.landmark_list[i].y_f};
     //  //landmarks.push_back( = {map_landmarks.landmark_list[i].id_i,map_landmarks.landmark_list[i].x_f, map_landmarks.landmark_list[i].y_f};
     // }

       // Map reducedLandmarks; //map instance that holds reduced number of landmarks based on sensor range.
       //
       // // determine max and min oberservation values
       //
       // double maxValX = 999;
       // double minValX = 999;
       // double maxValY = 999;
       // double minValY = 999;
       //
       // for(int w=0; w<transformObsVect.size(); w++) // loop through each observation
       // {
       //   if(transformObsVect.x < 999){
       //     minValX = transformObsVect.x
       //   }
       //
       //
       // }

       */

       /*

       // for(int i=0; i < map_landmarks.rows(); i++)
       // {
       //   if(abs(map_landmarks.x_f - ) > abs() || map_landmarks.y_f)
       //   // Declare single_landmark
       //   Map::single_landmark_s single_landmark_temp;
       //
       //   // Set values
       //   single_landmark_temp.id_i = map_landmarks.id_i;
       //   single_landmark_temp.x_f  = map_landmarks.x_f;
       //   single_landmark_temp.y_f  = map_landmarks.y_f;
       //
       //   // Add to landmark list of map
       //   reducedLandmarks.landmark_list.push_back(single_landmark_temp);
       // }
       */


       Map rangedLandmarks;

       LandmarkObs newLandMark;
       vector<LandmarkObs> newLandMarkVect;

       dataAssociation(map_landmarks)

       for(int g=0; g< map_landmarks.landmark_list.size(); g++)
       {

        if()
        {
          rangedLandmarks.landmark_list.push_back
        }
       }

       //std::cout << "Associated Landmarks for Particle # " << i+1 << " : " << std::endl;
       dataAssociation(map_landmarks,transformObsVect); // associate transformed observations with map landmarks for each particle

       // Insert landmark associations into the particle instance

       for(int w=0; w<transformObsVect.size(); w++)
       {
         particles[i].associations.push_back(transformObsVect[w].id);
         particles[i].sense_x.push_back(transformObsVect[w].x);
         particles[i].sense_y.push_back(transformObsVect[w].y);
         //std::cout << "Particle " << i << " association: " << particles[i].associations[w]<< " | x: " << particles[i].sense_x[w] << std::endl;
       }

       // Calculate weight for particle i

       double accumWeight=1;

       for(int z=0; z<particles[i].associations.size();z++)
       {
         int partLandMarkId = particles[i].associations[z]; //id of landmark associated to particle i
         double x = particles[i].sense_x[z]; // x position observation in map coordinates - what the sensor sensed
         double y = particles[i].sense_y[z]; // y position observation in map coordinates - what the sensor sensed
         double muX; //location of landmark in x
         double muY; //location of landmark in y
         for (int l = 0; l < map_landmarks.landmark_list.size(); l++)
         {
           int idLandmark = map_landmarks.landmark_list[l].id_i;
           double xLandmark = map_landmarks.landmark_list[l].x_f;
           double yLandmark = map_landmarks.landmark_list[l].y_f;
           if(idLandmark == partLandMarkId)
           {
             muX = xLandmark ; //location of landmark in x
             muY = yLandmark ; //location of landmark in y
           }
         }

         double expGauss = (pow(x - muX, 2) / (2 * pow(stdDevX, 2)))
                         + (pow(y - muY, 2) / (2 * pow(stdDevY, 2)));
         double gaussBase = 1 / (2 * M_PI * stdDevX * stdDevY);
         double Pxy = gaussBase*exp(-expGauss);
         std::cout << "x " << z << " value: " << x << std::endl;
         std::cout << "y " << z << " value: " << y << std::endl;
         std::cout << "muX " << z << " value: " << muX << std::endl;
         std::cout << "muY " << z << " value: " << muY << std::endl;
         std::cout << "Pxy " << z << " value: " << Pxy << std::endl;
         std::cout << "stdX " << z << " value: " << stdDevX << std::endl;
         std::cout << "stdY " << z << " value: " << stdDevY << std::endl;
         std::cout << " - - - - - - - - - - - - " << std::endl;
         accumWeight *= Pxy;

         // Update the weight vector

       }
       particles[i].weight = accumWeight;
       weights.push_back(particles[i].weight);
       std::cout << "Particle " << i << " weight: " << particles[i].weight << std::endl;


   }

}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional
   *   to their weight.
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */

   std::default_random_engine genIndex;

   vector<Particle> resampledP;

   // Create random noise values
   std::discrete_distribution<int> index(weights.begin(), weights.end()); //x

   for (int i = 0; i < particles.size(); i++){
     int weightedIndex = index(genIndex);
     resampledP.push_back(particles[weightedIndex]);
   }

   particles = resampledP;
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
  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
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
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}
