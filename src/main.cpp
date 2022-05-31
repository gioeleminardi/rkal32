/*
NAME
    rkal32.c - A third order Kalman filter for
    RDAS interpreted data files.
SYNOPSIS
    rkal32 [model] <infile >outfile
    model is an optional parameter to specify the ratio
    between the variance for the acceleration measurement
    (fixed) and the model (variable). Default value for
    this is that the model is about 100 times better
    than the measurement.
DESCRIPTION:
    Performs Kalman filtering on standard input
    data using a third order, or constant acceleration,
    propagation model.
    The standard input data is of the form:
    Column 1: Time of the measurement (seconds)
    Column 2: Acceleration Measurement (feet/sec/sec)
    Column 3: Pressure Measurement (altitude in feet)
    The standard output data is of the form:
    Column 1: Time of the estimate/measurement
    Column 2: Pressure Measurement
    Column 3: Acceleration Measurement
    Column 4: Position estimate
    Column 5: Rate estimate
    Column 6: Acceleration estimate
AUTHOR
    Hacked by David Schultz
    REPORTING BUGS AND REVISIONS
    mailto:david.schultz@earthlink.net
COMPILING
    cc -o rkal32.exe rkal32.c
*/
#include <cmath>
#include <cstdio>
#include <cstdlib>

#include "csv.hpp"

constexpr double ALTITUDESIGMA     = 4.572;
constexpr double ACCELERATIONSIGMA = 1.83;
constexpr double MODELSIGMA        = 0.183;

constexpr double altitude_variance     = ALTITUDESIGMA * ALTITUDESIGMA;
constexpr double acceleration_variance = ACCELERATIONSIGMA * ACCELERATIONSIGMA;
constexpr double model_variance        = MODELSIGMA * MODELSIGMA;

int main(int argc, char** argv)
{
  csv::CSVReader reader("../sim/out_matrix.csv");
  csv::CSVRow    row;

  char   buf[512];
  int    i, j, k, notdone;
  double alt_inovation, accel_inovation;
  double time, accel, pressure;
  double last_time, last_pressure;
  double det;
  double est[3]      = {0, 0, 0};
  double estp[3]     = {0, 0, 0};
  double pest[3][3]  = {2, 0, 0, 0, 9, 0, 0, 0, 9};
  double pestp[3][3] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
  double phi[3][3]   = {1, 0, 0, 0, 1, 0, 0, 0, 1.0};
  double phit[3][3]  = {1, 0, 0, 0, 1, 0, 0, 0, 1.0};
  double kgain[3][2] = {0.01, 0.01, 0.01, 0.01, 0.01, 0.01};
  double lastkgain[3][2], dt;
  double term[3][3];

  /* Initialize */

  reader.read_row(row);
  time  = row["TIME_s"].get<double>();
  accel = row["ACCELERATION_mss"].get<double>();
  if (accel < 0)
    accel *= -1.0f;
  pressure = row["ALTITUDE_m"].get<double>();

  est[0]    = pressure;
  last_time = time;

  reader.read_row(row);
  time     = row["TIME_s"].get<double>();
  accel    = row["ACCELERATION_mss"].get<double>();
  if (accel < 0)
    accel *= -1.0f;
  pressure = row["ALTITUDE_m"].get<double>();

  dt        = time - last_time;
  last_time = time;

  /*
   * Fill in state transition matrix and its transpose
   */
  phi[0][1]  = dt;
  phi[1][2]  = dt;
  phi[0][2]  = dt * dt / 2.0;
  phit[1][0] = dt;
  phit[2][1] = dt;
  phit[2][0] = dt * dt / 2.0;

  /* Compute the Kalman gain matrix. */
  for (i = 0; i <= 2; i++)
    for (j = 0; j <= 1; j++)
      lastkgain[i][j] = kgain[i][j];

  k = 0;
  while (true)
  {
    /* Propagate state covariance */
    term[0][0]  = phi[0][0] * pest[0][0] + phi[0][1] * pest[1][0] + phi[0][2] * pest[2][0];
    term[0][1]  = phi[0][0] * pest[0][1] + phi[0][1] * pest[1][1] + phi[0][2] * pest[2][1];
    term[0][2]  = phi[0][0] * pest[0][2] + phi[0][1] * pest[1][2] + phi[0][2] * pest[2][2];
    term[1][0]  = phi[1][0] * pest[0][0] + phi[1][1] * pest[1][0] + phi[1][2] * pest[2][0];
    term[1][1]  = phi[1][0] * pest[0][1] + phi[1][1] * pest[1][1] + phi[1][2] * pest[2][1];
    term[1][2]  = phi[1][0] * pest[0][2] + phi[1][1] * pest[1][2] + phi[1][2] * pest[2][2];
    term[2][0]  = phi[2][0] * pest[0][0] + phi[2][1] * pest[1][0] + phi[2][2] * pest[2][0];
    term[2][1]  = phi[2][0] * pest[0][1] + phi[2][1] * pest[1][1] + phi[2][2] * pest[2][1];
    term[2][2]  = phi[2][0] * pest[0][2] + phi[2][1] * pest[1][2] + phi[2][2] * pest[2][2];
    pestp[0][0] = term[0][0] * phit[0][0] + term[0][1] * phit[1][0] + term[0][2] * phit[2][0];
    pestp[0][1] = term[0][0] * phit[0][1] + term[0][1] * phit[1][1] + term[0][2] * phit[2][1];
    pestp[0][2] = term[0][0] * phit[0][2] + term[0][1] * phit[1][2] + term[0][2] * phit[2][2];
    pestp[1][0] = term[1][0] * phit[0][0] + term[1][1] * phit[1][0] + term[1][2] * phit[2][0];
    pestp[1][1] = term[1][0] * phit[0][1] + term[1][1] * phit[1][1] + term[1][2] * phit[2][1];
    pestp[1][2] = term[1][0] * phit[0][2] + term[1][1] * phit[1][2] + term[1][2] * phit[2][2];
    pestp[2][0] = term[2][0] * phit[0][0] + term[2][1] * phit[1][0] + term[2][2] * phit[2][0];
    pestp[2][1] = term[2][0] * phit[0][1] + term[2][1] * phit[1][1] + term[2][2] * phit[2][1];
    pestp[2][2] = term[2][0] * phit[0][2] + term[2][1] * phit[1][2] + term[2][2] * phit[2][2];

    pestp[2][2] = pestp[2][2] + model_variance;
    /*
    Calculate Kalman Gain
    */
    det         = (pestp[0][0] + altitude_variance) * (pestp[2][2] + acceleration_variance) - pestp[2][0] * pestp[0][2];
    kgain[0][0] = (pestp[0][0] * (pestp[2][2] + acceleration_variance) - pestp[0][2] * pestp[2][0]) / det;
    kgain[0][1] = (pestp[0][0] * (-pestp[0][2]) + pestp[0][2] * (pestp[0][0] + altitude_variance)) / det;
    kgain[1][0] = (pestp[1][0] * (pestp[2][2] + acceleration_variance) - pestp[1][2] * pestp[2][0]) / det;
    kgain[1][1] = (pestp[1][0] * (-pestp[0][2]) + pestp[1][2] * (pestp[0][0] + altitude_variance)) / det;
    kgain[2][0] = (pestp[2][0] * (pestp[2][2] + acceleration_variance) - pestp[2][2] * pestp[2][0]) / det;
    kgain[2][1] = (pestp[2][0] * (-pestp[0][2]) + pestp[2][2] * (pestp[0][0] + altitude_variance)) / det;
    pest[0][0]  = pestp[0][0] * (1.0 - kgain[0][0]) - kgain[0][1] * pestp[2][0];
    pest[0][1]  = pestp[0][1] * (1.0 - kgain[0][0]) - kgain[0][1] * pestp[2][1];
    pest[0][2]  = pestp[0][2] * (1.0 - kgain[0][0]) - kgain[0][1] * pestp[2][2];
    pest[1][0]  = pestp[0][0] * (-kgain[1][0]) + pestp[1][0] - kgain[1][1] * pestp[2][0];
    pest[1][1]  = pestp[0][1] * (-kgain[1][0]) + pestp[1][1] - kgain[1][1] * pestp[2][1];
    pest[1][2]  = pestp[0][2] * (-kgain[1][0]) + pestp[1][2] - kgain[1][1] * pestp[2][2];
    pest[2][0]  = (1.0 - kgain[2][1]) * pestp[2][0] - kgain[2][0] * pestp[2][0];
    pest[2][1]  = (1.0 - kgain[2][1]) * pestp[2][1] - kgain[2][0] * pestp[2][1];
    pest[2][2]  = (1.0 - kgain[2][1]) * pestp[2][2] - kgain[2][0] * pestp[2][2];

    /* Check for convergence. Criteria is less than .001% change from last
     * time through the mill.
     */
    notdone = 0;
    k++;
    for (i = 0; i <= 2; i++)
      for (j = 0; j <= 1; j++)
      {
        if ((kgain[i][j] - lastkgain[i][j]) / lastkgain[i][j] > 0.00001)
          notdone++;
        lastkgain[i][j] = kgain[i][j];
      }
    if (notdone)
      continue;
    else
      break;
  }

  printf("Input noise values used (standard deviation):\n");
  printf("# Altitude - %.4f m\n", sqrt(altitude_variance));
  printf("# Acceleration - %.4f m/s/s\n", sqrt(acceleration_variance));
  printf("# Model noise - %.4f m/s/s\n", sqrt(model_variance));
  printf("# Kalman gains converged after %d iterations.\n", k);
  for (i = 0; i <= 2; i++)
    for (j = 0; j <= 1; j++)
      printf("%15f ", kgain[i][j]);
  printf("\n");
  printf("# Estimated output first order statistics (standard deviation):\n");
  printf("# Altitude - %.4f m\n", sqrt(pest[0][0]));
  printf("# Velocity - %.4f m/s\n", sqrt(pest[1][1]));
  printf("# Acceleration - %.4f m/s/s\n", sqrt(pest[2][2]));

  /* Now run the Kalman filter on the data using previously
   * determined gains.
   */
  /* Output header for data. */
  printf("#\n# Output from rkal32:\n# A third order Kalman filter using acceleration and pressure measurements\n");
  printf("Time;Alt;Acc;EstPos;EstRate;EstAcc\n");
  while (reader.read_row(row))
  {
    time     = row["TIME_s"].get<double>();
    accel    = row["ACCELERATION_mss"].get<double>();
    if (accel < 0)
      accel *= -1.0f;
    pressure = row["ALTITUDE_m"].get<double>();

    /* remove offset and convert from G's to ft/sec/sec */
    //    accel = (accel - 1.0) * 32.17417;
    /* sanity check on time */
    if (last_time >= time)
    {
      fprintf(stderr, "Time does not increase.\n");
      exit(1);
    }
    /* Compute the innovations */
    alt_inovation   = pressure - estp[0];
    accel_inovation = accel - estp[2];

    /* Experimental code to modify Mach transition pressure
     * disturbances.
     */
    if (std::abs(alt_inovation) > 100)
    {
      /* We have a large error in altitude. Now see how fast we are
       * going.
       */
      if (estp[1] > 900 && estp[1] < 1200)
      {
        /* Somewhere in the neighborhood of Mach 1. Now check to
         * see if we are slowing down.
         */
        //
        if (estp[2] < 0)
        {
          /*
           * OK, now what do we do? Assume that velocity and
           * acceleration estimates are accurate. Adjust current
           * altitude estimate to be the same as the measured
           * altitude.
           */
          est[0]        = pressure;
          alt_inovation = 0;
        }
      }
    }
    /* Simple check for over-range on pressure measurement.
     * This is just hacked in based on a single data set. Actual
     * flight software needs something more sophisticated.
     */
    if (pressure > 11000)
      alt_inovation = 0;
    /* Propagate state */
    estp[0] = phi[0][0] * est[0] + phi[0][1] * est[1] + phi[0][2] * est[2];
    estp[1] = phi[1][0] * est[0] + phi[1][1] * est[1] + phi[1][2] * est[2];
    estp[2] = phi[2][0] * est[0] + phi[2][1] * est[1] + phi[2][2] * est[2];

    /*
     * Update state
     */
    est[0] = estp[0] + kgain[0][0] * alt_inovation + kgain[0][1] * accel_inovation;
    est[1] = estp[1] + kgain[1][0] * alt_inovation + kgain[1][1] * accel_inovation;
    est[2] = estp[2] + kgain[2][0] * alt_inovation + kgain[2][1] * accel_inovation;
    /*
     * Output
     */
    printf("%.2f;%.4f;%.4f;%.4f;%.4f;%.4f\n", time, pressure, accel, est[0], est[1], est[2]);
    last_time = time;
  }
}