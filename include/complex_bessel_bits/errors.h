/*! \file errors.h
 *
 * \author Joey Dumont <joey.dumont@gmail.com>
 * \author Denis Gagnon <gagnon88@gmail.com>
 *
 * \brief Definition of the functions computing the Bessel functions.
 *
 * Defines structs to hold error codes and error messages in the computation
 * of error functions.
 *
 * \copyright LGPL
 */

#ifndef ERRORS_H
#define ERRORS_H

#include <array>
#include <string>
#include <unordered_map>

// Error handling main code inspired by this:
// https://stackoverflow.com/questions/47841783/is-there-any-advantage-in-using-stdoptional-to-avoid-default-arguments-in-a-fu
enum ErrorCode
{
  Success,
  InputError,
  Overflow,
  PartialLossOfSignificance,
  FullLossOfSignificance,
  AlgorithmTermination
};

typedef struct BesselErrors
{
  ErrorCode   errorCode;    // Equivalent of IERR in the original FORTRAN.
  std::string errorMessage; // Error message associated with the error.
} BesselErrors;

static const std::unordered_map<std::string, std::array<BesselErrors, 6>> errorMessages = {
  { "besselJ",
    std::array<BesselErrors, 6>{
      BesselErrors{ Success, "Normal return        -- Computation completed." },
      BesselErrors{ InputError, "Input error          -- No computation" },
      BesselErrors{ Overflow,
                    "Overflow             -- No computation, Im(z) too large for scale=false" },
      BesselErrors{ PartialLossOfSignificance,
                    "Loss of significance -- abs(z) or order large \n computation done but losses of "
                    "significance by argument reduction produce less than half of machine accuracy" },
      BesselErrors{ FullLossOfSignificance,
                    "Loss of significance -- abs(z) or order too large\n no computation because of "
                    "complete loss of significance by argument reduction" },
      BesselErrors{
        AlgorithmTermination,
        "Error                -- no computation, algorithm termination condition not met" } } },
  { "besselY",
    std::array<BesselErrors, 6>{
      BesselErrors{ Success, "Normal return        -- Computation completed." },
      BesselErrors{ InputError, "Input error          -- No computation" },
      BesselErrors{
        Overflow,
        "Overflow             -- No computation, order is too large or abs(z) is too small or both" },
      BesselErrors{ PartialLossOfSignificance,
                    "Loss of significance -- abs(z) or order large \n computation done but losses of "
                    "significance by argument reduction produce less than half of machine accuracy" },
      BesselErrors{ FullLossOfSignificance,
                    "Loss of significance -- abs(z) or order too large\n no computation because of "
                    "complete loss of significance by argument reduction" },
      BesselErrors{
        AlgorithmTermination,
        "Error                -- no computation, algorithm termination condition not met" } } },
  { "besselI",
    std::array<BesselErrors, 6>{
      BesselErrors{ Success, "Normal return        -- Computation completed." },
      BesselErrors{ InputError, "Input error          -- No computation" },
      BesselErrors{ Overflow,
                    "Overflow             -- No computation, Re(z) too large for scale=false" },
      BesselErrors{ PartialLossOfSignificance,
                    "Loss of significance -- abs(z) or order large \n computation done but losses of "
                    "significance by argument reduction produce less than half of machine accuracy" },
      BesselErrors{ FullLossOfSignificance,
                    "Loss of significance -- abs(z) or order too large\n no computation because of "
                    "complete loss of significance by argument reduction" },
      BesselErrors{
        AlgorithmTermination,
        "Error                -- no computation, algorithm termination condition not met" } } },
  { "besselK",
    std::array<BesselErrors, 6>{
      BesselErrors{ Success, "Normal return        -- Computation completed." },
      BesselErrors{ InputError, "Input error          -- No computation" },
      BesselErrors{
        Overflow,
        "Overflow             -- No computation, order is too large or abs(z) is too small or both" },
      BesselErrors{ PartialLossOfSignificance,
                    "Loss of significance -- abs(z) or order large \n computation done but losses of "
                    "significance by argument reduction produce less than half of machine accuracy" },
      BesselErrors{ FullLossOfSignificance,
                    "Loss of significance -- abs(z) or order too large\n no computation because of "
                    "complete loss of significance by argument reduction" },
      BesselErrors{
        AlgorithmTermination,
        "Error                -- no computation, algorithm termination condition not met" } } },
  { "hankelH",
    std::array<BesselErrors, 6>{
      BesselErrors{ Success, "Normal return        -- Computation completed." },
      BesselErrors{ InputError, "Input error          -- No computation" },
      BesselErrors{
        Overflow,
        "Overflow             -- No computation, order is too large or abs(z) too small or both" },
      BesselErrors{ PartialLossOfSignificance,
                    "Loss of significance -- abs(z) or order large \n computation done but losses of "
                    "significance by argument reduction produce less than half of machine accuracy" },
      BesselErrors{ FullLossOfSignificance,
                    "Loss of significance -- abs(z) or order too large\n no computation because of "
                    "complete loss of significance by argument reduction" },
      BesselErrors{
        AlgorithmTermination,
        "Error                -- no computation, algorithm termination condition not met" } } },
  { "airy",
    std::array<BesselErrors, 6>{
      BesselErrors{ Success, "Normal return        -- Computation completed." },
      BesselErrors{ InputError, "Input error          -- No computation" },
      BesselErrors{
        Overflow,
        "Overflow             -- No computation, Re(2/3*z*sqrt(z)) too large for scale=false" },
      BesselErrors{ PartialLossOfSignificance,
                    "Loss of significance -- abs(z) large \n computation done but losses of "
                    "significance by argument reduction produce less than half of machine accuracy" },
      BesselErrors{ FullLossOfSignificance,
                    "Loss of significance -- abs(z) too large\n no computation because of complete "
                    "loss of significance by argument reduction" },
      BesselErrors{
        AlgorithmTermination,
        "Error                -- no computation, algorithm termination condition not met" } } },
  { "biry",
    std::array<BesselErrors, 6>{
      BesselErrors{ Success, "Normal return        -- Computation completed." },
      BesselErrors{ InputError, "Input error          -- No computation" },
      BesselErrors{ Overflow,
                    "Overflow             -- No computation, Re(z) too large for scale=false" },
      BesselErrors{ PartialLossOfSignificance,
                    "Loss of significance -- abs(z) large \n computation done but losses of "
                    "significance by argument reduction produce less than half of machine accuracy" },
      BesselErrors{ FullLossOfSignificance,
                    "Loss of significance -- abs(z) too large\n no computation because of complete "
                    "loss of significance by argument reduction" },
      BesselErrors{
        AlgorithmTermination,
        "Error                -- no computation, algorithm termination condition not met" } } }
};

#endif // ERRORS_H
