#include <complex_bessel.h>
#include <gtest/gtest.h>
#include <iostream>
#include <vector>

using namespace std::complex_literals;
using namespace sp_bessel;

class BesselErrorsTest : public ::testing::Test
{
public:
  BesselErrorsTest() {}

  BesselErrors              besselErrors;
  std::vector<BesselErrors> besselErrors_vec;

protected:
  void SetUp() override {}
  void TearDown() override {}
};

TEST_F(BesselErrorsTest, AllFunctionsSuccess)
{
  auto z = std::complex<double>(1.0, 0.0);

  besselJ(0, z, false, &besselErrors);
  EXPECT_EQ(0, besselErrors.errorCode);
  EXPECT_EQ("Normal return        -- Computation completed.", besselErrors.errorMessage);

  besselY(0, z, false, &besselErrors);
  EXPECT_EQ(0, besselErrors.errorCode);
  EXPECT_EQ("Normal return        -- Computation completed.", besselErrors.errorMessage);

  besselK(0, z, false, &besselErrors);
  EXPECT_EQ(0, besselErrors.errorCode);
  EXPECT_EQ("Normal return        -- Computation completed.", besselErrors.errorMessage);

  besselI(0, z, false, &besselErrors);
  EXPECT_EQ(0, besselErrors.errorCode);
  EXPECT_EQ("Normal return        -- Computation completed.", besselErrors.errorMessage);

  hankelH1(0, z, false, &besselErrors);
  EXPECT_EQ(0, besselErrors.errorCode);
  EXPECT_EQ("Normal return        -- Computation completed.", besselErrors.errorMessage);

  hankelH2(0, z, false, &besselErrors);
  EXPECT_EQ(0, besselErrors.errorCode);
  EXPECT_EQ("Normal return        -- Computation completed.", besselErrors.errorMessage);

  airy(z, 0, false, &besselErrors);
  EXPECT_EQ(0, besselErrors.errorCode);
  EXPECT_EQ("Normal return        -- Computation completed.", besselErrors.errorMessage);

  biry(z, 0, false, &besselErrors);
  EXPECT_EQ(0, besselErrors.errorCode);
  EXPECT_EQ("Normal return        -- Computation completed.", besselErrors.errorMessage);

  besselJp(0, z, 1, &besselErrors_vec);
  for (auto element : besselErrors_vec) {
    EXPECT_EQ(0, element.errorCode);
  }
}

TEST_F(BesselErrorsTest, besselJOverflow)
{
  auto z = std::complex<double>(1000.0, 1000.0);

  besselJ(0, z, false, &besselErrors);
  EXPECT_EQ(2, besselErrors.errorCode);
  EXPECT_EQ(
    "Overflow             -- No computation, Im(z) too large for scale=false",
    besselErrors.errorMessage);
}

TEST_F(BesselErrorsTest, besselJpOverflow)
{
  auto z = std::complex<double>(1000.0, 1000.0);

  besselJp(0, z, 1, &besselErrors_vec);
  for (auto element : besselErrors_vec) {
    EXPECT_EQ(2, element.errorCode);
    EXPECT_EQ(
      "Overflow             -- No computation, Im(z) too large for scale=false", element.errorMessage);
  }
}

TEST_F(BesselErrorsTest, besselJPartialLossOfSignificance)
{
  auto z = std::complex<double>(35000.0, 0.0);

  besselJ(0, z, false, &besselErrors);
  EXPECT_EQ(3, besselErrors.errorCode);
  EXPECT_EQ(
    "Loss of significance -- abs(z) or order large \n computation done but losses of significance by "
    "argument reduction produce less than half of machine accuracy",
    besselErrors.errorMessage);
}

TEST_F(BesselErrorsTest, besselJpPartialLossOfSignificance)
{
  auto z = std::complex<double>(35000.0, 0.0);

  besselJp(0, z, 1, &besselErrors_vec);

  for (auto element : besselErrors_vec) {
    EXPECT_EQ(3, element.errorCode);
    EXPECT_EQ(
      "Loss of significance -- abs(z) or order large \n computation done but losses of significance by "
      "argument reduction produce less than half of machine accuracy",
      element.errorMessage);
  }
}

TEST_F(BesselErrorsTest, besselJFulllLossOfSignificance)
{
  auto z = std::complex<double>(1'500'000'000.0, 0.0);

  besselJ(0, z, false, &besselErrors);
  EXPECT_EQ(4, besselErrors.errorCode);
  EXPECT_EQ(
    "Loss of significance -- abs(z) or order too large\n no computation because of complete loss of "
    "significance by argument reduction",
    besselErrors.errorMessage);
}

TEST_F(BesselErrorsTest, besselJpFulllLossOfSignificance)
{
  auto z = std::complex<double>(1'500'000'000.0, 0.0);

  besselJp(0, z, 1, &besselErrors_vec);

  for (auto element : besselErrors_vec) {
    EXPECT_EQ(4, element.errorCode);
    EXPECT_EQ(
      "Loss of significance -- abs(z) or order too large\n no computation because of complete loss of "
      "significance by argument reduction",
      element.errorMessage);
  }
}

// Not sure how to trigger this error...
// TEST_F(BesselErrorsTest, besselJAlgorithmTermination) {
//   auto z = std::complex<double>(1'500'000'000.0,0.0);
//
//   besselJ(0, z, false, &besselErrors);
//   EXPECT_EQ(5, besselErrors.errorCode);
//   EXPECT_EQ("Error                -- no computation, algorithm termination condition not met",
//   besselErrors.errorMessage);
// }

TEST_F(BesselErrorsTest, besselYOverflow)
{
  auto z = std::complex<double>(10.0, 10.0);
  besselY(300, z, false, &besselErrors);
  EXPECT_EQ(2, besselErrors.errorCode);
  EXPECT_EQ(
    "Overflow             -- No computation, order is too large or abs(z) is too small or both",
    besselErrors.errorMessage);
}

TEST_F(BesselErrorsTest, besselYpOverflow)
{
  auto z = std::complex<double>(10.0, 10.0);
  besselYp(300, z, 1, &besselErrors_vec);

  for (auto element : besselErrors_vec) {
    EXPECT_EQ(2, element.errorCode);
    EXPECT_EQ(
      "Overflow             -- No computation, order is too large or abs(z) is too small or both",
      element.errorMessage);
  }
}

TEST_F(BesselErrorsTest, besselYPartialLossOfSignificance)
{
  auto z = std::complex<double>(35'000.0, 0.0);
  besselY(0, z, false, &besselErrors);
  EXPECT_EQ(3, besselErrors.errorCode);
  EXPECT_EQ(
    "Loss of significance -- abs(z) or order large \n computation done but losses of significance by "
    "argument reduction produce less than half of machine accuracy",
    besselErrors.errorMessage);
}

TEST_F(BesselErrorsTest, besselYpPartialLossOfSignificance)
{
  auto z = std::complex<double>(35'000.0, 0.0);
  besselYp(0, z, 1, &besselErrors_vec);

  for (auto element : besselErrors_vec) {
    EXPECT_EQ(3, besselErrors.errorCode);
    EXPECT_EQ(
      "Loss of significance -- abs(z) or order large \n computation done but losses of significance by "
      "argument reduction produce less than half of machine accuracy",
      element.errorMessage);
  }
}

TEST_F(BesselErrorsTest, besselYFullLossOfSignificance)
{
  auto z = std::complex<double>(1'500'000'000.0, 0.0);
  besselY(0, z, false, &besselErrors);
  EXPECT_EQ(4, besselErrors.errorCode);
  EXPECT_EQ(
    "Loss of significance -- abs(z) or order too large\n no computation because of complete loss of "
    "significance by argument reduction",
    besselErrors.errorMessage);
}

TEST_F(BesselErrorsTest, besselYpFullLossOfSignificance)
{
  auto z = std::complex<double>(1'500'000'000.0, 0.0);
  besselYp(0, z, 1, &besselErrors_vec);

  for (auto element : besselErrors_vec) {
    EXPECT_EQ(4, element.errorCode);
    EXPECT_EQ(
      "Loss of significance -- abs(z) or order too large\n no computation because of complete loss of "
      "significance by argument reduction",
      element.errorMessage);
  }
}

// Not sure how to trigger this error...
// TEST_F(BesselErrorsTest, besselYAlgorithmTermination) {
//   auto z = std::complex<double>(1'500'000'000.0,0.0);
//   besselY(0, z, false, &besselErrors);
//   EXPECT_EQ(5, besselErrors.errorCode);
//   EXPECT_EQ("Error                -- no computation, algorithm termination condition not met",
//   besselErrors.errorMessage);
//
// }

TEST_F(BesselErrorsTest, besselIOverflow)
{
  auto z = std::complex<double>(1'000.0, 10.0);
  besselI(0, z, false, &besselErrors);
  EXPECT_EQ(2, besselErrors.errorCode);
  EXPECT_EQ(
    "Overflow             -- No computation, Re(z) too large for scale=false",
    besselErrors.errorMessage);
}

TEST_F(BesselErrorsTest, besselIpOverflow)
{
  auto z = std::complex<double>(1'000.0, 10.0);
  besselIp(0, z, 1, &besselErrors_vec);

  for (auto element : besselErrors_vec) {
    EXPECT_EQ(2, element.errorCode);
    EXPECT_EQ(
      "Overflow             -- No computation, Re(z) too large for scale=false", element.errorMessage);
  }
}

TEST_F(BesselErrorsTest, besselIPartialLossOfSignificance)
{
  auto z = std::complex<double>(0.0, 35'000.0);
  besselI(0, z, false, &besselErrors);
  EXPECT_EQ(3, besselErrors.errorCode);
  EXPECT_EQ(
    "Loss of significance -- abs(z) or order large \n computation done but losses of significance by "
    "argument reduction produce less than half of machine accuracy",
    besselErrors.errorMessage);
}

TEST_F(BesselErrorsTest, besselIpPartialLossOfSignificance)
{
  auto z = std::complex<double>(0.0, 35'000.0);
  besselIp(0, z, 1, &besselErrors_vec);

  for (auto element : besselErrors_vec) {
    EXPECT_EQ(3, element.errorCode);
    EXPECT_EQ(
      "Loss of significance -- abs(z) or order large \n computation done but losses of significance by "
      "argument reduction produce less than half of machine accuracy",
      element.errorMessage);
  }
}

TEST_F(BesselErrorsTest, besselIFullLossOfSignificance)
{
  auto z = std::complex<double>(0.0, 1'500'000'000.0);
  besselI(0, z, false, &besselErrors);
  EXPECT_EQ(4, besselErrors.errorCode);
  EXPECT_EQ(
    "Loss of significance -- abs(z) or order too large\n no computation because of complete loss of "
    "significance by argument reduction",
    besselErrors.errorMessage);
}

TEST_F(BesselErrorsTest, besselIpFullLossOfSignificance)
{
  auto z = std::complex<double>(0.0, 1'500'000'000.0);
  besselIp(0, z, 1, &besselErrors_vec);

  for (auto element : besselErrors_vec) {
    EXPECT_EQ(4, element.errorCode);
    EXPECT_EQ(
      "Loss of significance -- abs(z) or order too large\n no computation because of complete loss of "
      "significance by argument reduction",
      element.errorMessage);
  }
}

// Not sure how to trigger this error...
// TEST_F(BesselErrorsTest, besselIAlgorithmTermination) {
//   auto z = std::complex<double>(1'500'000'000.0,0.0);
//   besselI(0, z, false, &besselErrors);
//   EXPECT_EQ(5, besselErrors.errorCode);
//   EXPECT_EQ("Error                -- no computation, algorithm termination condition not met",
//   besselErrors.errorMessage);
//
// }

TEST_F(BesselErrorsTest, besselKOverflow)
{
  auto z = std::complex<double>(10.0, 10.0);
  besselK(300, z, false, &besselErrors);
  EXPECT_EQ(2, besselErrors.errorCode);
  EXPECT_EQ(
    "Overflow             -- No computation, order is too large or abs(z) is too small or both",
    besselErrors.errorMessage);
}

TEST_F(BesselErrorsTest, besselKpOverflow)
{
  auto z = std::complex<double>(10.0, 10.0);
  besselKp(300, z, false, &besselErrors_vec);

  for (auto element : besselErrors_vec) {
    EXPECT_EQ(2, element.errorCode);
    EXPECT_EQ(
      "Overflow             -- No computation, order is too large or abs(z) is too small or both",
      element.errorMessage);
  }
}

TEST_F(BesselErrorsTest, besselKPartialLossOfSignificance)
{
  auto z = std::complex<double>(35'000.0, 0.0);
  besselK(0, z, false, &besselErrors);
  EXPECT_EQ(3, besselErrors.errorCode);
  EXPECT_EQ(
    "Loss of significance -- abs(z) or order large \n computation done but losses of significance by "
    "argument reduction produce less than half of machine accuracy",
    besselErrors.errorMessage);
}

TEST_F(BesselErrorsTest, besselKpPartialLossOfSignificance)
{
  auto z = std::complex<double>(35'000.0, 0.0);
  besselKp(0, z, false, &besselErrors_vec);

  for (auto element : besselErrors_vec) {
    EXPECT_EQ(3, element.errorCode);
    EXPECT_EQ(
      "Loss of significance -- abs(z) or order large \n computation done but losses of significance by "
      "argument reduction produce less than half of machine accuracy",
      element.errorMessage);
  }
}

TEST_F(BesselErrorsTest, besselKFullLossOfSignificance)
{
  auto z = std::complex<double>(1'500'000'000.0, 0.0);
  besselK(0, z, false, &besselErrors);
  EXPECT_EQ(4, besselErrors.errorCode);
  EXPECT_EQ(
    "Loss of significance -- abs(z) or order too large\n no computation because of complete loss of "
    "significance by argument reduction",
    besselErrors.errorMessage);
}

TEST_F(BesselErrorsTest, besselKpFullLossOfSignificance)
{
  auto z = std::complex<double>(1'500'000'000.0, 0.0);
  besselKp(0, z, false, &besselErrors_vec);

  for (auto element : besselErrors_vec) {
    EXPECT_EQ(4, element.errorCode);
    EXPECT_EQ(
      "Loss of significance -- abs(z) or order too large\n no computation because of complete loss of "
      "significance by argument reduction",
      element.errorMessage);
  }
}

// Not sure how to trigger this error...
// TEST_F(BesselErrorsTest, besselKAlgorithmTermination) {
//   auto z = std::complex<double>(1'500'000'000.0,0.0);
//   besselK(0, z, false, &besselErrors);
//   EXPECT_EQ(5, besselErrors.errorCode);
//   EXPECT_EQ("Error                -- no computation, algorithm termination condition not met",
//   besselErrors.errorMessage);
//
// }

TEST_F(BesselErrorsTest, hankelH1Overflow)
{
  auto z = std::complex<double>(10.0, 10.0);
  hankelH1(300, z, false, &besselErrors);
  EXPECT_EQ(2, besselErrors.errorCode);
  EXPECT_EQ(
    "Overflow             -- No computation, order is too large or abs(z) too small or both",
    besselErrors.errorMessage);
}

TEST_F(BesselErrorsTest, hankelH1pOverflow)
{
  auto z = std::complex<double>(10.0, 10.0);
  hankelH1p(300, z, false, &besselErrors_vec);

  for (auto element : besselErrors_vec) {
    EXPECT_EQ(2, element.errorCode);
    EXPECT_EQ(
      "Overflow             -- No computation, order is too large or abs(z) too small or both",
      element.errorMessage);
  }
}

TEST_F(BesselErrorsTest, hankelH1PartialLossOfSignificance)
{
  auto z = std::complex<double>(35'000.0, 0.0);
  hankelH1(0, z, false, &besselErrors);
  EXPECT_EQ(3, besselErrors.errorCode);
  EXPECT_EQ(
    "Loss of significance -- abs(z) or order large \n computation done but losses of significance by "
    "argument reduction produce less than half of machine accuracy",
    besselErrors.errorMessage);
}

TEST_F(BesselErrorsTest, hankelH1pPartialLossOfSignificance)
{
  auto z = std::complex<double>(35'000.0, 0.0);
  hankelH1p(0, z, false, &besselErrors_vec);

  for (auto element : besselErrors_vec) {
    EXPECT_EQ(3, element.errorCode);
    EXPECT_EQ(
      "Loss of significance -- abs(z) or order large \n computation done but losses of significance by "
      "argument reduction produce less than half of machine accuracy",
      element.errorMessage);
  }
}

TEST_F(BesselErrorsTest, hankelH1FullLossOfSignificance)
{
  auto z = std::complex<double>(1'500'000'000.0, 0.0);
  hankelH1(0, z, false, &besselErrors);
  EXPECT_EQ(4, besselErrors.errorCode);
  EXPECT_EQ(
    "Loss of significance -- abs(z) or order too large\n no computation because of complete loss of "
    "significance by argument reduction",
    besselErrors.errorMessage);
}

TEST_F(BesselErrorsTest, hankelH1pFullLossOfSignificance)
{
  auto z = std::complex<double>(1'500'000'000.0, 0.0);
  hankelH1p(0, z, false, &besselErrors_vec);

  for (auto element : besselErrors_vec) {
    EXPECT_EQ(4, element.errorCode);
    EXPECT_EQ(
      "Loss of significance -- abs(z) or order too large\n no computation because of complete loss of "
      "significance by argument reduction",
      element.errorMessage);
  }
}

// Not sure how to trigger this error...
// TEST_F(BesselErrorsTest, hankelH1AlgorithmTermination) {
//   auto z = std::complex<double>(1'500'000'000.0,0.0);
//   hankelH1(0, z, false, &besselErrors);
//   EXPECT_EQ(5, besselErrors.errorCode);
//   EXPECT_EQ("Error                -- no computation, algorithm termination condition not met",
//   besselErrors.errorMessage);
//
// }

TEST_F(BesselErrorsTest, hankelH2Overflow)
{
  auto z = std::complex<double>(10.0, 10.0);
  hankelH2(300, z, false, &besselErrors);
  EXPECT_EQ(2, besselErrors.errorCode);
  EXPECT_EQ(
    "Overflow             -- No computation, order is too large or abs(z) too small or both",
    besselErrors.errorMessage);
}

TEST_F(BesselErrorsTest, hankelH2pOverflow)
{
  auto z = std::complex<double>(10.0, 10.0);
  hankelH2p(300, z, false, &besselErrors_vec);

  for (auto element : besselErrors_vec) {
    EXPECT_EQ(2, element.errorCode);
    EXPECT_EQ(
      "Overflow             -- No computation, order is too large or abs(z) too small or both",
      element.errorMessage);
  }
}

TEST_F(BesselErrorsTest, hankelH2PartialLossOfSignificance)
{
  auto z = std::complex<double>(35'000.0, 0.0);
  hankelH2(0, z, false, &besselErrors);
  EXPECT_EQ(3, besselErrors.errorCode);
  EXPECT_EQ(
    "Loss of significance -- abs(z) or order large \n computation done but losses of significance by "
    "argument reduction produce less than half of machine accuracy",
    besselErrors.errorMessage);
}

TEST_F(BesselErrorsTest, hankelH2pPartialLossOfSignificance)
{
  auto z = std::complex<double>(35'000.0, 0.0);
  hankelH2p(0, z, false, &besselErrors_vec);

  for (auto element : besselErrors_vec) {
    EXPECT_EQ(3, element.errorCode);
    EXPECT_EQ(
      "Loss of significance -- abs(z) or order large \n computation done but losses of significance by "
      "argument reduction produce less than half of machine accuracy",
      element.errorMessage);
  }
}

TEST_F(BesselErrorsTest, hankelH2FullLossOfSignificance)
{
  auto z = std::complex<double>(1'500'000'000.0, 0.0);
  hankelH2(0, z, false, &besselErrors);
  EXPECT_EQ(4, besselErrors.errorCode);
  EXPECT_EQ(
    "Loss of significance -- abs(z) or order too large\n no computation because of complete loss of "
    "significance by argument reduction",
    besselErrors.errorMessage);
}

TEST_F(BesselErrorsTest, hankelH2pFullLossOfSignificance)
{
  auto z = std::complex<double>(1'500'000'000.0, 0.0);
  hankelH2p(0, z, false, &besselErrors_vec);

  for (auto element : besselErrors_vec) {
    EXPECT_EQ(4, element.errorCode);
    EXPECT_EQ(
      "Loss of significance -- abs(z) or order too large\n no computation because of complete loss of "
      "significance by argument reduction",
      element.errorMessage);
  }
}

// Not sure how to trigger this error...
// TEST_F(BesselErrorsTest, hankelH2AlgorithmTermination) {
//   auto z = std::complex<double>(1'500'000'000.0,0.0);
//   hankelH2(0, z, false, &besselErrors);
//   EXPECT_EQ(5, besselErrors.errorCode);
//   EXPECT_EQ("Error                -- no computation, algorithm termination condition not met",
//   besselErrors.errorMessage);
//
// }

TEST_F(BesselErrorsTest, airyOverflow)
{
  auto z = std::complex<double>(-1'000.0, 1'000.0);
  airy(z, 0, false, &besselErrors);
  EXPECT_EQ(2, besselErrors.errorCode);
  EXPECT_EQ(
    "Overflow             -- No computation, Re(2/3*z*sqrt(z)) too large for scale=false",
    besselErrors.errorMessage);
}

TEST_F(BesselErrorsTest, airyPartialLossOfSignificance)
{
  auto z = std::complex<double>(35'000.0, 0.0);
  airy(z, 0, false, &besselErrors);
  EXPECT_EQ(3, besselErrors.errorCode);
  EXPECT_EQ(
    "Loss of significance -- abs(z) large \n computation done but losses of significance by argument "
    "reduction produce less than half of machine accuracy",
    besselErrors.errorMessage);
}

TEST_F(BesselErrorsTest, airyFullLossOfSignificance)
{
  auto z = std::complex<double>(1'500'000'000.0, 0.0);
  airy(z, 0, false, &besselErrors);
  EXPECT_EQ(4, besselErrors.errorCode);
  EXPECT_EQ(
    "Loss of significance -- abs(z) too large\n no computation because of complete loss of "
    "significance by argument reduction",
    besselErrors.errorMessage);
}

// Not sure how to trigger this error...
// TEST_F(BesselErrorsTest, airyAlgorithmTermination) {
//   auto z = std::complex<double>(1'500'000'000.0,0.0);
//   besselK(0, z, false, &besselErrors);
//   EXPECT_EQ(5, besselErrors.errorCode);
//   EXPECT_EQ("Error                -- no computation, algorithm termination condition not met",
//   besselErrors.errorMessage);
//
// }

TEST_F(BesselErrorsTest, biryOverflow)
{
  auto z = std::complex<double>(10'000.0, 0.0);
  biry(z, 0, false, &besselErrors);
  EXPECT_EQ(2, besselErrors.errorCode);
  EXPECT_EQ(
    "Overflow             -- No computation, Re(z) too large for scale=false",
    besselErrors.errorMessage);
}

// Not sure how to trigger this error. It returns the overflow error, then switches
// to full loss of significance. Should test against more precise library.
// TEST_F(BesselErrorsTest, biryPartialLossOfSignificance) {
//  auto z = std::complex<double>(1'000'000,100.0);
//  biry(z, 0, false, &besselErrors);
//  EXPECT_EQ(3, besselErrors.errorCode);
//  EXPECT_EQ("Loss of significance -- abs(z) large \n computation done but losses of significance by
//  argument reduction produce less than half of machine accuracy", besselErrors.errorMessage);
//}

TEST_F(BesselErrorsTest, biryFullLossOfSignificance)
{
  auto z = std::complex<double>(1'500'000'000.0, 0.0);
  biry(z, 0, false, &besselErrors);
  EXPECT_EQ(4, besselErrors.errorCode);
  EXPECT_EQ(
    "Loss of significance -- abs(z) too large\n no computation because of complete loss of "
    "significance by argument reduction",
    besselErrors.errorMessage);
}

// Not sure how to trigger this error...
// TEST_F(BesselErrorsTest, biryAlgorithmTermination) {
//   auto z = std::complex<double>(1'500'000'000.0,0.0);
//   besselK(0, z, false, &besselErrors);
//   EXPECT_EQ(5, besselErrors.errorCode);
//   EXPECT_EQ("Error                -- no computation, algorithm termination condition not met",
//   besselErrors.errorMessage);
//
// }

int
main(int argc, char* argv[])
{
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
