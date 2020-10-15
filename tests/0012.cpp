#include <complex_bessel.h>
#include <gtest/gtest.h>
#include <random>

using namespace std::complex_literals;
using namespace sp_bessel;

class BesselScalingTest : public ::testing::Test
{
public:
  BesselScalingTest()
  {
    std::random_device               rd;
    std::mt19937                     gen(rd());
    std::uniform_int_distribution<>  order_dist(-100, 100);
    std::uniform_real_distribution<> z_dist(-500.0, 500.0);

    order     = order_dist(gen);
    double zr = z_dist(gen);
    double zi = z_dist(gen);
    z         = std::complex<double>(zr, zi);
  }

  int                  order;
  std::complex<double> z;

protected:
  void SetUp() override {}
  void TearDown() override {}
};

TEST_F(BesselScalingTest, besselJ)
{
  // For besselJ, the expected scaling behaviour is EXP(-ABS(IMAG(Z)))*besselJ(nu,z)
  std::complex<double> expected = std::exp(-std::abs(std::imag(z))) * besselJ(order, z);
  std::complex<double> actual   = besselJ(order, z, true);

  // Compare real and imag parts separately.
  EXPECT_NEAR(std::real(expected), std::real(actual), 1.0e-3)
    << "BESSELJ scaling failed with z " << z << "\n";
  EXPECT_NEAR(std::imag(expected), std::imag(actual), 1.0e-3)
    << "BESSELJ scaling failed with z " << z << "\n";
}

TEST_F(BesselScalingTest, besselY)
{
  // For besselY, the expected scaling behaviour is EXP(-ABS(IMAG(Z)))*besselI(nu,z)
  std::complex<double> expected = std::exp(-std::abs(std::imag(z))) * besselY(order, z);
  std::complex<double> actual   = besselY(order, z, true);

  // Compare real and imag parts separately.
  EXPECT_NEAR(std::real(expected), std::real(actual), 1.0e-3)
    << "BESSELY scaling failed with z " << z << "\n";
  EXPECT_NEAR(std::imag(expected), std::imag(actual), 1.0e-3)
    << "BESSELY scaling failed with z " << z << "\n";
}

TEST_F(BesselScalingTest, besselI)
{
  // For besselI, the expected scaling behaviour is EXP(-ABS(REAL(Z)))*besseI(nu,z)
  order                         = std::abs(order);
  std::complex<double> expected = std::exp(-std::abs(std::real(z))) * besselI(order, z);
  std::complex<double> actual   = besselI(order, z, true);

  // Compare real and imag parts separately.
  EXPECT_NEAR(std::real(expected), std::real(actual), 1.0e-3)
    << "BESSELI scaling failed with z " << z << "\n";
  EXPECT_NEAR(std::imag(expected), std::imag(actual), 1.0e-3)
    << "BESSELI scaling failed with z " << z << "\n";
}

TEST_F(BesselScalingTest, besselK)
{
  // For besselK, the expected scaling behaviour is EXP(Z))*besselJ(nu,z)
  std::complex<double> expected = std::exp(z) * besselK(order, z);
  std::complex<double> actual   = besselK(order, z, true);

  // Compare real and imag parts separately.
  EXPECT_NEAR(std::real(expected), std::real(actual), 1.0e-3)
    << "BESSELK scaling failed with z " << z << "\n";
  EXPECT_NEAR(std::imag(expected), std::imag(actual), 1.0e-3)
    << "BESSELK scaling failed with z " << z << "\n";
}

TEST_F(BesselScalingTest, hankelH1)
{
  // For hankelH1, the expected scaling behaviour is EXP(-i*z))*hankelH1(nu,z)
  std::complex<double> expected = std::exp(-1i * z) * hankelH1(order, z);
  std::complex<double> actual   = hankelH1(order, z, true);

  // Compare real and imag parts separately.
  EXPECT_NEAR(std::real(expected), std::real(actual), 1.0e-3)
    << "HANKELH1 scaling failed with z " << z << "\n";
  EXPECT_NEAR(std::imag(expected), std::imag(actual), 1.0e-3)
    << "HANKELH1 scaling failed with z " << z << "\n";
}

TEST_F(BesselScalingTest, hankelH2)
{
  // For hankelH1, the expected scaling behaviour is EXP(i*z))*hankelH1(nu,z)
  std::complex<double> expected = std::exp(1i * z) * hankelH2(order, z);
  std::complex<double> actual   = hankelH2(order, z, true);

  // Compare real and imag parts separately.
  EXPECT_NEAR(std::real(expected), std::real(actual), 1.0e-3)
    << "HANKELH2 scaling failed with z " << z << "\n";
  EXPECT_NEAR(std::imag(expected), std::imag(actual), 1.0e-3)
    << "HANKELH2 scaling failed with z " << z << "\n";
}

int
main(int argc, char* argv[])
{
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
