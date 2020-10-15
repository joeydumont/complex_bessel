#include <array>
#include <boost/multi_array.hpp>
#include <complex_bessel.h>
#include <hdf5.h>
#include <iostream>

using namespace std::complex_literals;
using namespace sp_bessel;

static constexpr int    size    = 1000;
static constexpr int    nu_max  = 500;
static constexpr double exp_min = -15;
static constexpr double exp_max = 0;
static constexpr double dx      = (exp_max - exp_min) / size;

template<typename T>
struct hdf5_complex
{
  T real;
  T imag;
};

hid_t
get_hdf5_complex_type()
{

  hid_t type = H5Tcreate(H5T_COMPOUND, sizeof(hdf5_complex<double>));
  H5Tinsert(type, "r", HOFFSET(hdf5_complex<double>, real), H5T_NATIVE_DOUBLE);
  H5Tinsert(type, "i", HOFFSET(hdf5_complex<double>, imag), H5T_NATIVE_DOUBLE);
  return type;
}

int
main(void)
{

  // Prepare the HDF5 file
  hid_t  file, data_group, metadata_group, dataspace, dataset, dcpl, complex_type;
  hid_t  sumDset, orderDset, zDset;
  herr_t status;

  status = H5open();

  hsize_t dims[2] = { nu_max, size };

  file           = H5Fcreate("hankelOverflow.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  data_group     = H5Gcreate(file, "/BesselFunctionValues", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  metadata_group = H5Gcreate(file, "/Coordinates", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  // Create the datasets for the data group.
  dataspace    = H5Screate_simple(2, dims, NULL);
  dcpl         = H5Pcreate(H5P_DATASET_CREATE);
  complex_type = get_hdf5_complex_type();
  dataset = H5Dcreate(data_group, "hankelH2", complex_type, dataspace, H5P_DEFAULT, dcpl, H5P_DEFAULT);
  sumDset = H5Dcreate(data_group, "sumBessel", complex_type, dataspace, H5P_DEFAULT, dcpl, H5P_DEFAULT);

  // Create datasets for metadata group.
  hsize_t size_array[1] = { size };
  hid_t   dataspace_z   = H5Screate_simple(1, size_array, NULL);
  zDset =
    H5Dcreate(metadata_group, "z", H5T_NATIVE_DOUBLE, dataspace_z, H5P_DEFAULT, dcpl, H5P_DEFAULT);

  //
  hsize_t nu_max_array[1] = { nu_max };
  hid_t   dataspace_nu    = H5Screate_simple(1, nu_max_array, NULL);
  orderDset =
    H5Dcreate(metadata_group, "order", H5T_NATIVE_INT, dataspace_nu, H5P_DEFAULT, dcpl, H5P_DEFAULT);

  // Prepare arrays to store data and metadata.
  typedef boost::multi_array<int, 1>                  order_array_type;
  typedef boost::multi_array<double, 1>               metadata_array_type;
  typedef boost::multi_array<std::complex<double>, 2> data_array_type;

  metadata_array_type exponent(boost::extents[size]);
  metadata_array_type z(boost::extents[size]);
  order_array_type    nu(boost::extents[nu_max]);
  data_array_type     hankelH2Values(boost::extents[nu_max][size]);
  data_array_type     besselSumValues(boost::extents[nu_max][size]);
  data_array_type     difference(boost::extents[nu_max][size]);

  for (int i = 0; i < nu_max; i++) {
    nu[i] = i;
  }

  for (int j = 0; j < size; j++) {
    exponent[j] = exp_min + j * dx;
    z[j]        = std::pow(10, exponent[j]);
    std::cout << z[j] << "\n";
  }
  for (int i = 0; i < nu_max; i++) {
    for (int j = 0; j < size; j++) {
      hankelH2Values[i][j]  = hankelH2(i, z[j]);
      besselSumValues[i][j] = besselJ(i, z[j]) - 1i * besselY(i, z[j]);
      difference[i][j]      = hankelH2Values[i][j] - besselSumValues[i][j];
    }
  }

  // Write the data.
  status = H5Dwrite(zDset, H5T_NATIVE_DOUBLE, H5S_ALL, dataspace_z, H5P_DEFAULT, z.data());
  status = H5Dwrite(orderDset, H5T_NATIVE_INT, H5S_ALL, dataspace_nu, H5P_DEFAULT, nu.data());
  status = H5Dwrite(dataset, complex_type, H5S_ALL, dataspace, H5P_DEFAULT, hankelH2Values.data());
  status = H5Dwrite(sumDset, complex_type, H5S_ALL, dataspace, H5P_DEFAULT, besselSumValues.data());

  if (status) {
    std::cout << "There was an error." << std::endl;
  }

  return 0;
}