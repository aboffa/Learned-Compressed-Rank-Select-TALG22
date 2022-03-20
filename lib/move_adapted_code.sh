# include s18_vector implementation
curl https://raw.githubusercontent.com/mudetz/s18_vector/master/s18/head/s18_vector.hpp > lib/la_vector/lib/sdsl-lite/include/sdsl/s18_vector.hpp
curl https://raw.githubusercontent.com/mudetz/s18_vector/master/s18/head/constants.hpp > lib/la_vector/lib/sdsl-lite/include/sdsl/constants.hpp

# include sd_vector and rle_vector from https://github.com/vgteam/sdsl-lite
curl https://raw.githubusercontent.com/vgteam/sdsl-lite/master/lib/simple_sds.cpp > lib/la_vector/lib/sdsl-lite/lib/simple_sds.cpp
curl https://raw.githubusercontent.com/vgteam/sdsl-lite/master/include/sdsl/simple_sds.hpp > lib/la_vector/lib/sdsl-lite/include/sdsl/simple_sds.hpp
curl https://raw.githubusercontent.com/vgteam/sdsl-lite/master/include/sdsl/int_vector.hpp > lib/la_vector/lib/sdsl-lite/include/sdsl/int_vector.hpp
curl https://raw.githubusercontent.com/vgteam/sdsl-lite/master/lib/sd_vector.cpp > lib/la_vector/lib/sdsl-lite/lib/sd_vector.cpp
curl https://raw.githubusercontent.com/vgteam/sdsl-lite/master/include/sdsl/sd_vector.hpp > lib/la_vector/lib/sdsl-lite/include/sdsl/sd_vector.hpp
curl https://raw.githubusercontent.com/vgteam/sdsl-lite/master/include/sdsl/rle_vector.hpp > lib/la_vector/lib/sdsl-lite/include/sdsl/rle_vector.hpp
curl https://raw.githubusercontent.com/vgteam/sdsl-lite/master/include/sdsl/bits.hpp > lib/la_vector/lib/sdsl-lite/include/sdsl/bits.hpp

# move adapted code
cp lib/adapted_code/csa_wt.hpp lib/la_vector/lib/sdsl-lite/include/sdsl
cp lib/adapted_code/la_vector.hpp lib/la_vector/include
cp lib/adapted_code/piecewise_linear_model.hpp lib/la_vector/include

