#pragma once

#include <boost/preprocessor/seq/for_each.hpp>
#include <boost/preprocessor/stringize.hpp>
#include <boost/preprocessor/cat.hpp>

//#include "freq_index.hpp"
//#include "positive_sequence.hpp"
#include "partitioned_sequence.hpp"
#include "uniform_partitioned_sequence.hpp"
//#include "binary_freq_collection.hpp"
//#include "block_freq_index.hpp"
//#include "block_codecs.hpp"
//#include "mixed_block.hpp"


#define DS2I_INDEX_TYPES (ef)(single)(uniform)(opt)(block_optpfor)(block_varint)(block_interpolative)(block_mixed)(block_qmx)
#define DS2I_BLOCK_INDEX_TYPES (block_optpfor)(block_varint)(block_interpolative)(block_qmx)(block_mixed)
