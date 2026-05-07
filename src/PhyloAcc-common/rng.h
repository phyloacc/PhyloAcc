#ifndef PHYLOACC_COMMON_RNG_H
#define PHYLOACC_COMMON_RNG_H

#include "run.h"

#include <random>

namespace phyloacc {

enum class RngStream {
    WorkerGsl = 1,
    GtSiteShuffle = 2,
    GtGeneTreeShuffle = 3
};

unsigned long DeriveSeed(unsigned long base_seed,
                         ProgramKind program_kind,
                         int chain_index,
                         int element_index,
                         int model_or_block_index,
                         RngStream stream);

std::mt19937 MakeTwister(unsigned long seed);

}  // namespace phyloacc

#endif
