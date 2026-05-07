#include "rng.h"

#include <climits>
#include <cstdint>

namespace phyloacc {
namespace {

uint64_t SplitMix64(uint64_t value)
{
    value += 0x9e3779b97f4a7c15ULL;
    value = (value ^ (value >> 30)) * 0xbf58476d1ce4e5b9ULL;
    value = (value ^ (value >> 27)) * 0x94d049bb133111ebULL;
    return value ^ (value >> 31);
}

uint64_t AsUnsigned(int value)
{
    return static_cast<uint64_t>(static_cast<int64_t>(value));
}

void MixField(uint64_t& state, uint64_t value)
{
    state = SplitMix64(state ^ SplitMix64(value));
}

}  // namespace

unsigned long DeriveSeed(unsigned long base_seed,
                         ProgramKind program_kind,
                         int chain_index,
                         int element_index,
                         int model_or_block_index,
                         RngStream stream)
{
    uint64_t state = SplitMix64(0x7068796c6f616363ULL);
    MixField(state, static_cast<uint64_t>(base_seed));
    MixField(state, program_kind == ProgramKind::ST ? 0x5354ULL : 0x4754ULL);
    MixField(state, AsUnsigned(chain_index));
    MixField(state, AsUnsigned(element_index));
    MixField(state, AsUnsigned(model_or_block_index));
    MixField(state, static_cast<uint64_t>(stream));

    const unsigned long max_seed = ULONG_MAX;
    unsigned long seed = static_cast<unsigned long>((state % max_seed) + 1UL);
    return seed;
}

std::mt19937 MakeTwister(unsigned long seed)
{
    uint64_t state = SplitMix64(static_cast<uint64_t>(seed));
    std::seed_seq seq{
        static_cast<uint32_t>(state),
        static_cast<uint32_t>(state >> 32),
        static_cast<uint32_t>(SplitMix64(state)),
        static_cast<uint32_t>(SplitMix64(state) >> 32)
    };
    return std::mt19937(seq);
}

}  // namespace phyloacc
