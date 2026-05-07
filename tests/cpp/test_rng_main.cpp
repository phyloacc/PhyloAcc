#include "../../src/PhyloAcc-common/rng.h"

#include <cassert>
#include <iostream>
#include <set>

int main()
{
    static_assert(sizeof(unsigned long) >= 8, "Seed tests expect 64-bit unsigned long.");

    unsigned long st_seed = phyloacc::DeriveSeed(
        1, phyloacc::ProgramKind::ST, 0, 0, 0, phyloacc::RngStream::WorkerGsl);
    assert(st_seed == 15330747418239297920UL);
    assert(st_seed == phyloacc::DeriveSeed(
        1, phyloacc::ProgramKind::ST, 0, 0, 0, phyloacc::RngStream::WorkerGsl));

    std::set<unsigned long> seeds;
    seeds.insert(st_seed);
    seeds.insert(phyloacc::DeriveSeed(
        1, phyloacc::ProgramKind::ST, 2, 0, 0, phyloacc::RngStream::WorkerGsl));
    seeds.insert(phyloacc::DeriveSeed(
        1, phyloacc::ProgramKind::GT, 0, 3, 20, phyloacc::RngStream::WorkerGsl));
    seeds.insert(phyloacc::DeriveSeed(
        1, phyloacc::ProgramKind::GT, 0, 3, 20, phyloacc::RngStream::GtSiteShuffle));
    seeds.insert(phyloacc::DeriveSeed(
        1, phyloacc::ProgramKind::GT, 0, 3, 20, phyloacc::RngStream::GtGeneTreeShuffle));
    seeds.insert(phyloacc::DeriveSeed(
        12345, phyloacc::ProgramKind::GT, 1, 4, 7, phyloacc::RngStream::GtGeneTreeShuffle));
    assert(seeds.size() == 6);

    assert(phyloacc::DeriveSeed(
        1, phyloacc::ProgramKind::ST, 2, 0, 0, phyloacc::RngStream::WorkerGsl)
        == 3334740700155033886UL);
    assert(phyloacc::DeriveSeed(
        1, phyloacc::ProgramKind::GT, 0, 3, 20, phyloacc::RngStream::WorkerGsl)
        == 13519517523178014377UL);
    assert(phyloacc::DeriveSeed(
        1, phyloacc::ProgramKind::GT, 0, 3, 20, phyloacc::RngStream::GtSiteShuffle)
        == 5837018277079064570UL);
    assert(phyloacc::DeriveSeed(
        1, phyloacc::ProgramKind::GT, 0, 3, 20, phyloacc::RngStream::GtGeneTreeShuffle)
        == 10230629784440138387UL);
    assert(phyloacc::DeriveSeed(
        12345, phyloacc::ProgramKind::GT, 1, 4, 7, phyloacc::RngStream::GtGeneTreeShuffle)
        == 76732796917337790UL);

    std::mt19937 first = phyloacc::MakeTwister(st_seed);
    std::mt19937 second = phyloacc::MakeTwister(st_seed);
    assert(first() == second());

    std::cout << "RNG C++ unit tests passed." << std::endl;
    return 0;
}
