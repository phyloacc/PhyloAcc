#include "../../src/PhyloAcc-common/rng.h"

#include <cassert>
#include <iostream>
#include <set>

int main()
{
    static_assert(sizeof(unsigned long) >= 8, "Seed tests expect 64-bit unsigned long.");

    unsigned long st_seed = phyloacc::DeriveSeed(
        1, phyloacc::ProgramKind::ST, 0, 0, 0, phyloacc::RngStream::WorkerGsl);
    unsigned long st_run_seed = phyloacc::DeriveSeed(
        1, phyloacc::ProgramKind::ST, 0, -1, 0, phyloacc::RngStream::RunGsl);
    assert(st_seed == 15330747418239297920UL);
    assert(st_seed == phyloacc::DeriveSeed(
        1, phyloacc::ProgramKind::ST, 0, 0, 0, phyloacc::RngStream::WorkerGsl));
    assert(st_run_seed == phyloacc::DeriveSeed(
        1, phyloacc::ProgramKind::ST, 0, -1, 0, phyloacc::RngStream::RunGsl));

    std::set<unsigned long> seeds;
    seeds.insert(st_run_seed);
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
    assert(seeds.size() == 7);

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

    assert(phyloacc::DeriveSeed(
        1, phyloacc::ProgramKind::ST, 0, -1, 0, phyloacc::RngStream::RunGsl)
        == 3458786590926809558UL);

    std::mt19937 first = phyloacc::MakeTwister(st_seed);
    std::mt19937 second = phyloacc::MakeTwister(st_seed);
    std::mt19937 different = phyloacc::MakeTwister(st_run_seed);
    assert(first() == second());
    assert(first() != different());

    std::cout << "RNG C++ unit tests passed." << std::endl;
    return 0;
}
