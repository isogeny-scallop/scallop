
#include <iostream>
#include <iomanip>
#include <sstream>
#include <cassert>

////////////////////////////////////////////////////////////////

#include "params.hpp"   // from params.sage, includes bigint.hpp/fp.hpp/fq.hpp

#include "ec.hpp"
#include "ecp.hpp"
#include "isog.hpp"

#include "fact.hpp"
#include "ecgrp.hpp"

#include "action.hpp"

////////////////////////////////////////////////////////////////

constexpr auto p1 = p + 1;
using p1fact = _fact<p1>;
using L1fact = _fact<L1>;
using L2fact = _fact<L2>;
using L3fact = _fact<L3>;

auto const E0 = std::make_shared<ec<fq>>(a4, a6);
oriented_curve<fq> const _EPQ0 {E0, E0->lift_x(xP).value(), E0->lift_x(xQ).value()};

////////////////////////////////////////////////////////////////

static __inline__ uint64_t rdtsc(void)
{
    uint32_t hi, lo;
    __asm__ __volatile__ ("rdtsc" : "=a"(lo), "=d"(hi));
    return lo | (uint64_t) hi << 32;
}

#include <random>

int main()
{
    auto seed = std::random_device()() + time(nullptr);
    std::cerr << "\x1b[35mrandom seed: " << seed << "\x1b[0m" << std::endl;
    boost::random::minstd_rand rng(seed);

    using scalar_t = std::remove_cvref<decltype(hh)>::type;
    boost::random::uniform_int_distribution<scalar_t> const dist(0, hh-1);

    oriented_curve<fq> const EPQ0 = canonicalize<p1fact,L1fact,L2fact,L3fact>(_EPQ0);
#ifndef NDEBUG
    {
        auto again = canonicalize<p1fact,L1fact,L2fact,L3fact>(EPQ0);
        assert(again == EPQ0);  // idempotent
    }
#endif

    std::cout << "\nE0a = " << EPQ0.E->a();
    std::cout << "\nE0b = " << EPQ0.E->b();
    std::cout << "\nP0x = " << EPQ0.P.affine_x();
    std::cout << "\nQ0x = " << EPQ0.Q.affine_x();
    std::cout << std::endl;

    scalar_t v = dist(rng);
    std::cout << "\n v  = " << v << std::endl;

    auto EPQ1 = scallop<p1fact,L1fact,L2fact,L3fact,relations>(EPQ0, v);
    std::cout << "\nE1a = " << EPQ1.E->a();
    std::cout << "\nE1b = " << EPQ1.E->b();
    std::cout << "\nP1x = " << EPQ1.P.affine_x();
    std::cout << "\nQ1x = " << EPQ1.Q.affine_x();
    std::cout << std::endl;

#ifndef NDEBUG
    {
        auto again = canonicalize<p1fact,L1fact,L2fact,L3fact>(EPQ1);
        assert(again == EPQ1);  // idempotent
    }
#endif

    scalar_t w = dist(rng);
    std::cout << "\n w  = " << w << std::endl;

    auto EPQ2 = scallop<p1fact,L1fact,L2fact,L3fact,relations>(EPQ1, w);
    std::cout << "\nE2a = " << EPQ2.E->a();
    std::cout << "\nE2b = " << EPQ2.E->b();
    std::cout << "\nP2x = " << EPQ2.P.affine_x();
    std::cout << "\nQ2x = " << EPQ2.Q.affine_x();
    std::cout << std::endl;
#ifndef NDEBUG
    {
        auto again = canonicalize<p1fact,L1fact,L2fact,L3fact>(EPQ2);
        assert(again == EPQ2);  // idempotent
    }
#endif

    scalar_t vw = (v + w) % hh;
    std::cout << "\nv+w = " << vw << std::endl;

    auto EPQ3 = scallop<p1fact,L1fact,L2fact,L3fact,relations>(EPQ0, vw);
    std::cout << "\nE3a = " << EPQ3.E->a();
    std::cout << "\nE3b = " << EPQ3.E->b();
    std::cout << "\nP3x = " << EPQ3.P.affine_x();
    std::cout << "\nQ3x = " << EPQ3.Q.affine_x();
    std::cout << std::endl;
#ifndef NDEBUG
    {
        auto again = canonicalize<p1fact,L1fact,L2fact,L3fact>(EPQ3);
        assert(again == EPQ3);  // idempotent
    }
#endif

    if (EPQ2 == EPQ3) {
        std::cout << "\n\x1b[32mequal! :)\x1b[0m\n" << std::endl;
        return 0;
    }
    else {
        std::cout << "\n\x1b[31mNOT EQUAL!! :(\x1b[0m\n" << std::endl;
        return 1;
    }
}

