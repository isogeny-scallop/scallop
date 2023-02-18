
#pragma once

#include "fact.hpp"

#include <tuple>
#include <optional>

template <auto const &NORM, auto const &TRACE>
constexpr std::optional<std::pair<uint64_t,uint64_t>> eigenvalues(uint64_t ell)
{
    assert(ell < ((uint64_t) 1 << 27));  // or something

    uint64_t const norm  = static_cast<uint64_t>((NORM  % ell + ell) % ell),
                   trace = static_cast<uint64_t>((TRACE % ell + ell) % ell);

    std::pair<uint64_t,uint64_t> ret {-1,-1};
    size_t cnt = 0;
    for (uint64_t i = 0; cnt < 2 && i < ell; ++i) {
        if (((i + ell-trace)*i + norm) % ell == 0)
            (cnt++ ? ret.second : ret.first) = i;
    }
    if (cnt == 2)
        return ret;
    return {};
}

template <auto const &NORM, auto const &TRACE, auto const &ELLS>
consteval std::array<std::pair<uint64_t,uint64_t>, ELLS.size()> eigenvalues()
{
    constexpr auto ells = ELLS;
    std::array<std::pair<uint64_t,uint64_t>, ells.size()> ret;
    std::transform(ells.begin(), ells.end(), ret.begin(),
            [](uint64_t ell) {
                auto const r = eigenvalues<NORM, TRACE>(ell);
                if (!r)
                    throw std::logic_error("not split");
                return *r;
            });
    return ret;
}

////////////////////////////////////////////////////////////////

#include "isog.hpp"
#include "ecgrp.hpp"

template <class EXP, class O, class F>
std::pair<isog_chain<F>, isog_chain<F>> isogeny_with_dual(ecp<F> const &K, O const &ord, F const &u = F::one)
{
    auto phis = factored_isogeny(K, ord, u);

    assert(order<EXP>(K) == ord);
    auto P = complete_basis(K, EXP(ord));
    assert(order<EXP>(P) == ord);
//    assert(independence_checker<EXP>(K, ord)(P));

    auto psis = factored_isogeny(chain_eval(phis, P), ord);

    auto v = static_cast<F>(ord.value()) * u.inv();
    psis.back()._scale(v, phis.front().domain());

#ifndef NDEBUG
auto T = K.curve()->random_point();
assert(chain_eval(psis, chain_eval(phis, T)) == ord.value() * T);
#endif

    return {phis, psis};
}

////////////////////////////////////////////////////////////////

// P,Q basis of l-torsion, l prime
// wP,wQ images of P,Q under omega
// u,v eigenvalues of omega modulo l
// --> order-l point lying in u-eigenspace of omega
template <class F>
ecp<F> kernel_from_omega(ecp<F> const & P, ecp<F> const & Q,
                         ecp<F> const &wP, ecp<F> const &wQ,
                         uint64_t v)
{
    if (auto R = wP - v*P)
        return R;
    return wQ - v*Q;
}

template <class F>
struct oriented_curve
{
    std::shared_ptr<ec<F> const> E;
    ecp<F> P, Q;

    // NB: simply compares members, no canonicalization
    bool operator==(oriented_curve const &other) const
    { return *E == *other.E && P == other.P && Q == other.Q; }
};

template <class EXP, class L1, class L2, class L3, class F>
oriented_curve<F> scallop_steps(oriented_curve<F> const &EPQ,
                               L1 const &B1, L1 const &C1,
                               L2 const &B2, L2 const &C2,
                               L3 const &B3, L3 const &C3)
{
#ifndef NDEBUG
    for (size_t i = 0; i < L1::count; ++i)
        assert(L1::top[i] == 1 && B1[i] <= 1 && C1[i] <= 1 && !(B1[i] && C1[i]));
    for (size_t i = 0; i < L2::count; ++i)
        assert(L2::top[i] == 1 && B2[i] <= 1 && C2[i] <= 1 && !(B2[i] && C2[i]));
    for (size_t i = 0; i < L3::count; ++i)
        assert(L3::top[i] == 1 && B3[i] <= 1 && C3[i] <= 1 && !(B3[i] && C3[i]));
#endif

    constexpr EXP LL1(L1::top), LL2(L2::top), LL3(L3::top);
    constexpr static auto LL12val = (LL1 * LL2).value();
    constexpr static auto LL23val = (LL2 * LL3).value();
    using L12 = _fact<LL12val>;
    using L23 = _fact<LL23val>;
    auto const B12 = L12(B1)*L12(B2);

    assert(order<EXP>(EPQ.P) == LL1*LL2);
    assert(order<EXP>(EPQ.Q) == LL1);

    auto KC1 = (L1::top / C1).value() * EPQ.Q;
    assert(order<EXP>(KC1) == EXP(C1));
    auto [phi_E_C1, phi_E_C1_hat] = isogeny_with_dual<EXP>(KC1, C1);
    auto P_EC1 = chain_eval(phi_E_C1, EPQ.P),
         Q_EC1 = chain_eval(phi_E_C1, EPQ.Q);
    auto EC1 = chain_codomain(phi_E_C1);

    auto [P3,Q3] = basis(*EC1, LL3);
    auto P2 = L1::topval * P_EC1;
    assert(order<EXP>(P2) == LL2);
    auto Q2 = complete_basis(P2, LL2);
    assert(order<EXP>(Q2) == LL2);
//    assert(independence_checker<L2>(P2, L2::top)(Q2));

    auto wQ2 = chain_eval(phi_E_C1_hat, Q2);
    auto wP3 = chain_eval(phi_E_C1_hat, P3);
    auto wQ3 = chain_eval(phi_E_C1_hat, Q3);

    auto KB1B2 = (L12::top/B12).value() * EPQ.P;
    assert(order<EXP>(KB1B2) == EXP(B1)*EXP(B2));
    auto phi_E_B1B2 = factored_isogeny(KB1B2, B12);
    auto P_EB1B2 = chain_eval(phi_E_B1B2, EPQ.P);
    auto Q_EB1B2 = chain_eval(phi_E_B1B2, EPQ.Q);

    wQ2 = chain_eval(phi_E_B1B2, wQ2);
    wP3 = chain_eval(phi_E_B1B2, wP3);
    wQ3 = chain_eval(phi_E_B1B2, wQ3);

    auto const &KL1L2_B1B2 = P_EB1B2;
    assert(order<EXP>(KL1L2_B1B2) == EXP(L12::top/B12));
    auto [phi_E_L1L2, phi_E_L1L2_hat] = isogeny_with_dual<EXP>(KL1L2_B1B2, L12::top/B12);
    auto Q_EL1L2 = chain_eval(phi_E_L1L2, Q_EB1B2);
    auto EL1L2 = chain_codomain(phi_E_L1L2);

    wQ2 = chain_eval(phi_E_L1L2, wQ2);
    wP3 = chain_eval(phi_E_L1L2, wP3);
    wQ3 = chain_eval(phi_E_L1L2, wQ3);

    auto const &KL1_C1 = Q_EC1;
    assert(order<EXP>(KL1_C1) == EXP(L1::top / C1));
    auto const scaling_factor =
        ( //FIXME this is ugly
            (ALPHA1 >= 0 ?  F::make(static_cast<decltype(p)>( ALPHA1))
                         : -F::make(static_cast<decltype(p)>(-ALPHA1)))
            + F::gen *
            (ALPHA2 >= 0 ?  F::make(static_cast<decltype(p)>( ALPHA2))
                         : -F::make(static_cast<decltype(p)>(-ALPHA2)))
        ) * F::make(L12::topval).inv();
    auto [phi_E_L1_, phi_E_L1__hat] = isogeny_with_dual<EXP>(KL1_C1, L1::top/C1, scaling_factor);
    auto EL1_ = chain_codomain(phi_E_L1_);
    assert(*EL1_ == *EL1L2);
    auto P_EL1_ = chain_eval(phi_E_L1_, P_EC1);

    wQ2 = chain_eval(phi_E_L1__hat, wQ2);
    wP3 = chain_eval(phi_E_L1__hat, wP3);
    wQ3 = chain_eval(phi_E_L1__hat, wQ3);

    auto Pfull_EB1B2 = P_EB1B2 + chain_eval(phi_E_L1L2_hat, (L12::top/B12).value() * P_EL1_);
    auto Qfull_EC1 = Q_EC1 + chain_eval(phi_E_L1__hat, (L1::top/C1).value() * Q_EL1L2);

    auto KC1_C2 = (L2::top/C2).value() * wQ2;
    assert(order<EXP>(KC1_C2) == EXP(C2));

    auto KC1_A3 = EC1->point();
    {
        constexpr auto evs = eigenvalues<NORM_ALPHA, TRACE_ALPHA, L3::base>();
        auto const &smul = scalar_mul<typename L3::num_t, ecp<F>>;
        product_tree const tree(L3::top.prime_powers());
        auto const leaves_P3 = tree.template reduce<smul>(P3);
        auto const leaves_Q3 = tree.template reduce<smul>(Q3);
        auto const leaves_wP3 = tree.template reduce<smul>(wP3);
        auto const leaves_wQ3 = tree.template reduce<smul>(wQ3);
        for (size_t i = 0; i < L3::count; ++i) {
            assert(order<EXP>(leaves_P3[i]).value() == L3::base[i]);
            assert(order<EXP>(leaves_Q3[i]).value() == L3::base[i]);
            assert(order<EXP>(leaves_wP3[i]).value() == L3::base[i]);
            assert(order<EXP>(leaves_wQ3[i]).value() == L3::base[i]);
//            assert(!bsgs(leaves_P3[i], leaves_Q3[i], L3::base[i]));
//            assert(!bsgs(leaves_wP3[i], leaves_wQ3[i], L3::base[i]));
            assert(!(B3[i] && C3[i]));
            if (!B3[i] && !C3[i])
                continue;
            auto v = B3[i] ? evs[i].second : evs[i].first;
            auto ker = kernel_from_omega(leaves_P3[i], leaves_Q3[i], leaves_wP3[i], leaves_wQ3[i], v);
#ifndef NDEBUG
{
    auto u = B3[i] ? evs[i].first : evs[i].second;
    auto wker = ker;
    wker = chain_eval(phi_E_C1_hat, wker);
    wker = chain_eval(phi_E_B1B2, wker);
    wker = chain_eval(phi_E_L1L2, wker);
    wker = chain_eval(phi_E_L1__hat, wker);
    assert(*wker.curve() == *ker.curve());
    auto mer = kernel_from_omega(leaves_P3[i], leaves_Q3[i], -leaves_wP3[i], -leaves_wQ3[i], v);
    auto wmer = mer;
    wmer = chain_eval(phi_E_C1_hat, wmer);
    wmer = chain_eval(phi_E_B1B2, wmer);
    wmer = chain_eval(phi_E_L1L2, wmer);
    wmer = chain_eval(phi_E_L1__hat, wmer);
    assert(*wmer.curve() == *ker.curve());
    assert(wker == u*ker);
}
#endif
            KC1_A3 += ker;
        }
    }
    assert(order<EXP>(KC1_A3) == EXP(B3*C3));

    auto KC1_B1B2 = (L1::top/C1).value() * Q_EB1B2;
    assert(order<EXP>(KC1_B1B2) == EXP(C1));
    auto phi_EB1B2_C1 = factored_isogeny(KC1_B1B2, C1);
    auto Pfull_EA1C1 = chain_eval(phi_EB1B2_C1, Pfull_EB1B2);

    auto KB1B2_EC1 = (L12::top/B12).value() * P_EC1;
    assert(order<EXP>(KB1B2_EC1) == EXP(B12));
    auto phi_EC1_B1B2 = factored_isogeny(KB1B2_EC1, B12);
    assert(*chain_codomain(phi_EC1_B1B2) == *chain_codomain(phi_EB1B2_C1));
    auto Qfull_EA1C1 = chain_eval(phi_EC1_B1B2, Qfull_EC1);
    auto KEA1C1_C2A3 = chain_eval(phi_EC1_B1B2, KC1_C2 + KC1_A3);
    assert(order<EXP>(KEA1C1_C2A3) == EXP(C2)*EXP(B3*C3));

    auto phi_EA1C1_C2A3 = factored_isogeny(KEA1C1_C2A3, L23(C2)*L23(B3*C3));
    auto PEA = chain_eval(phi_EA1C1_C2A3, Pfull_EA1C1);
    auto QEA = chain_eval(phi_EA1C1_C2A3, Qfull_EA1C1);

    return {chain_codomain(phi_EA1C1_C2A3), PEA, QEA};
}

template <class L1, class L2, class L3>
constexpr size_t _count = L1::count + L2::count + L3::count;

template <class L1, class L2, class L3>
consteval std::array<size_t, _count<L1,L2,L3>> _mapping()
{
    std::array<size_t, _count<L1,L2,L3>> ret;
    for (size_t j1 = 0, j2 = 0, j3 = 0, k = 0; k < ret.size(); ++k) {
        uint64_t v1 = j1 < L1::count ? L1::base[j1] : static_cast<uint64_t>(-1);
        uint64_t v2 = j2 < L2::count ? L2::base[j2] : static_cast<uint64_t>(-1);
        uint64_t v3 = j3 < L3::count ? L3::base[j3] : static_cast<uint64_t>(-1);
        if (v1 < v2 && v1 < v3)
            ++j1, ret[k] = 1;
        else if (v2 < v1 && v2 < v3)
            ++j2, ret[k] = 2;
        else if (v3 < v1 && v3 < v2)
            ++j3, ret[k] = 3;
        else
            throw std::logic_error("not coprime?");
    }
    return ret;
}

template <class L1, class L2, class L3, typename T=int16_t>
using exponent_vec = std::array<T, L1::count + L2::count + L3::count>;

template <class L1, class L2, class L3>
constexpr std::tuple<L1,L1, L2,L2, L3,L3> _next_steps(exponent_vec<L1,L2,L3> &vec)
{
    constexpr auto _map = _mapping<L1,L2,L3>();
    typename L1::vec_t B1 {}, C1 {};
    typename L2::vec_t B2 {}, C2 {};
    typename L3::vec_t B3 {}, C3 {};
    for (size_t k = 0, j1 = 0, j2 = 0, j3 = 0; k < vec.size(); ++k) {
        assert(j1 <= L1::count && j2 <= L2::count && j3 <= L3::count);
        switch (_map[k]) {
        case 1:
            if      (vec[k] > 0) { ++B1[j1]; --vec[k]; }
            else if (vec[k] < 0) { ++C1[j1]; ++vec[k]; }
            ++j1;
            break;
        case 2:
            if      (vec[k] > 0) { ++B2[j2]; --vec[k]; }
            else if (vec[k] < 0) { ++C2[j2]; ++vec[k]; }
            ++j2;
            break;
        case 3:
            if      (vec[k] > 0) { ++B3[j3]; --vec[k]; }
            else if (vec[k] < 0) { ++C3[j3]; ++vec[k]; }
            ++j3;
            break;
        default:
            __builtin_unreachable();
        }
        if (k == vec.size() - 1)
            assert(j1 == L1::count && j2 == L2::count && j3 == L3::count);
    }
    return {B1,C1, B2,C2, B3,C3};
}

template <class L, class F>
ecp<F> canonicalize_point(ecp<F> const &R, ecp<F> const &S, F const &e, ecp<F> P)
{
    using num_t = signify_t<typename L::num_t>;
    constexpr num_t mod = L::topval;

    num_t u = dlp<L, num_t>(e, weil_pairing(P, S, L::topval), L::top).value();
    num_t v = dlp<L, num_t>(e, weil_pairing(R, P, L::topval), L::top).value();
    assert(u*R + v*S == P);
    if (u > (L::topval >> 1)) {
        u = L::topval - u;
        v = L::topval - v;
    }

    num_t mu = boost::integer::gcd(u, mod);
    num_t mv = mod / mu;
    num_t f;
    if (mu == 1)
        f = boost::integer::mod_inverse(u, mod);
    else if (mu == mod)
        f = boost::integer::mod_inverse(v, mod);
    else {
        auto fu = boost::integer::mod_inverse(v   , mu);
        auto fv = boost::integer::mod_inverse(u/mu, mv);
        f = crt<num_t>({fu, mu}, {fv, mv}).first;
    }
    assert(boost::integer::gcd(f, mod) == 1);
    assert(ext(u) * f % mod == ext(mu));
    assert(mu == 1 || ext(v) * f % mu == 1);
    P *= f;

    if (-P.affine_y() < P.affine_y())
        P = -P;
    return P;
}

template <class EXP, class L1, class L2, class L3, class F>
oriented_curve<F> canonicalize(oriented_curve<F> EPQ)
{
    auto const j = EPQ.E->j_invariant();

    std::shared_ptr<ec<F>> EE;
    if (!j)
        EE = std::make_shared<ec<F>>(F::one, F::zero);
    else if (j == F(1728))
        EE = std::make_shared<ec<F>>(F::zero, F::one);
    else
        EE = ec<F>::from_j(j);
    assert(EE->j_invariant() == j);

    {
        auto const &a = EPQ.E->a(), &b = EPQ.E->b();
        auto const &A = EE->a(), &B = EE->b();
        auto u2 = A * b * (a * B).inv();
        auto uu = u2.sqrt();
        if (!uu) {
            auto nonsq = F::gen;
            while (nonsq.is_square())
                nonsq += F::one;
            auto nonsq2 = nonsq * nonsq;
            EE = std::make_shared<ec<F>>(nonsq2*A, nonsq2*nonsq*B);
            assert(EE->j_invariant() == j);
            u2 *= nonsq.inv();
            uu = u2.sqrt();
        }
        auto const &u = uu.value();
        assert(a * u.pow(4).inv() == EE->a());
        assert(b * u.pow(6).inv() == EE->b());
        isog iso(EPQ.E->point());   // identity morphism
        iso._scale(u, EE);
        EPQ.E = EE;
        EPQ.P = iso(EPQ.P);
        EPQ.Q = iso(EPQ.Q);
    }

    constexpr EXP LL1 = L1::top, LL2 = L2::top, LL12 = LL1*LL2;
    constexpr static auto LL12val = LL12.value();
    F e;
    auto R = EE->point(), S = EE->point();
    {
        //TODO could be point_of_order(canonical=true) or something
        F x = F::gen - F::one;
        while (true) {
            x += F::one;
            auto pt = EE->lift_x(x);
            if (!pt)
                continue;
            auto ord = order<EXP>(*pt);
            if (!(ord >= LL12))
                continue;
            R = (ord / LL12).value() * *pt;
            assert(order<EXP>(R) == LL12);
            break;
        }
        x = F::gen+F::gen - F::one;
        while (true) {
            x += F::one;
            auto pt = EE->lift_x(x);
            if (!pt)
                continue;
            auto ord = order<EXP>(*pt);
            if (!(ord >= LL12))
                continue;
            S = (ord / LL12).value() * *pt;
            assert(order<EXP>(S) == LL12);
            e = weil_pairing(R, S, LL12val);
            if (multiplicative_order(e, EXP::top) == LL12val)
                break;
        }
    }

    using L12 = _fact<LL12val>;
    EPQ.P = canonicalize_point<L12>(R, S, e, EPQ.P);

    R *= L2::topval;
    S *= L2::topval;
    e = e.pow(L2::topval);
    assert(order<EXP>(R) == LL1);
    assert(order<EXP>(S) == LL1);
    assert(multiplicative_order(e, EXP::top) == LL1);
    EPQ.Q = canonicalize_point<L1>(R, S, e, EPQ.Q);

    assert(order<EXP>(EPQ.P) == LL12);
    assert(order<EXP>(EPQ.Q) == LL1);

    return EPQ;
}

template <class EXP, class L1, class L2, class L3, class F>
oriented_curve<F> scallop_vec(oriented_curve<F> EPQ,
                             exponent_vec<L1,L2,L3> vec)
{
    while (std::any_of(vec.begin(), vec.end(), [](int64_t e) { return e; })) {
        auto const tup = _next_steps<L1,L2,L3>(vec);
        EPQ = std::apply(scallop_steps<EXP,L1,L2,L3,F>, std::tuple_cat(std::make_tuple(EPQ), tup));
    }
    return canonicalize<EXP,L1,L2,L3>(EPQ);
}

#include "babai.hpp"

template <class EXP, class L1, class L2, class L3, auto const &relations, class I, class F>
oriented_curve<F> scallop(oriented_curve<F> const &EPQ, I const &x)
{
    auto y = static_cast<signify_t<I>>(x);
    auto const GS = gram_schmidt(relations);
    exponent_vec<L1,L2,L3> v = babai(y, relations, GS);
    return scallop_vec<EXP,L1,L2,L3>(EPQ, v);
}

