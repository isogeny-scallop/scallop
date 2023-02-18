
#pragma once

#include "ecp.hpp"
#include "fact.hpp"

template <class I>
class product_tree
{
    private:
        std::vector<std::vector<I>> layers;

    public:
        constexpr product_tree(std::vector<I> leaves)
        {
            if (leaves.empty())
                throw std::logic_error("no leaves");
            this->layers.push_back(leaves);
            while (this->layers.back().size() > 1) {
                std::vector<I> next;
                {
                    auto const &prev = this->layers.back();
                    for (size_t i = 0; i < prev.size()-1; i += 2)
                        next.push_back(prev[i] * prev[i+1]);
                    if (prev.size() % 2)
                        next.push_back(prev.back());
                }
                this->layers.push_back(next);
            }
        }

        template <size_t SZ>
        constexpr product_tree(std::array<I,SZ> const &leaves)
            : product_tree(std::vector(leaves.begin(), leaves.end())) { }

        constexpr std::vector<I> const &operator[](size_t i) const { return this->layers.at(i); }

        template <auto const &M, class G>
        std::vector<G> reduce(G const &el) const
        {
            std::vector<G> r {el};
            for (auto it = this->layers.rbegin()+1; it != this->layers.rend(); ++it) {
                decltype(r) next;
                assert(it[-1].size() == r.size());
                for (size_t i = 0; i < it->size()-1; i += 2) {
                    next.push_back(M((*it)[i+1], r[i/2]));
                    next.push_back(M((*it)[i+0], r[i/2]));
                }
                if (it->size() % 2)
                    next.push_back(r[it->size()/2]);
                r = std::move(next);
                assert(it->size() == r.size());
            }
            return r;
        }

        //TODO could precompute CRT basis stuff
        I crt(std::vector<size_t> const &rems) const
        {
            if (rems.size() != this->layers[0].size())
                throw std::logic_error("wrong size");

            std::vector<std::pair<I,I>> cur;
            cur.reserve(rems.size());
            for (size_t i = 0; i < rems.size(); ++i)
                cur.push_back({rems[i], this->layers[0][i]});

            while (cur.size() > 1) {
                decltype(cur) next;
                for (size_t i = 0; i < cur.size()-1; i += 2)
                    next.push_back(::crt(cur[i+0], cur[i+1]));
                if (cur.size() % 2)
                    next.push_back(cur.back());
                cur = std::move(next);
            }
            auto const &res = cur[0].first;

            return res;
        }
};


template <class I, class G>
std::function const scalar_mul = [](I const &k, G const &P) { return k * P; };

template <class EXP, class F>
EXP order(ecp<F> const &P)
{
    assert(!(EXP::top.value() * P));
    auto const &ps = EXP::base;
    constexpr auto qs = EXP::top.prime_powers();
    assert(std::equal(ps.begin(), ps.end(), qs.begin(), qs.end(),
                      [](uint64_t p, EXP::num_t const &q) { return q % p == 0; }));

    static product_tree const tree(qs);
    auto Qs = tree.template reduce<scalar_mul<typename EXP::num_t, ecp<F>>>(P);
    typename EXP::vec_t ord {};
    for (size_t i = 0; i < ps.size(); ++i) {
        while (Qs[i]) {
            ++ord[i];
            assert(ord[i] <= EXP::top[i]);
            Qs[i] *= ps[i];
        }
        assert(!Qs[i]);
    }
    return ord;
}

template <class I, class F>
std::function const exp_fun = [](I const &k, F const &x) { return x.pow(k); };

template <class EXP, class F>
EXP multiplicative_order(F const &x, EXP const &top)
{
    assert(x.pow(top.value()) == F::one);
    auto const &ps = EXP::base;
    auto const qs = top.prime_powers();
    assert(std::equal(ps.begin(), ps.end(), qs.begin(), qs.end(),
                      [](uint64_t p, EXP::num_t const &q) { return q == 1 || q % p == 0; }));

    static product_tree const tree(qs);
    auto ys = tree.template reduce<exp_fun<typename EXP::num_t, F>>(x);
    typename EXP::vec_t ord {};
    for (size_t i = 0; i < ps.size(); ++i) {
        while (ys[i] != F::one) {
            ++ord[i];
            assert(ord[i] <= top[i]);
            ys[i] = ys[i].pow(ps[i]);
        }
        assert(ys[i] == F::one);
    }
    return ord;
}

template <class EXP, class F>
bool orders_present(F const &x, EXP const &top)
{
    auto const &ps = EXP::base;
    auto const qs = top.prime_powers();
    assert(std::equal(ps.begin(), ps.end(), qs.begin(), qs.end(),
                      [](uint64_t p, EXP::num_t const &q) { return q == 1 || q % p == 0; }));

    static product_tree const tree(qs);
    auto ys = tree.template reduce<exp_fun<typename EXP::num_t, F>>(x);
    typename EXP::vec_t ord {};
    for (size_t i = 0; i < ps.size(); ++i) {
        if (!ord[i])
            continue;
        if (ys[i] == F::one)
            return false;
    }
    return true;
}

template <class EXP, class F>
ecp<F> point_of_order(ec<F> const &E, EXP const &n)
{
    while (true) {
        auto P = E.random_point();
        auto ord = order<EXP>(P);
        if (ord >= n) {
            typename EXP::num_t cof = ord / n;
            auto Q = cof * P;
            assert(order<EXP>(Q) == n);
            return Q;
        }
    }
}

template <class F>
std::optional<size_t> bsgs(F const &x, F const &y, size_t n)
{
    size_t m = std::sqrt(n-1) + 1;

    std::unordered_map<F, size_t> tab;
    tab.reserve(m);

    auto t = F::one;
    for (size_t i = 0; i < m; ++i) {
        tab[t] = i;
        t *= x;
    }

    auto const invmx = t.inv();
    t = y;
    for (size_t j = 0; j < m; ++j) {
        auto it = tab.find(t);
        if (it != tab.end())
            return it->second + m*j;
        t *= invmx;
    }

    return {};
}

template <class I, class F>
std::optional<I> dlp_primepower(F const &x, F const &y, I const &p, size_t k)
{
    F const x0 = x.pow(boost::multiprecision::pow(p, k-1));
    assert(x0 != F::one);
    assert(x0.pow(p) == F::one);
    assert(y.pow(boost::multiprecision::pow(p, k)) == F::one);

    I res = 0;

    for (size_t i = 0; i < k; ++i) {
        auto z = (x.pow(res).inv() * y).pow(boost::multiprecision::pow(p, k-1-i));
        auto s = bsgs(x0, z, static_cast<uint64_t>(p));
        if (!s)
            return {};
        res += boost::multiprecision::pow(p, i) * *s;
    }
    assert(x.pow(res) == y);
    return res;
}

template <class EXP, class I, class F>
std::optional<I> dlp(F const &x, F const &y, EXP const &ord)
{
    assert(multiplicative_order(x, EXP::top) == ord);
    assert(y.pow(ord.value()) == F::one);
    auto const &ps = EXP::base;
    constexpr auto qs = EXP::top.prime_powers();
    assert(std::equal(ps.begin(), ps.end(), qs.begin(), qs.end(),
                      [](uint64_t p, EXP::num_t const &q) { return q == 1 || q % p == 0; }));

    using num_t = signify_t<typename EXP::num_t>;

    std::array<num_t, EXP::count> myqs;
    std::copy(qs.begin(), qs.end(), myqs.begin());

    static product_tree const tree(myqs);
    auto xs = tree.template reduce<exp_fun<num_t, F>>(x);
    auto ys = tree.template reduce<exp_fun<num_t, F>>(y);
    std::vector<size_t> res;
    for (size_t i = 0; i < EXP::count; ++i) {
        auto r = dlp_primepower<num_t>(xs[i], ys[i], ps[i], ord[i]);
        if (!r)
            return {};
        res.push_back(static_cast<size_t>(*r));
    }
    auto r = tree.crt(res);
    assert(x.pow(r) == y);
    return static_cast<I>(r);
}

template <class EXP, class F>
ecp<F> complete_basis(ecp<F> const &P, EXP const &ord)
{
    assert(order<EXP>(P) == ord);
    auto const &E = P.curve();
    while (true) {
        auto Q = point_of_order(*E, ord);  //FIXME should do pairing test on large order first?
        auto const e = weil_pairing(P, Q, ord.value());
        assert(e.pow(ord.value()) == F::one);
        if (multiplicative_order(e, EXP::top) == ord)   //FIXME can EXP::top -> ord?
            return Q;
    }
}

template <class EXP, class F>
std::pair<ecp<F>,ecp<F>> basis(ec<F> const &E, EXP const &ord)
{
    auto const P = point_of_order<EXP>(E, ord);
    auto const Q = complete_basis<EXP>(P, ord);
    return {P, Q};
}

