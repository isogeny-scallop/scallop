
#pragma once

#include "bigint.hpp"

#define _STR(x) #x
#define STR(x) _STR(x)
#define STEP(I) \
        "mulxq " STR(8*(I)) "(%4), %0, %1 \n" \
        "mov " STR(8*(I)) "(%3), %2 \n" \
        "adcxq %0, %2 \n" \
        "mov %2, " STR(8*(I)) "(%3) \n" \
        "mov " STR(8*(I)+8) "(%3), %2 \n" \
        "adoxq %1, %2 \n" \
        "mov %2, " STR(8*(I)+8) "(%3) \n"

#include <boost/preprocessor/repetition/for.hpp>
#define     PRED(r, state)      state
#define     OP(r, state)        BOOST_PP_DEC(state)
#define     MACRO(r, state)     STEP(r-1)
#define     GENERATE(n)         BOOST_PP_FOR(n, PRED, OP, MACRO)

#define _DEFINE_MULSTEP(N) \
    template <> \
    inline void _mulstep<N>(uint64_t *acc, uint64_t const *ys, uint64_t x) \
    { \
        uint64_t lo, hi, tmp; \
        \
        asm volatile( \
            "xor %0, %0 \n" \
            GENERATE(N) \
            "adcq $0, " STR(8*N) "(%3) \n" \
            : "=&r"(lo), "=&r"(hi), "=&r"(tmp) \
            : "r"(acc), "r"(ys), "d"(x) \
            : "cc", "memory" \
        ); \
    }

template <size_t N> static inline void _mulstep(uint64_t *acc, uint64_t const *ys, uint64_t x);
_DEFINE_MULSTEP(1)
_DEFINE_MULSTEP(2)
_DEFINE_MULSTEP(3)
_DEFINE_MULSTEP(4)
_DEFINE_MULSTEP(5)
_DEFINE_MULSTEP(6)
_DEFINE_MULSTEP(7)
_DEFINE_MULSTEP(8)
_DEFINE_MULSTEP(9)
_DEFINE_MULSTEP(10)
_DEFINE_MULSTEP(11)
_DEFINE_MULSTEP(12)
_DEFINE_MULSTEP(13)
_DEFINE_MULSTEP(14)
_DEFINE_MULSTEP(15)
_DEFINE_MULSTEP(16)

#define LIMBS ((std::numeric_limits<decltype(p)>::digits + 63) / 64)


template <size_t N>
constexpr std::array<uint64_t,LIMBS> _chop(uint_t<N> const &v)
{
    std::array<uint64_t,LIMBS> r;
    for (size_t i = 0; i < LIMBS; ++i)
        r[i] = static_cast<uint64_t>(v >> 64*i);
    return r;
};

consteval std::array<uint64_t,LIMBS> _1()
{
    auto v = powm(static_cast<decltype(p)>(2), 64*LIMBS, p);
    return _chop(v);
}

consteval std::array<uint64_t,LIMBS> _R()
{
    auto v = powm(static_cast<decltype(p)>(2), 2*64*LIMBS, p);
    return _chop(v);
};

constexpr std::array<uint64_t,LIMBS> the_p = _chop(p);
constexpr std::array<uint64_t,LIMBS> the_1 = _1();
constexpr std::array<uint64_t,LIMBS> the_R = _R();
constexpr std::array<uint64_t,LIMBS> the_Rinv = {1};

#include <boost/integer/mod_inverse.hpp>

consteval uint64_t _inv_min_p()
{
    uint64_t x = -static_cast<uint64_t>(p);
    uint64_t y = 1;
    for (size_t k = 1; k < 64; k <<= 1)
        y *= 2 - x * y;
    if (x * y != 1)
        throw;
    return y;
}

#include <boost/random.hpp>

constexpr uint64_t the_inv_min_p = _inv_min_p();

struct no_montgomerize { };

template <typename __>
class _fp
{
    using num = std::array<uint64_t, LIMBS>;

    public:
        constexpr _fp() : limbs{0} { }
        _fp(uint64_t v)
            : limbs{v}
            { *this *= the_R; }
        _fp(std::string const &s)
            : limbs{_chop(decltype(p)(s))}
            { *this *= the_R; }
        _fp(decltype(p) &v)
            : limbs{_chop(v)}
            { *this *= the_R; }
        constexpr _fp(no_montgomerize _, std::array<uint64_t,LIMBS> const &v)
            : limbs{v}
            { (void) _; }

        consteval static decltype(p) &cardinality() { return p; }

        constexpr static bool _addc(num &z, num const &x, num const &y)
        {
            bool carry = false;
            for (size_t i = 0; i < LIMBS; ++i) {
                uint64_t t = x[i] + y[i];
                bool wrap = t < x[i];
                z[i] = t + carry;
                carry = wrap || z[i] < t;
            }
            return carry;
        }

        constexpr static bool _subc(num &z, num const &x, num const &y)
        {
            bool carry = false;
            for (size_t i = 0; i < LIMBS; ++i) {
                uint64_t t = x[i] - y[i];
                bool wrap = t > x[i];
                z[i] = t - carry;
                carry = wrap || z[i] > t;
            }
            return carry;
        }

        constexpr void _reduce_once()
        {
            num r;
            if (!_fp::_subc(r, this->limbs, the_p))
                this->limbs = r;
        }

        constexpr _fp &operator+=(_fp const &other)
        {
            _fp::_addc(this->limbs, this->limbs, other.limbs);
            this->_reduce_once();
            return *this;
        }
        constexpr _fp operator+(_fp const &other) const { _fp r = *this; return r += other; }

        constexpr _fp &operator-=(_fp const &other)
        {
            if (_fp::_subc(this->limbs, this->limbs, other.limbs))
                _fp::_addc(this->limbs, this->limbs, the_p);
            return *this;
        }
        constexpr _fp operator-(_fp const &other) const { _fp r = *this; return r -= other; }
        constexpr _fp operator-() const { _fp r; return r -= *this; }

        _fp operator*(num const &other) const
        {
            uint64_t rs[2*LIMBS] {0};

            for (size_t i = 0; i < LIMBS; ++i) {

                // add this[i] * other
                _mulstep<LIMBS>(rs + i, other.data(), this->limbs[i]);

                uint64_t f = rs[i] * the_inv_min_p;

                // add f*p
                _mulstep<LIMBS>(rs + i, the_p.data(), f);
            }

            _fp ret;
            for (size_t i = 0; i < LIMBS; ++i)
                ret.limbs[i] = rs[LIMBS+i];
            ret._reduce_once();
            return ret;
        }
        _fp &operator*=(num const &other) { return *this = *this * other; }
        _fp operator*(_fp const &other) const { return *this * other.limbs; }
        _fp &operator*=(_fp const &other) { return *this *= other.limbs; }

        _fp inv() const
        {
            static_assert(the_p[0] >= 2);
            assert(*this);
            _fp t = *this, r(1);
            for (size_t i = 0; i < the_p.size(); ++i) {
                uint64_t limb = the_p[i] - (!i ? 2 : 0);
                for (size_t j = 0; j < 64; ++j) {
                    if ((limb >> j) & 1)
                        r *= t;
                    t *= t;
                }
            }
            assert((*this*r).limbs == the_1);
            return r;
        }

        constexpr operator bool() const { return std::any_of(this->limbs.begin(), this->limbs.end(), [](uint64_t v) { return v; }); }
        constexpr bool operator==(_fp const &other) const { return std::equal(this->limbs.begin(), this->limbs.end(), other.limbs.begin()); }
        constexpr bool operator!=(_fp const &other) const { return !(*this == other); }

        constexpr bool operator<(_fp const &other) const { num dummy; return _fp::_subc(dummy, this->limbs, other.limbs); }  //XXX

        operator std::remove_const<decltype(p)>::type() const
        {
            _fp el = *this * the_Rinv;
            std::remove_const<decltype(p)>::type r = 0;
            for (auto it = el.limbs.rbegin(); it != el.limbs.rend(); ++it) {
                r <<= 64;
                r += *it;
            }
            return r;
        }

        friend std::ostream &operator<<(std::ostream &o, _fp const &el)
        {
            return o << static_cast<decltype(p)>(el);
        }

        static _fp random()
        {
            static boost::random::minstd_rand rng;
            static const boost::random::uniform_int_distribution<decltype(p)> dist(0, p-1);
            return {dist(rng)};
        }

        static constexpr _fp zero {no_montgomerize(), {0}};
        static constexpr _fp one {no_montgomerize(), the_1};

        friend size_t std::hash<_fp>::operator()(_fp const &) const;

    private:
        num limbs;
};

namespace std {
    template <typename __>
    struct hash<_fp<__>>
    {
        size_t operator()(_fp<__> const &el) const {
            char const *ptr = reinterpret_cast<char const *>(el.limbs.data());
            size_t const sz = el.limbs.size() * sizeof(uint64_t);
            std::string_view const s(ptr, sz);
            return std::hash<std::string_view>()(s);
        }
    };
}

