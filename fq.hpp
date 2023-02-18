
#pragma once

#include <optional>

#include "bigint.hpp"

// R[x] / (x^2 + 1)

template <class F>
class _fq
{
    public:
        constexpr _fq() : re{F::zero}, im{F::zero} { }
        constexpr _fq(F const &re, F const &im = F::zero) : re{re}, im{im} { }
        constexpr _fq(unsigned long v) : re{F(v)}, im{F::zero} { }

        template <size_t N>
        constexpr _fq(uint_t<N> const &v) : re{F(v)}, im{F::zero} { }

        using cardinality_t = dbl_t<typename std::remove_cvref<decltype(F::cardinality())>::type>;
        consteval static cardinality_t cardinality()
        {
            cardinality_t base_q = F::cardinality();
            return base_q * base_q;
        }

        constexpr _fq &operator+=(_fq const &other) {
            this->re += other.re;
            this->im += other.im;
            return *this;
        }
        constexpr _fq operator+(_fq const &other) const { _fq r = *this; return r += other; }

        constexpr _fq &operator-=(_fq const &other) {
            this->re -= other.re;
            this->im -= other.im;
            return *this;
        }
        constexpr _fq operator-(_fq const &other) const { _fq r = *this; return r -= other; }
        constexpr _fq operator-() const { return _fq(-this->re, -this->im); }

        _fq operator*(_fq const &other) const
        {
            auto const &a = this->re, &b = this->im, &c = other.re, &d = other.im;
            auto ad = a*d, bc = b*c, abcd = (a-b)*(c+d);
            auto re = abcd - ad + bc;
            auto im = ad + bc;
            return {re, im};
        }
        _fq &operator*=(_fq const &other) { return *this = *this * other; }

        template <class I>
        _fq pow(I e) const
        {
            _fq r = _fq::one;
            _fq x = *this;
            while (e) {
                if (e & 1)
                    r *= x;
                x *= x;
                e >>= 1;
            }
            return r;
        }

        _fq inv() const
        {
            assert(*this);
            auto d = (this->re*this->re + this->im*this->im).inv();
            _fq r(this->re * d, -this->im * d);
//            assert(r * (*this) == _fq::one);
            return r;
        }

        constexpr bool is_square() const
        {
            constexpr auto q = _fq::cardinality();
            static_assert((q & 1) == 1);
            return this->pow(q/2) == F::one;
        }

        constexpr _fq conjugate() const
        {
            return {this->re, -this->im};
        }

        std::optional<_fq> sqrt() const
        {
            if (!*this) return {F::zero};
//            if (!this->is_square()) return {};

            // https://eprint.iacr.org/2012/685.pdf Algorithm 9
            if (F::cardinality() % 4 == 3) {
                auto const &a = *this, &i = _fq::gen;
                constexpr auto q = F::cardinality();
                auto a1 = a.pow(q >> 2);
                auto alpha = a1 * a1 * a;
                auto a0 = alpha.conjugate() * alpha;
                if (a0 == -_fq::one) {
                    assert(!this->is_square());
                    return {};
                }
                auto x0 = a1 * a;
                _fq x;
                if (alpha == -_fq::one) {
                    x = i * x0;
                }
                else {
                    auto b = (_fq::one + alpha).pow(q >> 1);
                    x = b * x0;
                }
                assert(x * x == a);
                return x;
            }

            /* multiplication of two linear polynomials modulo x^2 - *this */
            auto mm = [this](std::array<_fq,2> const& f, std::array<_fq,2> const& g)
                {
                    std::array<_fq,2> ret = {
                        f[0]*g[0] + *this*f[1]*g[1],
                        f[0]*g[1] + f[1]*g[0],
                    };
                    return ret;
                }
            ;

            for (size_t tries = 0; tries < 64; ++tries) {
                std::array<_fq,2> a {_fq::random(), _fq::one};
                std::array<_fq,2> b = {_fq::one, _fq::zero};

                constexpr auto p_ = F::cardinality();   //FIXME could be reference
                static_assert(p_);
                constexpr size_t l = std::numeric_limits<decltype(p_)>::digits;
                uint_t<2*l> const p = p_;
                auto k = (p*p - 1) / 2;
                /* a^k % (x^2-this) */
                while (k) {
                    if (k & 1)
                        b = mm(b, a);
                    a = mm(a, a);
                    k >>= 1;
                }

                b[0] -= F::one;

                if (b[1]) {
                    _fq r = b[0] * b[1].inv();
                    if (r*r == *this)
                        return r;
                }
            }

            /* probably simply not a square */
            return {};
        }

        constexpr operator bool() const { return this->re || this->im; }
        constexpr bool operator==(_fq const &other) const { return this->re == other.re && this->im == other.im; }
        constexpr bool operator!=(_fq const &other) const { return this->re != other.re || this->im != other.im; }

        constexpr bool operator<(_fq const &other) const { return this->re < other.re || (this->re == other.re && this->im < other.im); }  //XXX

        friend std::ostream &operator<<(std::ostream &o, _fq const &el)
        {
            std::ostringstream ss;
//XXX            ss << el.re << "+" << el.im << "*i";
ss << el.im << "*i+" << el.re;
            return o << ss.str();
        }

        friend std::istream &operator>>(std::istream &i, _fq &el)
        {
            std::string re, im;
            i >> re >> im;
            el.re = F(re);
            el.im = F(im);
            return i;
        }

        F const &real() const { return this->re; }
        F const &imag() const { return this->im; }

        static _fq random() { return {F::random(), F::random()}; }

        static constexpr _fq zero {F::zero, F::zero};
        static constexpr _fq one {F::one, F::zero};
        static constexpr _fq gen {F::zero, F::one};

        template <class U>
        static _fq make(U const &u) { return {F(u), F::zero}; }

        template <class U, class V>
        static _fq make(U const &u, V const &v) { return {F(u), F(v)}; }

    private:
        F re, im;
};

namespace std {
    template <class F>
    struct hash<_fq<F>>
    {
        size_t operator()(_fq<F> const &el) const {
            return (0xff + 1111111*std::hash<F>()(el.real())) ^ el.imag();
        }
    };
}


#include "fp.hpp"

using fq = _fq<fp>;

