
#pragma once

#include "bigint.hpp"

template <auto const &M>
consteval size_t prime_divisor_count()
{
    if (M <= 0)
        throw std::domain_error("not positive");
    auto m = M;
    size_t r = 0;
    for (uint64_t i = 2; m != 1; ++i) {
        if (m % i == 0) {
            ++r;
            while (m % i == 0)
                m /= i;
        }
        if (i > (uint64_t) 1 << 42)
            throw std::domain_error("not smooth enough");
    }
    return r;
}

template <auto const &M>
consteval std::array<uint64_t, prime_divisor_count<M>()> prime_divisors()
{
    if (M <= 0)
        throw std::domain_error("not positive");
    auto m = M;
    std::array<uint64_t, prime_divisor_count<M>()> r;
    size_t idx = 0;
    for (uint64_t i = 2; m != 1; ++i) {
        if (m % i == 0) {
            r[idx++] = i;
            while (m % i == 0)
                m /= i;
        }
        if (i > (uint64_t) 1 << 42)
            throw std::domain_error("not smooth enough");
    }
    return r;
}

template <auto const &M, auto const &B, typename T=uint8_t>
consteval std::array<T, B.size()> factor()
{
    if (M <= 0)
        throw std::domain_error("not positive");
    auto m = M;
    std::array<T, B.size()> e {};
    for (size_t i = 0; i < B.size(); ++i) {
        while (m % B[i] == 0) {
            m /= B[i];
            ++e[i];
        }
    }
    if (m != 1)
        throw std::domain_error("does not factor over base");
    return e;
}


template <auto const &M, typename T=uint8_t>
class _fact
{
    public:
        using num_t = std::decay<decltype(M)>::type;

        constexpr static size_t count = prime_divisor_count<M>();
        using vec_t = std::array<T, count>;

    private:
        vec_t e;

    public:
        constexpr _fact() : e{} { }
        constexpr _fact(vec_t ee) : e{ee} { }

        template <auto const &M2>
        constexpr _fact(_fact<M2> const &other)
        {
            if (other.count > this->count)
                throw std::logic_error("impossible conversion");
            size_t j = 0;
            for (size_t i = 0; i < this->count; ++i) {
                if (j < other.count && this->base[i] == other.base[j])
                    this->e[i] = other[j++];
                else
                    this->e[i] = 0;
            }
            if (j != other.count)
                throw std::logic_error("impossible conversion");
        }

        constexpr static std::array<uint64_t, prime_divisor_count<M>()> base = prime_divisors<M>();

        constexpr static num_t topval = M;
        constexpr static _fact top = _fact(factor<M, base>());

        template <auto const &V>
        consteval static _fact vec()
        {
            return _fact(factor<V, base>());
        }

        constexpr size_t size() const { return count; }
        constexpr T const &operator[](size_t i) const { return this->e[i]; }

        //FIXME this should probably omit ones and return a std::vector
        constexpr std::array<num_t, count> prime_powers() const
        {
            std::array<num_t, count> ret;
            ret.fill(1);
            for (size_t i = 0; i < this->e.size(); ++i) {
                num_t x = base[i];
                for (auto c = this->e[i]; c; ) {
                    if (c & 1)
                        ret[i] *= x;
                    if (c >>= 1)
                        x *= x;
                }
            }
            return ret;
        }

        constexpr num_t value() const
        {
            if (!count)
                return 1;
            auto qs = this->prime_powers();
            std::vector t(std::make_move_iterator(qs.begin()), std::make_move_iterator(qs.end()));
            while (t.size() > 1) {
                for (size_t i = 0; i < t.size()-1; i += 2)
                    t[i] *= t[i+1];
                for (auto it = t.begin(); it != t.end() && ++it != t.end(); it = t.erase(it))
                    ;
            }
            return t[0];
        }
        constexpr operator num_t() const { return this->value(); }

        constexpr _fact &operator*=(_fact const &other)
        {
            std::transform(this->e.begin(), this->e.end(), other.e.begin(), this->e.begin(), std::plus<>());
            return *this;
        }

        constexpr _fact operator/=(_fact const &other)
        {
            if (!(*this >= other))
                throw std::logic_error("not divisible");
            std::transform(this->e.begin(), this->e.end(), other.e.begin(), this->e.begin(), std::minus<>());
            return *this;
        }

        constexpr _fact operator*(_fact const &other) const { auto f = *this; return f *= other; }
        constexpr _fact operator/(_fact const &other) const { auto f = *this; return f /= other; }

        constexpr bool operator==(_fact const &other) const { return this->e == other.e; }
        constexpr bool operator!=(_fact const &other) const { return !(*this == other); }
        constexpr bool operator>=(_fact const &other) const { return std::equal(this->e.begin(), this->e.end(), other.e.begin(), std::greater_equal<>()); }
        constexpr bool operator<=(_fact const &other) const { return other >= *this; }

        constexpr _fact inc(size_t i) const
        {
            vec_t ee = this->e;
            ++ee[i];
            return {ee};
        }

        constexpr _fact dec(size_t i) const
        {
            vec_t ee = this->e;
            --ee[i];
            return {ee};
        }

        friend std::ostream& operator<<(std::ostream& o, _fact const &f)
        {
            bool first = true;
            for (size_t i = 0; i < f.e.size(); ++i) {
                if (!f.e[i])
                    continue;
                o << (first ? "{" : ",");
                first = false;
                o << f.base[i] << ":" << (uint64_t) f.e[i];
            }
            return (first ? o << "{" : o) << "}";
        }
};

