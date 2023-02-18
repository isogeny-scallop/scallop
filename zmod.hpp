
#pragma once

#include <boost/random.hpp>

#include "bigint.hpp"

template <size_t N, uint_t<N> const &M>
class zmod
{
    public:
        constexpr zmod() : rep{0} { }
        constexpr zmod(std::string const &s) : rep{s} { }
        constexpr zmod(uint_t<N> const &v) : rep{v} { }

        consteval static uint_t<N> const &cardinality() { return M; };

        zmod &operator+=(zmod const &other) { this->rep += other.rep; this->rep %= M; return *this; }
        zmod operator+(zmod const &other) const { zmod r = *this; return r += other; }

        zmod &operator-=(zmod const &other)
        {
            if ((this->rep -= other.rep) >> (N-1))
                this->rep += M;
            return *this;
        }
        zmod operator-(zmod const &other) const { zmod r = *this; return r -= other; }
        zmod operator-() const { zmod r; return r -= *this; }

        zmod &operator*=(zmod const &other)
        {
            uint_t<2*N> r = this->rep;
            r *= other.rep;
            r %= M;
            this->rep = static_cast<uint_t<N>>(r);
            return *this;
        }
        zmod operator*(zmod const &other) const { zmod r = *this; return r *= other; }

        template <size_t L>
        zmod pow(uint_t<L> const &e) const { return powm(this->rep, e, M); }
        zmod inv() const { return this->pow(M-2); }  //XXX assumes prime

        operator bool() const { return !!this->rep; }
        bool operator==(zmod const &other) const { return this->rep == other.rep; }
        bool operator!=(zmod const &other) const { return this->rep != other.rep; }

        bool operator<(zmod const &other) const { return this->rep < other.rep; }  //XXX

        friend std::ostream &operator<<(std::ostream &o, zmod const &el) { return o << el.rep; }

        static zmod random()
        {
            static boost::random::minstd_rand rng;
            static const boost::random::uniform_int_distribution<uint_t<N>> dist(0, M-1);
            return {dist(rng)};
        };

        static constexpr zmod zero {};
        static constexpr zmod one {1};

    private:
        uint_t<N> rep;
};

