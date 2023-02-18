
#pragma once

#include <optional>
#include <memory>

#include <boost/random.hpp>

template <class F> class ecp;

template <class F>
class ec : public std::enable_shared_from_this<ec<F>>
{
    private:
        F _a, _b;

    public:
        using point_t = ecp<F>;

        ec(F const &a, F const &b) : _a{a}, _b{b}
        {
            if (!(4*this->_a*this->_a*this->_a + 27*this->_b*this->_b))
                throw std::logic_error("curve is singular");
        }

        F const &a() const { return this->_a; }
        F const &b() const { return this->_b; }

        point_t point() const { return point_t(this->shared_from_this()); }
        point_t point(F const &x, F const &y) { return point_t(this->shared_from_this(), x, y); }
        point_t point(F const &x, F const &y, F const &z) { return point_t(this->shared_from_this(), x, y, z); }

        std::optional<point_t> lift_x(F const &x) const
        {
            F rhs = (x*x + this->_a)*x + this->_b;
            std::optional<F> y = rhs.sqrt();
            if (y)
                return point_t(this->shared_from_this(), x, *y);
            return {};
        }

        // NB: not uniform
        point_t random_point() const
        {
            while (true) {
                auto pt = this->lift_x(F::random());
                if (pt) {
                    static boost::random::minstd_rand rng;
                    return rng() & 1 ? -*pt : *pt;
                }
            }
        }

    bool operator==(ec const &other) const { return this == &other || (this->_a == other._a && this->_b == other._b); }
    bool operator!=(ec const &other) const { return !(*this == other); }

    constexpr F j_invariant() const
    {
        auto const &a = _a, &b = _b;
        auto a3 = F(4)*a*a*a;
        auto b2 = F(27)*b*b;
        return F(1728) * a3 * (a3 + b2).inv();
    }

    constexpr static std::shared_ptr<ec> from_j(F const &j)
    {
        auto j2 = j*j;
        auto j3 = j2*j;
        F const f(1728);
        constexpr F two = F::one+F::one, three = two+F::one;
        F a = three * (f*j - j2);
        F b = two * (-j3 + two*f*j2 - f*f*j);
        auto E = std::make_shared<ec>(a, b);
        assert(E->j_invariant() == j);
        return E;
    }

    friend std::ostream& operator<<(std::ostream& o, ec const &E) { return o << "{y^2 = x^3 + (" << E._a << ")*x + (" << E._b << ")}"; }
};

namespace std {
    template <class F>
    struct hash<ec<F>>
    {
        size_t operator()(ec<F> const &E) const {
            return (0xec + 1111111*std::hash<F>()(E.a())) ^ std::hash<F>()(E.b());
        }
    };
}

