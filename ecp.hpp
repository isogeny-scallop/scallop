
#pragma once

#include <optional>
#include <tuple>
#include <memory>

#include "ec.hpp"

template <class F> F operator*(F const &x, std::optional<F> const &y) { return y ? x * y.value() : x; }
template <class F> std::optional<F> operator*(std::optional<F> const &x, std::optional<F> const &y) { return x ? x.value() * y : y; }

template <class F> class miller_fun;

template <class F>
class ecp
{
    private:
        std::shared_ptr<ec<F> const> E;

        /* z empty means z=1 */
        mutable F x,y;
        mutable std::optional<F> z;

    public:
        ecp(std::shared_ptr<ec<F> const> const &E) : E{E}, x{0}, y{1}, z{0} { }
        ecp(std::shared_ptr<ec<F> const> const &E, F const &xx, F const &yy) : E{E}, x{xx}, y{yy}, z{}
        {
            assert(!this->z);
            F const &x = this->x, &y = this->y;
            (void) x, (void) y, assert(y*y == x*x*x + this->E->a()*x + this->E->b());
        }
        ecp(std::shared_ptr<ec<F> const> const &E, F const &xx, F const &yy, F const &zz) : E{E}, x{xx}, y{yy}, z{zz}
        {
            F const &x = this->x, &y = this->y, &z = this->z.value();
            (void) x, (void) y, (void) z, assert(x || y || z);
            assert(y*y*z == x*x*x + this->E->a()*x*z*z + this->E->b()*z*z*z);
        }

        std::shared_ptr<ec<F> const> const &curve() const { return this->E; }

        F const &get_x() const { return this->x; }
        F const &get_y() const { return this->y; }
        std::optional<F> const &get_z() const { return this->z; }

        F const &affine_x() const {
            this->_normalize();
            if (this->z)
                throw std::logic_error("point at infinity has no affine coordinates");
            return this->x;
        }
        F const &affine_y() const {
            this->_normalize();
            if (this->z)
                throw std::logic_error("point at infinity has no affine coordinates");
            return this->y;
        }

        std::pair<F,F> _lambda(ecp const &other) const
        {
            auto const &P = *this, &Q = other;
            F const &x1 = P.x, &y1 = P.y;
            F const &x2 = Q.x, &y2 = Q.y;
            std::optional<F> const &z1 = P.z, &z2 = Q.z;

            if (x1*z2 == x2*z1) {
                F x1sq = x1 * x1;
                std::optional<F> z1sq = z1 * z1;
                return {x1sq+x1sq+x1sq + this->E->a()*z1sq, (y1+y1) * z1};
            }

            return {y2*z1 - y1*z2, x2*z1 - x1*z2};
        }

        ecp operator+(ecp const &other) const
        {
            if (*this->curve() != *other.curve())
                throw std::logic_error("points not on the same curve");

            if (!other) return *this;
            if (!*this) return other;

            F const &x1 = this->x, &y1 = this->y;
            std::optional<F> const &z1 = this->z;
            F const &x2 = other.x, &y2 = other.y;
            std::optional<F> const &z2 = other.z;

            if (x1*z2 == x2*z1 && y1*z2 != y2*z1) {
                assert(y2*z1 == -y1*z2);
                return {this->E};
            }

            auto lam = this->_lambda(other);
            F const &u = lam.first, &v = lam.second;

            F uu = u*u, vv = v*v;
            std::optional<F> zz = z1 * z2;

            F z3_ = vv*zz;
            F x3_ = uu*zz - vv*(x1*z2 + x2*z1);

            F vz1 = v*z1;
            F y3 = u*(x1*z3_ - x3_*z1) - y1*v*z3_;
            F x3 = x3_*vz1;
            F z3 = z3_*vz1;

            return ecp(this->curve(), x3,y3,z3);
        }
        ecp &operator+=(ecp const &other) { return *this = *this + other; }
        ecp operator-() const { ecp Q = *this; Q.y = -Q.y; return Q; }
        ecp operator-(ecp const &other) const { return *this + (-other); }
        ecp &operator-=(ecp const &other) { return *this = *this - other; }

        template <class Integer>
        friend ecp operator*(Integer const &k, ecp const &P) {
            if (!k)
                return P.E->point();
            ecp Q(P.E);
            for (ssize_t b = boost::multiprecision::msb(k); b >= 0; --b) {
                Q += Q;
                if (boost::multiprecision::bit_test(k, b))
                    Q += P;
            }
            return Q;
        }
        template <class Integer> friend ecp operator*(ecp const &P, Integer const &k) { return k * P; }
        template <class Integer> ecp &operator*=(Integer const &k) { return *this = k * *this; }

        operator bool() const { return !this->z || *this->z; }
        bool operator==(ecp const &other) const
        {
            return this->x * other.z == other.x * this->z
                && this->y * other.z == other.y * this->z;
        }
        bool operator!=(ecp const &other) const { return !(*this == other); }

        friend std::ostream &operator<<(std::ostream &o, ecp const &pt)
        {
            if (pt)
                return o << "(" << pt.affine_x() << "," << pt.affine_y() << ")";
            return o << "∞";
        }

        friend class miller_fun<F>;
        template <class FF, size_t N> friend FF weil_pairing(ecp<FF> const &, ecp<FF> const &, uint_t<N> const &);

    private:
        void _normalize() const
        {
            if (!this->z)
                return;
            if (*this->z) {
                F invz = this->z->inv();
                this->x *= invz;
                this->y *= invz;
                this->z.reset();
            }
            else {
                assert(!this->x && this->y);
                this->y = F::one;
            }
        }
};

namespace std {
    template <class F>
    struct hash<ecp<F>>
    {
        size_t operator()(ecp<F> const &P) const {
            size_t h = 0xec77;
            h += std::hash<ec<F>>()(*P.curve());
            h *= 1111111;
            if (P) {
                h ^= std::hash<F>()(P.affine_x());
                h *= 1111111;
                h ^= std::hash<F>()(P.affine_y());
            }
            return h;
        }
    };
}


template <class F>
class miller_fun
{
    private:
        class proj
        {
            public:
                proj(F const &a) : num{a} {}
                proj(std::pair<F,F> const &tup) : num{tup.first}, den{tup.second} { assert(*this->den); }
                proj(F const &a, F const &b) : num{a}, den{b} { assert(*this->den); }

                proj& operator*=(proj const &other)
                {
                    this->num *= other.num;
                    if (other.den) {
                        if (this->den)
                            *this->den *= *other.den;
                        else
                            this->den = *other.den;
                    }
                    return *this;
                }
                proj operator*(proj const &other) const { proj r = *this; return r *= other; }

                proj inv() const
                {
                    assert(this->num);
                    if (this->den)
                        return {*this->den, this->num};
                    return {F::one, this->num};
                }

                F const &x() const { return this->num; }
                std::optional<F> const &z() const { return this->den; }

                operator F() const
                {
                    this->_normalize();
                    return this->num;
                }

                void _normalize() const
                {
                    if (this->den) {
                        this->num *= this->den->inv();
                        this->den.reset();
                    }
                }

            private:
                mutable F num;
                mutable std::optional<F> den;
        };

        using line_fun = std::function<proj(ecp<F> const &)>;

        static line_fun _line_from_lambda(proj const &lam, ecp<F> const &P)
        {
            return [lam,P](ecp<F> const &T) -> proj {
                F r = lam.x() * (P.x - T.affine_x()*P.z) - lam.z().value() * (P.y - T.affine_y()*P.z);
                F s = lam.z().value() * P.z;
                return {r, s};
            };
        }

        static line_fun _line_(ecp<F> const &P, ecp<F> const &Q)
        {
            /* generic case */
            if (F zz = Q.x*P.z - P.x*Q.z) {
                F xx = Q.y*P.z - P.y*Q.z;
                return _line_from_lambda({xx, zz}, P);
            }
            if (!P) {
                /* line at infinity */
                if (!Q)
                    return [](ecp<F> const &T) -> proj { (void) T; return {F::one}; };
                /* vertical at Q */
                return _line_vert(Q);
            }
            /* vertical at P */
            if (!Q)
                return _line_vert(P);
            /* vertical at P, equivalently Q */
            if (P == -Q)
                return _line_vert(P);
            /* tangent at P=Q */
            if (P == Q)
                return _line_tang(P);
            __builtin_unreachable();
        }

        static line_fun _line_tang(ecp<F> const &P)
        {
            F xsq = P.x*P.x;
            std::optional<F> zsq = P.z*P.z;
            F u = xsq+xsq+xsq + P.curve()->a()*zsq;
            F v = (P.y+P.y) * P.z;
            return _line_from_lambda({u, v}, P);
        };

        static line_fun _line_vert(ecp<F> const &P)
        {
            return [P](ecp<F> const &T) -> proj {
                F r = T.affine_x()*P.z - P.x;
                std::optional<F> s = P.z;
                return s ? proj(r,*s) : proj(r);
            };
        };

    public:
        template <size_t N>
        miller_fun(ecp<F> P0, uint_t<N> const &n)
        {
            std::swap(this->P, P0);
            auto const &P = this->P;
            assert(n);
            assert(!(n * P));
            ecp<F> V = P;

            for (ssize_t i = boost::multiprecision::msb(n)-1; i >= 0; --i) {
                this->_data.push_back([    ](proj const &f, ecp<F> const &Q) { (void) Q; return f * f; });

                this->_data.push_back([V   ](proj const &f, ecp<F> const &Q) { return f * _line_(V, V)(Q); });
                V += V;
                this->_data.push_back([V   ](proj const &f, ecp<F> const &Q) -> std::optional<proj> { auto tmp = _line_(V,-V)(Q); if (tmp.x()) return f * tmp.inv(); return {}; });

                if (!boost::multiprecision::bit_test(n, i))
                    continue;

                this->_data.push_back([V,&P](proj const &f, ecp<F> const &Q) { return f * _line_(V, P)(Q); });
                V += P;
                this->_data.push_back([V   ](proj const &f, ecp<F> const &Q) -> std::optional<proj> { auto tmp = _line_(V,-V)(Q); if (tmp.x()) return f * tmp.inv(); return {}; });
            }
            assert(!V);
        }

        std::optional<proj> operator()(ecp<F> const &Q) const
        {
            proj f = F::one;
            if (!Q)
                return f;
            for (auto const &step: this->_data) {
                auto g = step(f, Q);
                if (!g)
                    return {};
                f = *g;
            }
            return f;
        }

    private:
        ecp<F> P{nullptr};
        std::vector<std::function<std::optional<proj>(proj const &, ecp<F> const &)>> _data;
};

template <class F, size_t N>
F weil_pairing(ecp<F> const &P, ecp<F> const &Q, uint_t<N> const &n)
{
    if (!P || !Q)
        return F::one;
    P._normalize();
    Q._normalize();
    auto f1 = miller_fun(P, n)(Q);
    if (!f1)
        return F::one;
    auto f2 = miller_fun(Q, n)(P);
    if (!f2)
        return F::one;
    F res = *f1 * f2->inv();
    return n & 1 ? -res : res;
}

