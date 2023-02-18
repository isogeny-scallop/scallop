
#pragma once

#include <optional>
#include <memory>

#include "ec.hpp"
#include "ecp.hpp"

template <class F>
class isog
{
    private:
        std::shared_ptr<ec<F> const> _domain;
        ecp<F> const _kernel;
        std::shared_ptr<ec<F> const> _codomain;
        size_t _degree;

        std::optional<std::pair<F,F>> _scaling;

        struct _tup { F X, Y, Z, GQxyz, TQ, UQ, VQ; };
        std::vector<_tup> _data;

    public:
        isog(ecp<F> const &K) : _domain{K.curve()}, _kernel{K}
        {
            auto const &a = this->_domain->a(), &b = this->_domain->b();

            this->_degree = 1;
            for (ecp<F> Q = this->_kernel, prevQ {this->_domain};
                    Q && Q != -prevQ;
                    Q += this->_kernel, this->_degree += 2) {
                auto &XQ = Q.get_x();
                auto &YQ = Q.get_y();
                auto &ZQ = Q.get_z();
                auto XQ2 = XQ*XQ;
                auto ZQ2 = ZQ*ZQ;
                auto GQx = XQ2+XQ2+XQ2 + a*ZQ2;
                auto GQy = (-YQ-YQ) * ZQ;
                auto GQxyz = GQx * GQy * ZQ;
                auto TQ = (YQ ? GQx+GQx : GQx) * ZQ2;
                auto UQ = GQy * GQy * ZQ;
                auto VQ = ZQ2 * ZQ2;
                this->_data.emplace_back(XQ, YQ, ZQ.value_or(F::one),
                                         GQxyz,
                                         TQ, UQ, VQ.value_or(F::one));
                if (Q == -Q) {
                    ++this->_degree;
                    break;
                }
                prevQ = Q;
            }

            F t = F::zero, w = F::zero, d = F::one;
            for (auto const &[XQ, YQ, ZQ, GQxyz, TQ, UQ, VQ]: this->_data) {
                t = (t*VQ + TQ*d) * ZQ;
                w = w*VQ*ZQ + UQ*d + XQ*TQ*d;
                d *= VQ * ZQ;
            }
            d = d.inv();
            t *= d;
            w *= d;

            auto A = a - F(5)*t;
            auto B = b - F(7)*w;

            this->_codomain = std::make_shared<ec<F> const>(A, B);
        }

        std::shared_ptr<ec<F> const> const &domain() const { return this->_domain; }
        std::shared_ptr<ec<F> const> const &codomain() const { return this->_codomain; }
        size_t degree() const { return this->_degree; }

        ecp<F> operator()(ecp<F> const &P) const
        {
            assert(*P.curve() == *this->_domain);
            if (!P)
                return {this->_codomain};
            F const &x = P.affine_x();
            F const &y = P.affine_y();
            F X = x, Y = y, Z = F::one;
            for (auto const &[XQ, YQ, ZQ, GQxyz, TQ, UQ, VQ]: this->_data) {
                auto const Dx = x*ZQ - XQ,
                           Dy = y*ZQ - YQ;
                auto xx = (TQ*Dx + UQ) * Dx*ZQ;
                auto yy = (UQ*ZQ*(y+y) + (TQ*Dy - GQxyz) * Dx) * ZQ;
                auto zz = Dx*Dx*Dx * VQ;
                X = X*zz + xx*Z;
                Y = Y*zz - yy*Z;
                Z *= zz;
            }
            if (this->_scaling) {
                X *= this->_scaling->first;
                Y *= this->_scaling->second;
            }
            return {this->_codomain, X, Y, Z};
        }

        void _scale(F const &u, std::shared_ptr<ec<F> const> new_codomain={})
        {
            if (!this->_scaling)
                this->_scaling = {F::one, F::one};
            auto ui = u.inv(), ui2 = ui*ui, ui3 = ui2*ui;
            this->_scaling->first *= ui2;
            this->_scaling->second *= ui3;
#ifdef NDEBUG
            if (!new_codomain)
#endif
            {
                auto A = this->_codomain->a() * ui2 * ui2,
                     B = this->_codomain->b() * ui3 * ui3;
#ifndef NDEBUG
                if (new_codomain)
                    assert(new_codomain->a() == A && new_codomain->b() == B);
                else
#endif
                new_codomain = std::make_shared<ec<F> const>(A, B);
            }
            std::swap(this->_codomain, new_codomain);
        }

        F _scaling_factor() const
        {
            if (!this->_scaling)
                return F::one;
            return this->_scaling->first * this->_scaling->second.inv();
        }

        friend std::ostream& operator<<(std::ostream& o, isog const &phi)
        {
            return o << "isogeny of degree " << phi._degree
                     << " from " << *phi._domain
                     << " to " << *phi._codomain
                     << " with kernel <" << phi._kernel << ">";
        }
};


template <class F>
using isog_chain = std::vector<isog<F>>;

template <class F>
ecp<F> chain_eval(isog_chain<F> const &chain, ecp<F> Q)
{
    for (auto const &phi: chain)
        Q = phi(Q);
    return Q;
}

//TODO optimize
template <class O, class F>
isog_chain<F> factored_isogeny(ecp<F> ker, O ord, F const &u = F::one)
{
    isog_chain<F> ret;
    for (size_t i = 0; i < O::base.size(); ++i) {
        while (ord[i]) {
            O const cof = ord.dec(i);
            auto K = cof.value() * ker;
            ord = std::move(cof);
            isog phi(K);
            assert(phi.degree() == O::base[i]);
            ker = phi(ker);
            ret.push_back(std::move(phi));
        }
    }
    assert(!ker);
    assert(ord.value() == 1);
    if (ret.empty())
        ret.push_back(isog(ker));

    if (&u != &F::one)
        ret.back()._scale(u);

    return ret;
}

template <class F>
std::shared_ptr<ec<F> const> const &chain_domain(isog_chain<F> const &chain)
{
    assert(!chain.empty());
    return chain.front().domain();
}

template <class F>
std::shared_ptr<ec<F> const> const &chain_codomain(isog_chain<F> const &chain)
{
    assert(!chain.empty());
    return chain.back().codomain();
}

