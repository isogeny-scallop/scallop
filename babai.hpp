
#pragma once

#include "bigint.hpp"
#include "params.hpp"

template <class T1, class T2, size_t N>
constexpr T1 inner_product(std::array<T1, N> const &a, std::array<T2, N> const &b) {
    T1 res;
    for (size_t i = 0; i < N; i++)
        res += a[i] * T1(b[i]);
    return res;
}

/* assumes matrix is non-singular */
template <class T, size_t M, size_t N, size_t P=200>
constexpr matrix<real_t<P>, M, N> gram_schmidt(matrix<T, M, N> const &mat) {
    matrix<real_t<P>, M, N> B;

    for (size_t i = 0; i < M; i++) {
        for (size_t k = 0; k < N; k++)
            B[i][k] = mat[i][k];
    }

    for (size_t i = 0; i < M; i++) {
        auto ip = inner_product(B[i], B[i]);
        for (size_t j = i + 1; j < M; j++) {
            auto ip2 = inner_product(B[i], mat[j]) / ip;
            for (size_t k = 0; k < N; k++)
                B[j][k] -= B[i][k] * ip2;
        }
    }

    return B;
}

/*
   Babai nearest-plane algorithm.

   Finds a vector close to `v` in the lattice spanned by
   `relations`. Overwrites `v`.

   `GS` is the Gram-Schmidt matrix of `relations`. Both `relations`
   and `GS` could be computed at compile time, if I knew how to do
   it...
*/
template <size_t N, size_t B, class R>
std::array<int_t<B>, N> babai(std::array<int_t<B>, N> &v,
                       matrix<int16_t, N, N> const &relations,
                       matrix<R, N, N> const &GS) {
    using I = int_t<B>;
    for (int i = N - 1; i >= 0; i--) {
        // wtf is this syntax?!
        I r = round( inner_product(GS[i], v) / inner_product(GS[i], GS[i]) ).template convert_to<I>();
        for (size_t j = 0; j < N; j++) {
            v[j] -= r * relations[i][j];
        }
    }
    return v;
}

template <size_t N, size_t B, class R>
std::array<int16_t, N> babai(int_t<B> const &exp,
                            matrix<int16_t, N, N> const &relations,
                            matrix<R, N, N> const &GS) {
    std::array<int_t<B>, N> v {0};
    v[geni] = exp;
    babai(v, relations, GS);
    std::array<int16_t, N> res;
    for (size_t i = 0; i < N; i++) {
        assert(v[i] >= std::numeric_limits<int16_t>::min() && v[i] <= std::numeric_limits<int16_t>::max());
        res[i] = static_cast<int16_t>(v[i]);
    }
    return res;
}
