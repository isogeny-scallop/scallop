
#pragma once

#include <boost/multiprecision/cpp_int.hpp>

using namespace boost::multiprecision::literals;

template <size_t N>
using uint_t = boost::multiprecision::number<
                   boost::multiprecision::cpp_int_backend<
                   N, N,
                   boost::multiprecision::unsigned_magnitude,
                   boost::multiprecision::unchecked,
                   void
                   >
               >;

template <size_t N>
using int_t = boost::multiprecision::number<
                   boost::multiprecision::cpp_int_backend<
                   N, N,
                   boost::multiprecision::signed_magnitude,
                   boost::multiprecision::unchecked,
                   void
                   >
               >;


template <class I>
using dbl_t = std::conditional<
        std::numeric_limits<I>::is_signed,
         int_t<2*std::numeric_limits<I>::digits>,
        uint_t<2*std::numeric_limits<I>::digits>
    >::type;
static_assert(std::is_same<dbl_t< int_t<333>>,  int_t<666>>::value);  // check
static_assert(std::is_same<dbl_t<uint_t<333>>, uint_t<666>>::value);  // check

template <class I>
constexpr dbl_t<I> ext(I const &v)
{ return static_cast<dbl_t<I>>(v); }

template <class I>
using signify_t = std::conditional<
        std::numeric_limits<I>::is_signed,
        int_t<std::numeric_limits<I>::digits+0>,
        int_t<std::numeric_limits<I>::digits+1>
    >::type;
static_assert(std::is_same<signify_t< int_t<666>>, int_t<666>>::value);     // check
static_assert(std::is_same<signify_t<uint_t<666>>, int_t<667>>::value);     // check

template <class I>
using unsignify_t = uint_t<std::numeric_limits<I>::digits>;
static_assert(std::is_same<unsignify_t< int_t<666>>, uint_t<666>>::value);  // check
static_assert(std::is_same<unsignify_t<uint_t<666>>, uint_t<666>>::value);  // check

/************** Bigfloats too ! *******************/

#include <boost/multiprecision/cpp_bin_float.hpp>

using namespace boost::multiprecision;

template <size_t N>
using real_t = number<cpp_bin_float<N>>;

template <class T, size_t M, size_t N>
using matrix = std::array<std::array<T, N>, M>;

/************** ...and CRT! *******************/

#include <boost/integer/mod_inverse.hpp>

// assumes m,n coprime and > 1
template <class I>
std::pair<I,I> crt(std::pair<I,I> const &xm, std::pair<I,I> const &yn)
{
    auto const &[x,m] = xm;
    auto const &[y,n] = yn;
    assert(m >= 2);
    assert(n >= 2);
    auto const f = boost::integer::mod_inverse(m, n);
    assert(f * m % n == 1);
    auto const g = boost::integer::mod_inverse(n, m);
    assert(g * n % m == 1);
    auto mn = m * n;
    auto z = static_cast<I>((ext(x) * g % mn * n
                           + ext(y) * f % mn * m) % mn);
    assert(z % m == x);
    assert(z % n == y);
    return {z, mn};
}
