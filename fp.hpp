
#pragma once

#include "bigint.hpp"

#if 0
    #include "zmod.hpp"
    using fp = zmod<std::numeric_limits<decltype(p)>::digits, p>;
#else
    #include "skylake.hpp"
    using fp = _fp<void>;
#endif

