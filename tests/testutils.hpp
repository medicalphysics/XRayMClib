/*This file is part of XRayMClib.

XRayMClib is free software : you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

XRayMClib is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with XRayMClib. If not, see < https://www.gnu.org/licenses/>.

Copyright 2026 Erlend Andersen
*/

#pragma once

// Minimal test harness shared by all test executables.
//
// Every check prints a single "[ PASS ]" or "[ FAIL ]" line and failures are
// counted globally. A test executable ends with `return xraymctest::summary("name");`
// which prints a summary and returns EXIT_SUCCESS only if no check failed.

#include "xraymc/xraymcrandom.hpp"

#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <string_view>

namespace xraymctest {

inline int& failureCount()
{
    static int n = 0;
    return n;
}

inline int& checkCount()
{
    static int n = 0;
    return n;
}

inline bool report(bool ok, std::string_view name, const std::string& detail = {})
{
    ++checkCount();
    if (!ok)
        ++failureCount();
    std::cout << (ok ? "[ PASS ] " : "[ FAIL ] ") << name;
    if (!detail.empty())
        std::cout << " (" << detail << ")";
    std::cout << '\n';
    return ok;
}

inline bool check(bool ok, std::string_view name)
{
    return report(ok, name);
}

// |value - expected| <= tolerance
inline bool checkNear(double value, double expected, double tolerance, std::string_view name)
{
    const bool ok = std::abs(value - expected) <= tolerance;
    std::ostringstream s;
    s.precision(10);
    s << "got " << value << ", expected " << expected << " +/- " << tolerance;
    return report(ok, name, s.str());
}

// |value / expected - 1| <= relativeTolerance
inline bool checkRelative(double value, double expected, double relativeTolerance, std::string_view name)
{
    const bool ok = std::abs(value - expected) <= std::abs(expected) * relativeTolerance;
    std::ostringstream s;
    s.precision(8);
    s << "got " << value << ", expected " << expected << " +/- " << relativeTolerance * 100 << "%";
    return report(ok, name, s.str());
}

inline bool checkVectorNear(const std::array<double, 3>& value, const std::array<double, 3>& expected, double tolerance, std::string_view name)
{
    bool ok = true;
    for (std::size_t i = 0; i < 3; ++i)
        ok = ok && std::abs(value[i] - expected[i]) <= tolerance;
    std::ostringstream s;
    s.precision(10);
    s << "got {" << value[0] << ", " << value[1] << ", " << value[2] << "}, expected {"
      << expected[0] << ", " << expected[1] << ", " << expected[2] << "}";
    return report(ok, name, s.str());
}

// Two independent estimates a +/- sa and b +/- sb agree within z standard deviations.
// Monte Carlo tests use this with z around 4 so that random failures are very unlikely.
inline bool checkStatisticallyEqual(double a, double sa, double b, double sb, double z, std::string_view name)
{
    const auto s = std::sqrt(sa * sa + sb * sb);
    const auto t = s > 0 ? std::abs(a - b) / s : (a == b ? 0.0 : std::numeric_limits<double>::infinity());
    std::ostringstream d;
    d.precision(6);
    d << a << " +/- " << sa << " vs " << b << " +/- " << sb << ", z = " << t << " (limit " << z << ")";
    return report(t <= z, name, d.str());
}

// Deterministic random state so that sampling tests are reproducible.
inline xraymc::RandomState seededState(std::uint64_t seed = 0x853c49e6748fea9bULL)
{
    std::uint64_t s[2] = { seed, 0xda3e39cb94b95bdbULL };
    return xraymc::RandomState(s);
}

// Path for temporary files produced by tests.
inline std::filesystem::path tempPath(const std::string& filename)
{
    return std::filesystem::temp_directory_path() / ("xraymclib_test_" + filename);
}

inline int summary(std::string_view suite)
{
    std::cout << '\n'
              << suite << ": " << checkCount() - failureCount() << " of " << checkCount() << " checks passed";
    if (failureCount() > 0)
        std::cout << ", " << failureCount() << " FAILED";
    std::cout << std::endl;
    return failureCount() == 0 && checkCount() > 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}

}
