/*This file is part of DXMClib.

DXMClib is free software : you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

DXMClib is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with DXMClib. If not, see < https://www.gnu.org/licenses/>.

Copyright 2023 Erlend Andersen
*/

#pragma once

#include "dxmc/particle.hpp"
#include "dxmc/vectormath.hpp"
#include "dxmc/world/visualizationintersectionresult.hpp"
#include "dxmc/world/worldintersectionresult.hpp"

#include <array>
#include <optional>

namespace dxmc {
namespace basicshape {
    namespace cylinderZ {

        static constexpr bool pointInside(const std::array<double, 3>& pos, const std::array<double, 3>& center, const double radii, const double half_height)
        {
            const std::array<double, 2> dp = { pos[0] - center[0], pos[1] - center[1] };
            return center[0] - radii <= pos[0] && pos[0] <= center[0] + radii && center[1] - radii <= pos[1] && pos[1] <= center[1] + radii
                && (center[2] - half_height < pos[2]) && (pos[2] < center[2] + half_height) && ((dp[0] * dp[0] + dp[1] * dp[1]) < radii * radii);
        }

        static constexpr std::optional<double> intersectCylinderWall(const Particle& p, const std::array<double, 3>& center, const double radii)
        {
            // nummeric stable ray sphere intersection in 2D
            const auto a = p.dir[0] * p.dir[0] + p.dir[1] * p.dir[1];
            if (a < std::numeric_limits<double>::epsilon())
                return std::nullopt;

            const auto r2 = radii * radii;
            const std::array f = { p.pos[0] - center[0], p.pos[1] - center[1] };

            // positive b mean center of sphere is in front of ray
            const auto b = -f[0] * p.dir[0] - f[1] * p.dir[1];

            // positive c means ray starts outside of sphere
            const auto c = f[0] * f[0] + f[1] * f[1] - r2;

            if ((c > 0) && (b < 0)) {
                // if ray starts outside sphere and center is begind ray
                // we exit early
                return std::nullopt;
            }

            const std::array delta1 = { f[0] + b * p.dir[0] / a, f[1] + b * p.dir[1] / a };

            const auto delta = r2 - (delta1[0] * delta1[0] + delta1[1] * delta1[1]);

            if (delta < 0) {
                // no solution to the quadratic equation (we miss)
                return std::nullopt;
            }

            const int sign = b > 0 ? 1 : -1;
            const auto q = b + sign * std::sqrt(a * delta);

            if (c < 0 && b > 0) {
                // inside sphere
                return std::make_optional(q / a);
            } else {
                return std::make_optional(c / q);
            }
        }

        static constexpr std::optional<std::array<double, 2>> intersectCylinderWallInterval(const Particle& p, const std::array<double, 3>& center, const double radii)
        {
            // nummeric stable ray sphere intersection in 2D
            const auto a = p.dir[0] * p.dir[0] + p.dir[1] * p.dir[1];
            if (a < std::numeric_limits<double>::epsilon())
                return std::nullopt;

            const auto r2 = radii * radii;
            const std::array f = { p.pos[0] - center[0], p.pos[1] - center[1] };

            // positive b mean center of sphere is in front of ray
            const auto b = -f[0] * p.dir[0] - f[1] * p.dir[1];

            // positive c means ray starts outside of sphere
            const auto c = f[0] * f[0] + f[1] * f[1] - r2;

            if ((c > 0) && (b < 0)) {
                // if ray starts outside sphere and center is begind ray
                // we exit early
                return std::nullopt;
            }

            const std::array delta1 = { f[0] + b * p.dir[0] / a, f[1] + b * p.dir[1] / a };

            const auto delta = r2 - (delta1[0] * delta1[0] + delta1[1] * delta1[1]);

            if (delta < 0) {
                // no solution to the quadratic equation (we miss)
                return std::nullopt;
            }

            const int sign = b > 0 ? 1 : -1;
            const auto q = b + sign * std::sqrt(a * delta);

            if (c < 0 && b > 0) {
                // inside sphere
                std::array t { c / q, q / a };
                return std::make_optional(t);
            } else {
                std::array t { c / q, q / a };
                return std::make_optional(t);
            }
        }

        static std::optional<double> intersectCylinderDiscIntervalZ(const Particle& p, const std::array<double, 3>& center, const double radii)
        {
            if (std::abs(p.dir[2]) <= std::numeric_limits<double>::epsilon())
                return std::nullopt;

            const auto tz = (center[2] - p.pos[2]) / p.dir[2];
            const auto xz = p.pos[0] + p.dir[0] * tz - center[0];
            const auto yz = p.pos[1] + p.dir[1] * tz - center[1];
            if (xz * xz + yz * yz <= radii * radii)
                return std::make_optional(tz);
            return std::nullopt;
        }

        static std::optional<double> intersectCylinderDiscZ(const Particle& p, const std::array<double, 3>& center, const double radii)
        {
            if (std::abs(p.dir[2]) <= std::numeric_limits<double>::epsilon())
                return std::nullopt;

            if (p.pos[2] > center[2] && p.dir[2] >= 0)
                return std::nullopt;
            if (p.pos[2] < center[2] && p.dir[2] <= 0)
                return std::nullopt;

            const auto tz = (center[2] - p.pos[2]) / p.dir[2];
            const auto xz = p.pos[0] + p.dir[0] * tz - center[0];
            const auto yz = p.pos[1] + p.dir[1] * tz - center[1];
            if (xz * xz + yz * yz <= radii * radii)
                return std::make_optional(tz);
            return std::nullopt;
        }

        static WorldIntersectionResult intersect(const Particle& p, const std::array<double, 3>& center, const double radii, const double half_height)
        {
            auto t_cand = intersectCylinderWall(p, center, radii);
            if (t_cand) {
                // we need to be conservative else we will miss intersections on plane cylinder intersection
                const auto tz = std::nextafter(p.pos[2] + p.dir[2] * *t_cand, 0.0);
                if (!(center[2] - half_height <= tz && tz <= center[2] + half_height))
                    t_cand.reset();
            }

            std::array centerDisc = { center[0], center[1], center[2] - half_height };
            auto tdisc1_cand = intersectCylinderDiscZ(p, centerDisc, radii);

            centerDisc[2] = center[2] + half_height;
            auto tdisc2_cand = intersectCylinderDiscZ(p, centerDisc, radii);

            const auto t_cand_min = std::min({ t_cand, tdisc1_cand, tdisc2_cand }, [](const auto& lh, const auto& rh) -> bool { return lh.value_or(std::numeric_limits<double>::max()) < rh.value_or(std::numeric_limits<double>::max()); });
            WorldIntersectionResult res;
            if (t_cand_min) {
                res.intersection = *t_cand_min;
                res.intersectionValid = true;
                res.rayOriginIsInsideItem = pointInside(p.pos, center, radii, half_height);
            }

            return res;
        }

        static constexpr std::optional<std::array<double, 2>> intersectForwardInterval(const Particle& p, const std::array<double, 3>& center, const double radii, const double half_height)
        {
            auto t_cand = intersectCylinderWallInterval(p, center, radii);
            if (!t_cand) {
                if (std::abs(p.dir[2]) == 1) {
                    const auto dx = p.pos[0] - center[0];
                    const auto dy = p.pos[1] - center[1];
                    if (dx * dx + dy * dy < radii * radii) {
                        const auto t1 = std::max(center[2] + half_height - p.pos[2], 0.0);
                        const auto t2 = std::max(center[2] - half_height - p.pos[2], 0.0);
                        std::array<double, 2> tz;
                        if (t1 < t2) {
                            tz = { t1, t2 };
                        } else {
                            tz = { t2, t1 };
                        }
                        if (tz[0] < tz[1]) {
                            return std::make_optional(tz);
                        }
                    }
                }
                return std::nullopt;
            }

            const auto tz1 = (center[2] + half_height - p.pos[2]) / p.dir[2];
            const auto tz2 = (center[2] - half_height - p.pos[2]) / p.dir[2];

            const auto tz_min = tz1 < tz2 ? tz1 : tz2;
            const auto tz_max = tz1 < tz2 ? tz2 : tz1;
            auto& t = t_cand.value();

            t[0] = std::max(t[0], tz_min);
            t[1] = std::min(t[1], tz_max);
            if (t[0] < 0)
                t[0] = 0;
            if (t[0] < t[1])
                return t_cand;
            return std::nullopt;
        }

        template <typename U>
        VisualizationIntersectionResult<U> intersectVisualization(const Particle& p, const std::array<double, 3>& center, const double radii, const double half_height)
        {
            auto t_cand = intersectCylinderWall(p, center, radii);
            if (t_cand) {
                // we need to be conservative else we will miss intersections on plane cylinder intersection
                const auto tz = std::nextafter(p.pos[2] + p.dir[2] * *t_cand, 0.0);
                if (!(center[2] - half_height <= tz && tz <= center[2] + half_height))
                    t_cand.reset();
            }

            std::array centerDisc = { center[0], center[1], center[2] - half_height };
            auto tdisc1_cand = intersectCylinderDiscZ(p, centerDisc, radii);

            centerDisc[2] = center[2] + half_height;
            auto tdisc2_cand = intersectCylinderDiscZ(p, centerDisc, radii);

            constexpr auto m = std::numeric_limits<double>::max();

            VisualizationIntersectionResult<U> res;

            if (t_cand.value_or(m) < tdisc1_cand.value_or(m)) {
                if (t_cand.value_or(m) < tdisc2_cand.value_or(m)) {
                    if (t_cand) {
                        res.intersection = *t_cand;
                        res.intersectionValid = true;
                        const auto x = center[0] - (p.pos[0] + res.intersection * p.dir[0]);
                        const auto y = center[1] - (p.pos[1] + res.intersection * p.dir[1]);
                        const auto ll = 1 / std::sqrt(x * x + y * y);
                        res.normal[0] = x * ll;
                        res.normal[1] = y * ll;
                    }
                } else {
                    if (tdisc2_cand) {
                        res.intersection = *tdisc2_cand;
                        res.intersectionValid = true;
                        res.normal[2] = -1;
                    }
                }
            } else {
                if (tdisc2_cand.value_or(m) < tdisc1_cand.value_or(m)) {
                    if (tdisc2_cand) {
                        res.intersection = *tdisc2_cand;
                        res.intersectionValid = true;
                        res.normal[2] = -1;
                    }
                } else {
                    if (tdisc1_cand) {
                        res.intersection = *tdisc1_cand;
                        res.intersectionValid = true;
                        res.normal[2] = 1;
                    }
                }
            }

            if (res.intersectionValid) {
                res.rayOriginIsInsideItem = pointInside(p.pos, center, radii, half_height);
            }
            if (res.rayOriginIsInsideItem) {
                res.normal = vectormath::scale(res.normal, -1.0);
            }
            return res;
        }
    }
}
}