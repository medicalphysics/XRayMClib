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

Copyright 2019 Erlend Andersen
*/

#pragma once // include guard

#include "dxmc/beamfilters.hpp"
#include "dxmc/constants.hpp"
#include "dxmc/dxmcrandom.hpp"
#include "dxmc/exposure.hpp"
#include "dxmc/floating.hpp"
#include "dxmc/lowenergycorrectionmodel.hpp"
#include "dxmc/progressbar.hpp"
#include "dxmc/transport.hpp"
#include "dxmc/tube.hpp"
#include "dxmc/world.hpp"

#include <array>
#include <concepts>
#include <memory>
#include <type_traits>
#include <vector>

namespace dxmc {

/**
 * @brief Abstract class for Sources
 * This class is intended to be a base class for implementation of arbitrary sources.
 * The main purpose for this class and it's derived classes is to specify geometry, x-ray source
 * and calibration to absolute dose.
 */
template <Floating T = double>
class Source {
public:
    /**
     * @brief This enumeration is not used by dxmc but is provided for
     * convenience in case of down casting of Source types.
     */
    enum class Type {
        None,
        CTSpiral,
        CTAxial,
        DX,
        CTDual,
        Pencil,
        Isotropic,
        IsotropicCT,
        CTTopogram,
        CBCT,
        Other
    };

private:
    /**
     * @brief Normalization of both cosine direction vectors to lenght of one
     * Note that direction cosines must be orthogonal.
     * @param
     */
    void normalizeDirectionCosines(void)
    {
        vectormath::normalize(&m_directionCosines[0]);
        vectormath::normalize(&m_directionCosines[3]);
    }

protected:
    std::array<T, 3> m_position = { 0, 0, 0 };
    std::array<T, 6> m_directionCosines = { 1, 0, 0, 0, 1, 0 };
    std::uint64_t m_historiesPerExposure = 1E6;
    Type m_type = Type::None;

public:
    /**
     * @brief Default constructor whos only job is to set source type to Type::None.
     */
    Source()
        : m_type(Type::None)
    {
    }
    virtual ~Source() = default;

    /**
     * @brief Generate an exposure
     * Generate an exposure that specifies emitting numberOfHistoriesPerExposure
     * from a specific point and with a specific direction and collimation.
     * @param i Exposure number, should be less that numberOfExposures.
     * @return Exposure object
     */
    virtual Exposure<T> getExposure(std::uint64_t i) const = 0;
    /**
     * @brief Return max photon energy produced by this source
     * This is used to generate attenuation look up tables
     * @return max photon energy with units keV
     */
    virtual T maxPhotonEnergyProduced() const { return Tube<T>::maxVoltage(); }
    /**
     * @brief Set source position
     * @param position
     */
    void setPosition(const std::array<T, 3>& position) { m_position = position; }
    /**
     * @brief Set source position
     * @param x
     * @param y
     * @param z
     */
    void setPosition(T x, T y, T z)
    {
        m_position[0] = x;
        m_position[1] = y;
        m_position[2] = z;
    }
    /**
     * @brief Get source position
     * @param
     * @return
     */
    std::array<T, 3>& position(void) { return m_position; }
    /**
     * @brief Get source position
     * @param
     * @return
     */
    const std::array<T, 3>& position(void) const { return m_position; }
    /**
     * @brief Set direction cosines
     * Direction cosines are two orthogonal vectors used to identify orientation of the source
     * the vectors are given as an array of size 6 with the x vector first and y vector last.
     * The cross product (x * y) gives the beam direction vector.
     * @param cosines
     */
    void setDirectionCosines(const std::array<T, 6>& cosines)
    {
        m_directionCosines = cosines;
        normalizeDirectionCosines();
    }
    /**
     * @brief Get direction cosines
     * Direction cosines are two orthogonal vectors used to identify orientation of the source
     * the vectors are given as an array of size 6 with the x vector first and y vector last.
     * The cross product (x * y) gives the beam direction vector.
     * @return
     */
    const std::array<T, 6>& directionCosines(void) const { return m_directionCosines; }
    /**
     * @brief Get direction cosines
     * Direction cosines are two orthogonal vectors used to identify orientation of the source
     * the vectors are given as an array of size 6 with the x vector first and y vector last.
     * The cross product (x * y) gives the beam direction vector.
     * @return
     */
    std::array<T, 6>& directionCosines(void) { return m_directionCosines; }
    /**
     * @brief Set histories per exposure
     * Histories per exposure denotes number of photon histories emitted per exposure
     * @param histories
     */
    void setHistoriesPerExposure(std::uint64_t histories) { m_historiesPerExposure = histories; }
    /**
     * @brief Get histories per exposure
     * Histories per exposure denotes number of photon histories emitted per exposure
     * @return
     */
    std::uint64_t historiesPerExposure(void) const { return m_historiesPerExposure; }
    /**
     * @brief Total exposures for this source
     * @return number of exposures for this source
     */
    virtual std::uint64_t totalExposures(void) const = 0;
    /**
     * @brief Get source type
     * @return
     */
    Source::Type type() const { return m_type; }
    /**
     * @brief Get calibration value for absolute dose
     * For sources this function must be implemented for calculation of absolute doses
     * See implementation examples for derived classes DXSource and CTAxialSource.
     * @param progress Progress bar for feedback if calculations are time
     * consuming.
     * @return calibration value such as dose = calibrationvalue * relative dose
     */
    virtual T getCalibrationValue(LOWENERGYCORRECTION model = LOWENERGYCORRECTION::NONE, ProgressBar<T>* progress = nullptr) const = 0;

    virtual bool isValid(void) const = 0;
    virtual bool validate(void) = 0;
    virtual void updateFromWorld(const World<T>& world) { }
};

template <Floating T = double>
class PencilSource final : public Source<T> {
public:
    PencilSource()
    {
        this->m_type = Source<T>::Type::Pencil;
    }
    virtual ~PencilSource() = default;
    Exposure<T> getExposure(std::uint64_t i) const override
    {
        constexpr std::array<T, 2> collimationAngles { 0, 0 };
        Exposure<T> exposure(this->m_position,
            this->m_directionCosines,
            collimationAngles,
            this->m_historiesPerExposure);

        exposure.setMonoenergeticPhotonEnergy(m_photonEnergy);
        return exposure;
    }

    void setPhotonEnergy(T energy)
    {
        if (energy > ELECTRON_REST_MASS<T>())
            m_photonEnergy = ELECTRON_REST_MASS<T>();
        else if (energy < 1)
            m_photonEnergy = 1;
        else
            m_photonEnergy = energy;
    }
    T photonEnergy() const { return m_photonEnergy; }
    T maxPhotonEnergyProduced() const override { return m_photonEnergy; }
    std::uint64_t totalExposures(void) const override { return m_totalExposures; }
    void setTotalExposures(std::uint64_t exposures)
    {
        if (exposures > 0)
            m_totalExposures = exposures;
    }
    void setAirDose(T Gycm2)
    {
        if (Gycm2 > 0.0)
            m_airDose = Gycm2;
    }
    T airDose(void) const { return m_airDose; } // Gycm2

    T getCalibrationValue(LOWENERGYCORRECTION model, ProgressBar<T>* = nullptr) const override
    {
        Material airMaterial("Air, Dry (near sea level)");
        const T nHistories = totalExposures() * this->historiesPerExposure();
        const T mea = static_cast<T>(airMaterial.getMassEnergyAbsorbtion(m_photonEnergy));
        const T calcOutput = nHistories * m_photonEnergy * mea * KEV_TO_MJ<T>();
        const T factor = m_airDose / calcOutput;
        return factor;
    }
    bool isValid(void) const override { return true; }
    bool validate(void) override { return true; }

protected:
    T m_photonEnergy = 100;
    T m_airDose = 1;
    std::uint64_t m_totalExposures = 10;
};

template <Floating T = double>
class IsotropicSource : public Source<T> {
public:
    IsotropicSource()
        : Source<T>()
    {
        this->m_type = Source<T>::Type::Isotropic;
        // initializing a specterdistribution
        const std::vector<T> energies = { 60.0 };
        const std::vector<T> weights = { 1.0 };
        m_maxPhotonEnergy = 60.0;
        m_specterDistribution = SpecterDistribution<T>(weights, energies);
    }
    virtual ~IsotropicSource() = default;
    Exposure<T> getExposure(std::uint64_t exposureNumber) const override
    {
        constexpr T weight { 1 };
        Exposure<T> exposure(this->m_position,
            this->m_directionCosines,
            this->m_collimationAngles,
            this->m_historiesPerExposure,
            weight,
            &m_specterDistribution);

        return exposure;
    }
    T maxPhotonEnergyProduced() const override
    {
        return m_maxPhotonEnergy;
    }
    void setTotalExposures(std::uint64_t nExposures)
    {
        m_totalExposures = nExposures;
    }
    std::uint64_t totalExposures() const override
    {
        return m_totalExposures;
    }

    T getCalibrationValue(LOWENERGYCORRECTION model, ProgressBar<T>* = nullptr) const override
    {
        return T { 1 };
    }

    bool isValid() const override { return true; }
    bool validate() override { return true; }

    void setSpecter(const std::vector<T>& weights, const std::vector<T>& energies)
    {
        const auto maxEnergyElement = std::max_element(energies.cbegin(), energies.cend());
        m_maxPhotonEnergy = *maxEnergyElement;
        m_specterDistribution = SpecterDistribution<T>(weights, energies);
    }

    void setCollimationAngles(T x0, T x1, T y0, T y1)
    {
        m_collimationAngles[0] = std::clamp(x0, -PI_VAL<T>() / 2, PI_VAL<T>() / 2);
        m_collimationAngles[1] = std::clamp(x1, -PI_VAL<T>() / 2, PI_VAL<T>() / 2);
        m_collimationAngles[2] = std::clamp(y0, -PI_VAL<T>(), PI_VAL<T>());
        m_collimationAngles[3] = std::clamp(y1, -PI_VAL<T>(), PI_VAL<T>());
    }

    void setCollimationAngles(T xRad, T yRad)
    {
        const T x0 = -xRad / 2;
        const T x1 = xRad / 2;
        const T y0 = -yRad / 2;
        const T y1 = yRad / 2;
        setCollimationAngles(x0, x1, y0, y1);
    }
    const std::array<T, 4>& collimationAngles() const
    {
        return m_collimationAngles;
    }

protected:
    std::uint64_t m_totalExposures = 1;
    std::array<T, 4> m_collimationAngles = { 0, 0, 0, 0 };
    SpecterDistribution<T> m_specterDistribution;
    T m_maxPhotonEnergy = 1.0;
};
template <Floating T = double>
class IsotropicCTSource final : public IsotropicSource<T> {
public:
    IsotropicCTSource()
        : IsotropicSource<T>()
    {
        this->m_type = Source<T>::Type::IsotropicCT;
    }
    virtual ~IsotropicCTSource() = default;
    Exposure<T> getExposure(std::uint64_t exposureNumber) const override
    {
        std::array<T, 3> rotAxis = { 0, 0, 1 };

        const auto angle = (exposureNumber * 2 * PI_VAL<T>()) / this->m_totalExposures;
        auto cosines = this->m_directionCosines;
        auto pos = this->m_position;

        vectormath::rotate(&pos[0], &rotAxis[0], angle);
        vectormath::rotate(&cosines[0], &rotAxis[0], angle);
        vectormath::rotate(&cosines[3], &rotAxis[0], angle);

        constexpr T weight { 1 };
        Exposure<T> exposure(pos,
            cosines,
            this->m_collimationAngles,
            this->m_historiesPerExposure,
            weight,
            &(this->m_specterDistribution));

        return exposure;
    }
};

template <Floating T = double>
class DAPSource : public Source<T> {
public:
    DAPSource()
        : Source<T>()
    {
        m_sdd = 1000.0;
        m_fieldSize[0] = 100.0;
        m_fieldSize[1] = 100.0;
        updateFieldSize(m_fieldSize);
        m_tube.setAlFiltration(2.0);
        this->setDirectionCosines(zeroDirectionCosines());
    }
    virtual ~DAPSource() = default;

    Tube<T>& tube(void)
    {
        m_specterValid = false;
        return m_tube;
    }

    const Tube<T>& tube(void) const { return m_tube; }

    T maxPhotonEnergyProduced() const override { return m_tube.voltage(); }

    void setCollimationAngles(const std::array<T, 2>& angles)
    {
        std::array<T, 2> ang = { std::abs(angles[0]), std::abs(angles[1]) };
        updateCollimationAngles(ang);
    }
    const std::array<T, 2>& collimationAngles(void) const
    {
        return m_collimationAngles;
    }

    void setCollimationAnglesDeg(const std::array<T, 2>& angles)
    {
        std::array<T, 2> ang = { std::abs(angles[0]) * DEG_TO_RAD<T>(), std::abs(angles[1]) * DEG_TO_RAD<T>() };
        updateCollimationAngles(ang);
    }
    const std::array<T, 2> collimationAnglesDeg(void) const
    {
        std::array<T, 2> a { m_collimationAngles[0] * RAD_TO_DEG<T>(),
            m_collimationAngles[1] * RAD_TO_DEG<T>() };
        return a;
    }

    void setFieldSize(const std::array<T, 2>& mm)
    {
        std::array<T, 2> amm = { std::abs(mm[0]), std::abs(mm[1]) };
        updateFieldSize(amm);
    }

    const std::array<T, 2>& fieldSize(void) const { return m_fieldSize; }

    void setSourceDetectorDistance(T mm)
    {
        m_sdd = std::abs(mm);
        updateFieldSize(m_fieldSize);
    }
    T sourceDetectorDistance() const { return m_sdd; }

    void setSourceAngles(T primaryAngle, T secondaryAngle)
    {
        constexpr T ANGLE_ERRF = 1E-6;
        if (secondaryAngle > PI_VAL<T>() / 2 - ANGLE_ERRF)
            secondaryAngle = PI_VAL<T>() / 2 - ANGLE_ERRF;
        if (secondaryAngle < -PI_VAL<T>() / 2 + ANGLE_ERRF)
            secondaryAngle = -PI_VAL<T>() / 2 + ANGLE_ERRF;
        while (primaryAngle > PI_VAL<T>())
            primaryAngle -= PI_VAL<T>();
        while (primaryAngle < -PI_VAL<T>())
            primaryAngle += PI_VAL<T>();

        std::array<T, 6> cos = zeroDirectionCosines();
        std::array<T, 3> z = { .0, .0, 1.0 };
        vectormath::rotate(cos.data(), z.data(), primaryAngle);
        vectormath::rotate(&cos[3], z.data(), primaryAngle);
        std::array<T, 3> x = { 1.0, .0, .0 };
        vectormath::rotate(cos.data(), x.data(), -secondaryAngle);
        vectormath::rotate(&cos[3], x.data(), -secondaryAngle);
        // handling tube rotation
        std::array<T, 3> beam_dir;
        vectormath::cross(cos.data(), beam_dir.data());
        vectormath::rotate(cos.data(), beam_dir.data(), m_tubeRotationAngle);
        vectormath::rotate(&cos[3], beam_dir.data(), m_tubeRotationAngle);

        this->setDirectionCosines(cos);
    }
    void setSourceAngles(const std::array<T, 2>& angles)
    {
        setSourceAngles(angles[0], angles[1]);
    }
    std::array<T, 2> sourceAngles() const
    {
        constexpr T ANGLE_ERRF = 1E-6;
        auto cos = this->directionCosines();
        std::array<T, 3> beam_direction;
        vectormath::cross(cos.data(), beam_direction.data());
        vectormath::rotate(cos.data(), beam_direction.data(), -m_tubeRotationAngle);
        vectormath::rotate(&cos[3], beam_direction.data(), -m_tubeRotationAngle);
        vectormath::cross(cos.data(), beam_direction.data());

        const T xy_lenght = std::sqrt(beam_direction[0] * beam_direction[0] + beam_direction[1] * beam_direction[1]);
        if (std::abs(xy_lenght) < ANGLE_ERRF) {
            if (beam_direction[2] > 0) {
                std::array<T, 2> angles = { 0, -PI_VAL<T>() / 2 };
                return angles;
            } else {
                std::array<T, 2> angles = { 0, PI_VAL<T>() / 2 };
                return angles;
            }
        }

        const T primAng = std::asin(-beam_direction[0]);

        const T zy_lenght = std::sqrt(beam_direction[2] * beam_direction[2] + beam_direction[1] * beam_direction[1]);
        if (std::abs(zy_lenght) < ANGLE_ERRF) {
            if (beam_direction[2] > 0) {
                std::array<T, 2> angles = { primAng, 0 };
                return angles;
            } else {
                std::array<T, 2> angles = { primAng, 0 };
                return angles;
            }
        }

        const T secAng = -std::asin(beam_direction[2] / zy_lenght);

        const std::array<T, 2> angles = { primAng, secAng };

        return angles;
    }
    void setSourceAnglesDeg(T primaryAngle, T secondaryAngle)
    {
        setSourceAngles(primaryAngle * DEG_TO_RAD<T>(), secondaryAngle * DEG_TO_RAD<T>());
    }

    void setSourceAnglesDeg(const std::array<T, 2>& angles)
    {
        setSourceAnglesDeg(angles[0], angles[1]);
    }

    std::array<T, 2> sourceAnglesDeg() const
    {
        auto angles = sourceAngles();
        angles[0] *= RAD_TO_DEG<T>();
        angles[1] *= RAD_TO_DEG<T>();
        return angles;
    }

    void setTubeRotation(T angle)
    {
        const T diff_angle = angle - m_tubeRotationAngle;
        std::array<T, 3> beam_direction;
        auto cos = this->directionCosines();
        vectormath::cross(cos.data(), beam_direction.data());
        vectormath::rotate(cos.data(), beam_direction.data(), diff_angle);
        vectormath::rotate(&cos[3], beam_direction.data(), diff_angle);
        this->setDirectionCosines(cos);
        m_tubeRotationAngle = angle;
    }
    T tubeRotation() const { return m_tubeRotationAngle; }
    void setTubeRotationDeg(T angle) { setTubeRotation(angle * DEG_TO_RAD<T>()); }
    T tubeRotationDeg() const { return tubeRotation() * RAD_TO_DEG<T>(); }

    void setDap(T Gycm2)
    {
        if (Gycm2 > 0.0)
            m_dap = Gycm2;
    }
    T dap(void) const { return m_dap; } // Gycm2

    T getCalibrationValue(LOWENERGYCORRECTION model, ProgressBar<T>* = nullptr) const override
    {
        auto specter = tube().getSpecter(true);
        std::vector<T> massAbsorb(specter.size(), T { 0 });
        Material airMaterial("Air, Dry (near sea level)");
        std::transform(specter.begin(), specter.end(), massAbsorb.begin(), [&](const auto& el) -> T { return airMaterial.getMassEnergyAbsorbtion(el.first); });

        T calcOutput = 0.0;
        for (std::size_t i = 0; i < specter.size(); ++i) {
            const auto& [keV, weight] = specter[i];
            calcOutput += keV * weight * massAbsorb[i]; // kev/g
        }

        calcOutput *= this->totalExposures() * this->historiesPerExposure();

        const T factor = m_dap / calcOutput;
        return factor;
    }

    bool isValid(void) const override { return m_specterValid; }
    bool validate(void) override
    {
        updateSpecterDistribution();
        return m_specterValid;
    }

    void setModelHeelEffect(bool on) { m_modelHeelEffect = on; }
    bool modelHeelEffect() const { return m_modelHeelEffect; }

protected:
    void updateFieldSize(const std::array<T, 2>& fieldSize)
    {
        for (std::size_t i = 0; i < 2; ++i) {
            m_fieldSize[i] = fieldSize[i];
            m_collimationAngles[i] = std::atan(m_fieldSize[i] * T { 0.5 } / m_sdd) * T { 2 };
        }
        m_specterValid = false;
    }

    void updateCollimationAngles(const std::array<T, 2>& collimationAngles)
    {
        for (std::size_t i = 0; i < 2; ++i) {
            m_collimationAngles[i] = collimationAngles[i];
            m_fieldSize[i] = std::tan(m_collimationAngles[i] / 2) * m_sdd * 2;
        }
        m_specterValid = false;
    }

    void updateSpecterDistribution()
    {
        if (!m_specterValid) {
            auto energies = m_tube.getEnergy();
            auto n_obs = m_tube.getSpecter(energies);
            m_specterDistribution = std::make_unique<SpecterDistribution<T>>(n_obs, energies);
            if (m_modelHeelEffect)
                m_heelFilter = std::make_unique<HeelFilter<T>>(m_tube, m_collimationAngles[1]);
            else
                m_heelFilter = nullptr;
            m_specterValid = true;
        }
    }
    std::array<T, 6> zeroDirectionCosines() const
    {
        // changing this will break setSourceAngles
        constexpr std::array<T, 6> cos = { -1.0, .0, .0, .0, .0, 1.0 };
        return cos;
    }

protected:
    T m_sdd = 1000.0;
    T m_dap = 1.0; // Gycm2
    std::array<T, 2> m_fieldSize;
    std::array<T, 2> m_collimationAngles;

    Tube<T> m_tube;
    T m_tubeRotationAngle = 0.0;
    std::shared_ptr<SpecterDistribution<T>> m_specterDistribution = nullptr;
    std::shared_ptr<HeelFilter<T>> m_heelFilter = nullptr;
    bool m_modelHeelEffect = true;
    bool m_specterValid = false;
};

template <Floating T = double>
class DXSource final : public DAPSource<T> {
public:
    DXSource()
        : DAPSource<T>()
    {
        this->m_type = Source<T>::Type::DX;
    }
    virtual ~DXSource() = default;
    Exposure<T> getExposure(std::uint64_t i) const override
    {
        constexpr T weight { 1 };
        Exposure<T> exposure(this->tubePosition(),
            this->m_directionCosines,
            this->m_collimationAngles,
            this->m_historiesPerExposure,
            weight,
            this->m_specterDistribution.get(),
            this->m_heelFilter.get());
        return exposure;
    }
    std::uint64_t totalExposures(void) const override
    {
        return m_totalExposures;
    }
    void setTotalExposures(std::uint64_t exposures)
    {
        m_totalExposures = std::max(exposures, std::uint64_t { 1 });
    }

    const std::array<T, 3> tubePosition(void) const
    {
        std::array<T, 3> beamDirection;
        vectormath::cross(this->m_directionCosines.data(), beamDirection.data());
        std::array<T, 3> pos;
        for (std::size_t i = 0; i < 3; ++i) {
            pos[i] = this->m_position[i] - beamDirection[i] * this->m_sdd;
        }
        return pos;
    }

private:
    std::uint64_t m_totalExposures = 1000;
};

template <Floating T = double>
class CBCTSource final : public DAPSource<T> {
public:
    CBCTSource()
        : DAPSource<T>()
    {
        this->m_type = Source<T>::Type::CBCT;
        this->setSourceDetectorDistance(500.0);
    }
    virtual ~CBCTSource() = default;

    const std::array<T, 3> rotationAxis() const
    {

        const auto& dc = this->m_directionCosines;
        std::array<T, 3> rot = { dc[3], dc[4], dc[5] };
        return rot;
    }

    void setSpanAngle(const T spanAngle)
    {
        m_angleSpan = std::max(spanAngle, m_angleStep);
        const auto steps = m_angleSpan / m_angleStep;
        const auto n_steps = static_cast<std::size_t>(steps);
        m_totalExposures = std::max(n_steps, std::size_t { 2 });
    }
    void setSpanAngleDeg(const T spanAngle)
    {
        setSpanAngle(spanAngle * DEG_TO_RAD<T>());
    }
    const T spanAngle() const { return m_angleSpan; }
    const T spanAngleDeg() const
    {
        return m_angleSpan * RAD_TO_DEG<T>();
    }

    void setStepAngle(const T stepAngle)
    {
        constexpr T minStep = PI_VAL<T>() / T { 360 };
        m_angleStep = std::max(stepAngle, minStep);
        const auto steps = m_angleSpan / m_angleStep;
        const auto n_steps = static_cast<std::size_t>(steps);
        m_totalExposures = std::max(n_steps, std::size_t { 2 });
    }
    void setStepAngleDeg(const T stepAngle)
    {
        setStepAngle(stepAngle * DEG_TO_RAD<T>());
    }
    const T stepAngle() const { return m_angleStep; }
    const T stepAngleDeg() const
    {
        return m_angleStep * RAD_TO_DEG<T>();
    }

    Exposure<T> getExposure(std::uint64_t i) const override
    {
        constexpr T weight { 1 };
        const auto rotAngle = i * m_angleStep;

        const auto tubePos = this->tubePosition();
        const auto& dPos = this->position();
        std::array<T, 3> normIso;
        for (std::size_t i = 0; i < 3; ++i) {
            normIso[i] = (tubePos[i] - dPos[i]);
        }
        const auto rotAxis = rotationAxis();
        vectormath::rotate(normIso.data(), rotAxis.data(), rotAngle);
        for (std::size_t i = 0; i < 3; ++i) {
            normIso[i] += dPos[i];
        }
        auto dirCosines = this->m_directionCosines;
        vectormath::rotate(dirCosines.data(), rotAxis.data(), rotAngle);
        vectormath::rotate(&dirCosines[3], rotAxis.data(), rotAngle);

        Exposure<T> exposure(normIso,
            dirCosines,
            this->m_collimationAngles,
            this->m_historiesPerExposure,
            weight,
            this->m_specterDistribution.get(),
            this->m_heelFilter.get());
        return exposure;
    }
    std::uint64_t totalExposures(void) const override
    {
        return m_totalExposures;
    }

    const std::array<T, 3> tubePosition(void) const
    {
        std::array<T, 3> beamDirection;
        vectormath::cross(this->m_directionCosines.data(), beamDirection.data());
        std::array<T, 3> pos;
        for (std::size_t i = 0; i < 3; ++i) {
            pos[i] = this->m_position[i] - beamDirection[i] * this->m_sdd * T { 0.5 };
        }
        return pos;
    }

private:
    std::size_t m_totalExposures = 180;
    T m_angleSpan = PI_VAL<T>();
    T m_angleStep = PI_VAL<T>() / T { 180 };
};

// forward decl
template <Floating T>
class CTSpiralSource;
template <Floating T>
class CTSpiralDualSource;
template <Floating T>
class CTAxialDualSource;
template <Floating T>
class CTAxialSource;

template <Floating T>
class CTBaseSource : public Source<T> {
public:
    CTBaseSource()
        : Source<T>()
    {
        this->m_type = Source<T>::Type::None;
        m_sdd = 1190.0;
        m_collimation = 38.4;
        m_fov = 500.0;
        m_startAngle = 0.0;
        m_scanLenght = 100.0;
        auto& t = tube();
        t.setAlFiltration(7.0);
        const std::array<T, 6> ct_cosines { -1, 0, 0, 0, 0, 1 };
        this->setDirectionCosines(ct_cosines);
    }
    virtual ~CTBaseSource() = default;
    Tube<T>& tube(void)
    {
        m_specterValid = false;
        return m_tube;
    }
    const Tube<T>& tube(void) const { return m_tube; }

    virtual T maxPhotonEnergyProduced() const override { return m_tube.voltage(); }

    void setBowTieFilter(std::shared_ptr<BowTieFilter<T>> filter) { m_bowTieFilter = filter; }
    std::shared_ptr<BowTieFilter<T>> bowTieFilter(void) { return m_bowTieFilter; }
    const std::shared_ptr<BowTieFilter<T>> bowTieFilter() const { return m_bowTieFilter; }

    void setSourceDetectorDistance(T sdd)
    {
        m_sdd = std::abs(sdd);
        m_specterValid = false;
    }
    T sourceDetectorDistance(void) const { return m_sdd; }

    void setCollimation(T collimation)
    {
        m_collimation = std::abs(collimation);
        m_specterValid = false;
    }
    T collimation(void) const { return m_collimation; }

    void setFieldOfView(T fov)
    {
        m_fov = std::abs(fov);
    }
    T fieldOfView(void) const
    {
        return m_fov;
    }

    void setGantryTiltAngle(T angle)
    {
        if (angle < -PI_VAL<T>())
            m_gantryTiltAngle = -PI_VAL<T>();
        else if (angle > PI_VAL<T>())
            m_gantryTiltAngle = PI_VAL<T>();
        else
            m_gantryTiltAngle = angle;
    }
    T gantryTiltAngle() const
    {
        return m_gantryTiltAngle;
    }

    void setGantryTiltAngleDeg(T angle)
    {
        setGantryTiltAngle(angle * DEG_TO_RAD<T>());
    }

    T gantryTiltAngleDeg() const
    {
        return m_gantryTiltAngle * RAD_TO_DEG<T>();
    }

    void setStartAngle(T angle)
    {
        m_startAngle = angle;
    }

    T startAngle(void) const
    {
        return m_startAngle;
    }
    T startAngleDeg(void) const
    {
        return RAD_TO_DEG<T>() * m_startAngle;
    }
    void setStartAngleDeg(T angle)
    {
        m_startAngle = DEG_TO_RAD<T>() * angle;
    }
    virtual void setScanLenght(T scanLenght)
    {
        m_scanLenght = std::abs(scanLenght);
    }
    T scanLenght(void) const
    {
        return m_scanLenght;
    }

    void setCtdiVol(T ctdivol)
    {
        if (ctdivol > 0.0)
            m_ctdivol = ctdivol;
    }
    T ctdiVol(void) const { return m_ctdivol; }

    virtual std::uint64_t totalExposures(void) const override = 0;

    void setCtdiPhantomDiameter(std::uint64_t mm)
    {
        constexpr std::uint64_t minValue = 160;
        m_ctdiPhantomDiameter = std::max(mm, minValue);
    }
    std::uint64_t ctdiPhantomDiameter(void) const { return m_ctdiPhantomDiameter; }

    void setModelHeelEffect(bool on) { m_modelHeelEffect = on; }
    bool modelHeelEffect() const { return m_modelHeelEffect; }

    bool isValid(void) const override { return m_specterValid; }
    virtual bool validate(void) override
    {
        updateSpecterDistribution();
        return m_specterValid;
    }

protected:
    template <typename U>
        requires std::is_same<CTAxialSource<T>, U>::value || std::is_same<CTAxialDualSource<T>, U>::value
    static T ctCalibration(U& sourceCopy, LOWENERGYCORRECTION model, ProgressBar<T>* progressBar = nullptr)
    {
        T meanWeight = 0;
        for (std::uint64_t i = 0; i < sourceCopy.totalExposures(); ++i) {
            auto dummy = sourceCopy.getExposure(i);
            meanWeight += dummy.beamIntensityWeight();
        }
        meanWeight /= sourceCopy.totalExposures();

        const std::array<T, 6> cosines({ -1, 0, 0, 0, 0, 1 });
        sourceCopy.setDirectionCosines(cosines);
        const std::array<T, 3> pos = { 0, 0, 0 };
        sourceCopy.setPosition(pos);
        sourceCopy.setScanLenght(sourceCopy.collimation());

        sourceCopy.setUseXCareFilter(false); // we need to disable organ aec for ctdi statistics, this should be ok
        std::size_t statCounter = CTDIPhantom<T>::ctdiMinHistories() / (sourceCopy.exposuresPerRotatition() * sourceCopy.historiesPerExposure());
        if (statCounter < 1) {
            statCounter = 1;
        }

        CTDIPhantom<T> world(sourceCopy.ctdiPhantomDiameter());

        sourceCopy.updateFromWorld(world);

        const auto histories = sourceCopy.historiesPerExposure();
        sourceCopy.setHistoriesPerExposure(histories * statCounter); // ensuring enough histories for ctdi measurement
        sourceCopy.validate();

        if (progressBar) {
            progressBar->setPlaneNormal(ProgressBar<T>::Axis::Z);
            progressBar->setPrefixMessage("CTDI calibration ");
        }

        Transport<T> transport;
        transport.setLowEnergyCorrectionModel(model);
        auto result = transport(world, &sourceCopy, progressBar, false);

        using holePosition = typename CTDIPhantom<T>::HolePosition;
        std::array<holePosition, 5> position = { holePosition::Center, holePosition::West, holePosition::East, holePosition::South, holePosition::North };

        const auto spacing = world.spacing();
        const auto voxelVolume = (spacing[0] * spacing[1] * spacing[2]) / 1000; // cm3
        const auto voxelMass = world.airDensity() * voxelVolume / 1000; // kg

        std::array<T, 5> measureDose;
        measureDose.fill(T { 0 });
        for (std::size_t i = 0; i < 5; ++i) {
            const auto& holeIndices = world.holeIndices(position[i]);
            for (const auto& idx : holeIndices) {
                measureDose[i] += result.dose[idx];
            }
            measureDose[i] /= holeIndices.size(); // n_hist * keV/kg
        }

        const T ctdiPher = (measureDose[1] + measureDose[2] + measureDose[3] + measureDose[4]) / T { 4 };
        const T ctdiCent = measureDose[0];
        T ctdiw = (ctdiCent + 2 * ctdiPher) / 3;
        ctdiw *= T { 100 } / sourceCopy.collimation();
        const T factor = sourceCopy.ctdiVol() / ctdiw / meanWeight;
        return factor;
    }
    virtual void updateSpecterDistribution()
    {
        if (!m_specterValid) {
            auto energies = m_tube.getEnergy();
            auto n_obs = m_tube.getSpecter(energies);
            m_specterDistribution = std::make_shared<SpecterDistribution<T>>(n_obs, energies);
            const T heel_span_angle = std::atan(m_collimation * T { 0.5 } / m_sdd) * T { 2.0 };
            if (m_modelHeelEffect)
                m_heelFilter = std::make_shared<HeelFilter<T>>(m_tube, heel_span_angle);
            else
                m_heelFilter = nullptr;
            m_specterValid = true;
        }
    }

protected:
    T m_sdd;
    T m_collimation;
    T m_fov;
    T m_startAngle;
    T m_scanLenght;
    T m_ctdivol = 1;
    T m_gantryTiltAngle = 0;
    std::uint64_t m_ctdiPhantomDiameter = 320;
    std::shared_ptr<BowTieFilter<T>> m_bowTieFilter = nullptr;
    Tube<T> m_tube;
    std::shared_ptr<SpecterDistribution<T>> m_specterDistribution = nullptr;
    std::shared_ptr<HeelFilter<T>> m_heelFilter = nullptr;
    bool m_modelHeelEffect = true;
    bool m_specterValid = false;
};

template <Floating T>
class CTSource : public CTBaseSource<T> {
public:
    CTSource()
        : CTBaseSource<T>()
    {
        m_exposureAngleStep = DEG_TO_RAD<T>();
    }
    virtual ~CTSource() = default;
    virtual Exposure<T> getExposure(std::uint64_t i) const override = 0;

    void setExposureAngleStep(T angleStep)
    {
        const auto absAngle = std::abs(angleStep);
        m_exposureAngleStep = std::clamp(absAngle, DEG_TO_RAD<T>() / 10, PI_VAL<T>() / 2);
    }

    T exposureAngleStep(void) const
    {
        return m_exposureAngleStep;
    }

    void setExposureAngleStepDeg(T angleStep)
    {
        setExposureAngleStep(angleStep * DEG_TO_RAD<T>());
    }

    T exposureAngleStepDeg(void) const
    {
        return m_exposureAngleStep * RAD_TO_DEG<T>();
    }

    void setAecFilter(std::shared_ptr<AECFilter<T>> filter) { m_aecFilter = filter; }
    std::shared_ptr<AECFilter<T>> aecFilter(void) { return m_aecFilter; }

    bool useXCareFilter() const { return m_useXCareFilter; }
    void setUseXCareFilter(bool use) { m_useXCareFilter = use; }
    XCareFilter<T>& xcareFilter() { return m_xcareFilter; }
    const XCareFilter<T>& xcareFilter() const { return m_xcareFilter; }

    virtual void updateFromWorld(const World<T>& world) override
    {
        if (m_aecFilter)
            m_aecFilter->updateFromWorld(world);
    }
    virtual std::uint64_t exposuresPerRotatition() const
    {
        constexpr T pi_2 = 2 * PI_VAL<T>();
        return static_cast<std::size_t>(pi_2 / m_exposureAngleStep);
    }

    T m_exposureAngleStep = RAD_TO_DEG<T>();
    std::shared_ptr<AECFilter<T>> m_aecFilter = nullptr;
    XCareFilter<T> m_xcareFilter;
    bool m_useXCareFilter = false;
};

template <Floating T = double>
class CTDualSource : public CTSource<T> {
public:
    CTDualSource()
        : CTSource<T>()
    {
        this->m_type = Source<T>::Type::None;
        m_sddB = this->m_sdd;
        m_fovB = this->m_fov;
        m_startAngleB = this->m_startAngle + PI_VAL<T>() * T { 0.5 };
        m_tubeB.setAlFiltration(this->m_tube.AlFiltration());
    }
    virtual ~CTDualSource() = default;

    T tubeAmas() const { return m_tubeAmas; }
    T tubeBmas() const { return m_tubeBmas; }
    void setTubeAmas(T mas)
    {
        this->m_specterValid = false;
        m_tubeAmas = std::max(T { 0.0 }, mas);
    }
    void setTubeBmas(T mas)
    {
        this->m_specterValid = false;
        m_tubeBmas = std::max(T { 0.0 }, mas);
    }

    Tube<T>& tubeB(void)
    {
        this->m_specterValid = false;
        return m_tubeB;
    }
    const Tube<T>& tubeB(void) const { return m_tubeB; }

    T maxPhotonEnergyProduced() const override { return std::max(this->m_tube.voltage(), m_tubeB.voltage()); }
    std::uint64_t exposuresPerRotatition() const override
    {
        return 2 * static_cast<std::size_t>((2 * PI_VAL<T>()) / this->m_exposureAngleStep);
    }
    void setBowTieFilterB(std::shared_ptr<BowTieFilter<T>> filter) { m_bowTieFilterB = filter; }
    std::shared_ptr<BowTieFilter<T>> bowTieFilterB(void) { return m_bowTieFilterB; }
    const std::shared_ptr<BowTieFilter<T>> bowTieFilterB() const { return m_bowTieFilterB; }

    void setSourceDetectorDistanceB(T sdd)
    {
        this->m_specterValid = false;
        m_sddB = std::abs(sdd);
    }
    T sourceDetectorDistanceB(void) const { return m_sddB; }

    void setFieldOfViewB(T fov) { m_fovB = std::abs(fov); }
    T fieldOfViewB(void) const { return m_fovB; }
    void setStartAngleB(T angle) { m_startAngleB = angle; }
    T startAngleB(void) const { return m_startAngleB; }
    void setStartAngleDegB(T angle)
    {
        m_startAngleB = DEG_TO_RAD<T>() * angle;
    }

    T startAngleDegB(void) const
    {
        return RAD_TO_DEG<T>() * m_startAngleB;
    }

    bool validate(void) override
    {
        updateSpecterDistribution();
        return this->m_specterValid;
    }

protected:
    void updateSpecterDistribution() override
    {
        if (!this->m_specterValid) {
            const auto energyA = this->m_tube.getEnergy();
            const auto energyB = m_tubeB.getEnergy();
            auto specterA = this->m_tube.getSpecter(energyA, false);
            auto specterB = m_tubeB.getSpecter(energyB, false);

            const auto sumA = std::accumulate(specterA.cbegin(), specterA.cend(), T { 0.0 });
            const auto sumB = std::accumulate(specterB.cbegin(), specterB.cend(), T { 0.0 });
            const auto weightA = m_tubeAmas * sumA;
            const auto weightB = m_tubeBmas * sumB;

            std::transform(std::execution::par_unseq, specterA.cbegin(), specterA.cend(), specterA.begin(), [=](auto i) { return i / sumA; });
            std::transform(std::execution::par_unseq, specterB.cbegin(), specterB.cend(), specterB.begin(), [=](auto i) { return i / sumB; });

            m_tubeAweight = weightA * T { 2 } / (weightA + weightB);
            m_tubeBweight = weightB * T { 2 } / (weightA + weightB);

            this->m_specterDistribution = std::make_shared<SpecterDistribution<T>>(specterA, energyA);
            m_specterDistributionB = std::make_shared<SpecterDistribution<T>>(specterB, energyB);

            const auto heel_span_angle = std::atan(this->m_collimation * T { 0.5 } / this->m_sdd) * T { 2 };
            this->m_heelFilter = std::make_shared<HeelFilter<T>>(this->m_tube, heel_span_angle);
            m_heelFilterB = std::make_shared<HeelFilter<T>>(m_tubeB, heel_span_angle);

            this->m_specterValid = true;
        }
    }

    T m_sddB;
    T m_fovB;
    T m_startAngleB;
    T m_tubeAmas = 100.0;
    T m_tubeBmas = 100.0;
    T m_tubeBweight = -1.0;
    T m_tubeAweight = -1.0;
    std::shared_ptr<BowTieFilter<T>> m_bowTieFilterB = nullptr;
    Tube<T> m_tubeB;
    std::shared_ptr<SpecterDistribution<T>> m_specterDistributionB = nullptr;
    std::shared_ptr<HeelFilter<T>> m_heelFilterB = nullptr;
};

template <Floating T = double>
class CTAxialSource final : public CTSource<T> {
public:
    CTAxialSource()
        : CTSource<T>()
    {
        this->m_type = Source<T>::Type::CTAxial;
        m_step = this->m_collimation;
        this->m_scanLenght = m_step;
    }

    CTAxialSource(const CTSpiralSource<T>& other);
    virtual ~CTAxialSource() = default;

    Exposure<T> getExposure(std::uint64_t exposureIndex) const override
    {
        // calculating position
        std::array<T, 3> pos = { 0, -this->m_sdd / T { 2 }, 0 };

        constexpr T PI_2 = 2 * PI_VAL<T>();

        const std::uint64_t anglesPerRotation = static_cast<std::uint64_t>(PI_2 / this->m_exposureAngleStep);
        const std::uint64_t rotationNumber = exposureIndex / anglesPerRotation;

        const auto angle = this->m_startAngle + this->m_exposureAngleStep * (exposureIndex - (rotationNumber * anglesPerRotation));

        auto directionCosines = this->m_directionCosines;
        T* rotationAxis = &directionCosines[3];
        T* otherAxis = &directionCosines[0];
        std::array<T, 3> tiltAxis = { 1, 0, 0 };
        auto tiltCorrection = pos;
        vectormath::rotate(tiltCorrection.data(), tiltAxis.data(), this->m_gantryTiltAngle);
        vectormath::rotate(rotationAxis, tiltAxis.data(), this->m_gantryTiltAngle);
        vectormath::rotate(otherAxis, tiltAxis.data(), this->m_gantryTiltAngle);
        vectormath::rotate(pos.data(), rotationAxis, angle);
        pos[2] += m_step * rotationNumber + tiltCorrection[2];

        vectormath::rotate(otherAxis, rotationAxis, angle);
        for (std::size_t i = 0; i < 3; ++i) {
            pos[i] += this->m_position[i];
        }
        const std::array<T, 2> collimationAngles = {
            std::atan(this->m_fov / this->m_sdd) * T { 2 },
            std::atan(this->m_collimation / this->m_sdd) * T { 2 }
        };

        T weight { 1 };
        if (this->m_aecFilter)
            weight *= this->m_aecFilter->sampleIntensityWeight(pos);
        if (this->m_useXCareFilter)
            weight *= this->m_xcareFilter.sampleIntensityWeight(angle);

        Exposure<T> exposure(pos,
            directionCosines,
            collimationAngles,
            this->m_historiesPerExposure,
            weight,
            this->m_specterDistribution.get(),
            this->m_heelFilter.get(),
            this->m_bowTieFilter.get());

        return exposure;
    }

    void setStep(T step)
    {
        const auto absStep = std::abs(step);
        const auto n_steps = this->m_scanLenght / m_step;
        m_step = absStep > 0.01 ? absStep : 0.01;
        setScanLenght(m_step * n_steps);
    }

    T step(void) const { return m_step; }
    void setScanLenght(T scanLenght) override
    {
        this->m_scanLenght = std::max(m_step * std::ceil(std::abs(scanLenght) / m_step), m_step);
    }

    std::uint64_t totalExposures(void) const override
    {
        const std::uint64_t anglesPerRotation = static_cast<std::uint64_t>(PI_VAL<T>() * T { 2 } / this->m_exposureAngleStep);
        const std::uint64_t rotationNumbers = static_cast<std::uint64_t>(std::round(this->m_scanLenght / m_step));
        return anglesPerRotation * rotationNumbers;
    }

    T getCalibrationValue(LOWENERGYCORRECTION model, ProgressBar<T>* progressBar = nullptr) const override
    {
        auto copy = *this;
        return CTSource<T>::ctCalibration(copy, model, progressBar);
    }

private:
    T m_step;
};
template <Floating T = double>
class CTSpiralSource final : public CTSource<T> {
public:
    CTSpiralSource()
        : CTSource<T>()
    {
        this->m_type = Source<T>::Type::CTSpiral;
        m_pitch = 1.0;
    }
    virtual ~CTSpiralSource() = default;
    Exposure<T> getExposure(std::uint64_t exposureIndex) const override
    {
        std::array<T, 3> pos = { 0, -this->m_sdd / 2, 0 };

        const auto angle = this->m_startAngle + this->m_exposureAngleStep * exposureIndex;

        auto directionCosines = this->m_directionCosines;
        T* rotationAxis = &directionCosines[3];
        T* otherAxis = &directionCosines[0];
        std::array<T, 3> tiltAxis = { 1, 0, 0 };
        auto tiltCorrection = pos;
        vectormath::rotate(tiltCorrection.data(), tiltAxis.data(), this->m_gantryTiltAngle);
        vectormath::rotate(rotationAxis, tiltAxis.data(), this->m_gantryTiltAngle);
        vectormath::rotate(otherAxis, tiltAxis.data(), this->m_gantryTiltAngle);
        vectormath::rotate(pos.data(), rotationAxis, angle);
        constexpr T PI_2 = T { 2 } * PI_VAL<T>();
        pos[2] += (exposureIndex * this->m_exposureAngleStep) * this->m_collimation * m_pitch / PI_2 + tiltCorrection[2];

        vectormath::rotate(otherAxis, rotationAxis, angle);
        for (std::size_t i = 0; i < 3; ++i) {
            pos[i] += this->m_position[i];
        }

        const std::array<T, 2> collimationAngles = {
            std::atan(this->m_fov / this->m_sdd) * T { 2.0 },
            std::atan(this->m_collimation / this->m_sdd) * T { 2.0 }
        };

        T weight { 1.0 };
        if (this->m_aecFilter)
            weight *= this->m_aecFilter->sampleIntensityWeight(pos);
        if (this->m_useXCareFilter)
            weight *= this->m_xcareFilter.sampleIntensityWeight(angle);

        Exposure<T> exposure(pos,
            directionCosines,
            collimationAngles,
            this->m_historiesPerExposure,
            weight,
            this->m_specterDistribution.get(),
            this->m_heelFilter.get(),
            this->m_bowTieFilter.get());

        return exposure;
    }

    void setPitch(T pitch)
    {
        m_pitch = std::max(T { 0.01 }, pitch);
    }

    T pitch(void) const { return m_pitch; }
    void setScanLenght(T scanLenght) override
    {
        this->m_scanLenght = std::max(std::abs(scanLenght), this->m_collimation * m_pitch * T { 0.5 });
    }
    std::uint64_t totalExposures(void) const override
    {
        constexpr T PI_2 = 2 * PI_VAL<T>();
        return static_cast<std::uint64_t>(this->m_scanLenght * PI_2 / (this->m_collimation * m_pitch * this->m_exposureAngleStep));
    }

    T getCalibrationValue(LOWENERGYCORRECTION model, ProgressBar<T>* progressBar = nullptr) const override
    {
        CTAxialSource<T> copy = *this;

        return this->ctCalibration(copy, model, progressBar) * m_pitch;
    }

private:
    T m_pitch;
};

template <Floating T = double>
class CTAxialDualSource final : public CTDualSource<T> {
public:
    CTAxialDualSource()
        : CTDualSource<T>()
    {
        this->m_type = Source<T>::Type::CTDual;
        m_step = this->m_collimation;
        this->m_scanLenght = m_step;
    }

    CTAxialDualSource(const CTSpiralDualSource<T>& other);

    virtual ~CTAxialDualSource() = default;

    Exposure<T> getExposure(std::uint64_t exposureIndexTotal) const override
    {
        T sdd, startAngle, fov;
        uint64_t exposureIndex = exposureIndexTotal / 2;
        BeamFilter<T>* bowTieFilter = nullptr;
        SpecterDistribution<T>* specterDistribution = nullptr;
        HeelFilter<T>* heelFilter = nullptr;
        T weight { 1 };
        if (exposureIndexTotal % 2 == 0) {
            sdd = this->m_sdd;
            startAngle = this->m_startAngle;
            fov = this->m_fov;
            bowTieFilter = this->m_bowTieFilter.get();
            specterDistribution = this->m_specterDistribution.get();
            heelFilter = this->m_heelFilter.get();
            weight = this->m_tubeAweight;
        } else {
            sdd = this->m_sddB;
            startAngle = this->m_startAngleB;
            fov = this->m_fovB;
            bowTieFilter = this->m_bowTieFilterB.get();
            specterDistribution = this->m_specterDistributionB.get();
            heelFilter = this->m_heelFilterB.get();
            weight = this->m_tubeBweight;
        }
        // calculating position
        std::array<T, 3> pos = { 0, -this->m_sdd / T { 2.0 }, 0 };

        constexpr T PI_2 = 2 * PI_VAL<T>();

        const std::uint64_t anglesPerRotation = static_cast<std::uint64_t>(PI_2 / this->m_exposureAngleStep);
        const std::uint64_t rotationNumber = exposureIndex / anglesPerRotation;

        const auto angle = startAngle + this->m_exposureAngleStep * (exposureIndex - (rotationNumber * anglesPerRotation));

        auto directionCosines = this->m_directionCosines;
        T* rotationAxis = &directionCosines[3];
        T* otherAxis = &directionCosines[0];
        std::array<T, 3> tiltAxis = { 1, 0, 0 };
        auto tiltCorrection = pos;
        vectormath::rotate(tiltCorrection.data(), tiltAxis.data(), this->m_gantryTiltAngle);
        vectormath::rotate(rotationAxis, tiltAxis.data(), this->m_gantryTiltAngle);
        vectormath::rotate(otherAxis, tiltAxis.data(), this->m_gantryTiltAngle);
        vectormath::rotate(pos.data(), rotationAxis, angle);
        pos[2] += m_step * rotationNumber + tiltCorrection[2];

        vectormath::rotate(otherAxis, rotationAxis, angle);
        for (std::size_t i = 0; i < 3; ++i) {
            pos[i] += this->m_position[i];
        }
        const std::array<T, 2> collimationAngles = {
            std::atan(fov / sdd) * T { 2.0 },
            std::atan(this->m_collimation / sdd) * T { 2.0 }
        };

        if (this->m_aecFilter)
            weight *= this->m_aecFilter->sampleIntensityWeight(pos);
        if (this->m_useXCareFilter)
            weight *= this->m_xcareFilter.sampleIntensityWeight(angle);

        Exposure<T> exposure(pos,
            directionCosines,
            collimationAngles,
            this->m_historiesPerExposure,
            weight,
            specterDistribution,
            heelFilter,
            bowTieFilter);

        return exposure;
    }

    void setStep(T step)
    {
        const auto absStep = std::abs(step);
        const auto n_steps = this->m_scanLenght / m_step;
        m_step = absStep > 0.01 ? absStep : 0.01;
        setScanLenght(m_step * n_steps);
    }

    T step(void) const { return m_step; }
    void setScanLenght(T scanLenght) override
    {
        this->m_scanLenght = std::max(m_step * std::ceil(std::abs(scanLenght) / m_step), m_step);
    }

    std::uint64_t totalExposures(void) const override
    {
        const std::uint64_t anglesPerRotation = static_cast<std::uint64_t>(PI_VAL<T>() * T { 2 } / this->m_exposureAngleStep);
        const std::uint64_t rotationNumbers = static_cast<std::uint64_t>(std::round(this->m_scanLenght / m_step));
        return anglesPerRotation * rotationNumbers * 2;
    }

    T getCalibrationValue(LOWENERGYCORRECTION model, ProgressBar<T>* progressBar = nullptr) const override
    {
        auto copy = *this;
        return this->ctCalibration(copy, model, progressBar);
    }

private:
    T m_step;
};

template <Floating T = double>
class CTSpiralDualSource final : public CTDualSource<T> {
public:
    // Source overrides
    CTSpiralDualSource()
        : CTDualSource<T>()
    {
        this->m_type = Source<T>::Type::CTDual;
        m_pitch = 1.0;
    }
    virtual ~CTSpiralDualSource() = default;
    Exposure<T> getExposure(std::uint64_t exposureIndexTotal) const override
    {
        T sdd, startAngle, fov;
        uint64_t exposureIndex = exposureIndexTotal / 2;
        BeamFilter<T>* bowTie = nullptr;
        SpecterDistribution<T>* specterDistribution = nullptr;
        HeelFilter<T>* heelFilter = nullptr;
        T weight { 1 };
        if (exposureIndexTotal % 2 == 0) {
            sdd = this->m_sdd;
            startAngle = this->m_startAngle;
            fov = this->m_fov;
            bowTie = this->m_bowTieFilter.get();
            specterDistribution = this->m_specterDistribution.get();
            heelFilter = this->m_heelFilter.get();
            weight = this->m_tubeAweight;
        } else {
            sdd = this->m_sddB;
            startAngle = this->m_startAngleB;
            fov = this->m_fovB;
            bowTie = this->m_bowTieFilterB.get();
            specterDistribution = this->m_specterDistributionB.get();
            heelFilter = this->m_heelFilterB.get();
            weight = this->m_tubeBweight;
        }

        std::array<T, 3> pos = { 0, -this->m_sdd / 2, 0 };
        const auto angle = startAngle + this->m_exposureAngleStep * exposureIndex;

        auto directionCosines = this->m_directionCosines;
        T* rotationAxis = &directionCosines[3];
        T* otherAxis = &directionCosines[0];
        std::array<T, 3> tiltAxis = { 1, 0, 0 };
        auto tiltCorrection = pos;
        vectormath::rotate(tiltCorrection.data(), tiltAxis.data(), this->m_gantryTiltAngle);
        vectormath::rotate(rotationAxis, tiltAxis.data(), this->m_gantryTiltAngle);
        vectormath::rotate(otherAxis, tiltAxis.data(), this->m_gantryTiltAngle);
        vectormath::rotate(pos.data(), rotationAxis, angle);
        constexpr T PI_2 = T { 2 } * PI_VAL<T>();
        pos[2] += (exposureIndex * this->m_exposureAngleStep) * this->m_collimation * m_pitch / PI_2 + tiltCorrection[2];

        vectormath::rotate(otherAxis, rotationAxis, angle);
        for (std::size_t i = 0; i < 3; ++i) {
            pos[i] += this->m_position[i];
        }

        const std::array<T, 2> collimationAngles = {
            std::atan(fov / sdd) * T { 2.0 },
            std::atan(this->m_collimation / sdd) * T { 2.0 }
        };

        if (this->m_aecFilter)
            weight *= this->m_aecFilter->sampleIntensityWeight(pos);
        if (this->m_useXCareFilter)
            weight *= this->m_xcareFilter.sampleIntensityWeight(angle);

        Exposure<T> exposure(pos,
            directionCosines,
            collimationAngles,
            this->m_historiesPerExposure,
            weight,
            specterDistribution,
            heelFilter,
            bowTie);

        return exposure;
    }

    std::uint64_t totalExposures(void) const override
    {
        const auto coll = this->collimation();
        const auto p = this->pitch();
        const auto angStep = this->exposureAngleStep();
        const auto scanL = this->scanLenght();
        auto singleSourceExposure = static_cast<std::uint64_t>(scanL * 2 * PI_VAL<T>() / (coll * p * angStep));
        return singleSourceExposure * 2;
    }

    T getCalibrationValue(LOWENERGYCORRECTION model, ProgressBar<T>* progressBar = nullptr) const override
    {
        CTAxialDualSource<T> copy = *this;

        return CTSource<T>::ctCalibration(copy, model, progressBar) * m_pitch;
    }

    T pitch(void) const { return m_pitch; }
    void setPitch(T pitch)
    {
        m_pitch = std::max(T { 0.01 }, pitch);
    }

    void setScanLenght(T scanLenght) override
    {
        const auto coll = this->collimation();
        const auto p = this->pitch();
        const auto scanLenght_val = std::max(std::abs(scanLenght), coll * p * T { 0.5 });
        CTBaseSource<T>::setScanLenght(scanLenght_val);
    }

private:
    T m_pitch = 1.0;
};

// Implementation that requires declared classes
template <Floating T>
CTAxialSource<T>::CTAxialSource(const CTSpiralSource<T>& other)
    : CTSource<T>(other)
{
    this->m_step = this->m_collimation;
    setScanLenght(other.scanLenght());
}

template <Floating T>
CTAxialDualSource<T>::CTAxialDualSource(const CTSpiralDualSource<T>& other)
    : CTDualSource<T>(other)
{
    m_step = this->m_collimation;
    setScanLenght(other.scanLenght());
}

template <Floating T>
class CTTopogramSource : public CTBaseSource<T> {
public:
    CTTopogramSource()
        : CTBaseSource<T>()
    {
        this->m_type = Source<T>::Type::CTTopogram;
    }
    virtual ~CTTopogramSource() = default;
    Exposure<T> getExposure(std::uint64_t i) const override
    {

        // calculating position
        std::array<T, 3> pos = { 0, -this->m_sdd / T { 2 }, 0 };

        const auto totalExposures = this->totalExposures();
        const auto angle = this->m_startAngle;
        const auto step = this->scanLenght() / (totalExposures - 1);

        auto directionCosines = this->m_directionCosines;
        T* rotationAxis = &directionCosines[3];
        T* otherAxis = &directionCosines[0];
        std::array<T, 3> tiltAxis = { 1, 0, 0 };
        auto tiltCorrection = pos;
        vectormath::rotate(tiltCorrection.data(), tiltAxis.data(), this->m_gantryTiltAngle);
        vectormath::rotate(rotationAxis, tiltAxis.data(), this->m_gantryTiltAngle);
        vectormath::rotate(otherAxis, tiltAxis.data(), this->m_gantryTiltAngle);
        vectormath::rotate(pos.data(), rotationAxis, angle);
        pos[2] += step * i + tiltCorrection[2];

        vectormath::rotate(otherAxis, rotationAxis, angle);
        for (std::size_t i = 0; i < 3; ++i) {
            pos[i] += this->m_position[i];
        }
        const std::array<T, 2> collimationAngles = {
            std::atan(this->m_fov / this->m_sdd) * 2,
            std::atan(this->m_collimation / this->m_sdd) * 2
        };

        const T weight = 1;

        Exposure<T> exposure(pos,
            directionCosines,
            collimationAngles,
            this->m_historiesPerExposure,
            weight,
            this->m_specterDistribution.get(),
            this->m_heelFilter.get(),
            this->m_bowTieFilter.get());

        return exposure;
    }

    std::uint64_t totalExposures() const
    {
        const auto n = this->scanLenght();
        const auto un = static_cast<std::uint64_t>(std::ceil(n));
        constexpr std::uint64_t min = 1;
        return std::max(un, min);
    }

    T getCalibrationValue(LOWENERGYCORRECTION model, ProgressBar<T>* progressBar = nullptr) const override
    {
        // This should be prettier
        CTAxialSource<T> copy;
        CTSource<T>& ctsource = copy;
        CTBaseSource<T>& base = ctsource;
        base = *this;

        const auto ctdivol_axial = this->ctdiVol() * this->scanLenght() / this->collimation();
        copy.setCtdiVol(ctdivol_axial);
        copy.setScanLenght(0);
        copy.setStep(this->m_collimation);

        // matching angle step to number of exposures
        constexpr auto maxStep = (2 * PI_VAL<T>()) / 72; // min step
        const auto step = std::min(2 * PI_VAL<T>() / totalExposures(), maxStep);
        copy.setExposureAngleStep(step);

        const auto totalExposures = this->totalExposures();

        const auto factor = CTSource<T>::ctCalibration(copy, model, progressBar);
        return (factor * totalExposures) / copy.totalExposures();
    }
};

}
