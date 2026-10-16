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

Copyright 2022 Erlend Andersen
*/

#include "testutils.hpp"

#include "xraymc/beams/cbctbeam.hpp"
#include "xraymc/beams/ctsequentialbeam.hpp"
#include "xraymc/beams/ctspiralbeam.hpp"
#include "xraymc/beams/ctspiraldualenergybeam.hpp"
#include "xraymc/beams/dxbeam.hpp"
#include "xraymc/beams/flatmonoenergyfield.hpp"
#include "xraymc/beams/isotropicbeam.hpp"
#include "xraymc/beams/isotropicbeamcircle.hpp"
#include "xraymc/beams/isotropiccircularbeam.hpp"
#include "xraymc/beams/isotropiccircularmonoenergybeam.hpp"
#include "xraymc/beams/isotropicmonoenergybeam.hpp"
#include "xraymc/beams/isotropicmonoenergybeamcircle.hpp"
#include "xraymc/beams/pencilbeam.hpp"
#include "xraymc/beams/tube/tube.hpp"
#include "xraymc/serializer.hpp"
#include "xraymc/world/world.hpp"
#include "xraymc/world/worlditems/depthdose.hpp"
#include "xraymc/world/worlditems/enclosedroom.hpp"
#include "xraymc/world/worlditems/worldbox.hpp"
#include "xraymc/world/worlditems/worldsphere.hpp"

#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

using namespace xraymc;
using xraymctest::check;
using xraymctest::checkNear;

void testPrimitives()
{
    auto buffer = Serializer::getEmptyBuffer();
    const double d = 3.14159;
    const std::uint64_t u = 1234567890123ULL;
    const std::uint8_t b = 200;
    const std::array<double, 3> arr = { 1, -2, 3.5 };
    const std::vector<double> vec = { 0.1, 0.2, 0.3, 0.4 };
    const std::vector<std::uint64_t> uvec = { 7, 8, 9 };
    const std::vector<std::string> strs = { "Test 1", "", "Testing 2 ", "a much longer string with spaces" };
    Serializer::serialize(d, buffer);
    Serializer::serialize(u, buffer);
    Serializer::serialize(b, buffer);
    Serializer::serialize(arr, buffer);
    Serializer::serialize(vec, buffer);
    Serializer::serialize(uvec, buffer);
    Serializer::serialize(strs, buffer);

    double d2;
    std::uint64_t u2;
    std::uint8_t b2;
    std::array<double, 3> arr2;
    std::vector<double> vec2;
    std::vector<std::uint64_t> uvec2;
    std::vector<std::string> strs2;
    auto span = std::span<const char>(buffer);
    span = Serializer::deserialize(d2, span);
    span = Serializer::deserialize(u2, span);
    span = Serializer::deserialize(b2, span);
    span = Serializer::deserialize(arr2, span);
    span = Serializer::deserialize(vec2, span);
    span = Serializer::deserialize(uvec2, span);
    span = Serializer::deserialize(strs2, span);
    check(d2 == d && u2 == u && b2 == b, "Serializer: scalar round trip");
    check(arr2 == arr, "Serializer: array round trip");
    check(vec2 == vec && uvec2 == uvec, "Serializer: vector round trip");
    check(strs2 == strs, "Serializer: string vector round trip");
    check(span.empty(), "Serializer: entire buffer consumed");

    bool threw = false;
    try {
        double x;
        Serializer::deserialize(x, std::span<const char>(buffer.data(), 3));
    } catch (const std::length_error&) {
        threw = true;
    }
    check(threw, "Serializer: deserialize from too short buffer throws length_error");
}

void testCompositeTypes()
{
    auto buffer = Serializer::getEmptyBuffer();
    const std::map<std::uint8_t, double> w = { { 1, 0.111894 }, { 8, 0.888106 } };
    const std::vector<std::map<std::uint8_t, double>> ws = { w, { { 82, 1.0 } }, {} };
    DoseScore ds;
    ds.set(1.5, 0.25, 42);
    std::vector<DoseScore> dsv(3);
    dsv[1].set(2, 3, 4);
    Serializer::serializeMaterialWeights(w, buffer);
    Serializer::serializeMaterialWeights(ws, buffer);
    Serializer::serializeDoseScore(ds, buffer);
    Serializer::serializeDoseScore(std::span<const DoseScore>(dsv), buffer);
    Serializer::serializeDoseScore(std::span<const DoseScore>(dsv), buffer);

    std::map<std::uint8_t, double> w2;
    std::vector<std::map<std::uint8_t, double>> ws2;
    DoseScore ds2;
    std::vector<DoseScore> dsv2;
    std::array<DoseScore, 3> dsa;
    auto span = std::span<const char>(buffer);
    span = Serializer::deserializeMaterialWeights(w2, span);
    span = Serializer::deserializeMaterialWeights(ws2, span);
    span = Serializer::deserializeDoseScore(ds2, span);
    span = Serializer::deserializeDoseScore(dsv2, span);
    span = Serializer::deserializeDoseScore(dsa, span);
    check(w2 == w, "Serializer: material weights round trip");
    check(ws2 == ws, "Serializer: vector of material weights round trip");
    check(ds2 == ds, "Serializer: dose score round trip");
    check(dsv2 == dsv, "Serializer: dose score vector round trip");
    check(dsa[1] == dsv[1] && dsa[0] == dsv[0], "Serializer: dose score array round trip");
    check(span.empty(), "Serializer: entire composite buffer consumed");
}

void testFiles()
{
    std::vector<std::string> ss;
    for (int i = 0; i < 1000; ++i)
        ss.push_back("Repeating test string number " + std::to_string(i % 10));
    auto buffer = Serializer::getEmptyBuffer();
    Serializer::serialize(ss, buffer);

    for (bool compress : { false, true }) {
        const auto path = xraymctest::tempPath(compress ? "compressed.xr" : "plain.xr").string();
        check(Serializer::write(path, buffer, compress), std::string("Serializer: write ") + (compress ? "compressed" : "uncompressed") + " file");
        const auto r = Serializer::read(path);
        check(r.has_value() && r.value() == buffer, std::string("Serializer: read back ") + (compress ? "compressed" : "uncompressed") + " file");
        if (compress) {
            const auto size = std::filesystem::file_size(path);
            check(size < buffer.size() / 2, "Serializer: compressed file is smaller (" + std::to_string(size) + " vs " + std::to_string(buffer.size()) + " bytes)");
        }
        std::filesystem::remove(path);
    }

    Serializer s(xraymctest::tempPath("member.xr").string());
    check(s.write(buffer), "Serializer: member write with filename");
    check(Serializer::read(xraymctest::tempPath("member.xr").string()).value() == buffer, "Serializer: member write can be read back");
    std::filesystem::remove(xraymctest::tempPath("member.xr"));
    check(!Serializer().write(buffer), "Serializer: member write without filename fails");

    const auto missing = Serializer::read(xraymctest::tempPath("does_not_exist.xr").string());
    check(!missing && missing.error() == Serializer::parse_error::buffer_file_error, "Serializer: reading missing file gives file error");

    const auto badPath = xraymctest::tempPath("bad.xr");
    {
        std::ofstream f(badPath, std::ios::binary);
        f << "this is not a valid xraymc file header at all";
    }
    const auto bad = Serializer::read(badPath.string());
    check(!bad && bad.error() == Serializer::parse_error::buffer_version_mismatch, "Serializer: wrong header gives version mismatch");
    {
        std::ofstream f(badPath, std::ios::binary);
        f << "short";
    }
    const auto shortFile = Serializer::read(badPath.string());
    check(!shortFile && shortFile.error() == Serializer::parse_error::buffer_size_short, "Serializer: short file gives size error");
    std::filesystem::remove(badPath);
}

void testItems()
{
    // Several items in one buffer can be read back sequentially by name
    WorldSphere<5, 1> sphere(4, { 1, 2, 3 });
    WorldBox<5, 1> box({ -1, -2, -3, 3, 2, 1 });
    auto buffer = Serializer::getEmptyBuffer();
    Serializer::serializeItem(sphere, buffer);
    Serializer::serializeItem(box, buffer);

    auto span = std::span<const char>(buffer);
    check(Serializer::getCurrentItemName(span) == WorldSphere<5, 1>::magicID(), "Serializer: first item name");
    std::array<char, 32> name;
    std::vector<char> item;
    span = Serializer::deserializeItem(name, item, span);
    check(WorldSphere<5, 1>::validMagicID(name) && !WorldBox<5, 1>::validMagicID(name), "Serializer: item magic ID identifies type");
    const auto s2 = WorldSphere<5, 1>::deserialize(item);
    check(s2.has_value() && s2.value() == sphere, "Serializer: sphere round trip");
    span = Serializer::deserializeItem(name, item, span);
    const auto b2 = WorldBox<5, 1>::deserialize(item);
    check(WorldBox<5, 1>::validMagicID(name) && b2.has_value() && b2.value() == box, "Serializer: box round trip");
    check(span.empty(), "Serializer: all items consumed");

    bool threw = false;
    try {
        Serializer::deserializeItem(name, item, std::span<const char>(buffer.data(), 40));
    } catch (const std::length_error&) {
        threw = true;
    }
    check(threw, "Serializer: truncated item throws length_error");

    Tube tube(95, 11);
    tube.setAlFiltration(3);
    auto tbuf = Serializer::getEmptyBuffer();
    Serializer::serializeItem(tube, tbuf);
    std::vector<char> tdata;
    Serializer::deserializeItem(name, tdata, tbuf);
    const auto t2 = Tube::deserialize(tdata);
    bool tubeSame = t2.has_value();
    if (t2) {
        // spectra are normalized with parallel reductions, allow for rounding differences
        const auto a = t2.value().getSpecter();
        const auto b = tube.getSpecter();
        tubeSame = a.size() == b.size();
        for (std::size_t i = 0; i < a.size() && tubeSame; ++i)
            tubeSame = a[i].first == b[i].first && std::abs(a[i].second - b[i].second) <= 1E-12;
    }
    check(tubeSame, "Serializer: tube round trip");
}

// Beams are equal after a round trip if they produce identical particles from the same random state
template <typename B>
void testBeam(const B& beam, const std::string& name, bool compareCalibration = true)
{
    auto buffer = Serializer::getEmptyBuffer();
    Serializer::serializeItem(beam, buffer);
    std::array<char, 32> id;
    std::vector<char> data;
    Serializer::deserializeItem(id, data, buffer);
    check(B::validMagicID(id), name + ": magic ID round trip");
    const auto copy_opt = B::deserialize(data);
    check(copy_opt.has_value(), name + ": deserialize succeeds");
    if (!copy_opt)
        return;
    const auto& copy = copy_opt.value();
    check(copy.numberOfExposures() == beam.numberOfExposures() && copy.numberOfParticles() == beam.numberOfParticles(), name + ": number of exposures and particles preserved");
    bool same = true;
    for (std::size_t e : { std::size_t { 0 }, beam.numberOfExposures() / 2, beam.numberOfExposures() - 1 }) {
        auto s1 = xraymctest::seededState(e + 1);
        auto s2 = xraymctest::seededState(e + 1);
        const auto ea = beam.exposure(e);
        const auto eb = copy.exposure(e);
        same = same && ea.position() == eb.position() && ea.numberOfParticles() == eb.numberOfParticles();
        for (int i = 0; i < 100; ++i) {
            const auto pa = ea.sampleParticle(s1);
            const auto pb = eb.sampleParticle(s2);
            // directions may be renormalized when deserialized, allow for rounding differences
            for (std::size_t k = 0; k < 3; ++k)
                same = same && std::abs(pa.pos[k] - pb.pos[k]) < 1E-12 && std::abs(pa.dir[k] - pb.dir[k]) < 1E-12;
            same = same && std::abs(pa.energy - pb.energy) < 1E-12 && std::abs(pa.weight - pb.weight) < 1E-12;
        }
    }
    check(same, name + ": deserialized beam samples identical particles");
    if (compareCalibration)
        checkNear(copy.calibrationFactor(), beam.calibrationFactor(), 1E-12 * std::abs(beam.calibrationFactor()), name + ": calibration factor preserved");
}

void testBeams()
{
    {
        PencilBeam<> b({ 1, 2, 3 }, { 0, 1, 1 }, 77);
        b.setAirKerma(3);
        b.setNumberOfExposures(7);
        testBeam(b, "PencilBeam");
    }
    {
        FlatMonoEnergyField<> b({ 0, 1, 0 }, 66);
        b.setDirectionCosines({ 0, 1, 0 }, { 0, 0, 1 });
        b.setLenght(4, 3);
        b.setAirKerma(2);
        testBeam(b, "FlatMonoEnergyField");
    }
    {
        IsotropicMonoEnergyBeam<> b({ 1, 0, 0 }, { { { 0, 1, 0 }, { 0, 0, 1 } } }, 45);
        b.setCollimationHalfAngles(-0.1, -0.2, 0.3, 0.1);
        testBeam(b, "IsotropicMonoEnergyBeam");
    }
    {
        IsotropicBeam<> b({ 1, 0, 0 });
        b.setEnergySpecter(Tube(90).getSpecter());
        b.setCollimationHalfAngles(0.1);
        testBeam(b, "IsotropicBeam");
    }
    {
        IsotropicBeamCircle<> b({ 0, 0, 1 }, { 1, 1, 0 });
        b.setEnergySpecter(Tube(80).getSpecter());
        b.setCollimationHalfAngle(0.2);
        testBeam(b, "IsotropicBeamCircle");
    }
    {
        IsotropicMonoEnergyBeamCircle<> b({ 0, 0, 1 }, { 1, 1, 0 });
        b.setEnergy(70);
        b.setCollimationHalfAngle(0.2);
        testBeam(b, "IsotropicMonoEnergyBeamCircle");
    }
    {
        IsotropicCircularBeam<> b({ 0, 0, 1 });
        b.setRadius(30);
        b.setEnergySpecter(Tube(110).getSpecter());
        b.setCollimationHalfAngles(-0.1, -0.1, 0.1, 0.1);
        testBeam(b, "IsotropicCircularBeam");
    }
    {
        IsotropicCircularMonoEnergyBeam<> b({ 0, 0, 1 });
        b.setRadius(30);
        b.setEnergy(50);
        b.setCollimationHalfAngles(-0.1, -0.1, 0.1, 0.1);
        testBeam(b, "IsotropicCircularMonoEnergyBeam");
    }
    {
        DXBeam<> b({ 0, 0, -50 });
        b.setTubeVoltage(85);
        b.addTubeFiltrationMaterial(13, 2);
        b.setCollimationHalfAnglesDeg(4, 6);
        b.setDAPvalue(5);
        testBeam(b, "DXBeam");
    }
    {
        CBCTBeam<> b;
        b.setTubeVoltage(100);
        b.setStartAngleDeg(-30);
        b.setStopAngleDeg(40);
        b.setStepAngleDeg(2);
        b.setCollimationHalfAnglesDeg(8, 6);
        b.setDAPvalue(10);
        testBeam(b, "CBCTBeam");
    }
    // CT beams calibrate with a CTDI phantom simulation, calibration is not compared
    {
        CTSpiralBeam<> b;
        b.setStartStopPosition({ 0, 0, -2 }, { 0, 0, 2 });
        b.setCollimation(1);
        b.setPitch(1.2);
        b.setStepAngleDeg(10);
        b.setTubeVoltage(100);
        testBeam(b, "CTSpiralBeam", false);
    }
    {
        CTSequentialBeam<> b;
        b.setNumberOfSlices(2);
        b.setSliceSpacing(3);
        b.setStepAngleDeg(10);
        b.setTubeVoltage(110);
        testBeam(b, "CTSequentialBeam", false);
    }
    {
        CTSpiralDualEnergyBeam<> b;
        b.setStartStopPosition({ 0, 0, -2 }, { 0, 0, 2 });
        b.setTubeAVoltage(140);
        b.setTubeBVoltage(80);
        b.setStepAngleDeg(10);
        testBeam(b, "CTSpiralDualEnergyBeam", false);
    }
}

void testWorld()
{
    using Cylinder = DepthDose<5, 1>;
    using Room = EnclosedRoom<5, 1>;
    World<Cylinder, Room> world;
    world.reserveNumberOfItems(2);
    auto& cylinder = world.addItem<Cylinder>({ 1, 10 }, "Cylinder");
    cylinder.setMaterial(Material<5>::byZ(13).value(), 2.27);
    auto& room = world.addItem<Room>({ 2, 200 }, "Room");
    room.setMaterial(Material<5>::byNistName("Concrete, Ordinary").value(), NISTMaterials::density("Concrete, Ordinary"));
    world.setMaterialDensity(0.002);
    world.build();

    PencilBeam<> beam({ 0, 0, -10 }, { 0, 0, 1 }, 60);
    beam.setNumberOfExposures(4);
    beam.setNumberOfParticlesPerExposure(5000);
    Transport::run(world, beam, 2);

    // save world and beam to one file and load both
    auto save = Serializer::getEmptyBuffer();
    Serializer::serializeItem(world.magicID(), world.serialize(), save);
    Serializer::serializeItem(beam, save);
    const auto path = xraymctest::tempPath("world.xrl").string();
    check(Serializer::write(path, save, true), "World: write world and beam to file");
    const auto load = Serializer::read(path);
    std::filesystem::remove(path);
    check(load.has_value(), "World: read world and beam from file");
    if (!load)
        return;
    auto next = std::span<const char> { load.value() };
    std::array<char, 32> name;
    std::vector<char> buffer;
    next = Serializer::deserializeItem(name, buffer, next);
    check(World<Cylinder, Room>::validMagicID(name), "World: magic ID");
    auto wl = World<Cylinder, Room>::deserialize(buffer);
    next = Serializer::deserializeItem(name, buffer, next);
    auto bl = PencilBeam<>::deserialize(buffer);
    check(wl.has_value() && bl.has_value(), "World: world and beam deserialized");
    if (!wl || !bl)
        return;
    auto& w2 = wl.value();
    check(w2.items().size() == 2 && w2.itemNames() == world.itemNames(), "World: items and names preserved");
    checkNear(w2.fillMaterialDensity(), 0.002, 1E-15, "World: fill material density preserved");
    check(w2.fillMaterial().composition() == world.fillMaterial().composition(), "World: fill material preserved");
    const auto* c2 = std::get_if<Cylinder>(w2.getItemPointerFromName("Cylinder"));
    check(c2 != nullptr, "World: item found by name after deserialization");
    if (c2) {
        bool doseEqual = true;
        for (std::size_t i = 0; i < 100; ++i)
            doseEqual = doseEqual && c2->doseScored(i) == cylinder.doseScored(i);
        check(doseEqual && c2->doseScored(0).dose() > 0, "World: scored dose preserved");
    }
    w2.build();
    check(w2.AABB() == world.AABB(), "World: rebuilt world has same AABB");
    Particle p { .pos = { 0, 0, -10 }, .dir = { 0, 0, 1 }, .energy = 60, .weight = 1 };
    const auto i1 = world.intersect(p);
    const auto i2 = w2.intersect(p);
    check(i1.valid() && i2.valid() && i1.intersection == i2.intersection, "World: rebuilt world has same intersections");
}

int main()
{
    std::cout << "Testing serialization\n";
    testPrimitives();
    testCompositeTypes();
    testFiles();
    testItems();
    testBeams();
    testWorld();
    return xraymctest::summary("testserializer");
}
