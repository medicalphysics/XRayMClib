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

// NIST material compositions obtained from https://physics.nist.gov/PhysRefData/XrayMassCoef/tab2.html

#pragma once

#include <algorithm>
#include <map>
#include <string>
#include <vector>

namespace xraymc {

class NISTMaterials {
public:
    static const std::vector<std::string> listNames()
    {
        const auto& instance = Instance();
        std::vector<std::string> names(instance.nistdata.size());
        std::transform(instance.nistdata.cbegin(), instance.nistdata.cend(), names.begin(),
            [](const auto& n) { return n.first; });
        return names;
    }
    static std::map<std::size_t, double> Composition(const std::string& name)
    {
        const auto& instance = Instance();
        if (instance.nistdata.contains(name))
            return instance.nistdata.at(name).massFractions;
        std::map<std::size_t, double> empty;
        return empty;
    }
    static double density(const std::string& name)
    {
        const auto& instance = Instance();
        if (instance.nistdata.contains(name))
            return instance.nistdata.at(name).density;
        return 0;
    }

    NISTMaterials(const NISTMaterials&) = delete;
    void operator=(const NISTMaterials&) = delete;

protected:
    struct NISTdata {
        double density = 0;
        std::map<std::size_t, double> massFractions;
    };
    static const NISTMaterials& Instance()
    {
        static NISTMaterials instance;
        return instance;
    }

    NISTMaterials()
    {
        nistdata["A-150 Tissue-Equivalent Plastic"] = { .density = 1.127f, .massFractions = { { 1, 0.101327f }, { 6, 0.775501f }, { 7, 0.035057f }, { 8, 0.052316f }, { 9, 0.017422f }, { 20, 0.018378f } } };
        nistdata["Acetone"] = { .density = 0.7899f, .massFractions = { { 1, 0.104122f }, { 6, 0.620405f }, { 8, 0.275473f } } };
        nistdata["Acetylene"] = { .density = 0.001097f, .massFractions = { { 1, 0.077418f }, { 6, 0.922582f } } };
        nistdata["Adenine"] = { .density = 1.35f, .massFractions = { { 1, 0.037294f }, { 6, 0.44443f }, { 7, 0.518275f } } };
        nistdata["Adipose Tissue (ICRP)"] = { .density = 0.92f, .massFractions = { { 1, 0.119477f }, { 6, 0.63724f }, { 7, 0.00797f }, { 8, 0.232333f }, { 11, 0.0005f }, { 12, 2e-05f }, { 15, 0.00016f }, { 16, 0.00073f }, { 17, 0.00119f }, { 19, 0.00032f }, { 20, 2e-05f }, { 26, 2e-05f }, { 30, 2e-05f } } };
        nistdata["Air, Dry (near sea level)"] = { .density = 0.001205f, .massFractions = { { 6, 0.000124f }, { 7, 0.755267f }, { 8, 0.231781f }, { 18, 0.012827f } } };
        nistdata["Alanine"] = { .density = 1.42f, .massFractions = { { 1, 0.07919f }, { 6, 0.404439f }, { 7, 0.157213f }, { 8, 0.359159f } } };
        nistdata["Aluminum Oxide"] = { .density = 3.97f, .massFractions = { { 8, 0.470749f }, { 13, 0.529251f } } };
        nistdata["Amber"] = { .density = 1.1f, .massFractions = { { 1, 0.10593f }, { 6, 0.788973f }, { 8, 0.105096f } } };
        nistdata["Ammonia"] = { .density = 0.000826f, .massFractions = { { 1, 0.177547f }, { 7, 0.822453f } } };
        nistdata["Aniline"] = { .density = 1.0235f, .massFractions = { { 1, 0.075759f }, { 6, 0.773838f }, { 7, 0.150403f } } };
        nistdata["Anthracene"] = { .density = 1.283f, .massFractions = { { 1, 0.05655f }, { 6, 0.94345f } } };
        nistdata["B-100 Bone-Equivalent Plastic"] = { .density = 1.45f, .massFractions = { { 1, 0.065471f }, { 6, 0.536945f }, { 7, 0.0215f }, { 8, 0.032085f }, { 9, 0.167411f }, { 20, 0.176589f } } };
        nistdata["Bakelite"] = { .density = 1.25f, .massFractions = { { 1, 0.057441f }, { 6, 0.774591f }, { 8, 0.167968f } } };
        nistdata["Barium Fluoride"] = { .density = 4.89f, .massFractions = { { 9, 0.21672f }, { 56, 0.78328f } } };
        nistdata["Barium Sulfate"] = { .density = 4.5f, .massFractions = { { 8, 0.274212f }, { 16, 0.137368f }, { 56, 0.58842f } } };
        nistdata["Benzene"] = { .density = 0.87865f, .massFractions = { { 1, 0.077418f }, { 6, 0.922582f } } };
        nistdata["Beryllium oxide"] = { .density = 3.01f, .massFractions = { { 4, 0.36032f }, { 8, 0.63968f } } };
        nistdata["Bismuth Germanium oxide"] = { .density = 7.13f, .massFractions = { { 8, 0.154126f }, { 32, 0.17482f }, { 83, 0.671054f } } };
        nistdata["Blood (ICRP)"] = { .density = 1.06f, .massFractions = { { 1, 0.101866f }, { 6, 0.10002f }, { 7, 0.02964f }, { 8, 0.759414f }, { 11, 0.00185f }, { 12, 4e-05f }, { 14, 3e-05f }, { 15, 0.00035f }, { 16, 0.00185f }, { 17, 0.00278f }, { 19, 0.00163f }, { 20, 6e-05f }, { 26, 0.00046f }, { 30, 1e-05f } } };
        nistdata["Bone, Compact (ICRU)"] = { .density = 1.85f, .massFractions = { { 1, 0.063984f }, { 6, 0.278f }, { 7, 0.027f }, { 8, 0.410016f }, { 12, 0.002f }, { 15, 0.07f }, { 16, 0.002f }, { 20, 0.147f } } };
        nistdata["Bone, Cortical (ICRP)"] = { .density = 1.85f, .massFractions = { { 1, 0.047234f }, { 6, 0.14433f }, { 7, 0.04199f }, { 8, 0.446096f }, { 12, 0.0022f }, { 15, 0.10497f }, { 16, 0.00315f }, { 20, 0.20993f }, { 30, 0.0001f } } };
        nistdata["Boron Carbide"] = { .density = 2.52f, .massFractions = { { 5, 0.78261f }, { 6, 0.21739f } } };
        nistdata["Boron Oxide"] = { .density = 1.812f, .massFractions = { { 5, 0.310551f }, { 8, 0.689449f } } };
        nistdata["Brain (ICRP)"] = { .density = 1.03f, .massFractions = { { 1, 0.110667f }, { 6, 0.12542f }, { 7, 0.01328f }, { 8, 0.737723f }, { 11, 0.00184f }, { 12, 0.00015f }, { 15, 0.00354f }, { 16, 0.00177f }, { 17, 0.00236f }, { 19, 0.0031f }, { 20, 9e-05f }, { 26, 5e-05f }, { 30, 1e-05f } } };
        nistdata["Butane"] = { .density = 0.002493f, .massFractions = { { 1, 0.173408f }, { 6, 0.826592f } } };
        nistdata["N-Butyl Alcohol"] = { .density = 0.8098f, .massFractions = { { 1, 0.135978f }, { 6, 0.648171f }, { 8, 0.215851f } } };
        nistdata["C-552 Air-Equivalent Plastic"] = { .density = 1.76f, .massFractions = { { 1, 0.02468f }, { 6, 0.50161f }, { 8, 0.004527f }, { 9, 0.465209f }, { 14, 0.003973f } } };
        nistdata["Cadmium Telluride"] = { .density = 6.2f, .massFractions = { { 48, 0.468355f }, { 52, 0.531645f } } };
        nistdata["Cadmium Tungstate"] = { .density = 7.9f, .massFractions = { { 8, 0.177644f }, { 48, 0.312027f }, { 74, 0.510329f } } };
        nistdata["Calcium Carbonate"] = { .density = 2.8f, .massFractions = { { 6, 0.120003f }, { 8, 0.479554f }, { 20, 0.400443f } } };
        nistdata["Calcium Fluoride"] = { .density = 3.18f, .massFractions = { { 9, 0.486659f }, { 20, 0.513341f } } };
        nistdata["Calcium Oxide"] = { .density = 3.3f, .massFractions = { { 8, 0.285299f }, { 20, 0.714701f } } };
        nistdata["Calcium Sulfate"] = { .density = 2.96f, .massFractions = { { 8, 0.470095f }, { 16, 0.235497f }, { 20, 0.294408f } } };
        nistdata["Calcium Tungstate"] = { .density = 6.062f, .massFractions = { { 8, 0.22227f }, { 20, 0.139202f }, { 74, 0.638529f } } };
        nistdata["Carbon Dioxide"] = { .density = 0.001842f, .massFractions = { { 6, 0.272916f }, { 8, 0.727084f } } };
        nistdata["Carbon Tetrachloride"] = { .density = 1.594f, .massFractions = { { 6, 0.078083f }, { 17, 0.921917f } } };
        nistdata["Cellulose Acetate, Cellophane"] = { .density = 1.42f, .massFractions = { { 1, 0.062162f }, { 6, 0.444462f }, { 8, 0.493376f } } };
        nistdata["Cellulose Acetate Butyrate"] = { .density = 1.2f, .massFractions = { { 1, 0.067125f }, { 6, 0.545403f }, { 8, 0.387472f } } };
        nistdata["Cellulose Nitrate"] = { .density = 1.49f, .massFractions = { { 1, 0.029216f }, { 6, 0.271296f }, { 7, 0.121276f }, { 8, 0.578212f } } };
        nistdata["Ceric Sulfate Dosimeter Solution"] = { .density = 1.03f, .massFractions = { { 1, 0.107596f }, { 7, 0.0008f }, { 8, 0.874976f }, { 16, 0.014627f }, { 58, 0.002001f } } };
        nistdata["Cesium Fluoride"] = { .density = 4.115f, .massFractions = { { 9, 0.125069f }, { 55, 0.874931f } } };
        nistdata["Cesium Iodide"] = { .density = 4.51f, .massFractions = { { 53, 0.488451f }, { 55, 0.511549f } } };
        nistdata["Chlorobenzene"] = { .density = 1.1058f, .massFractions = { { 1, 0.044772f }, { 6, 0.640254f }, { 17, 0.314974f } } };
        nistdata["Chloroform"] = { .density = 1.4832f, .massFractions = { { 1, 0.008443f }, { 6, 0.100613f }, { 17, 0.890944f } } };
        nistdata["Concrete, Portland"] = { .density = 2.3f, .massFractions = { { 1, 0.01f }, { 6, 0.001f }, { 8, 0.529107f }, { 11, 0.016f }, { 12, 0.002f }, { 13, 0.033872f }, { 14, 0.337021f }, { 19, 0.013f }, { 20, 0.044f }, { 26, 0.014f } } };
        nistdata["Concrete, Ordinary"] = { .density = 2.3f, .massFractions = { { 1, 0.022100f }, { 6, 0.002484f }, { 8, 0.574930f }, { 11, 0.01520f }, { 12, 0.00126f }, { 13, 0.01995f }, { 14, 0.30462f }, { 19, 0.01004f }, { 20, 0.04295f }, { 26, 0.00643f } } };
        nistdata["Cyclohexane"] = { .density = 0.779f, .massFractions = { { 1, 0.143711f }, { 6, 0.856289f } } };
        nistdata["1,2-Ddihlorobenzene"] = { .density = 1.3048f, .massFractions = { { 1, 0.027425f }, { 6, 0.490233f }, { 17, 0.482342f } } };
        nistdata["Dichlorodiethyl Ether"] = { .density = 1.2199f, .massFractions = { { 1, 0.056381f }, { 6, 0.335942f }, { 8, 0.111874f }, { 17, 0.495802f } } };
        nistdata["1,2-Dichloroethane"] = { .density = 1.2351f, .massFractions = { { 1, 0.04074f }, { 6, 0.242746f }, { 17, 0.716515f } } };
        nistdata["Diethyl Ether"] = { .density = 0.71378f, .massFractions = { { 1, 0.135978f }, { 6, 0.648171f }, { 8, 0.215851f } } };
        nistdata["N,N-Dimethyl Formamide"] = { .density = 0.9487f, .massFractions = { { 1, 0.096523f }, { 6, 0.492965f }, { 7, 0.191625f }, { 8, 0.218887f } } };
        nistdata["Dimethyl Sulfoxide"] = { .density = 1.1014f, .massFractions = { { 1, 0.077403f }, { 6, 0.307467f }, { 8, 0.204782f }, { 16, 0.410348f } } };
        nistdata["Ethane"] = { .density = 0.001253f, .massFractions = { { 1, 0.201115f }, { 6, 0.798885f } } };
        nistdata["Ethyl Alcohol"] = { .density = 0.7893f, .massFractions = { { 1, 0.131269f }, { 6, 0.521438f }, { 8, 0.347294f } } };
        nistdata["Ethyl Cellulose"] = { .density = 1.13f, .massFractions = { { 1, 0.090027f }, { 6, 0.585182f }, { 8, 0.324791f } } };
        nistdata["Ethylene"] = { .density = 0.001175f, .massFractions = { { 1, 0.143711f }, { 6, 0.856289f } } };
        nistdata["Eye Lens (ICRP)"] = { .density = 1.1f, .massFractions = { { 1, 0.099269f }, { 6, 0.19371f }, { 7, 0.05327f }, { 8, 0.653751f } } };
        nistdata["Ferric Oxide"] = { .density = 5.2f, .massFractions = { { 8, 0.300567f }, { 26, 0.699433f } } };
        nistdata["Ferroboride"] = { .density = 7.15f, .massFractions = { { 5, 0.162174f }, { 26, 0.837826f } } };
        nistdata["Ferrous Oxide"] = { .density = 5.7f, .massFractions = { { 8, 0.222689f }, { 26, 0.777311f } } };
        nistdata["Ferrous Sulfate Dosimeter Solution"] = { .density = 1.024f, .massFractions = { { 1, 0.108259f }, { 7, 2.7e-05f }, { 8, 0.878636f }, { 11, 2.2e-05f }, { 16, 0.012968f }, { 17, 3.4e-05f }, { 26, 5.4e-05f } } };
        nistdata["Freon-12"] = { .density = 1.12f, .massFractions = { { 6, 0.099335f }, { 9, 0.314247f }, { 17, 0.586418f } } };
        nistdata["Freon-12B2"] = { .density = 1.8f, .massFractions = { { 6, 0.057245f }, { 9, 0.181096f }, { 35, 0.761659f } } };
        nistdata["Freon-13"] = { .density = 0.95f, .massFractions = { { 6, 0.114983f }, { 9, 0.545622f }, { 17, 0.339396f } } };
        nistdata["Freon-13B1"] = { .density = 1.5f, .massFractions = { { 6, 0.080659f }, { 9, 0.382749f }, { 35, 0.536592f } } };
        nistdata["Freon-13I1"] = { .density = 1.8f, .massFractions = { { 6, 0.061309f }, { 9, 0.290924f }, { 53, 0.647767f } } };
        nistdata["Gadolinium Oxysulfide"] = { .density = 7.44f, .massFractions = { { 8, 0.084528f }, { 16, 0.08469f }, { 64, 0.830782f } } };
        nistdata["Gallium Arsenide"] = { .density = 5.31f, .massFractions = { { 31, 0.482019f }, { 33, 0.517981f } } };
        nistdata["Gel in Photographic Emulsion"] = { .density = 1.2914f, .massFractions = { { 1, 0.08118f }, { 6, 0.41606f }, { 7, 0.11124f }, { 8, 0.38064f }, { 16, 0.01088f } } };
        nistdata["Glass, Pyrex"] = { .density = 2.23f, .massFractions = { { 5, 0.040064f }, { 8, 0.539562f }, { 11, 0.028191f }, { 13, 0.011644f }, { 14, 0.37722f }, { 19, 0.003321f } } };
        nistdata["Glass, Lead"] = { .density = 6.22f, .massFractions = { { 8, 0.156453f }, { 14, 0.080866f }, { 22, 0.008092f }, { 33, 0.002651f }, { 82, 0.751938f } } };
        nistdata["Glass, Plate"] = { .density = 2.4f, .massFractions = { { 8, 0.4598f }, { 11, 0.096441f }, { 14, 0.336553f }, { 20, 0.107205f } } };
        nistdata["Glucose"] = { .density = 1.54f, .massFractions = { { 1, 0.071204f }, { 6, 0.363652f }, { 8, 0.565144f } } };
        nistdata["Glutamine"] = { .density = 1.46f, .massFractions = { { 1, 0.068965f }, { 6, 0.410926f }, { 7, 0.191681f }, { 8, 0.328427f } } };
        nistdata["Glycerol"] = { .density = 1.2613f, .massFractions = { { 1, 0.087554f }, { 6, 0.391262f }, { 8, 0.521185f } } };
        nistdata["Guanine"] = { .density = 1.58f, .massFractions = { { 1, 0.033346f }, { 6, 0.39738f }, { 7, 0.463407f }, { 8, 0.105867f } } };
        nistdata["Gypsum, Plaster of Paris"] = { .density = 2.32f, .massFractions = { { 1, 0.023416f }, { 8, 0.557572f }, { 16, 0.186215f }, { 20, 0.232797f } } };
        nistdata["N-Heptane"] = { .density = 0.68376f, .massFractions = { { 1, 0.160937f }, { 6, 0.839063f } } };
        nistdata["N-Hexane"] = { .density = 0.6603f, .massFractions = { { 1, 0.163741f }, { 6, 0.836259f } } };
        nistdata["Kapton Polyimide Film"] = { .density = 1.42f, .massFractions = { { 1, 0.026362f }, { 6, 0.691133f }, { 7, 0.07327f }, { 8, 0.209235f } } };
        nistdata["Lanthanum Oxybromide"] = { .density = 6.28f, .massFractions = { { 8, 0.068138f }, { 35, 0.340294f }, { 57, 0.591568f } } };
        nistdata["Lanthanum Oxysulfide"] = { .density = 5.86f, .massFractions = { { 8, 0.0936f }, { 16, 0.093778f }, { 57, 0.812622f } } };
        nistdata["Lead Oxide"] = { .density = 9.53f, .massFractions = { { 8, 0.071682f }, { 82, 0.928318f } } };
        nistdata["Lithium Amide"] = { .density = 1.178f, .massFractions = { { 1, 0.087783f }, { 3, 0.302262f }, { 7, 0.609955f } } };
        nistdata["Lithium Carbonate"] = { .density = 2.11f, .massFractions = { { 3, 0.187871f }, { 6, 0.16255f }, { 8, 0.649579f } } };
        nistdata["Lithium Fluoride"] = { .density = 2.635f, .massFractions = { { 3, 0.267585f }, { 9, 0.732415f } } };
        nistdata["Lithium Hydride"] = { .density = 0.82f, .massFractions = { { 1, 0.126797f }, { 3, 0.873203f } } };
        nistdata["Lithium Iodide"] = { .density = 3.494f, .massFractions = { { 3, 0.051858f }, { 53, 0.948142f } } };
        nistdata["Lithium Oxide"] = { .density = 2.013f, .massFractions = { { 3, 0.46457f }, { 8, 0.53543f } } };
        nistdata["Lithium Tetraborate"] = { .density = 2.44f, .massFractions = { { 3, 0.082085f }, { 5, 0.25568f }, { 8, 0.662235f } } };
        nistdata["Lung (ICRP)"] = { .density = 1.05f, .massFractions = { { 1, 0.101278f }, { 6, 0.10231f }, { 7, 0.02865f }, { 8, 0.757072f }, { 11, 0.00184f }, { 12, 0.00073f }, { 15, 0.0008f }, { 16, 0.00225f }, { 17, 0.00266f }, { 19, 0.00194f }, { 20, 9e-05f }, { 26, 0.00037f }, { 30, 1e-05f } } };
        nistdata["M3 Wax"] = { .density = 1.05f, .massFractions = { { 1, 0.114318f }, { 6, 0.655823f }, { 8, 0.092183f }, { 12, 0.134792f }, { 20, 0.002883f } } };
        nistdata["Magnesium Carbonate"] = { .density = 2.958f, .massFractions = { { 6, 0.142455f }, { 8, 0.569278f }, { 12, 0.288267f } } };
        nistdata["Magnesium Fluoride"] = { .density = 3.0f, .massFractions = { { 9, 0.609883f }, { 12, 0.390117f } } };
        nistdata["Magnesium Oxide"] = { .density = 3.58f, .massFractions = { { 8, 0.396964f }, { 12, 0.603036f } } };
        nistdata["Magnesium Tetraborate"] = { .density = 2.53f, .massFractions = { { 5, 0.240837f }, { 8, 0.62379f }, { 12, 0.135373f } } };
        nistdata["Mercuric Iodide"] = { .density = 6.36f, .massFractions = { { 53, 0.55856f }, { 80, 0.44144f } } };
        nistdata["Methane"] = { .density = 0.000667f, .massFractions = { { 1, 0.251306f }, { 6, 0.748694f } } };
        nistdata["Methanol"] = { .density = 0.7914f, .massFractions = { { 1, 0.125822f }, { 6, 0.374852f }, { 8, 0.499326f } } };
        nistdata["Mix D Wax"] = { .density = 0.99f, .massFractions = { { 1, 0.13404f }, { 6, 0.77796f }, { 8, 0.03502f }, { 12, 0.038594f }, { 22, 0.014386f } } };
        nistdata["MS20 Tissue Substitute"] = { .density = 1.0f, .massFractions = { { 1, 0.081192f }, { 6, 0.583442f }, { 7, 0.017798f }, { 8, 0.186381f }, { 12, 0.130287f }, { 17, 0.0009f } } };
        nistdata["Muscle, Skeletal"] = { .density = 1.04f, .massFractions = { { 1, 0.100637f }, { 6, 0.10783f }, { 7, 0.02768f }, { 8, 0.754773f }, { 11, 0.00075f }, { 12, 0.00019f }, { 15, 0.0018f }, { 16, 0.00241f }, { 17, 0.00079f }, { 19, 0.00302f }, { 20, 3e-05f }, { 26, 4e-05f }, { 30, 5e-05f } } };
        nistdata["Muscle, Striated"] = { .density = 1.04f, .massFractions = { { 1, 0.101997f }, { 6, 0.123f }, { 7, 0.035f }, { 8, 0.729003f }, { 11, 0.0008f }, { 12, 0.0002f }, { 15, 0.002f }, { 16, 0.005f }, { 19, 0.003f } } };
        nistdata["Muscle-Equivalent Liquid, with Sucrose"] = { .density = 1.11f, .massFractions = { { 1, 0.098234f }, { 6, 0.156214f }, { 7, 0.035451f }, { 8, 0.7101f } } };
        nistdata["Muscle-Equivalent Liquid, without Sucrose"] = { .density = 1.07f, .massFractions = { { 1, 0.101969f }, { 6, 0.120058f }, { 7, 0.035451f }, { 8, 0.742522f } } };
        nistdata["Naphthalene"] = { .density = 1.145f, .massFractions = { { 1, 0.062909f }, { 6, 0.937091f } } };
        nistdata["Nitrobenzene"] = { .density = 1.19867f, .massFractions = { { 1, 0.040935f }, { 6, 0.585374f }, { 7, 0.113773f }, { 8, 0.259918f } } };
        nistdata["Nitrous Oxide"] = { .density = 0.001831f, .massFractions = { { 7, 0.636483f }, { 8, 0.363517f } } };
        nistdata["Nylon, Du Pont ELVAmide 8062"] = { .density = 1.08f, .massFractions = { { 1, 0.103509f }, { 6, 0.648415f }, { 7, 0.099536f }, { 8, 0.148539f } } };
        nistdata["Nylon, type 6 and type 6/6"] = { .density = 1.14f, .massFractions = { { 1, 0.097976f }, { 6, 0.636856f }, { 7, 0.123779f }, { 8, 0.141389f } } };
        nistdata["Nylon, type 6/10"] = { .density = 1.14f, .massFractions = { { 1, 0.107062f }, { 6, 0.680449f }, { 7, 0.099189f }, { 8, 0.1133f } } };
        nistdata["Nylon, type 11 (Rilsan)"] = { .density = 1.425f, .massFractions = { { 1, 0.115476f }, { 6, 0.720819f }, { 7, 0.076417f }, { 8, 0.087289f } } };
        nistdata["Octane, Liquid"] = { .density = 0.7026f, .massFractions = { { 1, 0.158821f }, { 6, 0.841179f } } };
        nistdata["Paraffin Wax"] = { .density = 0.93f, .massFractions = { { 1, 0.148605f }, { 6, 0.851395f } } };
        nistdata["N-Pentane"] = { .density = 0.6262f, .massFractions = { { 1, 0.167635f }, { 6, 0.832365f } } };
        nistdata["Photographic Emulsion"] = { .density = 3.815f, .massFractions = { { 1, 0.0141f }, { 6, 0.072261f }, { 7, 0.01932f }, { 8, 0.066101f }, { 16, 0.00189f }, { 35, 0.349103f }, { 47, 0.474105f }, { 53, 0.00312f } } };
        nistdata["Plastic Scintillator (Vinyltoluene based)"] = { .density = 1.032f, .massFractions = { { 1, 0.085f }, { 6, 0.915f } } };
        nistdata["Plutonium Dioxide"] = { .density = 11.46f, .massFractions = { { 8, 0.118055f }, { 94, 0.881945f } } };
        nistdata["Polyacrylonitrile"] = { .density = 1.17f, .massFractions = { { 1, 0.056983f }, { 6, 0.679056f }, { 7, 0.263962f } } };
        nistdata["Polycarbonate (Makrolon, Lexan)"] = { .density = 1.2f, .massFractions = { { 1, 0.055491f }, { 6, 0.755751f }, { 8, 0.188758f } } };
        nistdata["Polychlorostyrene"] = { .density = 1.3f, .massFractions = { { 1, 0.061869f }, { 6, 0.696325f }, { 17, 0.241806f } } };
        nistdata["Polyethylene"] = { .density = 0.94f, .massFractions = { { 1, 0.143711f }, { 6, 0.856289f } } };
        nistdata["Polyethylene Terephthalate (Mylar)"] = { .density = 1.4f, .massFractions = { { 1, 0.041959f }, { 6, 0.625017f }, { 8, 0.333025f } } };
        nistdata["Polymethyl Methacralate (Lucite, Perspex)"] = { .density = 1.19f, .massFractions = { { 1, 0.080538f }, { 6, 0.599848f }, { 8, 0.319614f } } };
        nistdata["Polyoxymethylene"] = { .density = 1.425f, .massFractions = { { 1, 0.067135f }, { 6, 0.400017f }, { 8, 0.532848f } } };
        nistdata["Polypropylene"] = { .density = 0.9f, .massFractions = { { 1, 0.143711f }, { 6, 0.856289f } } };
        nistdata["Polystyrene"] = { .density = 1.06f, .massFractions = { { 1, 0.077418f }, { 6, 0.922582f } } };
        nistdata["Polytetrafluoroethylene (Teflon)"] = { .density = 2.2f, .massFractions = { { 6, 0.240183f }, { 9, 0.759817f } } };
        nistdata["Polytrifluorochloroethylene"] = { .density = 2.1f, .massFractions = { { 6, 0.20625f }, { 9, 0.489354f }, { 17, 0.304395f } } };
        nistdata["Polyvinyl Acetate"] = { .density = 1.19f, .massFractions = { { 1, 0.070245f }, { 6, 0.558066f }, { 8, 0.371689f } } };
        nistdata["Polyvinyl Alcohol"] = { .density = 1.3f, .massFractions = { { 1, 0.091517f }, { 6, 0.545298f }, { 8, 0.363185f } } };
        nistdata["Polyvinyl Butyral"] = { .density = 1.12f, .massFractions = { { 1, 0.092802f }, { 6, 0.680561f }, { 8, 0.226637f } } };
        nistdata["Polyvinyl Chloride"] = { .density = 1.3f, .massFractions = { { 1, 0.04838f }, { 6, 0.38436f }, { 17, 0.56726f } } };
        nistdata["Polyvinylidene Chloride, Saran"] = { .density = 1.7f, .massFractions = { { 1, 0.020793f }, { 6, 0.247793f }, { 17, 0.731413f } } };
        nistdata["Polyvinylidene Fluoride"] = { .density = 1.76f, .massFractions = { { 1, 0.03148f }, { 6, 0.375141f }, { 9, 0.593379f } } };
        nistdata["Polyvinyl Pyrrolidone"] = { .density = 1.25f, .massFractions = { { 1, 0.081616f }, { 6, 0.648407f }, { 7, 0.126024f }, { 8, 0.143953f } } };
        nistdata["Potassium Iodide"] = { .density = 3.13f, .massFractions = { { 19, 0.235528f }, { 53, 0.764472f } } };
        nistdata["Potassium Oxide"] = { .density = 2.32f, .massFractions = { { 8, 0.169852f }, { 19, 0.830148f } } };
        nistdata["Propane"] = { .density = 0.001879f, .massFractions = { { 1, 0.182855f }, { 6, 0.817145f } } };
        nistdata["Propane, Liquid"] = { .density = 0.43f, .massFractions = { { 1, 0.182855f }, { 6, 0.817145f } } };
        nistdata["N-Propyl Alcohol"] = { .density = 0.8035f, .massFractions = { { 1, 0.134173f }, { 6, 0.599595f }, { 8, 0.266232f } } };
        nistdata["Pyridine"] = { .density = 0.9819f, .massFractions = { { 1, 0.06371f }, { 6, 0.759217f }, { 7, 0.177073f } } };
        nistdata["Rubber, Butyl"] = { .density = 0.92f, .massFractions = { { 1, 0.143711f }, { 6, 0.856289f } } };
        nistdata["Rubber, Natural"] = { .density = 0.92f, .massFractions = { { 1, 0.118371f }, { 6, 0.881629f } } };
        nistdata["Rubber, Neoprene"] = { .density = 1.23f, .massFractions = { { 1, 0.05692f }, { 6, 0.542646f }, { 17, 0.400434f } } };
        nistdata["Silicon Dioxide"] = { .density = 2.32f, .massFractions = { { 8, 0.532565f }, { 14, 0.467435f } } };
        nistdata["Silver Bromide"] = { .density = 6.473f, .massFractions = { { 35, 0.425537f }, { 47, 0.574463f } } };
        nistdata["Silver Chloride"] = { .density = 5.56f, .massFractions = { { 17, 0.247368f }, { 47, 0.752632f } } };
        nistdata["Silver Halides in Photographic Emulsion"] = { .density = 6.47f, .massFractions = { { 35, 0.422895f }, { 47, 0.573748f }, { 53, 0.003357f } } };
        nistdata["Silver Iodide"] = { .density = 6.01f, .massFractions = { { 47, 0.459458f }, { 53, 0.540542f } } };
        nistdata["Skin (ICRP)"] = { .density = 1.1f, .massFractions = { { 1, 0.100588f }, { 6, 0.22825f }, { 7, 0.04642f }, { 8, 0.619002f }, { 11, 7e-05f }, { 12, 6e-05f }, { 15, 0.00033f }, { 16, 0.00159f }, { 17, 0.00267f }, { 19, 0.00085f }, { 20, 0.00015f }, { 26, 1e-05f }, { 30, 1e-05f } } };
        nistdata["Sodium Carbonate"] = { .density = 2.532f, .massFractions = { { 6, 0.113323f }, { 8, 0.452861f }, { 11, 0.433815f } } };
        nistdata["Sodium Iodide"] = { .density = 3.667f, .massFractions = { { 11, 0.153373f }, { 53, 0.846627f } } };
        nistdata["Sodium Monoxide"] = { .density = 2.27f, .massFractions = { { 8, 0.258143f }, { 11, 0.741857f } } };
        nistdata["Sodium Nitrate"] = { .density = 2.261f, .massFractions = { { 7, 0.164795f }, { 8, 0.56472f }, { 11, 0.270485f } } };
        nistdata["Stilbene"] = { .density = 0.9707f, .massFractions = { { 1, 0.067101f }, { 6, 0.932899f } } };
        nistdata["Sucrose"] = { .density = 1.5805f, .massFractions = { { 1, 0.064779f }, { 6, 0.42107f }, { 8, 0.514151f } } };
        nistdata["Terphenyl"] = { .density = 1.234f, .massFractions = { { 1, 0.044543f }, { 6, 0.955457f } } };
        nistdata["Testes (ICRP)"] = { .density = 1.04f, .massFractions = { { 1, 0.104166f }, { 6, 0.09227f }, { 7, 0.01994f }, { 8, 0.773884f }, { 11, 0.00226f }, { 12, 0.00011f }, { 15, 0.00125f }, { 16, 0.00146f }, { 17, 0.00244f }, { 19, 0.00208f }, { 20, 0.0001f }, { 26, 2e-05f }, { 30, 2e-05f } } };
        nistdata["Tetrachloroethylene"] = { .density = 1.625f, .massFractions = { { 6, 0.144856f }, { 17, 0.855144f } } };
        nistdata["Thallium Chloride"] = { .density = 7.004f, .massFractions = { { 17, 0.147822f }, { 81, 0.852178f } } };
        nistdata["Tissue, Soft (ICRP)"] = { .density = 1.0f, .massFractions = { { 1, 0.104472f }, { 6, 0.23219f }, { 7, 0.02488f }, { 8, 0.630238f }, { 11, 0.00113f }, { 12, 0.00013f }, { 15, 0.00133f }, { 16, 0.00199f }, { 17, 0.00134f }, { 19, 0.00199f }, { 20, 0.00023f }, { 26, 5e-05f }, { 30, 3e-05f } } };
        nistdata["Tissue, Soft (ICRU four-component)"] = { .density = 1.0f, .massFractions = { { 1, 0.101172f }, { 6, 0.111f }, { 7, 0.026f }, { 8, 0.761828f } } };
        nistdata["Tissue-Equivalent GAS (Methane based)"] = { .density = 0.001064f, .massFractions = { { 1, 0.101869f }, { 6, 0.456179f }, { 7, 0.035172f }, { 8, 0.40678f } } };
        nistdata["Tissue-Equivalent GAS (Propane based)"] = { .density = 0.001826f, .massFractions = { { 1, 0.102672f }, { 6, 0.56894f }, { 7, 0.035022f }, { 8, 0.293366f } } };
        nistdata["Titanium Dioxide"] = { .density = 4.26f, .massFractions = { { 8, 0.400592f }, { 22, 0.599408f } } };
        nistdata["Toluene"] = { .density = 0.8669f, .massFractions = { { 1, 0.08751f }, { 6, 0.91249f } } };
        nistdata["Trichloroethylene"] = { .density = 1.46f, .massFractions = { { 1, 0.007671f }, { 6, 0.182831f }, { 17, 0.809498f } } };
        nistdata["Triethyl Phosphate"] = { .density = 1.07f, .massFractions = { { 1, 0.082998f }, { 6, 0.395628f }, { 8, 0.351334f }, { 15, 0.17004f } } };
        nistdata["Tungsten Hexafluoride"] = { .density = 2.4f, .massFractions = { { 9, 0.382723f }, { 74, 0.617277f } } };
        nistdata["Uranium Dicarbide"] = { .density = 11.28f, .massFractions = { { 6, 0.091669f }, { 92, 0.908331f } } };
        nistdata["Uranium Monocarbide"] = { .density = 13.63f, .massFractions = { { 6, 0.048036f }, { 92, 0.951964f } } };
        nistdata["Uranium Oxide"] = { .density = 10.96f, .massFractions = { { 8, 0.118502f }, { 92, 0.881498f } } };
        nistdata["Urea"] = { .density = 1.323f, .massFractions = { { 1, 0.067131f }, { 6, 0.199999f }, { 7, 0.466459f }, { 8, 0.266411f } } };
        nistdata["Valine"] = { .density = 1.23f, .massFractions = { { 1, 0.094641f }, { 6, 0.512645f }, { 7, 0.119565f }, { 8, 0.27315f } } };
        nistdata["Viton Fluoroelastomer"] = { .density = 1.8f, .massFractions = { { 1, 0.009417f }, { 6, 0.280555f }, { 9, 0.710028f } } };
        nistdata["Water, Liquid"] = { .density = 1.0f, .massFractions = { { 1, 0.111894f }, { 8, 0.888106f } } };
        nistdata["Water Vapor"] = { .density = 0.000756f, .massFractions = { { 1, 0.111894f }, { 8, 0.888106f } } };
        nistdata["Xylene"] = { .density = 0.87f, .massFractions = { { 1, 0.094935f }, { 6, 0.905065f } } };
    }

private:
    std::map<std::string, NISTdata> nistdata;
};
}
