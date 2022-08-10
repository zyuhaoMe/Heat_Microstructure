/* 
 * Author: zyuhao 
 * Date: 2022-04-08
 * Content: head file of MaterialClass
 * contain the physic parameters of a material relative to heat transfer
 */

//

#ifndef COMPARE_ANALYTICAL_FVM_MATERIALCLASS_H
#define COMPARE_ANALYTICAL_FVM_MATERIALCLASS_H

#include <iostream>
#include <string>

using namespace std;

class MaterialClass {
  public:
    string materialName; // material name
    double dens; // material density
    double absorp; // power absorptivity
    double specHeat; // specific Heat
    double thermalCond; // thermal conductivity
    double thermalDiff; // thermal diffusivity
    double tempLiquid; // liquid temperature
    double convecCoeff; // convective coefficient (used in convective boundary)
    double emissivity; // emissivity (used in radiation of boundary)

    // default ctor
    MaterialClass(string name);

    // set function of materialName
    void setMaterialName(const string inMaterialName);

    // set function of density
    void setDens(const double inDens);

    // set function of absorp
    void setAbsorp(const double inAbsorp);

    // set function of specHeat
    void setSpecHeat(const double inSpecHeat);

    // set function of thermalCond
    void setThermalCond(const double inThermalCond);

    // set function of thermalDiff
    void setThermalDiff(const double inThermalDiff);

    // set function of tempLiquid
    void setTempLiquid(const double inTempLiquid);

    // set function of convecCoeff
    void setConvecCoeff(const double inConvecCoeff);

    // set function of emissivity
    void setEmissivity(const double inEmissivity);

    // friend function to allow operator overloading
    friend ostream& operator<<(ostream& os, const MaterialClass& myMaterial);
};

// overloading << for MaterialClass
ostream& operator<<(ostream& os, const MaterialClass& myMaterial);

#endif //COMPARE_ANALYTICAL_FVM_MATERIALCLASS_H
