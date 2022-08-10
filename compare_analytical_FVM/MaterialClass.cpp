/* 
 * Author: zyuhao 
 * Date: 2022-04-08
 * Content: source code of material class
 * contain the physic parameters of a material relative to heat transfer
 */

//

#include "MaterialClass.h"

// default ctor
MaterialClass::MaterialClass(string name)
{
  materialName = name;
  dens = 7780;
  absorp = 0.3;
  specHeat = 523;
  thermalCond = 20.;
  thermalDiff = thermalCond / (dens * specHeat);
  tempLiquid = 3200;
  convecCoeff = 25;
  emissivity = 0.8;
}


// set function of materialName
void MaterialClass::setMaterialName(const string inMaterialName)
{
  materialName = inMaterialName;
}

// set function of density
void MaterialClass::setDens(const double inDens)
{
  dens = inDens;
}

// set function of absorp
void MaterialClass::setAbsorp(const double inAbsorp)
{
  absorp = inAbsorp;
}
// set function of specHeat
void MaterialClass::setSpecHeat(const double inSpecHeat)
{
  specHeat = inSpecHeat;
}

// set function of thermalCond
void MaterialClass::setThermalCond(const double inThermalCond)
{
  thermalCond = inThermalCond;
}

// set function of thermalDiff
void MaterialClass::setThermalDiff(const double inThermalDiff)
{
  thermalDiff = inThermalDiff;
}

// set function of tempLiquid
void MaterialClass::setTempLiquid(const double inTempLiquid)
{
  tempLiquid = inTempLiquid;
}

// set function of convecCoeff
void MaterialClass::setConvecCoeff(const double inConvecCoeff)
{
  convecCoeff = inConvecCoeff;
}

// set function of emissivity
void MaterialClass::setEmissivity(const double inEmissivity)
{
  emissivity = inEmissivity;
}


ostream& operator<<(ostream& os, const MaterialClass& myMaterial)
{
  os << endl << "Material parameters list: " << endl
       << "  Material name:        " << myMaterial.materialName << endl
       << "  Density:              " << myMaterial.dens << endl
       << "  Power absorptivity:   " << myMaterial.absorp << endl
       << "  Specific heat:        " << myMaterial.specHeat << endl
       << "  Thermal conductivity: " << myMaterial.thermalCond << endl
       << "  Thermal diffusivity:  " << myMaterial.thermalDiff << endl
       << "  Liquid temperature:   " << myMaterial.tempLiquid << endl
       << "  Convective coeff:     " << myMaterial.convecCoeff << endl
       << "  Emissivity:           " << myMaterial.emissivity << endl
       << endl;

  return os;
}