#include <CLHEP/Units/PhysicalConstants.h>
using CLHEP::pi; using CLHEP::twopi; using CLHEP::halfpi;

#include <RAT/DU/LightPathCalculator.hh>
#include <RAT/Log.hh>
#include <RAT/DB.hh>
#include <RAT/DU/PMTInfo.hh>
#include <RAT/DU/Utility.hh>
using namespace RAT;
using namespace RAT::DU;

#include <TVector3.h>
#include <TMath.h>

#include <algorithm>
#include <string>

#include <iostream>
using namespace std;

ClassImp( RAT::DU::LightPathCalculator )

//////////////////////////////////////
//////////////////////////////////////

namespace RAT {

namespace DU {

void 
LightPathCalculator::LoadRefractiveIndex( DBLinkPtr dbTable,
                                          TGraph& property )
{

  // Define the constants for the conversion of the wavelength values [nm] -> [MeV]

  // Planck's Constant * Speed of Light in Vacuum ( 197.3270 MeV fm )
  // equivalent to 197.3270e-6 MeV nm
  const double hbarc = 197.32705e-6;

  // Enumeration for database value type of the refractive indices
  enum EWavelength { eWavelength, eDyWavelength, eEnergy };
  int wavelengthOption = eEnergy;

  // Obtain the value type from the RATDB table for the refractive indices
  const string option = dbTable->GetS( "RINDEX_option" );
  if( option == string( "wavelength" ) )
    wavelengthOption = eWavelength;
  else if( option == string( "dy_dwavelength" ) )
    wavelengthOption = eDyWavelength;
  else if( option == string( "energy" ) )
    wavelengthOption = eEnergy;

  // Fill in the values into two vectors, the wavelengths (val1) and the corresponding refractive indices (val2)
  vector< double > val1 = dbTable->GetDArray( "RINDEX_value1" );
  vector< double > val2 = dbTable->GetDArray( "RINDEX_value2" );

  // Check that the vectors actually have the wavelength/refractive indices values in them
  Log::Assert( val1.size() == val2.size() && !val1.empty(), "LightPathCalculator::LoadPropertyVector: Value arrays are miss-sized for RINDEX in " + dbTable->GetName() + "[" + dbTable->GetIndex() + "]" );
  int start = 0, dir = +1;

  // Now, based on the value type (i.e. wavelength, energy, dy_dWavelength) declare the values in the TGraph

  if( wavelengthOption != eEnergy && val1.front() < val1.back() ) // Must invert the order as energy is 1/wavelength
    {
      start = val1.size() - 1;
      dir = -1;
    }
  else if( wavelengthOption == eEnergy && val1.front() > val1.back() ) // Must invert the energy order
    {
      start = val1.size() - 1;
      dir = -1;
    }
  int point = 0;
  for( int index = start; index >= 0 && index < static_cast<int>( val1.size() ); index += dir )
    {
      double energy = val1[index];
      double value = val2[index];
      if( wavelengthOption == eDyWavelength )
        {
          energy = twopi * hbarc / ( val1[index] );
          value *= val1[index] / energy; // Times value by the
        }
      else if( wavelengthOption == eWavelength )
        energy = twopi * hbarc / ( val1[index] );
      property.SetPoint( point, energy, value );
      point++;
    }
}

//////////////////////////////////////
//////////////////////////////////////

void
LightPathCalculator::BeginOfRun()
{
  DB *db = DB::Get();

  // Obtain the inner/outer radii of the acrylic vessel from the RAT database
  fAVInnerRadius = db->GetLink( "SOLID", "acrylic_vessel_inner" )->GetD( "r_sphere" );
  fAVOuterRadius = db->GetLink( "SOLID", "acrylic_vessel_outer" )->GetD( "r_sphere" ); 

  // Obtain the radius of the PMT bucket from the RAT database
  fPMTRadius = db->GetLink( "GREY_DISC_PARAMETERS", "DiscOptics0_0" )->GetD( "disc_radius" );

  // Obtain the inner/outer radii of the acrylic vessel neck from the RAT database
  fNeckInnerRadius = db->GetLink( "SOLID", "acrylic_vessel_inner" )->GetD( "r_neck" );
  fNeckOuterRadius = db->GetLink( "SOLID", "acrylic_vessel_outer" )->GetD( "r_neck" );

  // If the geometry is using a partial fill configuration, obtain the fill fraction (fFillZ) as well as the refractive indices
  // for the upper/lower targets
  try 
    { 
      fFillZ = db->GetLink( "GEO", "inner_av" )->GetD( "split_z" ); 
      LoadRefractiveIndex( db->GetLink( "OPTICS", db->GetLink( "GEO", "inner_av" )->GetS( "material_top" ) ), fUpperTargetRI );
      LoadRefractiveIndex( db->GetLink( "OPTICS", db->GetLink( "GEO", "inner_av" )->GetS( "material_bottom" ) ), fLowerTargetRI );
    }
  catch( DBNotFoundError& e ) 
    { 
      fFillZ = fAVInnerRadius; 
      warn << "LightPathCalculator::BeginOfRun: Homogeneous not partial filled scint assumed.\n"; 
    }

  // Now obtain the refractive indices of the three main detector media, the material in the 'scint' region, the acrylic of
  // the acrylic vessel region and the water in the water region surrounding the acrylic vessel
  try { LoadRefractiveIndex( db->GetLink( "OPTICS", db->GetLink( "GEO", "inner_av" )->GetS( "material" ) ), fScintRI ); }
  catch( DBNotFoundError& e ) { warn << "LightPathCalculator::BeginOfRun: scint material refractive index cannot be loaded.\n"; }
  try { LoadRefractiveIndex( db->GetLink( "OPTICS", db->GetLink( "GEO", "av" )->GetS( "material" ) ), fAVRI ); }
  catch( DBNotFoundError& e ) { warn << "LightPathCalculator::BeginOfRun: av material refractive index cannot be loaded.\n"; }
  try { LoadRefractiveIndex( db->GetLink( "OPTICS", db->GetLink( "GEO", "cavity" )->GetS( "material" ) ), fWaterRI ); }
  catch( DBNotFoundError& e ) { warn << "LightPathCalculator::BeginOfRun: cavity material refractive index cannot be loaded.\n"; }

  // Declare the enumerator values for the path types. 

  fLightPathTypeMap[ SAW ] = "Scint->AV->Water";
  fLightPathTypeMap[ AW ] = "AV->Water";
  fLightPathTypeMap[ ASAW ] = "AV->Scint->AV->Water";
  fLightPathTypeMap[ WASAW ] = "Water->AV->Scint->AV->Water";
  fLightPathTypeMap[ WAW ] = "Water->AV->Water";
  fLightPathTypeMap[ W ] = "Water";
  fLightPathTypeMap[ WRefl ] = "Water->Reflection->Water";
  fLightPathTypeMap[ Null ] = "Light Path Uninitialised";
  

  // Reset/Initialise the values of all the private member variables in this class to dummy values.
  // By 'dummy', I mean non-physical values, values for which there is no physical interpretation in the
  // context under the function of this class.
  ResetValues();

}

//////////////////////////////////////
//////////////////////////////////////

void
LightPathCalculator::ResetValues()
{
  
  // Reset the values of the private member variables to dummy values
  
  fFillFraction = -10.0;
  fLoopCeiling = -10.0;
  fFinalLoopSize = -10.0;
  fPathPrecision = -10.0;
  
  fScintRIVal = -10.0;
  fUpperTargetRIVal = -10.0;
  fLowerTargetRIVal = -10.0;
  fAVRIVal = -10.0;
  fWaterRIVal = -10.0;
  
  fIncidentVecOnPMT.SetXYZ( 0.0, 0.0, 0.0 );
  fInitialLightVec.SetXYZ( 0.0, 0.0, 0.0 );
  
  fStartPos.SetXYZ( -10.0e6, -10.0e6, -10.0e6 );
  fEndPos.SetXYZ( -10.0e6, -10.0e6, -10.0e6 );
  fLightPathEndPos.SetXYZ( -10.0e6, -10.0e6, -10.0e6 );
  
  fIsTIR = false;
  fResvHit = false;
  fXAVNeck = false;
  fELLIEReflect = false;
  fStraightLine = false;
  
  fPointOnAV1st.SetXYZ( -10.0e6, -10.0e6, -10.0e6 );
  fPointOnAV2nd.SetXYZ( -10.0e6, -10.0e6, -10.0e6 );
  fPointOnAV3rd.SetXYZ( -10.0e6, -10.0e6, -10.0e6 );
  fPointOnAV4th.SetXYZ( -10.0e6, -10.0e6, -10.0e6 );
  
  fPointOnNeck1st.SetXYZ( -10.0e6, -10.0e6, -10.0e6 );
  fPointOnNeck2nd.SetXYZ( -10.0e6, -10.0e6, -10.0e6 );
  
  fLightPathType = Null;

  fDistInScint = -100.0;
  fDistInUpperTarget = -100.0;
  fDistInLowerTarget = -100.0;
  fDistInAV = -100.0;
  fDistInWater = -100.0;
  fDistInNeckScint = 0.0;
  fDistInNeckAV = 0.0;
  fDistInNeckWater = 0.0;
  
  fEnergy = WavelengthToEnergy( 400.0e-6 );
  
  fSolidAngle = -100.0;
  fCosThetaAvg = -2.0;
  
  fFresnelTCoeff = -10.0;
  fFresnelRCoeff = -10.0;
  
}

//////////////////////////////////////
//////////////////////////////////////

Double_t
LightPathCalculator::EnergyToWavelength( const Double_t energy )
{
  
  // Define the constants for the conversion of the wavelength values [mm] -> [MeV]
  // 2 * Pi...
  const Double_t twopi = 2.0 * 3.14159265358979323846;
  
  // Planck's Constant * Speed of Light in Vacuum ( 197.3270 MeV fm )
  // equivalent to 197.3270e-12 MeV mm (MeV, mm being the RAT units)
  const Double_t hbarc = 197.32705e-12;
  
  // Scale 'hbarc' in MeV mm based on units of 'energy'
  Double_t eFactor = twopi * hbarc;
  
  // Check that the energy is non-zero
  if ( energy == 0.0 ){ 
    warn << "LightPathCalculator::EnergyToWavelength: Warning, zero energy provided!\n"; 
    return 0.0;
  }
  Double_t energyInmm = eFactor / energy;
  
  return energyInmm;
  
}


//////////////////////////////////////
//////////////////////////////////////

Double_t
LightPathCalculator::WavelengthToEnergy( const Double_t wavelength )
{
  
  // Define the constants for the conversion of the wavelength values [mm] -> [MeV]
  // 2 * Pi...
  const Double_t twopi = 2.0 * 3.14159265358979323846;
  
  // Planck's Constant * Speed of Light in Vacuum ( 197.3270 MeV fm )
  // equivalent to 197.3270e-12 MeV mm (MeV, mm being the RAT units)
  const Double_t hbarc = 197.32705e-12;
  
  // Scale 'hbarc' in MeV mm
  Double_t eFactor = twopi * hbarc;
  
  // Check that the wavelength is non-zero
  if ( wavelength == 0.0 ){ 
    warn << "LightPathCalculator::WavelengthToEnergy: Warning, zero wavelength provided!\n";
    return 0.0;
  }
  Double_t wavelengthInMeV = eFactor / wavelength;
  
  return wavelengthInMeV;
  
}

//////////////////////////////////////
//////////////////////////////////////

void 
LightPathCalculator::CalcByPosition( const TVector3& eventPos,
                                     const TVector3& pmtPos,
                                     const Double_t& energyMeV,
                                     const Double_t& localityVal )
{

  // Ensure all booleans are set to false (except fELLIEReflect) 
  fResvHit = false;
  fIsTIR = false;
  fStraightLine = false;

  // Set the start and end position requirements of the path
  fStartPos = eventPos;
  fEndPos = pmtPos;

  // Set the energy of the source 
  // and the tolerance of the final path position
  fEnergy = energyMeV;
  fPathPrecision = localityVal;

  // Check to see if the straight line approximation is required.
  if ( fPathPrecision == 0.0 ){ 
    // Set the refractive indices all to unity.
    fScintRIVal = 1.0;
    fAVRIVal = 1.0;
    fWaterRIVal = 1.0;
    
    // Set the straight line path boolean
    fStraightLine = true;
  }

  else{
    // Initalise the refractive indices of the scintillator, acrylic and water regions
    // based on the energy value provided
    fScintRIVal = GetScintRI( fEnergy );
    fAVRIVal = GetAVRI( fEnergy );
    fWaterRIVal = GetWaterRI( fEnergy ); 

    fStraightLine = false;
  }

  // Check that the refractive indices are Initialised Properly
  if ( std::isnan( fScintRIVal * fAVRIVal * fWaterRIVal ) ){
    warn << "LightPathCalculator::CalcByPosition: Error, one or all of the refractive indices is 'nan', check your wavelength value is within limits!\n";
    return;
  }

  // Check that the refractive indices are not zero
  if ( fScintRIVal * fAVRIVal * fWaterRIVal == 0.0 ){
    warn << "LightPathCalculator::CalcByPosition: Error, one or all of the refractive indices is zero, check your wavelength value is within limits!\n";
    return;
  }
  
  // Check the size of the loop (i.e. maximum number of iterations)
  // if not declared - set to 20
  if ( fLoopCeiling <= 0 ){ fLoopCeiling = 20; }
  
  // Check for the path precision - the proximity required to the PMT bucket
  // If not declared, set to 10 mm to within the PMT.
  if ( fPathPrecision <= 0.0 ){ fPathPrecision = 10.0; }
  
  // Begin with the initial light path direction as the straight line direction
  TVector3 initOffset = ( fEndPos - fStartPos ).Unit();
  
  // Beginning of algorithm loop
  for ( Int_t iVal = 0; iVal < fLoopCeiling; iVal++ ){
    
    // Calculate refracted path based on the initial direction 'initOffset'
    PathCalculation( initOffset );
    
    // Check for locality, loop ceiling or total internal reflection (TIR)
    // break loop if 'fLightPathEndPos' is close enough to 'fEndPos' or if
    // total internal reflection was encountered (TIR)
    if( ( LocalityCheck( iVal ) )  || fIsTIR ){

      if ( fIsTIR ){
        fScintRIVal = 1.0;
        fAVRIVal = 1.0;
        fWaterRIVal = 1.0;
        PathCalculation( ( fEndPos - fStartPos ).Unit() );
        fFinalLoopSize = 0.0;
        fIsTIR = true;
        fStraightLine = true;
      }

      else{ fFinalLoopSize = iVal; fIsTIR = false; }

      break; 

    }
    
    // Readjust initial photon vector to be used in the next iteration and
    // recalculate the path
    else{ ReadjustOffset( fLightPathEndPos, initOffset ); }

  }

  return;

}

//////////////////////////////////////
//////////////////////////////////////

TVector3
LightPathCalculator::PathRefraction( const TVector3& incidentVec,
                                     const TVector3& incidentSurfVec,
                                     const Double_t& incRIndex,
                                     const Double_t& refRIndex )
{
  
  const Double_t ratioRI = incRIndex / refRIndex;
  const Double_t cosTheta1 = incidentSurfVec.Dot( -1.0 * incidentVec ); // Incident angle [Snell's Law]
  const Double_t cosTheta2 = TMath::Sqrt( 1 - TMath::Power( ratioRI, 2 ) * (  1 - TMath::Power( cosTheta1, 2 ) ) ); // Refracted Angle [Snell's Law]

  // Initialise the refracted photon vector
  TVector3 refractedVec( 0.0, 0.0, 0.0 );
  
  // Check for Total Internal Reflection (TIR) (equivalent to a 'nan' radicand i.e. a negative squareroot)
  if ( std::isnan( cosTheta2 ) ){ 
    fIsTIR = true;
    // Set the refracted vec to the straight equivalent
    refractedVec = incidentVec;
  }
  
  // Define the refracted vector
  else if ( cosTheta1 >= 0.0 ){ refractedVec = ( ratioRI * incidentVec ) + ( ( ratioRI * cosTheta1 ) - cosTheta2 ) * incidentSurfVec; }

  else { refractedVec = ( ratioRI * incidentVec ) - ( ( ratioRI * cosTheta1 ) - cosTheta2 ) * incidentSurfVec; }
  
  // Ensure the refracted vector is unit normalised
  return refractedVec.Unit();

}

//////////////////////////////////////
//////////////////////////////////////

void 
LightPathCalculator::PathCalculation( const TVector3& initOffset )
{

  fInitialLightVec = initOffset;
  
  // Calculation for events that originate within the scintillator region
  if ( fStartPos.Mag() < fAVInnerRadius ){
    
    // Calculate the vector from the event position in the scintillator to the inner AV edge
    fPointOnAV1st = VectorToSphereEdge( fStartPos, initOffset, fAVInnerRadius, 0 );
    
    // Calculate the refraction between the scintillator and the acrylic
    TVector3 vec1 = PathRefraction( initOffset, ( -1.0 * fPointOnAV1st ).Unit(), fScintRIVal, fAVRIVal );
    
    // Calculate the vector from the scintillator/av intersection to 
    // the outer AV edge
    fPointOnAV2nd = VectorToSphereEdge( fPointOnAV1st, vec1, fAVOuterRadius, 0 );
    
    // Calculate the refraction between the acrylic and the water
    TVector3 vec2 = PathRefraction( vec1, ( -1.0 * fPointOnAV2nd ).Unit(), fAVRIVal, fWaterRIVal );
    
    // Check whether path enters the neck region
    SetAVNeckInformation( fPointOnAV1st, initOffset );

    // Set the incident vector on the PMT bucket
    fIncidentVecOnPMT = vec2;
    
    // Calculate the vector from the av/water intersection to the 
    // hypothesised PMT position
    fLightPathEndPos = VectorToSphereEdge( fPointOnAV2nd, vec2, fEndPos.Mag(), 0 );

    // Set the light path type
    // Scint -> AV -> Water -> PMT
    fLightPathType = SAW;

    fDistInScint = ( fStartPos - fPointOnAV1st ).Mag();
    fDistInAV = ( fPointOnAV1st - fPointOnAV2nd ).Mag();
    fDistInWater = ( fPointOnAV2nd - fLightPathEndPos ).Mag();

    // If the light path went through the neck, then the values need to account for this
    // and are adjusted based on information concerning the neck calculated by 'SetNeckInformation'
    if ( fXAVNeck ){
      
      fDistInScint += fDistInNeckScint;
      fDistInAV = fDistInNeckAV;
      fDistInWater = fDistInNeckWater;
      fIncidentVecOnPMT = ( fLightPathEndPos - fPointOnNeck2nd ).Unit();
      
    }
    
    return;
    
  }
  
  // Calculation for events that originate within the AV
  else if ( ( fStartPos.Mag() > fAVInnerRadius ) && ( fStartPos.Mag() < fAVOuterRadius  ) ){
    
    // First need to check if the path enters the scintillator region...

    // Calculate the maximum angle threshold between the event position
    // and the light direction for the path to enter the scintillator region
    Double_t approachAngle = ClosestAngle( fStartPos, fAVInnerRadius );

    // The angle between event position and path direction
    Double_t angle = ( -1.0 * fStartPos ).Angle( initOffset );

    // These are calculated to ensure a non-zero discriminant in the calculation
    // ( see 'VectorToSphereEdge' )
    Double_t sinAngle = TMath::Sin( angle );
    Double_t rRatio = fAVInnerRadius / fStartPos.Mag();
    
    // If event enters the scintillator region...
    if ( ( angle < approachAngle ) && ( sinAngle < rRatio ) ){
      
      // Calculate the vector from the event position to the inner AV edge
      fPointOnAV1st = VectorToSphereEdge( fStartPos, initOffset, fAVInnerRadius, 1 );
      
      // Calculate the refraction between the acrylic and the scintillator
      TVector3 vec1 = PathRefraction( initOffset, ( fPointOnAV1st ).Unit(), fAVRIVal, fScintRIVal );
      
      // Calculate the vector from the inner AV radius to the otherside
      // of the scintillator region
      fPointOnAV2nd = VectorToSphereEdge( fPointOnAV1st, vec1, fAVInnerRadius, 0 );
      
      // Calculate the refraction between the scintillator and the acrylic
      TVector3 vec2 = PathRefraction( vec1, ( -1.0 * fPointOnAV2nd ).Unit(), fScintRIVal, fAVRIVal );
      
      // Calculate the vector from the scintillator/av intersection to the
      // av/water intersection
      fPointOnAV3rd = VectorToSphereEdge( fPointOnAV2nd, vec2, fAVOuterRadius, 0 );
      
      // Calculate the refraction between the acrylic and the water
      TVector3 vec3 = PathRefraction( vec2, ( -1.0 * fPointOnAV2nd ).Unit(), fAVRIVal, fWaterRIVal );
      
      // Check to see if the path enters the neck region
      SetAVNeckInformation( fPointOnAV2nd, vec1 );

      // The incident vector on the PMT bucket
      fIncidentVecOnPMT = vec3;
      
      // Calculate the vector to the hypothesised PMT position
      fLightPathEndPos = VectorToSphereEdge( fPointOnAV3rd, vec3, fEndPos.Mag(), 0 );

      // Set the light path type
      // AV -> Scint -> AV -> Water -> PMT
      fLightPathType = ASAW;

      fDistInScint = ( fPointOnAV1st - fPointOnAV2nd ).Mag();
      fDistInAV = ( fStartPos - fPointOnAV1st ).Mag() + ( fPointOnAV2nd - fPointOnAV3rd ).Mag();
      fDistInWater = ( fPointOnAV3rd - fLightPathEndPos ).Mag();

      if ( fXAVNeck ){
        
        fDistInScint += fDistInNeckScint;
        fDistInAV = ( fPointOnAV1st - fStartPos ).Mag() + fDistInNeckAV;
        fDistInWater = fDistInNeckWater;
        fIncidentVecOnPMT = ( fLightPathEndPos - fPointOnNeck2nd ).Unit();
        
      }
      
      return;
      
    }
    
    // Event exits the AV and immediately into the water
    else{
      
      // Calculate the vector from the event position to the outer AV edge
      fPointOnAV1st = VectorToSphereEdge( fStartPos, initOffset, fAVOuterRadius, 0 );
      
      // Calculate the refraction between the acrylic and the water
      TVector3 vec1 = PathRefraction( initOffset, ( -1.0 * fPointOnAV1st ).Unit(), fAVRIVal, fWaterRIVal );
      
      // Check to see if the light path enters the neck region
      SetAVNeckInformation( fPointOnAV1st, initOffset );

      // The incident vector on the PMT bucket
      fIncidentVecOnPMT = vec1;
      
      // Calculate the vector from the outer AV radius to the hypothesised
      // PMT position
      fLightPathEndPos = VectorToSphereEdge( fPointOnAV1st, vec1, fEndPos.Mag(), 0 );

      // Set the light path type
      // AV -> Water -> PMT
      fLightPathType = AW;

      fDistInScint = 0.0;
      fDistInAV = ( fStartPos - fPointOnAV1st ).Mag();
      fDistInWater = ( fPointOnAV1st - fLightPathEndPos ).Mag();

      if ( fXAVNeck ){
        
        fDistInScint += fDistInNeckScint;
        fDistInAV = fDistInNeckAV;
        fDistInWater = fDistInNeckWater;
        fIncidentVecOnPMT = ( fLightPathEndPos - fPointOnNeck2nd ).Unit();
        
      }
      
      return;
    }
  }
  
  // Calculation for events that originate outside of the AV (in the water)
  else {

    Double_t maxReflectionAngle = ReflectionAngle( fStartPos, fAVOuterRadius );

    // Check to see if reflect path off of the AV are required (fELLIEReflect = true),
    // and if so, whether they are eligible for being reflected off of the AV
    if ( fELLIEReflect && ( fStartPos.Angle( fEndPos ) < maxReflectionAngle ) ){

      // Assume that light reflected off of the outer AV surface
      
      // Find the magnitude weighted average vector between the pmt position
      // and event position.
      TVector3 avVec = ( fEndPos + fStartPos ).Unit();
      // Note: This is different to (fEndPos.Unit() + fStartPos.Unit()).Unit()
      // which is the bisection of the two unit vectors.
      
      // Calculate the point on the AV where the light reflected
      TVector3 reflectedPoint = avVec * fAVOuterRadius;
      fPointOnAV1st = reflectedPoint;
      
      // Distance through the acrylic and scintillator is zero
      fDistInScint = 0.0;
      fDistInAV = 0.0;
      Double_t a = (fStartPos.Mag()+fAVOffset) /(fStartPos.Mag());
      // Distance through the water
      fDistInWater = ( (a*fStartPos) - reflectedPoint ).Mag() + ( reflectedPoint - fEndPos ).Mag();
      
      fIncidentVecOnPMT = ( fEndPos - reflectedPoint ).Unit();
      fInitialLightVec = ( reflectedPoint - fStartPos ).Unit();
      
      fLightPathType = WRefl;
      
      return;
      
    }

    // Calculate the path as normal
    else{
    
      // First need to check if the path enters the AV region...
      
      // Calculate the maximum angle threshold between the event position
      // and the light direction for the path to enter the AV region
      Double_t approachAngle = ClosestAngle( fStartPos, fAVOuterRadius );
      
      // The angle between event position and path direction
      Double_t angle = ( -1.0 * fStartPos ).Angle( initOffset );
      
      // These are calculated to ensure a non-zero discriminant in the calculation
      // ( see 'VectorToSphereEdge' )
      Double_t sinAngle = TMath::Sin( fStartPos.Angle( initOffset ) );
      Double_t rRatio = fAVOuterRadius / fStartPos.Mag();
      
      // If the path enters the AV Aregion 
      if ( ( angle < approachAngle ) && ( sinAngle < rRatio ) ){
        
        // Calculate the vector from the event position to the outer AV edge
        fPointOnAV1st = VectorToSphereEdge( fStartPos, initOffset, fAVOuterRadius, 1 );
        
        // Calculate the refraction between the water and the acrylic
        TVector3 vec1 = PathRefraction( initOffset, ( fPointOnAV1st ).Unit(), fWaterRIVal, fAVRIVal );
        
        // Now need to check if the path enters the scintillator region...
        approachAngle = ClosestAngle( fPointOnAV1st, fAVInnerRadius );
        angle = ( -1.0 * fPointOnAV1st ).Angle( vec1 );
        
        sinAngle = TMath::Sin( fPointOnAV1st.Angle( vec1 ) );
        rRatio = fAVInnerRadius / fPointOnAV1st.Mag();
        
        // If the path then enters the scintillator region
        if ( ( angle < approachAngle ) && ( sinAngle < rRatio ) ){
          
          // Calculate path to the scintillator through the AV region
          fPointOnAV2nd = VectorToSphereEdge( fPointOnAV1st, vec1,  fAVInnerRadius, 1 );
          
          // Calculate the refraction between the acrylic and the scintillator
          TVector3 vec2 = PathRefraction( vec1, ( fPointOnAV2nd ).Unit(), fAVRIVal, fScintRIVal );
          
          // Calculate the vector from the inner AV radius to the otherside
          // of the scintillator region
          fPointOnAV3rd = VectorToSphereEdge( fPointOnAV2nd, vec2, fAVInnerRadius, 0 );
          
          // Calculate the refraction between the acrylic and the scintillator
          TVector3 vec3 = PathRefraction( vec2, ( -1.0 * fPointOnAV3rd ).Unit(), fAVRIVal, fScintRIVal );
          
          // Calculate the vector from the scint/av intersection to the av/water
          // intersection
          fPointOnAV4th = VectorToSphereEdge( fPointOnAV3rd, vec3, fAVOuterRadius, 0 );
          
          // Calculate the refraction between the acrylic and the water
          TVector3 vec4 =  PathRefraction( vec3, ( -1.0 * fPointOnAV4th ).Unit(), fAVRIVal, fWaterRIVal );
          
          // Check to see if the path enters the neck region
          SetAVNeckInformation( fPointOnAV3rd, vec2 ); 
          
          // The incident vector on the PMT bucket
          fIncidentVecOnPMT = vec4;
          
          // Calculate the vector to the hypothesised PMT position
          fLightPathEndPos = VectorToSphereEdge( fPointOnAV4th, vec4, fEndPos.Mag(), 0 );
          
          // Set the light path type
          // Water -> AV -> Scint -> AV -> Water -> PMT
          fLightPathType = WASAW;
          
          fDistInScint = ( fPointOnAV2nd - fPointOnAV3rd ).Mag();
          fDistInAV = ( fPointOnAV1st - fPointOnAV2nd ).Mag() + ( fPointOnAV3rd - fPointOnAV4th ).Mag();
          fDistInWater = ( fStartPos - fPointOnAV1st ).Mag() + ( fPointOnAV4th - fLightPathEndPos ).Mag();
          
          if ( fXAVNeck ){
            
            fDistInScint += fDistInNeckScint;
            fDistInAV = ( fPointOnAV1st - fPointOnAV2nd ).Mag() + fDistInNeckAV;
            fDistInWater = ( fStartPos - fPointOnAV1st ).Mag() + fDistInNeckWater;
            fIncidentVecOnPMT = ( fLightPathEndPos - fPointOnNeck2nd ).Unit();
            
          }
          
          return;
          
        }
        
        else{
          
          
          // Calculate the vector from the outer AV radius to the otherside
          // of the AV region
          fPointOnAV2nd = VectorToSphereEdge( fPointOnAV1st, vec1, fAVOuterRadius, 0 );    
          
          // Calculate the refraction between the acrylic and the water
          TVector3 vec2 = PathRefraction( vec1, ( -1.0 * fPointOnAV2nd ).Unit(), fAVRIVal, fWaterRIVal );
          
          // Check to see if the path enters the neck region
          SetAVNeckInformation( fPointOnAV2nd, vec1 ); 
          
          // Set the incident vector on the PMT bucket
          fIncidentVecOnPMT = vec2;
          
          // Calculate the vector from the av/water intersection to the
          // hypothesised PMT position
          fLightPathEndPos = VectorToSphereEdge( fPointOnAV2nd, vec2, fEndPos.Mag(), 0 );
          
          // Set the light path type
          // Water -> AV -> Water -> PMT
          fLightPathType = WAW;
          
          fDistInScint = 0.0;
          fDistInAV = ( fPointOnAV1st - fPointOnAV2nd ).Mag();
          fDistInWater = ( fStartPos - fPointOnAV1st ).Mag() + ( fPointOnAV2nd - fLightPathEndPos ).Mag();
          
          if ( fXAVNeck ){
            
            fDistInScint += fDistInNeckScint;
            fDistInAV = ( fPointOnAV1st - fPointOnAV2nd ).Mag() + fDistInNeckAV;
            fDistInWater = ( fStartPos - fPointOnAV1st ).Mag() + fDistInNeckWater;
            fIncidentVecOnPMT = ( fLightPathEndPos - fPointOnNeck2nd ).Unit();
            
          }
          
          return;
          
        }
      }
      
      // If event does not enter AV or scintillator region.
      else{
        
        fLightPathEndPos = fEndPos;

        fDistInScint = 0.0;
        fDistInAV = 0.0;
        fDistInWater = ( fStartPos - fEndPos ).Mag();
        
        // The incident vector on the PMT bucket
        fIncidentVecOnPMT = ( fEndPos - fStartPos ).Unit();
        
        // Set the light path type
        // Water -> PMT
        fLightPathType = W;
        
        return;
        
      }
    }
  }
}

//////////////////////////////////////
//////////////////////////////////////

Bool_t LightPathCalculator::LocalityCheck( const Int_t iVal )
  
{
  
  // Check the locality condition
  if ( ( ( fLightPathEndPos - fEndPos ).Mag() < fPathPrecision ) || ( iVal == fLoopCeiling - 1 ) ){
    
    // Meets the conditions
    
    // If in-precise measurement (only occurs at final loop value i == fLoopCeiling-1)
    if ( ( fLightPathEndPos - fEndPos ).Mag() > fPathPrecision ){
      fResvHit = true;
    }    
    return true;
  }
  
  else{
    // Does not meet the conditions
    return false;
  }
}

//////////////////////////////////////
//////////////////////////////////////

void 
LightPathCalculator::ReadjustOffset( const TVector3& distWater,
                                     TVector3& initOffset )
  
{
  
  Double_t scale = fStartPos.Mag() / ( 10000.0 );
  
  if ( distWater.Theta() > fEndPos.Theta() ){
    Double_t thetaAdj =  ( initOffset.Theta() ) - scale * ( TMath::Abs( distWater.Theta() - fEndPos.Theta() ) );
    initOffset.SetTheta( thetaAdj );
  }
  
  else{
    Double_t thetaAdj =  ( initOffset.Theta() ) + scale * ( TMath::Abs( distWater.Theta() - fEndPos.Theta() ) );
    initOffset.SetTheta( thetaAdj );
  }
  
  if ( distWater.Phi() > fEndPos.Phi() ){
    Double_t phiAdj =  ( initOffset.Phi() ) - scale * ( TMath::Abs(distWater.Phi() - fEndPos.Phi() ) );
    initOffset.SetPhi( phiAdj );
  }
  
  else{
    Double_t phiAdj =  ( initOffset.Phi() ) + scale * ( TMath::Abs( distWater.Phi() - fEndPos.Phi() ) );
    initOffset.SetPhi( phiAdj );
  }
}

//////////////////////////////////////
//////////////////////////////////////

Double_t 
LightPathCalculator::ClosestAngle( const TVector3& eventPos,
                                   const Double_t& edgeRadius )
{
  
  // Calculate one of the angles of the Right-Angled triangle to obtain
  // the other
  Double_t Phi = TMath::ACos( edgeRadius / ( eventPos.Mag() ) );
  Double_t angle = halfpi - Phi;
  
  return angle;

}

//////////////////////////////////////
//////////////////////////////////////

Double_t 
LightPathCalculator::ReflectionAngle( const TVector3& eventPos,
                                      const Double_t& edgeRadius )
{
  
  Double_t angle = 2.0 * ( TMath::ACos( edgeRadius / eventPos.Mag() ) ); 
  return angle;

}

//////////////////////////////////////
//////////////////////////////////////

void 
LightPathCalculator::CalcByPositionPartial( const TVector3& eventPos,
                                            const TVector3& pmtPos,
                                            const Double_t& energyMeV,
                                            const Double_t& localityVal )
{

  // Ensure all booleans are set to false (except fELLIEReflect) 
  fResvHit = false;
  fIsTIR = false;
  fStraightLine = false;

  // Set the start and end position requirements of the path
  fStartPos = eventPos;
  fEndPos = pmtPos;

  // Set the energy of the source 
  // and the tolerance of the final path position
  fEnergy = energyMeV;
  fPathPrecision = localityVal;

  // Check to see if the straight line approximation is required.
  if ( fPathPrecision == 0.0 ){ 
    // Set the refractive indices all to unity.
    fUpperTargetRIVal = 1.0;
    fLowerTargetRIVal = 1.0;
    fAVRIVal = 1.0;
    fWaterRIVal = 1.0;
    
    // Set the straight line path boolean
    fStraightLine = true;
  }

  else{
    // Initalise the refractive indices of the upper/lower target, acrylic and water regions
    // based on the energy value provided
    fUpperTargetRIVal = GetUpperTargetRI( fEnergy );
    fLowerTargetRIVal = GetLowerTargetRI( fEnergy );
    fAVRIVal = GetAVRI( fEnergy );
    fWaterRIVal = GetWaterRI( fEnergy ); 

    fStraightLine = false;
  }

  // Check that the refractive indices are Initialised Properly
  if ( std::isnan( fUpperTargetRIVal * fLowerTargetRIVal * fAVRIVal * fWaterRIVal ) ){
    warn << "LightPathCalculator::CalcByPosition: Error, one or all of the refractive indices is 'nan', check your wavelength value is within limits!\n";
    return;
  }

  // Check that the refractive indices are not zero
  if ( fUpperTargetRIVal * fLowerTargetRIVal * fAVRIVal * fWaterRIVal == 0.0 ){
    warn << "LightPathCalculator::CalcByPosition: Error, one or all of the refractive indices is zero, check your wavelength value is within limits!\n";
    return;
  }

  // Check the size of the loop (i.e. maximum number of iterations)
  // if not declared - set to 20
  if ( fLoopCeiling <= 0 ){ fLoopCeiling = 20; }

  // Check for the path precision - the proximity required to the PMT bucket
  // If not declared, set to 10 mm to within the PMT.
  if ( fPathPrecision <= 0.0 ){ fPathPrecision = 10.0; }
 
  
  // Begin with the initial light path direction as the straight line direction
  TVector3 initOffset = ( fEndPos - fStartPos ).Unit();
  
  // Beginning of algorithm loop
  for ( Int_t iVal = 0; iVal < fLoopCeiling; iVal++ ){  
    
    // Calculate refracted path based on the initial direction 'initOffset'
    PathCalculationPartial( initOffset );
    
    // Check for locality or total internal reflection (TIR)
    if( ( LocalityCheck( iVal ) ) || fIsTIR ){

      if ( fIsTIR ){
        fAVRIVal = 1.0;
        fWaterRIVal = 1.0;
        fUpperTargetRIVal = 1.0;
        fLowerTargetRIVal = 1.0;
        PathCalculationPartial( ( fEndPos - fStartPos ).Unit() );
        fIsTIR = true;
        fStraightLine = true;
        fFinalLoopSize = 0.0;
      }

      else{ fFinalLoopSize = iVal; fIsTIR = false; }

      break;
    }

    // Readjust initial photon vector to be used in the next iteration and
    // recalculate the path
    else{ ReadjustOffset( fLightPathEndPos, initOffset ); }
    
  }

  return;
}

//////////////////////////////////////
//////////////////////////////////////

TVector3
LightPathCalculator::VectorToSphereEdge( const TVector3 startPos,
                                         const TVector3 startDir,
                                         const Double_t radiusFromCentre,
                                         const bool outside )
{
  
  // The a, b, and c coefficients of a typical quadratic equation
  const Double_t aCoeff = 1.0;
  const Double_t bCoeff = 2.0 * startPos.Mag() * TMath::Cos( startPos.Angle( startDir ) );
  const Double_t cCoeff = startPos.Mag2() - ( radiusFromCentre * radiusFromCentre );
  
  // Parameter (set to 0.0 to begin with)
  Double_t distParam = 0.0;
  
  // Discriminant for quadratic equation
  const Double_t discrim = TMath::Power( bCoeff, 2 ) - ( 4.0 * aCoeff * cCoeff );
  
  // Intersection Calculation
  if ( discrim >= 0.0 ){
    
    const Double_t discrimPlus = ( ( -1.0 * bCoeff ) + TMath::Sqrt(discrim) ) / ( 2.0 * aCoeff );
    const Double_t discrimMinus = ( ( -1.0 * bCoeff ) - TMath::Sqrt(discrim) ) / ( 2.0 * aCoeff );
    
    if ( discrimPlus > 0.0 && !outside ){ distParam = discrimPlus; }   
    else { distParam = std::min( discrimMinus, discrimPlus ); }

  }

  else {
    fIsTIR = true;
    warn << "LightPathCalculator::VectorToSphereEdge: Erroneous calculation - negative discriminant" << endl;
    warn << "Discriminant Value: " << discrim << endl;
    warn << "ACoeff: " << aCoeff << ", BCoeff: " << bCoeff << ", CCoeff: " << cCoeff << endl;
    warn << "StartPos.Mag(): " << startPos.Mag() << " and R: " << radiusFromCentre << endl;
  }
  
  // Implementation of Parameterisation to find intersection point
  TVector3 endPointVec = startPos + ( distParam * startDir );

  return endPointVec;
  
}

//////////////////////////////////////
//////////////////////////////////////

TVector3
LightPathCalculator::VectorToCylinderEdge( const TVector3& startPos,
                                           const TVector3& startDir,
                                           const TVector3& cylinderBaseOrigin,
                                           const Double_t cylinderRadius )
{
  
  TVector3 startDirXY = startDir;
  startDirXY.SetZ( 0.0 );

  TVector3 startPosXY = startPos;
  startPosXY.SetZ( 0.0 );

  TVector3 pointOnCylinderBase = VectorToSphereEdge( startPosXY, startDirXY, cylinderRadius, 0 );

  Double_t angleTmp = ( pointOnCylinderBase.Unit() ).Angle( startDir );
  Double_t heightInCylinder = 0.0;

  if ( angleTmp != halfpi ){
    heightInCylinder = ( ( pointOnCylinderBase - startPosXY ).Mag() ) * TMath::Tan( angleTmp );
  }

  TVector3 pointOnCylinder = pointOnCylinderBase + cylinderBaseOrigin;
  pointOnCylinder.SetZ( pointOnCylinder.Z() + heightInCylinder );

  return pointOnCylinder;
    
}


//////////////////////////////////////
//////////////////////////////////////

void
LightPathCalculator::PathThroughTarget( const TVector3& enterPos,
                                        const TVector3& enterDir,
                                        TVector3& exitPos,
                                        TVector3& exitDir )
{

  // The z-coordinate of the interface between the lower and upper target in the AV
  Double_t fillZ = fAVInnerRadius * ( 1.0 - ( 2.0 * fFillFraction ) );

  // Initialise two refractive indices
  Double_t firstRI = 1.0;
  Double_t secondRI = 1.0;

  // Check which target the start position begins in...

  // ...the upper target
  if ( enterPos.Z() > fillZ ){
    firstRI = fUpperTargetRIVal;
    secondRI = fLowerTargetRIVal;
  }

  // ...the lower target
  else{
    firstRI = fLowerTargetRIVal;
    secondRI = fUpperTargetRIVal;
  }
  
  Double_t d = ( fillZ - enterPos.Z() ) / enterDir.Z();

  // check if it enters both targets
  if ( enterDir.Z() != 0 && d > 0 
       && ( pow ( d * enterDir.X() + enterPos.X(), 2 ) + pow( d * enterDir.Y() + enterPos.Y(), 2 ) 
            <= pow( fAVInnerRadius, 2 ) + pow( fillZ, 2 ) ) ){
    
    // Calculate the vector from the event position to the lower target
    TVector3 firstTargetPos = d * enterDir + enterPos;
    
    // Calculate the refraction between the upper target and the lower target
    TVector3 boundaryDir;
    if ( enterDir.Z() > 0 ){ boundaryDir = TVector3( 0, 0, -1 ); }

    else{ boundaryDir = TVector3( 0, 0, 1 ); }


    TVector3 vec0 = PathRefraction( enterDir, boundaryDir, firstRI, secondRI );

    // calculate the vector from the last point to the inner AV edge
    exitPos = VectorToSphereEdge( firstTargetPos, vec0, fAVInnerRadius, 0 );
    // calculate the refraction between the lower target and the acrylic
    exitDir = PathRefraction( vec0, ( -1.0 * exitPos ).Unit(), secondRI, fAVRIVal );

    if ( enterPos.Z() > fillZ ){
      fDistInUpperTarget = d;
      fDistInLowerTarget = ( exitPos - enterPos ).Mag() - d;
    }

    else{
      fDistInLowerTarget = d;
      fDistInUpperTarget = ( exitPos - enterPos ).Mag() - d;
    }
  }

  else{

    // Calculate the vector from the event position to the inner AV edge
    exitPos = VectorToSphereEdge( enterPos, enterDir, fAVInnerRadius, 0 );

    // Calculate the refraction between the scintillator and the acrylic
    exitDir = PathRefraction( enterDir,( -1.0 * exitPos ).Unit(), firstRI, fAVRIVal );

    if (enterPos.Z() > fillZ){ fDistInUpperTarget = ( exitPos - enterPos ).Mag(); }
    else{ fDistInLowerTarget = ( exitPos - enterPos ).Mag(); }
  }
}

//////////////////////////////////////
//////////////////////////////////////

void
LightPathCalculator::PathCalculationPartial( const TVector3& initOffset )
{

  fInitialLightVec = initOffset;

  Double_t fillZ = fAVInnerRadius - 2.0 * fAVInnerRadius * fFillFraction;

  // Calculation for events that originate within the target volume
  if ( fStartPos.Mag() < fAVInnerRadius ){

    TVector3 exitTargetPos, exitTargetDir;

    // Calculate the path through the upper/lower targets
    PathThroughTarget( fStartPos, fInitialLightVec, fPointOnAV1st, exitTargetDir );
   
    // Calculate the vector from inner AV to outer AV edge
    fPointOnAV2nd = VectorToSphereEdge( fPointOnAV1st, exitTargetDir, fAVOuterRadius, 0 );

    // calculate the refraction between the acrylic and the water interface
    TVector3 exitAVDir = PathRefraction( exitTargetDir, ( -1.0 * fPointOnAV2nd ).Unit(), fAVRIVal, fWaterRIVal );

    // Check whether path enters the neck region
    SetAVNeckInformation( fPointOnAV1st, exitTargetDir );

    // Set the incident vector on the PMT bucket
    fIncidentVecOnPMT = exitAVDir;

    // calculate the vector from the av/water intersection to the hypothesised PMT position
    fLightPathEndPos = VectorToSphereEdge( fPointOnAV2nd, exitAVDir, fEndPos.Mag(), 0 );

    // Set the light path type
    fLightPathType = SAW;

    fDistInAV = ( fPointOnAV1st - fPointOnAV2nd ).Mag();
    fDistInWater = ( fPointOnAV2nd - fLightPathEndPos ).Mag();

    // If the light path went through the neck, then the values need to account for this
    // and are adjusted based on information concerning the neck calculated by 'SetNeckInformation'
    if ( fXAVNeck ){
      
      fDistInUpperTarget += fDistInNeckScint;
      fDistInAV = fDistInNeckAV;
      fDistInWater = fDistInNeckWater;
      fIncidentVecOnPMT = ( fEndPos - fPointOnNeck2nd ).Unit();
      
    }

    return;   

  }
  // Calculation for events that originate within the AV
  else if ( ( fStartPos.Mag() > fAVInnerRadius ) && ( fStartPos.Mag() < fAVOuterRadius ) ){

    // First need to check if the path enters inside the scintillator region (containing the uperr/lower targets)...

    // Calculate the maximum angle threshold between the event position
    // and the light direction for the path to enter the scintillator region
    Double_t approachAngle = ClosestAngle( fStartPos, fAVInnerRadius );

    // The angle between event position and path direction
    Double_t angle = ( -1.0 * fStartPos ).Angle( initOffset );

    // These are calculated to ensure a non-zero discriminant in the calculation
    // ( see 'VectorToSphereEdge' )
    Double_t sinAngle = TMath::Sin( fStartPos.Angle( initOffset ) );
    Double_t rRatio = fAVInnerRadius / fStartPos.Mag();

    // If event enters the scintillator (lower/upper target) regions...
    if ( ( angle < approachAngle ) && ( sinAngle < rRatio ) ){

      // Calculate the vector from the event position to the inner AV edge
      fPointOnAV1st = VectorToSphereEdge( fStartPos, initOffset, fAVInnerRadius, 1 );

      Double_t targetRI = 1.0;
      // Check which target we are intersecting
      if ( fPointOnAV1st.Z() > fillZ )
        targetRI = fUpperTargetRIVal;
      else
        targetRI = fLowerTargetRIVal;

      // Calculate the refraction between the acrylic and the scintillator
      TVector3 enterTargetDir = PathRefraction( initOffset, ( fPointOnAV1st ).Unit(), fAVRIVal, targetRI );

      // Calculate the vector from the inner AV radius to the otherside
      // of the scintillator region
      TVector3 exitTargetDir;
      PathThroughTarget( fPointOnAV1st, enterTargetDir, fPointOnAV2nd, exitTargetDir );

      // Calculate the vector from the scintillator/AV intersection to the
      // AV/water intersection
      fPointOnAV3rd = VectorToSphereEdge( fPointOnAV2nd, exitTargetDir, fAVOuterRadius, 0 );

      // Calculate the refraction between the acrylic and the water
      TVector3 exitAVDir = PathRefraction( exitTargetDir, ( -1.0 * fPointOnAV3rd ).Unit(), fAVRIVal, fWaterRIVal );

      // Check to see if the path enters the neck region
      SetAVNeckInformation( fPointOnAV2nd, exitTargetDir );

      // The incident vector on the PMT bucket
      fIncidentVecOnPMT = exitAVDir;

      // Calculate the vector to the hypothesised PMT position
      fLightPathEndPos = VectorToSphereEdge( fPointOnAV3rd, exitAVDir, ( fEndPos ).Mag(), 0 );

      // Set the light path type
      // AV -> Scint -> AV -> Water -> PMT
      fLightPathType = ASAW;
      
      fDistInAV = ( fStartPos - fPointOnAV1st ).Mag() + ( fPointOnAV2nd - fPointOnAV3rd ).Mag();
      fDistInWater = ( fPointOnAV3rd - fLightPathEndPos ).Mag();

      if ( fXAVNeck ){
        
        fDistInUpperTarget += fDistInNeckScint;
        fDistInAV = ( fPointOnAV1st - fStartPos ).Mag() + fDistInNeckAV;
        fDistInWater = fDistInNeckWater;
        fIncidentVecOnPMT = ( fLightPathEndPos - fPointOnNeck2nd ).Unit();
        
      }

    }

    // Event exits the AV and immediately into the water
    else{

      // Calculate the vector from the event position to the outer AV edge
      fPointOnAV1st = VectorToSphereEdge( fStartPos, initOffset,  fAVOuterRadius, 0 );

      // Calculate the refraction between the acrylic and the water
      TVector3 exitAVDir = PathRefraction( initOffset, ( -1.0 * fPointOnAV1st ).Unit(), fAVRIVal, fWaterRIVal );

      // Check to see if the light path enters the neck region
      SetAVNeckInformation( fPointOnAV1st, initOffset );

      // The incident vector on the PMT bucket
      fIncidentVecOnPMT = exitAVDir;

      // Calculate the vector from the outer AV radius to the hypothesised
      // PMT position
      fLightPathEndPos = VectorToSphereEdge( fPointOnAV1st, exitAVDir,  fEndPos.Mag(), 0 ); 

      // Set the light path type
      // AV -> Water -> PMT
      fLightPathType = AW;

      fDistInAV = ( fStartPos - fPointOnAV1st ).Mag();
      fDistInWater = ( fPointOnAV1st - fLightPathEndPos ).Mag();  

      if ( fXAVNeck ){
        
        fDistInUpperTarget += fDistInNeckScint;
        fDistInAV = fDistInNeckAV;
        fDistInWater = fDistInNeckWater;
        fIncidentVecOnPMT = ( fLightPathEndPos - fPointOnNeck2nd ).Unit();
        
      }

      return;
    }
  }

  // Calculation for events that originate outside of the AV (in the water)
  else {

    Double_t maxReflectionAngle = ReflectionAngle( fStartPos, fAVOuterRadius );

    // Check to see if reflect path off of the AV are required (fELLIEReflect = true),
    // and if so, whether they are eligible for being reflected off of the AV
    if ( fELLIEReflect && ( fStartPos.Angle( fEndPos ) < maxReflectionAngle ) ){

      // Assume that light reflected off of the outer AV surface
      
      // Find the magnitude weighted average vector between the pmt position
      // and event position.
      TVector3 avVec = ( fEndPos + fStartPos ).Unit();
      // Note: This is different to (fEndPos.Unit() + fStartPos.Unit()).Unit()
      // which is the bisection of the two unit vectors.
      
      // Calculate the point on the AV where the light reflected
      TVector3 reflectedPoint = avVec * fAVOuterRadius;
      fPointOnAV1st = reflectedPoint;
      
      // Distance through the acrylic and scintillator is zero
      fDistInScint = 0.0;
      fDistInAV = 0.0;
      
      // Distance through the water
      fDistInWater = ( fStartPos - reflectedPoint ).Mag() + ( reflectedPoint - fEndPos ).Mag();
      
      fIncidentVecOnPMT = ( fEndPos - reflectedPoint ).Unit();
      fInitialLightVec = ( reflectedPoint - fStartPos ).Unit();
      
      fLightPathType = WRefl;
      
      return;

    }

    // Calculate the path as normal
    else{

      // First need to check if the path enters the AV region...
      
      // Calculate the maximum angle threshold between the event position
      // and the light direction for the path to enter the AV region
      Double_t approachAngle = ClosestAngle( fStartPos, fAVOuterRadius );
      
      // The angle between event position and path direction
      Double_t angle = ( -1.0 * fStartPos ).Angle( initOffset );
      
      // These are calculated to ensure a non-zero discriminant in the calculation
      // ( see 'VectorToSphereEdge' )
      Double_t sinAngle = TMath::Sin( fStartPos.Angle( initOffset ) );
      Double_t rRatio = fAVOuterRadius / fStartPos.Mag();
      
      // If the path enters the AV Aregion 
      if ( ( angle < approachAngle ) && ( sinAngle < rRatio ) ){
        
        // Calculate the vector from the event position to the outer AV edge
        fPointOnAV1st = VectorToSphereEdge( fStartPos, initOffset, fAVOuterRadius, 1 );
        
        // Calculate the refraction between the water and the acrylic
        TVector3 enterAVDir = PathRefraction( initOffset, ( fPointOnAV1st ).Unit(), fWaterRIVal, fAVRIVal );
        
        // Now need to check if the path enters the scintillator region...
        approachAngle = ClosestAngle( fPointOnAV1st, fAVInnerRadius );
        angle = ( -1.0 * fPointOnAV1st ).Angle( enterAVDir );
        
        sinAngle = TMath::Sin( fPointOnAV1st.Angle( enterAVDir ) );
        rRatio = fAVInnerRadius / fPointOnAV1st.Mag();
        
        // If the path then enters the scintillator region
        if ( ( angle < approachAngle ) && ( sinAngle < rRatio ) ){
          
          // Calculate path to the scintillator through the AV region
          fPointOnAV2nd = VectorToSphereEdge( fPointOnAV1st, enterAVDir, fAVInnerRadius, 1 );
          
          Double_t targetRI = 1.0;
          // Check which target we are intersecting
          if ( fPointOnAV2nd.Z() > fillZ ){ targetRI = fUpperTargetRIVal; }
          
          else{ targetRI = fLowerTargetRIVal; }
          
          // Calculate the refraction between the acrylic and the scintillator
          TVector3 enterTargetDir = PathRefraction( enterAVDir, ( fPointOnAV2nd ).Unit(), fAVRIVal, targetRI );
          
          // Calculate the vector from the inner AV radius to the otherside
          // of the scintillator region
          TVector3 exitTargetDir;
          PathThroughTarget( fPointOnAV2nd, enterTargetDir, fPointOnAV3rd, exitTargetDir );
          
          // Calculate the vector from the scint/av intersection to the av/water
          // intersection
          fPointOnAV4th = VectorToSphereEdge( fPointOnAV3rd, exitTargetDir, fAVOuterRadius, 0 );
          
          // Calculate the refraction between the acrylic and the water
          TVector3 exitAVDir =  PathRefraction( exitTargetDir, ( -1.0 * fPointOnAV4th ).Unit(), fAVRIVal, fWaterRIVal );
          
          // Check to see if the path enters the neck region
          SetAVNeckInformation( fPointOnAV3rd, exitTargetDir ); 
          
          // The incident vector on the PMT bucket
          fIncidentVecOnPMT = exitAVDir;
          
          // Calculate the vector to the hypothesised PMT position
          fLightPathEndPos = VectorToSphereEdge( fPointOnAV4th, exitAVDir, fEndPos.Mag(), 0 );
          
          // Set the light path type
          // Water -> AV -> Scint -> AV -> Water -> PMT
          fLightPathType = WASAW;
          
          fDistInAV = ( fPointOnAV1st - fPointOnAV2nd ).Mag() + ( fPointOnAV3rd - fPointOnAV4th ).Mag();
          fDistInWater = ( fStartPos - fPointOnAV1st ).Mag() + ( fPointOnAV4th - fLightPathEndPos ).Mag();
          
          if ( fXAVNeck ){
            
            fDistInUpperTarget += fDistInNeckScint;
            fDistInAV = ( fPointOnAV1st - fPointOnAV2nd ).Mag() + fDistInNeckAV;
            fDistInWater = ( fStartPos - fPointOnAV1st ).Mag() + fDistInNeckWater;
            fIncidentVecOnPMT = ( fLightPathEndPos - fPointOnNeck2nd ).Unit();
            
          }  
          
          return;
          
        }
        
        else{
          
          // Calculate the vector from the outer AV radius to the otherside
          // of the AV region
          fPointOnAV2nd = VectorToSphereEdge( fPointOnAV1st, enterAVDir, fAVOuterRadius, 0 );    
          
          // Calculate the refraction between the acrylic and the water
          TVector3 vec2 = PathRefraction( enterAVDir, ( -1.0 * fPointOnAV2nd ).Unit(), fAVRIVal, fWaterRIVal );
          
          // Check to see if the path enters the neck region
          SetAVNeckInformation( fPointOnAV2nd, enterAVDir ); 
          
          // Set the incident vector on the PMT bucket
          fIncidentVecOnPMT = vec2;
          
          // Calculate the vector from the av/water intersection to the
          // hypothesised PMT position
          fLightPathEndPos = VectorToSphereEdge( fPointOnAV2nd, vec2, fEndPos.Mag(), 0 );
          
          // Set the light path type
          // Water -> AV -> Water -> PMT
          fLightPathType = WAW;
          
          fDistInAV = ( fPointOnAV1st - fPointOnAV2nd ).Mag();
          fDistInWater = ( fStartPos - fPointOnAV1st ).Mag() + ( fPointOnAV2nd - fLightPathEndPos ).Mag();
          
          if ( fXAVNeck ){
            
            fDistInUpperTarget += fDistInNeckScint;
            fDistInAV = ( fPointOnAV1st - fPointOnAV2nd ).Mag() + fDistInNeckAV;
            fDistInWater = ( fStartPos - fPointOnAV1st ).Mag() + fDistInNeckWater;
            fIncidentVecOnPMT = ( fLightPathEndPos - fPointOnNeck2nd ).Unit();
            
          }
          
          return;
          
        }
      }
      
      // If event does not enter AV or scintillator region.
      else{
        
        fLightPathEndPos = fEndPos;
        
        fDistInScint = 0.0;
        fDistInAV = 0.0;
        fDistInWater = ( fStartPos - fEndPos ).Mag();
          
        // The incident vector on the PMT bucket
        fIncidentVecOnPMT = ( fEndPos - fStartPos ).Unit();
        
        // Set the light path type
        // Water -> PMT
        fLightPathType = W;
        
        return;
        
      }
    }
  }
}

//////////////////////////////////////
//////////////////////////////////////

void 
LightPathCalculator::SetAVNeckInformation( const TVector3& pointOnAV,
                                           const TVector3& dirVec )
{
  
  fDistInNeckScint = 0.0;
  fDistInNeckAV = 0.0;
  fDistInNeckWater = 0.0;

  Double_t minHeight = TMath::Sqrt( TMath::Power( fAVInnerRadius, 2.0 ) - TMath::Power( fNeckInnerRadius, 2.0 ) );
  TVector3 pointOnAVXY = pointOnAV;
  pointOnAVXY.SetZ( 0.0 );

  TVector3 zAxis( 0.0, 0.0, 1.0 );

  if ( pointOnAV.Z() > minHeight 
       && pointOnAVXY.Mag() < fNeckInnerRadius
       && dirVec.Angle( zAxis ) < halfpi ){

    fXAVNeck = true;

    // The refractive indices of the scintillator, av and water
    Double_t scRI = fScintRIVal;
    Double_t avRI = fAVRIVal;
    Double_t waterRI = fWaterRIVal;

    TVector3 neckBase( 0, 0, fAVInnerRadius );
    TVector3 pointOnInnerNeck = VectorToCylinderEdge( pointOnAV, dirVec, neckBase, fNeckInnerRadius );
    fDistInNeckScint = ( pointOnInnerNeck - pointOnAV ).Mag();
    fPointOnNeck1st = pointOnInnerNeck;

    TVector3 innerNeckNorm = pointOnInnerNeck;
    innerNeckNorm.SetZ( 0.0 );
    TVector3 vec1 = PathRefraction( dirVec, ( -1.0 * innerNeckNorm ).Unit(), scRI, avRI );

    TVector3 pointOnOuterNeck = VectorToCylinderEdge( pointOnInnerNeck, vec1, neckBase, fNeckOuterRadius );
    fDistInNeckAV = ( pointOnOuterNeck - pointOnInnerNeck ).Mag();
    fPointOnNeck2nd = pointOnOuterNeck;

    TVector3 outerNeckNorm = pointOnOuterNeck;
    outerNeckNorm.SetZ( 0.0 );
    TVector3 vec2 = PathRefraction( vec1, ( -1.0 * outerNeckNorm ).Unit(), avRI, waterRI );

    TVector3 pointInWater = VectorToSphereEdge( fPointOnNeck2nd, vec2, fEndPos.Mag(), 0 ); 

    fDistInNeckWater = ( fEndPos - pointInWater ).Mag();

  }

  else { fXAVNeck = false; }
}

//////////////////////////////////////
//////////////////////////////////////

void 
LightPathCalculator::CalculateSolidAngle( const TVector3& pmtNorm,
                                          const Int_t nVal )
{
  
  if ( nVal != 0 ){
    CalculateSolidAnglePolygon( pmtNorm, nVal );
    return;
  }

  Double_t cosThetaAvg = 0.0;

  TVector3 av1(0.0, 0.0, 0.0);
  TVector3 av2(0.0, 0.0, 0.0);
  TVector3 av3(0.0, 0.0, 0.0);
  TVector3 av4(0.0, 0.0, 0.0);

  TVector3 hypEndPos( 0, 0, 0 );

  if ( fStraightLine ){ CalcByPosition( fStartPos, fEndPos ); }
  else{ CalcByPosition( fStartPos, fEndPos, fEnergy, 10.0 ); }
  LightPathCalculator tmpStore = *this;

  TVector3 av = fAVInnerRadius * fInitialLightVec;
  
  TVector3 zPrime = pmtNorm;
  TVector3 yPrime = ( zPrime.Orthogonal() ).Unit();
  TVector3 xPrime = ( yPrime.Cross( zPrime ) ).Unit();

  TVector3 pmtVec1 = fEndPos + ( fPMTRadius * xPrime );
  TVector3 pmtVec2 = fEndPos - ( fPMTRadius * xPrime );
  TVector3 pmtVec3 = fEndPos + ( fPMTRadius * yPrime );
  TVector3 pmtVec4 = fEndPos - ( fPMTRadius * yPrime );

  if( !fStraightLine ){

    if( !fIsTIR && !fResvHit ){
      
      CalcByPosition( fStartPos, pmtVec1, fEnergy, fPathPrecision );
      if ( !fIsTIR && !fResvHit ) {av1 = fPointOnAV1st; cosThetaAvg += fIncidentVecOnPMT.Dot( -1.0 * pmtNorm.Unit() );}
      else {av1 = av;}
  
      CalcByPosition( fStartPos, pmtVec2, fEnergy, fPathPrecision );    
      if ( !fIsTIR && !fResvHit ) {av2 = fPointOnAV1st; cosThetaAvg += fIncidentVecOnPMT.Dot( -1.0 * pmtNorm.Unit() );}
      else {av2 = av;}
 
      CalcByPosition( fStartPos, pmtVec3, fEnergy, fPathPrecision );     
      if ( !fIsTIR && !fResvHit ) {av3 = fPointOnAV1st; cosThetaAvg += fIncidentVecOnPMT.Dot( -1.0 * pmtNorm.Unit() );}
      else {av3 = av;}
  
      CalcByPosition( fStartPos, pmtVec4, fEnergy, fPathPrecision );    
      if ( !fIsTIR && !fResvHit ) {av4 = fPointOnAV1st; cosThetaAvg += fIncidentVecOnPMT.Dot( -1.0 * pmtNorm.Unit() );}
      else {av4 = av;}
      
    }
    
    else{
      
      av1 = pmtVec1;
      av2 = pmtVec2;
      av3 = pmtVec3;
      av4 = pmtVec4;
      cosThetaAvg = ( pmtVec1.Unit() * zPrime + pmtVec2.Unit() * zPrime +
                      pmtVec3.Unit() * zPrime + pmtVec4.Unit() * zPrime );
      
    }

  }

  else{

    CalcByPosition( fStartPos, pmtVec1 );
    av1 = fPointOnAV1st; cosThetaAvg += fIncidentVecOnPMT.Dot( -1.0 * pmtNorm.Unit() );
    
    
    CalcByPosition( fStartPos, pmtVec2 );
    av2 = fPointOnAV1st; cosThetaAvg += fIncidentVecOnPMT.Dot( -1.0 * pmtNorm.Unit() );
    
    
    CalcByPosition( fStartPos, pmtVec3 );
    av3 = fPointOnAV1st; cosThetaAvg += fIncidentVecOnPMT.Dot( -1.0 * pmtNorm.Unit() );
    
    
    CalcByPosition( fStartPos, pmtVec4 );
    av4 = fPointOnAV1st; cosThetaAvg += fIncidentVecOnPMT.Dot( -1.0 * pmtNorm.Unit() );
    
  }

  *this = tmpStore;

  fCosThetaAvg = cosThetaAvg;
  fCosThetaAvg /= 4.0;

  Double_t angAlpha = ( ( av1 - fStartPos ).Unit() ).Angle( ( ( av2 - fStartPos ).Unit() ) )/2.0;
  Double_t angBeta = ( ( av3 - fStartPos ).Unit() ).Angle( ( ( av4 - fStartPos ).Unit() ) )/2.0;

  fSolidAngle = pi * TMath::Sin( angAlpha ) * TMath::Sin( angBeta );


  if ( GetSolidAngle() < 0 || GetSolidAngle() > 4.0 * pi ){

    warn << "LightPathCalculator: Solid Angle is out of range!\n";
    fSolidAngle = -10.0;

  }

}

//////////////////////////////////////
//////////////////////////////////////

void 
LightPathCalculator::CalculateSolidAnglePolygon( const TVector3& pmtNorm,
                                                 const Int_t nVal )
{
  
  Double_t cosThetaAvg = 0.0;

  TVector3 hypEndPos( 0, 0, 0 );
  CalcByPosition( fStartPos, fEndPos, fEnergy, fPathPrecision );

  TVector3 avMid = fPointOnAV1st;

  LightPathCalculator tmpStore = *this;

  std::vector<TVector3> nVertices;
    
  TVector3 zPrime = pmtNorm;
  TVector3 yPrime = ( zPrime.Orthogonal() ).Unit();
  TVector3 xPrime = ( yPrime.Cross( zPrime ) ).Unit();

  for (Int_t j = 0; j < nVal; j++){

    TVector3 tmpVec = yPrime;
    TVector3 nVec;
    Double_t rotationStep = twopi/ nVal;
    tmpVec.Rotate( j * rotationStep, fEndPos );
    nVec = fEndPos + ( fPMTRadius * tmpVec );
    nVertices.push_back( nVec );

  }

  std::vector<TVector3> nAVPoints;

  if( !fIsTIR && !fResvHit && ( fDistInScint + fDistInAV > 1.0 ) ){

    for (Int_t k = 0; k < nVal; k++){

      CalcByPosition( fStartPos, nVertices[k], fEnergy, fPathPrecision );

      if ( !fIsTIR && !fResvHit ) { nAVPoints.push_back( fPointOnAV1st ); cosThetaAvg += fIncidentVecOnPMT.Dot( -1.0 * pmtNorm );}
      else { nAVPoints.push_back( avMid );}

    }
  }
  
  else{
    for ( Int_t k = 0; k < nVal; k++ ){
      nAVPoints.push_back( nVertices[k] );
    }
  }

  *this = tmpStore;
  fCosThetaAvg = cosThetaAvg;
  fCosThetaAvg /= (Double_t)nVal;

  TVector3 evToAVMid = ( avMid - fStartPos ).Unit();

  for ( Int_t k = 0; k < nVal; k++ ){
    TVector3 vecA, vecB, vecAMid, vecBMid;
    if ( k == ( nVal-1 ) ){
      vecA =  ( nAVPoints[k] - fStartPos ).Unit();
      vecAMid = vecA - evToAVMid;
      vecB =  ( nAVPoints[0] - fStartPos ).Unit();
      vecBMid = vecB - evToAVMid;
    }
    
    else{
      vecA =  ( nAVPoints[k] - fStartPos ).Unit();
      vecAMid = vecA - evToAVMid;
      vecB =  ( nAVPoints[k+1] - fStartPos ).Unit();
      vecBMid = vecB - evToAVMid;
    }
    fSolidAngle += ( 0.5 * ( ( vecBMid.Cross( vecAMid ) ).Mag() ) );
  }

}

//////////////////////////////////////
//////////////////////////////////////

void 
LightPathCalculator::FresnelTRCoeff( const TVector3& dir,
                                     const TVector3& norm,
                                     const Double_t n1,
                                     const Double_t n2,
                                     Double_t& T,
                                     Double_t& R )
{
  
  Double_t cos1 = dir.Dot( norm );
  if ( cos1 < 0.0 ){ cos1 = -1.0 * cos1; }

  const Double_t sin1_2 = 1.0 - ( cos1 * cos1 );
  const std::complex<Double_t> n12 = n1 / n2;
  const std::complex<Double_t> cos2 = std::sqrt( std::complex<Double_t>( 1.0, 0.0 ) - ( n12 * n12 ) * sin1_2 );
  
  const std::complex<Double_t> Ds = ( n2 * cos2 + n1 * cos1 );
  const std::complex<Double_t> Dp = ( n2 * cos1 + n1 * cos2 );

  const Double_t Rsc = std::abs( ( n1 * cos1 - n2 * cos2 ) / Ds );
  const Double_t Rpc = std::abs( ( n1 * cos2 - n2 * cos1 ) / Dp );
  const Double_t Rs = Rsc * Rsc;
  const Double_t Rp = Rpc * Rpc;

  const std::complex<Double_t> Nt = std::complex<Double_t>( 4.0, 0.0 ) * n12 * cos1 * cos2;
  const Double_t Ts = std::abs( n2 * n2 * Nt / ( Ds * Ds ) );
  const Double_t Tp = std::abs( n2 * n2 * Nt / ( Dp * Dp ) );

  const Double_t sFraction = 0.5; //GetSPolarisationFraction( norm, dir, pol );
  if ( sFraction < 0.0 || sFraction > 1.0 ){
    warn << "LightPathCalculator::CalculateFresnelTRCoeff incorrect polarisation calculates as : " << sFraction << "\n";
  }

  T = sFraction * Ts + ( 1.0 - sFraction ) * Tp;
  R = sFraction * Rs + ( 1.0 - sFraction ) * Rp;

  if ( T > 1.1 || R < 0.0 ){
    warn << "LightPathCalculator::CalculateFrsnelTRCoeff, fFresnelRCoeff is too high " << Rs << " " << Rp << " " << R << "\n";
  }

}

//////////////////////////////////////
//////////////////////////////////////

Double_t 
LightPathCalculator::CalculateCosThetaPMT( const Int_t pmtID )
{

  Double_t cosTheta = 0.0;

  const RAT::DU::PMTInfo pmtInfo = RAT::DU::Utility::Get()->GetPMTInfo();
  TVector3 pmtNorm = pmtInfo.GetDirection( pmtID );
  Double_t angle = fIncidentVecOnPMT.Angle( pmtNorm );

  cosTheta = TMath::Cos( angle );
  return cosTheta;

}

//////////////////////////////////////
//////////////////////////////////////

void 
LightPathCalculator::CalculateFresnelTRCoeff()
{

  // Type 'SAW' : Scint -> AV -> Water -> PMT
  if ( fLightPathType == SAW ){
    Double_t T1 = 0.0;
    Double_t R1 = 0.0;

    Double_t T2 = 0.0;
    Double_t R2 = 0.0;
    
    FresnelTRCoeff( GetIncidentVecOn1stSurf(),
                    -1.0 * ( GetPointOnAV1st().Unit() ),
                    fScintRIVal,
                    fAVRIVal,
                    T1, R1 );
    
    FresnelTRCoeff( GetIncidentVecOn2ndSurf(),
                    -1.0 * ( GetPointOnAV2nd().Unit() ),
                    fAVRIVal,
                    fWaterRIVal,
                    T2, R2 );
    
    fFresnelTCoeff = T1 * T2;
    fFresnelRCoeff = R1 * R2;
  }
  
  // Type 'AW' : AV -> Water -> PMT
  else if ( fLightPathType == AW ){
    Double_t T1 = 0.0;
    Double_t R1 = 0.0;
    
    FresnelTRCoeff( GetIncidentVecOn1stSurf(),
                    -1.0 * ( GetPointOnAV1st().Unit() ),
                    fAVRIVal,
                    fWaterRIVal,
                    T1, R1 );
    
    fFresnelTCoeff = T1;
    fFresnelRCoeff = R1;
  }
  
  // Type 'ASAW' : AV -> Scint -> AV -> Water -> PMT
  else if ( fLightPathType == ASAW ){
    Double_t T1 = 0.0;
    Double_t R1 = 0.0;
    
    Double_t T2 = 0.0;
    Double_t R2 = 0.0;
    
    Double_t T3 = 0.0;
    Double_t R3 = 0.0;
    
    FresnelTRCoeff( GetIncidentVecOn1stSurf(),
                    ( GetPointOnAV1st().Unit() ),
                    fAVRIVal,
                    fScintRIVal,
                    T1, R1 );
    
    FresnelTRCoeff( GetIncidentVecOn2ndSurf(),
                    -1.0 * ( GetPointOnAV2nd().Unit() ),
                    fScintRIVal,
                    fAVRIVal,
                    T2, R2 );
    
    FresnelTRCoeff( GetIncidentVecOn3rdSurf(),
                    -1.0 * ( GetPointOnAV3rd().Unit() ),
                    fAVRIVal,
                    fWaterRIVal,
                    T3, R3 );
    
    fFresnelTCoeff = T1 * T2 * T3;
    fFresnelRCoeff = R1 * R2 * R3;
  }
  
  // Type 'WASAW' : Water -> AV -> Scint -> AV -> Water -> PMT
  else if ( fLightPathType == WASAW ){
    Double_t T1 = 0.0;
    Double_t R1 = 0.0;
    
    Double_t T2 = 0.0;
    Double_t R2 = 0.0;
    
    Double_t T3 = 0.0;
    Double_t R3 = 0.0;
    
    Double_t T4 = 0.0;
    Double_t R4 = 0.0;
    
    FresnelTRCoeff( GetIncidentVecOn1stSurf(),
                    ( GetPointOnAV1st().Unit() ),
                    fWaterRIVal,
                    fAVRIVal,
                    T1, R1 );
    
    FresnelTRCoeff( GetIncidentVecOn2ndSurf(),
                    ( GetPointOnAV2nd().Unit() ),
                    fAVRIVal,
                    fScintRIVal,
                    T2, R2 );
    
    FresnelTRCoeff( GetIncidentVecOn3rdSurf(),
                    -1.0 * ( GetPointOnAV3rd().Unit() ),
                    fScintRIVal,
                    fAVRIVal,
                    T3, R3 );
    
    FresnelTRCoeff( GetIncidentVecOn4thSurf(),
                    -1.0 * ( GetPointOnAV4th().Unit() ),
                    fAVRIVal,
                    fWaterRIVal,
                    T4, R4 );
    
    fFresnelTCoeff = T1 * T2 * T3 * T4;
    fFresnelRCoeff = R1 * R2 * R3 * R4;
  }

  // Type 'WAW' : Water -> AV -> Water -> PMT
  else if ( fLightPathType == WAW ){
    Double_t T1 = 0.0;
    Double_t R1 = 0.0;
    
    Double_t T2 = 0.0;
    Double_t R2 = 0.0;
    
    FresnelTRCoeff( GetIncidentVecOn1stSurf(),
                    ( GetPointOnAV1st().Unit() ),
                    fWaterRIVal,
                    fAVRIVal,
                    T1, R1 );
    
    FresnelTRCoeff( GetIncidentVecOn2ndSurf(),
                    -1.0 * ( GetPointOnAV2nd().Unit() ),
                    fAVRIVal,
                    fWaterRIVal,
                    T2, R2 );
    
    fFresnelTCoeff = T1 * T2;
    fFresnelRCoeff = R1 * R2;
  }
  
  // Type 'W' : Water -> PMT
  else if ( fLightPathType == W ){
    fFresnelTCoeff = 1.0;
    fFresnelRCoeff = 0.0;
  }

  // Type 'WRefl' :  Water -> Reflection -> Water -> PMT
  else if ( fLightPathType == WRefl ){
    fFresnelTCoeff = 1.0;
    fFresnelRCoeff = 0.0;
  }

}

} // namespace DU

} // namespace RAT
