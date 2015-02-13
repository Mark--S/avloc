////////////////////////////////////////////////////////////////////
/// \class RAT::DU::LightPathCalculator
///
/// \brief Calculates the refracted (or straight) light path distances
///        between two points within the detector
///          
/// \author Rob Stainforth <r.stainforth@liv.ac.uk>
///
/// REVISION HISTORY:\n
///     02/2013 : Rob Stainforth - First Revision, new file. \n
///     12/2013 : Rob Stainforth - Updated to include ELLIE reflections \n
///   2014-03-01 : P Jones - refactor as part of ds review. \n
///   2014-04-28 : Rob Stainforth - rename class from LightPath -> LightPathCalculator \n
///
/// \details Returns the refracted path through the scintillator 
/// AV and water of the detector region. Currently requires single 
/// value of the wavelength in MeV - if not specified defaults to
/// 400 nm equivalent.
/// If total internal reflection 
/// is present at one of the material interfaces, then the 
/// straight line calculation of the path is performed instead.
///
/// It can be used in RAT, ROOT or external programs that link to
/// the libRATEvent library.  
///
////////////////////////////////////////////////////////////////////

#ifndef __RAT_DU_LightPathCalculator__
#define __RAT_DU_LightPathCalculator__

#include <RAT/DB.hh>

#include <TObject.h>
#include <TGraph.h>
#include <TVector3.h>

#include <map>
#include <iostream>

namespace RAT 
{
  namespace DU
  {
    
    // Light Path types
    // Type 'SAW' - Scint -> AV -> Water -> PMT
    // Type 'AW' - AV -> Water -> PMT
    // Type 'ASAW' - AV -> Scint -> AV -> Water -> PMT
    // Type 'WASAW' - Water -> AV -> Scint -> AV -> Water -> PMT
    // Type 'WAW' - Water -> AV -> Water -> PMT
    // Type 'W' - Water -> PMT
    // Type 'WRefl' - Water -> Reflection -> PMT (Reflection off of the AV)
    // Type 'Null' - Light Path Uninitialised
    enum eLightPathType { SAW, AW, ASAW, WASAW, WAW, W, WRefl, Null };
    
    class LightPathCalculator : public TObject
    {

    public:
      
      /////////////////////////////////
      ////////     METHODS     ////////
      /////////////////////////////////

      /// Default constructor
      LightPathCalculator() : TObject() { }

      /// Called at the start of a run, loads from the database
      /// Initialise the inner and outer AV radii, 
      /// inner and outer neck radii, PMT Bucket radii,
      /// fill fractions (partial fill) and the refractive indices
      void BeginOfRun();

      /// Initalise the refractive indices from the database
      ///
      /// @param[in] dbTable Link to the RATDB table
      /// @param{in,out] property The TGraph to store the refractive indices in
      void LoadRefractiveIndex( DBLinkPtr dbTable,
                                TGraph& property );

      /// Reset/Initalise all the values of the private member variables. Variables
      /// are set to values for which there is no physical interpretation.
      void ResetValues();

      /// Converts the energy to the equivalent wavelength, remember RAT units are MeV [Energy], mm [Length]
      ///
      /// @param[in] energy The energy to convert
      Double_t EnergyToWavelength( const Double_t energy );

      /// Converts the wavelengths to equivalent energies, remember RAT units are MeV [Energy], mm [Length]
      ///
      /// @param[in] wavelength The wavelength to convert
      Double_t WavelengthToEnergy( const Double_t wavelength );

      /// Use this to calculate the path of the light from 'eventPos' to 'pmtPos'.
      /// Refraction is modelled for values of 'localityVal' > 0.0. If 'localityVal' = 0.0,
      /// then the straight line approximation is used. For refractive calculations, the refraction 
      /// is modelled based on the energy of the light provided (default 3.103125 * 1e-6 MeV = 400nm).
      ///
      /// @param[in] eventPos The starting point of the light path (typically an event position)
      /// @param[in] pmtPos The end point of the PMT (typically a PMT position)
      /// @param[in] energy The photon energy in MeV (defaults to 400 nm -> 3.103125 * 1e-6 MeV)
      /// @param[in] localityVal The accepted tolerance [mm] for how close the path is calculated to the 'pmtPos' (defaults to 0.0 -> Straight Line Calculation)
      void CalcByPosition( const TVector3& eventPos,
                           const TVector3& pmtPos,
                           const Double_t& energyMeV = 3.103125 * 1e-6,
                           const Double_t& localityVal = 0.0 );

      /// Used for partial fill geometry.
      /// Use this to calculate the path of the light from 'eventPos' to 'pmtPos'.
      /// Refraction is modelled for values of 'localityVal' > 0.0. If 'localityVal' = 0.0,
      /// then the straight line approximation is used. For refractive calculations, the refraction 
      /// is modelled based on the energy of the light provided (default 3.103125 * 1e-6 MeV = 400nm).      
      ///
      /// @param[in] eventPos The starting point of the light path (typically an event position)
      /// @param[in] pmtPos The end point of the PMT (typically a PMT position)
      /// @param[in] energy The photon energy in MeV (defaults to 400 nm -> 3.103125 * 1e-6 MeV)
      /// @param[in] localityVal The accepted tolerance [mm] for how close the path is calculated to the 'pmtPos' (defaults to 0.0 -> Straight Line Calculation)
      void CalcByPositionPartial( const TVector3& eventPos,
                                  const TVector3& pmtPos, 
                                  const Double_t& energy = 3.103125 * 1e-6,
                                  const Double_t& localityVal = 0.0 );

      /// Calculate the solid angle for this light path, as subtended at the start point from the PMT
      /// at the end point of the light path.
      /// A call to one of the 'CalcByPosition' functions must be made first before calling this
      ///
      /// @param[in] pmtNorm The PMT-Bucket normal vector (pointing IN, towards the AV)
      /// @param[in] nVal The nVal-sided polygon superimposed onto the PMT bucket for the calculation,
      /// use nVal = 0 for original LOCAS ellipse approximation (quicker). nVal must be > 4.
      void CalculateSolidAngle( const TVector3& pmtNorm,
                                const Int_t nVal );

      /// Calculate the cos of the angle, theta, the light path makes with the bucket face.
      /// This must be called following a call to 'CalcByPosition'
      ///
      /// @param[in] pmtID The ID of the PMT for which the CosTheta value needs to be calculated
      /// @param[out] The value of the CosTheta for the light path incident on this PMT
      Double_t CalculateCosThetaPMT( const Int_t pmtID );

      /// Calculate the total Frensel transmission/reflection coefficients. 
      /// Use the respective 'getters' to return these values.
      void CalculateFresnelTRCoeff();

      /////////////////////////////////
      ////////     GETTERS     ////////
      /////////////////////////////////
      
      /// Return refractive index in scintillator for a given wavelength in MeV
      /// @param[in] energy The wavelength in MeV
      ///
      /// @return The refractive index in the scintillator for this wavelength
      Double_t GetScintRI( const Double_t energy ) const { return fScintRI.Eval( energy ); }

      /// Return refractive index in the upper target (partial fill) for a given wavelength in MeV
      /// @param[in] energy The wavelength in MeV
      ///
      /// @return The refractive index in the upper filled part of the detector for this wavelength
      Double_t GetUpperTargetRI( const Double_t energy ) const { return fUpperTargetRI.Eval( energy ); }

      /// Return refractive index in the lower target (partial fill) for a given wavelength in MeV
      /// @param[in] energy The wavelength in MeV
      ///
      /// @return The refractive index in the lower filled part of the detector for this wavelength
      Double_t GetLowerTargetRI( const Double_t energy ) const { return fLowerTargetRI.Eval( energy ); }
      
      /// Return refractive index in AV for a given wavelength in MeV
      /// @param[in] energy The wavelength in MeV
      ///
      /// @return The refractive index in the AV for this wavelength
      Double_t GetAVRI( const Double_t energy ) const { return fAVRI.Eval( energy ); }
      
      /// Return refractive index in water for a given wavelength in MeV
      /// @param[in] energy The wavelength
      ///
      /// @return The refractive index in the water for this wavelength
      Double_t GetWaterRI( const Double_t energy ) const { return fWaterRI.Eval( energy ); }
      
      /// Return the loop ceiling value 
      /// (i.e. max number of possible iterations to be made for the refracted path calculation)
      ///
      /// @return The maximum number of iterations for the light path calculation
      Double_t GetLoopCeiling() const { return fLoopCeiling; }
      
      /// Return final loop value used by iterative scheme used to calculate the refracted path
      /// (i.e. the value of 'i' for the final, i-th iteration)
      ///
      /// @return The number of iterations made for the light path calculation
      Double_t GetFinalLoopValue() const { return fFinalLoopSize; }
      
      /// Return whether total internal reflection was detected
      /// for the path. If so, then the values of the distances in the scintillator, AV and water are
      /// from straight line path approximations.
      ///
      /// @return TRUE: if the path encountered total internal reflection, FALSE: if not
      Bool_t GetTIR() const { return fIsTIR; }
      
      /// Return whether the given end point was not calculated to within 1000 mm of the calculated end point.
      /// i.e. final path location > localityVal away from PMT position
      ///
      /// @return TRUE: if the end position of the path is far from the required end position FALSE: if not
      Bool_t GetResvHit() const { return fResvHit; }

      /// Return whether the light path enters the neck of the AV
      /// If TRUE: then there are also distances calculated for the scintillator, AV and water from 
      /// passing through the neck, see 'GetDistInNeckScint', GetDistInNeckAV' and 'GetDistInNeckWater'
      ///
      /// @return TRUE: if the light path enters the neck of the AV FALSE: if not
      Bool_t GetXAVNeck() const { return fXAVNeck; }

      /// Return whether the light path calculation used the straight line method
      ///
      /// @return TRUE: if the light path is a straight line FALSE: if not
      Bool_t GetStraightLine() const { return fStraightLine; }

      /// Return whether the ELLIE reflected distances are required for start positions outside the AV (i.e. in the water)
      /// If TRUE: then it is assumed that PMTs surrounding the starting position in the water (most likely a fibre position)
      /// are hit from reflections off of the AV. PMTs on the far side are assumed to have light paths associated
      /// with direct light traveling through the entire detector. See LightPath on docDB for how the different regions of PMTs
      /// surrounding the fibre position are determined to have occurred from reflections or direct light.
      ///
      /// @return Whether the ELLIE reflected distances are required for start positions outside the AV
      Bool_t GetELLIEReflect() const { return fELLIEReflect; }
      
      /// Return the current precision for the path
      /// (the 'localityVal' variable used by CalcByPosition prototype 2). This is essentially the tolerance
      /// of how close the final calculated end position of the light path needs to be to the one passed
      /// as an argument to any of the 'CalcByPosition' calls.
      ///
      /// @return The locality value specified for the current light path calculation
      Double_t GetPathPrecision() const { return fPathPrecision; }

      /// @return The wavelength value used for the light path calculation in nm.
      Double_t GetEnergy() const { return fEnergy; }

      /// @return The fill fraction of the detector
      Double_t GetFillFraction() const { return fFillFraction; }

      /// @return The distance in Scintillator
      Double_t GetDistInScint() const { return fDistInScint; }

      /// @return The distance in the acrylic of the AV
      Double_t GetDistInAV() const { return fDistInAV; }

      /// @return The distance in the water
      Double_t GetDistInWater() const { return fDistInWater; }

      /// @return The distance in the upper target of the detector (above fill - for partial fill geometry)
      Double_t GetDistInUpperTarget() const { return fDistInUpperTarget; }

      /// @return The distance in the lower target of the detector (below fill - for partial fill geometry)
      Double_t GetDistInLowerTarget() const { return fDistInLowerTarget; }

      /// @return The distance in the scintillator in the neck region
      Double_t GetDistInNeckScint() const { return fDistInNeckScint; }

      /// @return The distance in the acrylic in the neck region
      Double_t GetDistInNeckAV() const { return fDistInNeckAV; }

      /// @return The distance in the water following the path going through the neck region
      Double_t GetDistInNeckWater() const { return fDistInNeckWater; }

      /// @return The total distance on the light path
      Double_t GetTotalDist() { return ( fDistInScint + fDistInAV + fDistInWater ); }

      /// @return The total distance of the light path for a partial fill geometry
      Double_t GetTotalDistPartial() { return ( fDistInScint + fDistInAV + fDistInWater ); }

      /// @return The solid angle as calculated by 'CalculateSolidAngle'
      Double_t GetSolidAngle() const { return fSolidAngle; }

      /// @return The average CosTheta of the incident path on the PMT bucket face following the 'CalculateSolidAngle' call
      Double_t GetCosThetaAvg() const { return fCosThetaAvg; }

      /// @return The total Fresnel transmission coefficient for the light path, following a call to 'CalculateFresnelTRCoeff'
      Double_t GetFresnelTCoeff() const { return fFresnelTCoeff; }

      /// @return The total Fresnel reflectivity coefficient for the light path, following a call to 'CalculateFresnelTRCoeff'
      Double_t GetFresnelRCoeff() const { return fFresnelRCoeff; }

      /// @return the light path start point
      TVector3 GetStartPos() const { return fStartPos; }

      /// @return The 'required' light path end point
      TVector3 GetEndPos() const { return fEndPos; }

      /// @return The light path end position as calculated by 'CalcByPosition'
      TVector3 GetLightPathEndPos() const { return fLightPathEndPos; }
      
      /// @return The (unit normalised) incident vector at the PMT bucket, going INTO the PMT
      TVector3 GetIncidentVecOnPMT() const { return fIncidentVecOnPMT; }

      /// @return The (unit normalised) initial light vector from the source position
      TVector3 GetInitialLightVec() const { return fInitialLightVec; }

      /// NOTE: For the below points on AV, depending on the Light Path
      /// type, it may not have 1st, 2nd, 3rd or 4th points on the AV
      /// where the path intersected. This is based on the LightPathType,
      /// see comments under 'fLightPathType' private member declaration
      /// for details.

      /// @return The first point on the AV where the light path intersects
      TVector3 GetPointOnAV1st() const { return fPointOnAV1st; }

      /// @return The second point on the AV where the light path intersects
      TVector3 GetPointOnAV2nd() const { return fPointOnAV2nd; }
      
      /// @return The third point on the AV where the light path intersects
      TVector3 GetPointOnAV3rd() const { return fPointOnAV3rd; }

      /// @return The fourth point on the AV where the light path intersects
      TVector3 GetPointOnAV4th() const { return fPointOnAV4th; }

      /// @return The first point on the AV where the light path intersects
      TVector3 GetPointOnNeck1st() const { return fPointOnNeck1st; }

      /// @return The second point on the AV where the light path intersects
      TVector3 GetPointOnNeck2nd() const { return fPointOnNeck2nd; }

      /// @return The light path type
      std::string GetLightPathType() { return fLightPathTypeMap[ fLightPathType ]; }

      /// @return The incident vector (locally) at the first incident surface on the AV
      TVector3 GetIncidentVecOn1stSurf() const { return ( fPointOnAV1st - fStartPos ).Unit(); }

      /// @return The incident vector (locally) at the second incident surface on the AV
      TVector3 GetIncidentVecOn2ndSurf() const { return ( fPointOnAV2nd - fPointOnAV1st ).Unit(); }

      /// @return The incident vector (locally) at the third incident surface on the AV
      TVector3 GetIncidentVecOn3rdSurf() const { return ( fPointOnAV3rd - fPointOnAV2nd ).Unit(); }

      /// @return The incident vector (locally) at the fourth incident surface on the AV
      TVector3 GetIncidentVecOn4thSurf() const { return ( fPointOnAV4th - fPointOnAV3rd ).Unit(); }

      /// @return The inner AV radius
      Double_t GetAVInnerRadius() const { return fAVInnerRadius; }
      
      /// @return The outer AV radius
      Double_t GetAVOuterRadius() const { return fAVOuterRadius; }

      /// @return The inner AV Neck radius
      Double_t GetNeckInnerRadius() const { return fNeckInnerRadius; }

      /// @return The outer AV neck radius
      Double_t GetNeckOuterRadius() const { return fNeckOuterRadius; }

      /// @return The PMT bucket radius
      Double_t GetPMTRadius() const { return fPMTRadius; }


      /////////////////////////////////
      ////////     SETTERS     ////////
      /////////////////////////////////

      /// For calculations where the event position is in the water region/PSUP
      /// Calculations of the light path distance can be calculated assuming it first 
      /// reflected off of the AV first. In which case SetELLIEReflect must be passed
      /// a TRUE bool. ( Default: False - No reflection off of AV )
      ///
      /// @param[in] reflect TRUE: ELLIE reflected distances required FALSE: Not required
      void SetELLIEReflect( const Bool_t reflect )
      {
        fELLIEReflect = reflect;
      }
      
      //Set the radius of the AV for AVloc
      void SetAVRadius(Double_t newRadius){fAVOuterRadius = newRadius;}
      
      void SetAVOffset(Double_t newOffset){fAVOffset = newOffset;}
      // This ROOT macro adds dictionary methods to this class.
      // The number is 0 as this class is never, and should never be written to disc.
      // It assumes this class has no virtual methods, use ClassDef if change this.
      ClassDefNV( LightPathCalculator, 0 );
    private:
      
      /// Calculate the refracted path. Performs most of the work required to
      /// obtain a refracted path
      ///
      /// @param[in] initOffset The initial unit vector pointing from the starting position for the light path direction
      void PathCalculation( const TVector3& initOffset );  
      
      /// Readjust the initial photon direction. Used if the previous path
      /// does not meet the locality conditions
      ///
      /// @param[in] eventPos The starting point of the light path (typically an event position)
      /// @param[in] pmtPos The end point of the PMT (typically a PMT position)
      /// @param[in] distWater distWater The point in the water where the light path EXITS OUT OF
      /// @param[in] initOffset The initial unit vector pointing from the starting position for the light path direction
      void ReadjustOffset( const TVector3& distWater,
                           TVector3& initOffset );
      
      /// Test the locality conditions for the hypothesised path end point
      /// against the actual 'required' light path end point (usually a PMT position )
      ///
      /// @param[in] pmtPos The end point of the PMT (typically a PMT position)
      /// @param[in] distWater distWater The point in the water where the light path EXITS OUT OF
      /// @param[in] i The final i-th value of the iteration that successfully calculated the light path
      ///
      /// @return TRUE: if the light path is close to the required position FALSE: if not
      Bool_t LocalityCheck( const Int_t i );
      
      /// Calculate the maximum angle between the event position
      /// and the path direction for the path to intersect with the 
      /// sphere of radius 'edgeRadius'
      ///
      /// @param[in] eventPos The starting point of the light path (typically an event position)
      /// @param[in] edgeRadius The radius of the nearest sphere to the event position
      ///
      /// @return Calculate the closest angular displacement of a path close to a surface interface
      Double_t ClosestAngle( const TVector3& eventPos,
                             const Double_t& edgeRadius );

      /// Calculate the maximum allowed angle between the event position
      /// and the PMT position for it to reflect off of the AV
      ///
      /// @param[in] eventPos The starting point of the light path (typically an event position)
      /// @param[in] edgeRadius The radius of the nearest sphere to the event position
      ///
      /// @return Calculate the maximum allowed angle between the event position and the PMT position for it to reflect off of the AV
      Double_t ReflectionAngle( const TVector3& eventPos,
                                const Double_t& edgeRadius );
      
      
      /// Calculate refracted photon vector (unit normalised)
      ///
      /// @param[in] incidentVec The unit vector incident on the surface surface
      /// @param[in] incientSurfVec The unit surface vector of the incident surface
      /// @param[in] incRIndex The incident refractive index
      /// @param[in] refRIndex The refractive index of the refractive media
      ///
      /// @return The refracted unit vector
      TVector3 PathRefraction( const TVector3& incidentVec,
                               const TVector3& incidentSurfVec,
                               const Double_t& incRIndex,
                               const Double_t& refRIndex );
      
      /// Calculate vector from some initial point ('startPos'), with an initial
      /// starting direction, 'startDir' to the edge of a sphere of given radius ('radiusFromCentre')
      ///
      /// @param[in] startPos starting position within (or outside of) a sphere
      /// @param[in] startDir Unit starting direction vector from startPos
      /// @param[in] radiusFromCentre The radius of the spherical edge required to be calculated
      /// @param[in] outside Whether the starting position is outside of the sphere in question
      ///
      /// @return The vector at the spheres edge
      TVector3 VectorToSphereEdge( const TVector3 startPos,
                                   const TVector3 startDir,
                                   const Double_t radiusFromCentre,
                                   const Bool_t outside );

      /// Calculate vector from some initial point ('startPos'), with an initial
      /// starting direction, 'startDir' to the edge of a cylinder of given radius ('cylinderRadius')
      ///
      /// @param[in] startPos starting position within (or outside of) a cylinder
      /// @param[in] startDir Unit starting direction vector from startPos
      /// @param[in] cylinderBaseOrigin The location of the origin of the base of the cylinder
      /// @param[in] cylinderRadius The radius of the cylinder
      ///
      /// @return The vector at the spheres edge
      TVector3 VectorToCylinderEdge( const TVector3& startPos,
                                     const TVector3& startDir,
                                     const TVector3& cylinderBaseOrigin,
                                     const Double_t cylinderRadius );

      /// Calculate the path through the upper and lower regions of a partial fill geometry
      ///
      /// @param[in] enterPos The entry point at the target
      /// @param[in] enterDir The entry direction at the target
      /// @param[in,out] exitPos The exit position of the target
      /// @param[in,out] exitDir The exit direction at the exit point of the target
      void PathThroughTarget( const TVector3& enterPos,
                              const TVector3& enterDir,
                              TVector3& exitPos,
                              TVector3& exitDir );

      /// Calculate the refracted path. Performs most of the work required to
      /// obtained a refracted path
      /// @param[in] initialDir The initial direction of the partial path
      ///
      /// @return The partial path
      void PathCalculationPartial( const TVector3& initialDir );

      /// Set the AV neck variables (i.e. those pertaining to whether the
      /// path entered the neck or not, and if it did, the distance it traveled in the scintillator, acrylic and water ).
      /// This is called for all light paths, which, at some point, leave the AV region
      ///
      /// @param[in] pointOnAV The point on the inner AV where the light path enters
      /// @param[in] dirVec The (unit normalised) direction vector of the light path at the above 'pointOnAV'
      void SetAVNeckInformation( const TVector3& pointOnAV,
                                 const TVector3& dirVec );

      /// Calculates the solid angle using a more rigorous method.
      /// This is reserved for calculating the solid angle for locations
      /// close to the AV (i.e. near AV region, large R ( e.g. R > 5000 mm ), but can be generally used by specifying
      /// an 'nVal' value to the 'CalculateSolidAngle' call.
      ///
      /// @param[in] pmtNorm The PMT bucket normal (unit normalised) pointing INTO the detector, towards the AV
      /// @param[in] nVal The nVal-sided polygon to be super-imposed on the PMT bucket and used as part of the solidangle calculation
      void CalculateSolidAnglePolygon( const TVector3& pmtNorm,
                                       const Int_t nVal );

      /// Calculates the Fresnel transmission/reflectivity coefficients for the light path, assumed 50:50 polarisation
      /// of the light (spherically and linearly)
      ///
      /// @param[in] dir The (unit normalised) direction vector incident on the material surface
      /// @param[in] norm The (unit normalised) surface normal vector pointing in the opposite direction from the incident vector above
      /// @param[in] n1 The incident material refractive index
      /// @param[in] n2 The refracted material refractive index
      /// @param[in,out] T The calculated transmission coefficient at this boundary
      /// @param[in,out] R The calculated reflectivity coefficient at this boundary
      void FresnelTRCoeff( const TVector3& dir,
                           const TVector3& norm,
                           const Double_t n1,
                           const Double_t n2,
                           Double_t& T,
                           Double_t& R );
      
      
      /// Set if Total Internal Reflection is detected
      ///
      /// @param[in] val TRUE: The path calculation encountered total internal reflection FALSE: It didn't
      void SetTIR( const Bool_t val ) { fIsTIR = val; }
      
      /// Set if the calculated path was difficult to resolve.
      /// (i.e. whether the path > 1000 mm away from the end path)
      ///
      /// @param[in] val TRUE: The end calculated position of the path is far from the required positions FALSE: It wasn't
      void SetResvHit( const Bool_t val ) { fResvHit = val; }

      ///////////////////////////////////////////////////
      /////////     PRIVATE MEMBER VARIABLES     ////////
      ///////////////////////////////////////////////////
           

      Double_t fNeckInnerRadius;                                         ///< Radius of the inner neck region
      Double_t fNeckOuterRadius;                                         ///< Radius of the outer neck
      Double_t fAVInnerRadius;                                           ///< Radius of the scint region
      Double_t fAVOuterRadius;                                           ///< Radius of the AV region
      Double_t fPMTRadius;                                               ///< Radius of the PMT bucket
      Double_t fFillZ;                                                   ///< z position of the partial fill
      
      TGraph fScintRI;                                                   ///< Scintillator refractive index TGraph
      TGraph fUpperTargetRI;                                             ///< The 'un-filled' region of the detector refractive index
      TGraph fLowerTargetRI;                                             ///< The 'filled' region of the detector refractive index
      TGraph fAVRI;                                                      ///< AV refractive index TGraph
      TGraph fWaterRI;                                                   ///< Water refractive index TGraph

      Double_t fFillFraction;                                            ///< The fill fraction of the detector (from the bottom of the detector)
      Double_t fLoopCeiling;                                             ///< Iteration Ceiling for algortithm loop
      Double_t fFinalLoopSize;                                           ///< Final loop value which meets locality conditions
      Double_t fPathPrecision;                                           ///< The accepted path proximity/tolerance to the PMT location [mm]

      Double_t fScintRIVal;                                              ///< The value of the scintillator refractive index used for this path
      Double_t fUpperTargetRIVal;                                        ///< The value of the upper target volume index used for this path (partial fill)
      Double_t fLowerTargetRIVal;                                        ///< The value of the lower target volume index used for this path (partial fill)
      Double_t fAVRIVal;                                                 ///< The value of the AV refractive index used for this path
      Double_t fWaterRIVal;                                              ///< The value of the water refractive index used for this path
      
      TVector3 fIncidentVecOnPMT;                                        ///< Final light path direction (unit normalised)
      TVector3 fInitialLightVec;                                         ///< Initial light path direction (unit normalised)

      TVector3 fStartPos;                                                ///< Start position of the light path
      TVector3 fEndPos;                                                  ///< Required end position of the light path
      TVector3 fLightPathEndPos;                                         ///< Calculated end position of the light path

      Double_t  fAVOffset;                                               ///< Offset of the AV from the origin (Used for AVloc)
      

      Bool_t fIsTIR;                                                     ///< TRUE: Total Internal Reflection encountered FALSE: It wasn't
      Bool_t fResvHit;                                                   ///< TRUE: Difficult path to resolve and calculate FALSE: It wasn't
      Bool_t fXAVNeck;                                                   ///< TRUE: Path entered neck region FALSE: It didn't
      Bool_t fELLIEReflect;                                              ///< TRUE: Reflected distances in water off of AV on PMTs near starting position required
                                                                         ///< FALSE: Reflected distances not required
      Bool_t fStraightLine;                                              ///< TRUE: Light Path is a straight line approximation FALSE: It isn't


      // Note: Depending on the light path type (see 'fLightPathType'), the path may
      // intersect the AV/Neck once, twice or three (times a lady?...sorry)
      // ...or four times

      TVector3 fPointOnAV1st;                                            ///< Point on AV where light path first hits the AV
      TVector3 fPointOnAV2nd;                                            ///< Point on AV where light path hits the AV a second time
      TVector3 fPointOnAV3rd;                                            ///< Point on AV where light path hits the AV a third time
      TVector3 fPointOnAV4th;                                            ///< Point on AV where light path hits the AV a fourth time

      TVector3 fPointOnNeck1st;                                          ///< Point on the Neck where the light path hits a first time
      TVector3 fPointOnNeck2nd;                                          ///< Point on the Neck where the light path hits a second time

      eLightPathType fLightPathType;                                     ///< Light path type, based on what regions of the detector the path enters

      std::map< eLightPathType, std::string > fLightPathTypeMap;         ///< Map containing a descriptor for the light path type

      Double_t fDistInScint;                                             ///< Distance in the scintillator region
      Double_t fDistInUpperTarget;                                       ///< Distance in the upper target region (partial fill geometry)
      Double_t fDistInLowerTarget;                                       ///< Distance in the lower target region (partial fill geometry)
      Double_t fDistInAV;                                                ///< Distance in the acrylic region of the AV
      Double_t fDistInWater;                                             ///< Distance in the water region
      Double_t fDistInNeckScint;                                         ///< Distance through the scintillator region in the neck (if GetXAVNeck() = TRUE)
      Double_t fDistInNeckAV;                                            ///< Distance through the acrylic of the AV region in the neck (if GetXAVNeck() = TRUE)
      Double_t fDistInNeckWater;                                         ///< Distance through the water region in the neck (if GetXAVNeck() = TRUE)

      Double_t fEnergy;                                                  ///< The value of the wavelength in MeV

      Double_t fSolidAngle;                                              ///< The solid angle subtended by the PMT for this light path
      Double_t fCosThetaAvg;                                             ///< Average incident angle on the PMT for this path.
                                                                         ///< This is only calculated after a call to CalculateSolidAngle
    
      Double_t fFresnelTCoeff;                                           ///< The combined Fresnel TRANSMISSION coefficient for this path
      Double_t fFresnelRCoeff;                                           ///< The combined Fresnel REFLECTIVITY coefficient for this path  

    };
    
  } // namespace DU
  
} // namespace RAT

#endif
