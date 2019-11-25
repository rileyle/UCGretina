
#include "TrackerGammaHit.hh"
#include "G4UnitsTable.hh"
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"

G4Allocator<TrackerGammaHit> TrackerGammaHitAllocator;


TrackerGammaHit::TrackerGammaHit() {}



TrackerGammaHit::~TrackerGammaHit() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TrackerGammaHit::TrackerGammaHit(const TrackerGammaHit& right)
  : G4VHit()
{
  trackID        = right.trackID;
  particleID     = right.particleID;
  process        = right.process;
  parentTrackID  = right.parentTrackID;
  creatorProcess = right.creatorProcess;
  detNumb        = right.detNumb;
  segNumb        = right.segNumb;
  edep           = right.edep;
  ke             = right.ke;
  pos            = right.pos;
  trackOrigin    = right.trackOrigin;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const TrackerGammaHit& TrackerGammaHit::operator=(const TrackerGammaHit& right)
{
  trackID        = right.trackID;
  particleID     = right.particleID;
  process        = right.process;
  parentTrackID  = right.parentTrackID;
  creatorProcess = right.creatorProcess;
  detNumb        = right.detNumb;
  segNumb        = right.segNumb;
  edep           = right.edep;
  ke             = right.ke;
  pos            = right.pos;
  trackOrigin    = right.trackOrigin;
  return *this;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int TrackerGammaHit::operator==(const TrackerGammaHit& right) const
{
  return (this==&right) ? 1 : 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackerGammaHit::Draw()
{
  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
  if(pVVisManager)
  {
    G4Circle circle(pos);
    circle.SetScreenSize(4);
    circle.SetFillStyle(G4Circle::filled);
    G4Colour colour(1.,1.,0.);
    G4VisAttributes attribs(colour);
    circle.SetVisAttributes(attribs);
    pVVisManager->Draw(circle);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackerGammaHit::Print()
{

  G4cout << std::setw(3)  << std::right << trackID 
	 << std::setw(12) << particleID
	 << std::setw(12) << process
	 << std::setw(3) << parentTrackID
	 << std::setw(20) << creatorProcess
	 << std::setw(4)  << detNumb
	 << std::setw(3)  << segNumb
	 << std::fixed << std::setprecision(2) << std::setw(10) << std::right
	 << edep/keV
 	 << std::setw(10) << std::right
	 << ke/keV
 	 << std::setw(10) << std::right
         << pos.getX()/mm
	 << std::setw(10) << std::right
         << pos.getY()/mm
	 << std::setw(10) << std::right
         << pos.getZ()/mm
 	 << std::setw(10) << std::right
         << trackOrigin.getX()/mm
	 << std::setw(10) << std::right
         << trackOrigin.getY()/mm
	 << std::setw(10) << std::right
         << trackOrigin.getZ()/mm
         << G4endl;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

