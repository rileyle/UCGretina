#ifndef cache_h
#define cache_h 1

typedef struct cacheTrajectoryPoint {
    G4double x;
    G4double y;
    G4double z;
    G4double bx;
    G4double by;
    G4double bz;
    G4double t;
} CTP;

typedef struct cacheS800Data {
  G4double ata;
  G4double bta;
  G4double dta;
  G4double yta;
} CS;

#endif
