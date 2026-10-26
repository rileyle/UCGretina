#ifndef UCG_COMPATIBILITY_HH
#define UCG_COMPATIBILITY_HH 1

// Small compatibility layer for building UCGretina against both Geant4 10 and
// Geant4 11.  Geant4 11 moved the free string helpers into G4StrUtil; older
// releases provide the equivalent member functions on G4String.

#include "G4Types.hh"
#include "G4String.hh"

#ifdef G4V11
// Older Geant4 11 releases ship G4StrUtil.hh; 11.4 declares G4StrUtil in
// G4String.hh instead.
#if __has_include("G4StrUtil.hh")
#include "G4StrUtil.hh"
#endif
#endif

namespace UCG {

inline G4bool contains(const G4String& text, const G4String& part)
{
#ifdef G4V11
  return G4StrUtil::contains(text, part);
#else
  return text.contains(part);
#endif
}

inline G4bool contains(const G4String& text, char c)
{
#ifdef G4V11
  return G4StrUtil::contains(text, c);
#else
  return text.contains(c);
#endif
}

}  // namespace UCG

#endif
