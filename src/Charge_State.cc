#include "Charge_State.hh"


Charge_State::Charge_State()
{
  Charge=0;
  UnReactedFraction=1.;
  ReactedFraction=1.;
  setKEu=30*MeV;
  useSetKEu=true;
}

Charge_State::~Charge_State()
{
  ;
}

