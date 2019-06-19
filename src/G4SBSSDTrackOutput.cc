#include "G4SBSSDTrackOutput.hh"


G4SBSSDTrackOutput::G4SBSSDTrackOutput(){
  Clear();
}

G4SBSSDTrackOutput::G4SBSSDTrackOutput(G4String name){
  sdname = name;
  Clear();
}

G4SBSSDTrackOutput::~G4SBSSDTrackOutput(){;}

void G4SBSSDTrackOutput::Clear(){
  notracks = 0;
  nptracks = 0;
  nsdtracks = 0;
  
  //"Original track" info:
  otrid.clear();
  omid.clear();
  opid.clear();

  oposx.clear();
  omomx.clear();
  opolx.clear();

  oposy.clear();
  omomy.clear();
  opoly.clear();

  oposz.clear();
  omomz.clear();
  opolz.clear();

  oenergy.clear();
  otime.clear();

  //"Primary track" info:
  ptrid.clear();
  // pmid.clear();
  ppid.clear();

  pposx.clear();
  pmomx.clear();
  ppolx.clear();

  pposy.clear();
  pmomy.clear();
  ppoly.clear();

  pposz.clear();
  pmomz.clear();
  ppolz.clear();
  
  penergy.clear();
  ptime.clear();

  //"SD boundary crossing t3rack info":
  sdtrid.clear();
  sdmid.clear();
  sdpid.clear();

  sdposx.clear();
  sdmomx.clear();
  sdpolx.clear();

  sdposy.clear();
  sdmomy.clear();
  sdpoly.clear();

  sdposz.clear();
  sdmomz.clear();
  sdpolz.clear();
  
  sdenergy.clear();
  sdtime.clear();

  //maps:
  otracklist.clear();
  ptracklist.clear();
  sdtracklist.clear();

  // otrack_hits.clear();
  // ptrack_hits.clear();
  // sdtrack_hits.clear();
}

void G4SBSSDTrackOutput::SetSDname(G4String name){
  sdname = name;
}

void G4SBSSDTrackOutput::ConvertToTreeUnits(){
  //Units in the ROOT Tree are meters, GeV, ns:
  for( int i=0; i<notracks; i++ ){
    oposx[i] /= CLHEP::m;
    oposy[i] /= CLHEP::m;
    oposz[i] /= CLHEP::m;

    omomx[i] /= CLHEP::GeV;
    omomy[i] /= CLHEP::GeV;
    omomz[i] /= CLHEP::GeV;

    oenergy[i] /= CLHEP::GeV;
    otime[i] /= CLHEP::ns;
  }

  for( int i=0; i<nptracks; i++ ){
    pposx[i] /= CLHEP::m;
    pposy[i] /= CLHEP::m;
    pposz[i] /= CLHEP::m;

    pmomx[i] /= CLHEP::GeV;
    pmomy[i] /= CLHEP::GeV;
    pmomz[i] /= CLHEP::GeV;

    penergy[i] /= CLHEP::GeV;
    ptime[i] /= CLHEP::ns;
  }

  for( int i=0; i<nsdtracks; i++ ){
    sdposx[i] /= CLHEP::m;
    sdposy[i] /= CLHEP::m;
    sdposz[i] /= CLHEP::m;

    sdmomx[i] /= CLHEP::GeV;
    sdmomy[i] /= CLHEP::GeV;
    sdmomz[i] /= CLHEP::GeV;

    sdenergy[i] /= CLHEP::GeV;
    sdtime[i] /= CLHEP::ns;
  }
}

G4int G4SBSSDTrackOutput::InsertOriginalTrackInformation( G4Track *aTrack ){

  G4SBSTrackInformation *aTrackInfo = (G4SBSTrackInformation*) aTrack->GetUserInformation();

  int tidtemp = aTrackInfo->GetOriginalTrackID();
  
  std::pair< map<int,int>::iterator, bool > newtrack = otracklist.insert( std::pair<int,int>(tidtemp,notracks) );

  if( newtrack.second ){ //new track found: add its info to the arrays:
    otrid.push_back( tidtemp );
    omid.push_back( aTrackInfo->GetOriginalParentID() );
    opid.push_back( aTrackInfo->GetOriginalDefinition()->GetPDGEncoding() );

    G4ThreeVector postemp = aTrackInfo->GetOriginalPosition();
    G4ThreeVector momtemp = aTrackInfo->GetOriginalMomentum();
    G4ThreeVector poltemp = aTrackInfo->GetOriginalPolarization();

    oposx.push_back( postemp.x() );
    oposy.push_back( postemp.y() );
    oposz.push_back( postemp.z() );

    omomx.push_back( momtemp.x() );
    omomy.push_back( momtemp.y() );
    omomz.push_back( momtemp.z() );

    opolx.push_back( poltemp.x() );
    opoly.push_back( poltemp.y() );
    opolz.push_back( poltemp.z() );

    oenergy.push_back( aTrackInfo->GetOriginalEnergy() );
    otime.push_back( aTrackInfo->GetOriginalTime() );

    notracks++;
  }

  return otracklist[tidtemp]; //this is the position/index in the (unsorted) otrack array 
}

G4int G4SBSSDTrackOutput::InsertPrimaryTrackInformation( G4Track *aTrack ){
  G4SBSTrackInformation *aTrackInfo = (G4SBSTrackInformation*) aTrack->GetUserInformation();

  int tidtemp = aTrackInfo->GetPrimaryTrackID();

  std::pair< map<int,int>::iterator, bool > newtrack = ptracklist.insert( std::pair<int,int>(tidtemp,nptracks) );

  if( newtrack.second ){ //new track found: add its info to the arrays:
    ptrid.push_back( tidtemp );
    //pmid.push_back( aTrackInfo->GetPrimaryParentID() );
    ppid.push_back( aTrackInfo->GetPrimaryDefinition()->GetPDGEncoding() );

    G4ThreeVector postemp = aTrackInfo->GetPrimaryPosition();
    G4ThreeVector momtemp = aTrackInfo->GetPrimaryMomentum();
    G4ThreeVector poltemp = aTrackInfo->GetPrimaryPolarization();
  

    pposx.push_back( postemp.x() );
    pposy.push_back( postemp.y() );
    pposz.push_back( postemp.z() );

    pmomx.push_back( momtemp.x() );
    pmomy.push_back( momtemp.y() );
    pmomz.push_back( momtemp.z() );

    ppolx.push_back( poltemp.x() );
    ppoly.push_back( poltemp.y() );
    ppolz.push_back( poltemp.z() );
    
    penergy.push_back( aTrackInfo->GetPrimaryEnergy() );
    ptime.push_back( aTrackInfo->GetPrimaryTime() );

    nptracks++;
  }

  return ptracklist[tidtemp]; //this is the position/index in the (unsorted) ptrack array 
}

G4int G4SBSSDTrackOutput::InsertSDTrackInformation( G4Track *aTrack ){
  G4SBSTrackInformation *aTrackInfo = (G4SBSTrackInformation*) aTrack->GetUserInformation();

  if( (aTrackInfo->fSDlist).find( sdname ) != (aTrackInfo->fSDlist).end() ){

    //G4cout << "found track information for SD " << sdname << G4endl;
    
    int tidtemp = (aTrackInfo->fSDTrackID)[sdname];
    
    std::pair< std::map<int,int>::iterator, bool > newtrack = sdtracklist.insert( std::pair<int,int>(tidtemp,nsdtracks) );

    if( newtrack.second ){ //new track found: add its info to the arrays:
      sdtrid.push_back( tidtemp );
      sdmid.push_back( (aTrackInfo->fSDParentID)[sdname] );
      sdpid.push_back( (aTrackInfo->fSDDefinition)[sdname]->GetPDGEncoding() );
      
      G4ThreeVector postemp = (aTrackInfo->fSDPosition)[sdname];
      G4ThreeVector momtemp = (aTrackInfo->fSDMomentum)[sdname];
      G4ThreeVector poltemp = (aTrackInfo->fSDPolarization)[sdname];
      

      sdposx.push_back( postemp.x() );
      sdposy.push_back( postemp.y() );
      sdposz.push_back( postemp.z() );

      sdmomx.push_back( momtemp.x() );
      sdmomy.push_back( momtemp.y() );
      sdmomz.push_back( momtemp.z() );

      sdpolx.push_back( poltemp.x() );
      sdpoly.push_back( poltemp.y() );
      sdpolz.push_back( poltemp.z() );
      
      sdenergy.push_back( (aTrackInfo->fSDEnergy)[sdname] );
      sdtime.push_back( (aTrackInfo->fSDTime)[sdname] );

      nsdtracks++;
    }
    return sdtracklist[tidtemp]; //this is the position/index in the (unsorted) sdtrack array 
  } else {
    return -1;
  }
}
