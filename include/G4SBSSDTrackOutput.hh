#ifndef G4SBSSDTrackOutput_h
#define G4SBSSDTrackOutput_h 1


#include "G4Track.hh"
#include "G4SBSTrackInformation.hh"
#include "G4String.hh"
#include "TVector3.h"

#include <vector>
#include <set>
#include <map>

using namespace std;

//Class to efficiently store "Original track", "Primary track", and "SD boundary track" information
//for each sensitive detector: In conjunction with this routine 

class G4SBSSDTrackOutput {

public:
  G4SBSSDTrackOutput();
  G4SBSSDTrackOutput(G4String name);
  ~G4SBSSDTrackOutput();

  void Clear();
  void SetSDname(G4String name);

  //Returns index in the track output arrays:
  G4int InsertOriginalTrackInformation( G4Track *aTrack );
  G4int InsertPrimaryTrackInformation( G4Track *aTrack );
  G4int InsertSDTrackInformation( G4Track *aTrack );

  void Merge( G4SBSSDTrackOutput & );
  
  void ConvertToTreeUnits();
  
  G4String sdname;

  //key value is track ID; mapped value is index in relevant track array:
  map<int,int> otracklist;
  map<int,int> ptracklist;

  //sdtracklist has to be done differently since the same SD track can cross multiple SD boundary volumes,
  //But at the same time we only want to record each G4 Track ID exactly once for each SD boundary volume it crosses:

  //So let's modify this structure to map<int,map<G4String,int> >,
  //where the int key is the G4 track ID, and G4String key is an SD name, and the mapped value is the index in the
  //sd track array:
  map<int,map<G4String,int> > sdtracklist;

  // //Hit lists mapped by track index:
  // map<int,vector<int> > otrack_hits; //list of SD "hits" associated with each otrack
  // map<int,vector<int> > ptrack_hits; //list of SD "hits" associated with each ptrack
  // map<int,vector<int> > sdtrack_hits; //list of SD "hits" associated with each sdtrack
  
  //"Original track" info:
  int notracks;
  vector<int> otrid, omid, opid, ompid;
  //vector<TVector3> opos, omom, opol;
  vector<double> oposx,oposy,oposz,omomx,omomy,omomz,opolx,opoly,opolz;
  vector<double> oenergy, otime;

  //"Primary track" info:
  int nptracks;
  vector<int> ptrid, ppid;
  //vector<TVector3> ppos, pmom, ppol;
  vector<double> pposx,pposy,pposz,pmomx,pmomy,pmomz,ppolx,ppoly,ppolz;
  vector<double> penergy, ptime;

  //SD boundary crossing info:
  int nsdtracks;
  vector<int> sdtrid, sdmid, sdpid, sdmpid;
  //vector<TVector3> sdpos, sdmom, sdpol;
  vector<double> sdposx,sdposy,sdposz,sdmomx,sdmomy,sdmomz,sdpolx,sdpoly,sdpolz;
  vector<double> sdenergy, sdtime;
  vector<double> sdvx,sdvy,sdvz,sdvnx,sdvny,sdvnz,sdEkin;

};

#endif
