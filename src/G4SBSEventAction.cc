#include "G4SBSEventAction.hh"
#include "G4SBSEventGen.hh"
#include "G4SBSCalHit.hh"
#include "G4SBSGEMHit.hh"

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4VHitsCollection.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4VVisManager.hh"
#include "G4SDManager.hh"
#include "G4UImanager.hh"
#include "G4ios.hh"

#include "CLHEP/Matrix/Matrix.h"
#include "CLHEP/Matrix/Vector.h"

#include "G4SBSIO.hh"

#define MAXHIT 20

G4SBSEventAction::G4SBSEventAction()
{
    fGEMres = 70.0*um;

    // Load up resolution file if it exists

    int idx;

    for( idx = 0; idx < __MAXGEM; idx++ ){
	fGEMsigma[idx] = 1.0;
    }


}

G4SBSEventAction::~G4SBSEventAction()
{;}


void G4SBSEventAction::LoadSigmas(const char filename[] ){
    printf("Loading sigmas file %s\n", filename);
    FILE *sigmafile = fopen(filename, "r");

    int idx;
    if( !sigmafile ){ 
	printf("WARNING:  file %s not found.  Initializing to 1.0\n", filename);

	for( idx = 0; idx < __MAXGEM; idx++ ){
	    fGEMsigma[idx] = 1.0;
	}
	return;
    }

    int nread = 2;
    int gemidx;

    idx = 0;
    while( nread == 2 ){
	nread = fscanf(sigmafile, "%d%lf", &gemidx, &fGEMsigma[idx]);
	if( nread==2 && idx+1 == gemidx){ idx++; }
    }

    return;
}


void G4SBSEventAction::BeginOfEventAction(const G4Event*ev) {
   if( (ev->GetEventID()%1000)==0 ){
	printf("Event %8d\n", ev->GetEventID());
    }

    return;
}

void G4SBSEventAction::EndOfEventAction(const G4Event* evt )
{
  G4SDManager * SDman = G4SDManager::GetSDMpointer();
  G4String colNam;
  gemCollID   = SDman->GetCollectionID(colNam="GEMcol");
  hcalCollID  = SDman->GetCollectionID(colNam="HCALcol");
  bbcalCollID = SDman->GetCollectionID(colNam="BBCalcol");
  
   //G4cout << ">>> Event " << evt->GetEventID() << G4endl;
  
  G4HCofThisEvent * HCE = evt->GetHCofThisEvent();
  G4SBSCalHitsCollection* bbcalHC    = 0;
  G4SBSCalHitsCollection* hcalHC = 0;
  G4SBSGEMHitsCollection* gemHC = 0;
  if(HCE)
    {
      gemHC   = (G4SBSGEMHitsCollection*)(HCE->GetHC(gemCollID));
      hcalHC  = (G4SBSCalHitsCollection*)(HCE->GetHC(hcalCollID));
      bbcalHC = (G4SBSCalHitsCollection*)(HCE->GetHC(bbcalCollID));
    }

  bool hasbb   = false;
  bool hashcal = false;

  tr_t trdata;

  if(bbcalHC) {
      if( bbcalHC->entries() > 0 ){
	  hasbb = true;
	  trdata.bcx = (*bbcalHC)[0]->GetPos().x()/cm;
	  trdata.bcy = (*bbcalHC)[0]->GetPos().y()/cm;
      }
  }
  
  if(hcalHC) {
      if( hcalHC->entries() > 0 ){
	  hashcal = true;
      	  trdata.hcx = (*hcalHC)[0]->GetPos().x()/cm;
	  trdata.hcy = (*hcalHC)[0]->GetPos().y()/cm;
	  trdata.hct = (*hcalHC)[0]->GetTime()/ns + CLHEP::RandGauss::shoot(0.0, fevgen->GetToFres());

	  trdata.hclx = (*hcalHC)[0]->GetLabPos().x()/cm;
	  trdata.hcly = (*hcalHC)[0]->GetLabPos().y()/cm;
	  trdata.hclz = (*hcalHC)[0]->GetLabPos().z()/cm;

	  // Calculate expected time of flight
	  G4ThreeVector q3m = fevgen->GetBeamP()-fevgen->GetElectronP();
	  G4ThreeVector path = (*hcalHC)[0]->GetLabPos()-fevgen->GetV();
	  double hcd = path.mag();

	  trdata.hctex = hcd/(q3m.mag()*(0.3*m/ns)/sqrt(q3m.mag()*q3m.mag()+proton_mass_c2*proton_mass_c2))/ns;
	  // Angular difference between q and reconstructed vector
	  double cosang = q3m.unit()*path.unit();
	  if( cosang > 1.0 ){ cosang = 1.0; } //  Apparent numerical problems in this dot product
	  trdata.hcdang = acos(cosang);
      }
  }

  // If we don't have something in both arms end
  // and don't fill
  if( !hasbb && !hashcal ){
      return;
  }


  trdata.hcal = hashcal;
  trdata.bb   = hasbb;

  trdata.x = trdata.y = trdata.xp = trdata.yp = -1e9;
  trdata.gemtr = 0;
  int idx, i, j, ierr, nhit, gid;

  int    map = 0;
  // Just use 4 GEMs for now
  double lx[MAXHIT], ly[MAXHIT], lz[MAXHIT];

  double txp, typ, tx, ty;

  hit_t hitdata;

  hitdata.ndata = 0;

  if(gemHC) {
      // Need at least three hits to draw a line

      nhit = 0;
      for( idx = 0; idx < gemHC->entries() && idx < MAXHIT; idx++ ){
	  gid = (*gemHC)[idx]->GetGEMID();

	  if( gid == 0 ) continue;


	  if( gid == 1 ){
	      //  Project back to z = 0 plane
	      tx  =  (*gemHC)[idx]->GetPos().getX() - (*gemHC)[idx]->GetXp()*(*gemHC)[idx]->GetPos().getZ();
	      ty  =  (*gemHC)[idx]->GetPos().getY() - (*gemHC)[idx]->GetYp()*(*gemHC)[idx]->GetPos().getZ();
	      txp =  (*gemHC)[idx]->GetXp();
	      typ =  (*gemHC)[idx]->GetYp();
	  };

	  map |= (1 << (*gemHC)[idx]->GetGEMID());

	  // Smear by resolution
	  lx[nhit] = (*gemHC)[idx]->GetPos().getX() + CLHEP::RandGauss::shoot(0.0, fGEMres);
	  ly[nhit] = (*gemHC)[idx]->GetPos().getY() + CLHEP::RandGauss::shoot(0.0, fGEMres);
	  lz[nhit] = (*gemHC)[idx]->GetPos().getZ();

	  hitdata.gid[nhit] = (*gemHC)[idx]->GetGEMID();
	  hitdata.x[nhit] =  ly[nhit]/m;
	  hitdata.y[nhit] = -lx[nhit]/m;
	  hitdata.z[nhit] = lz[nhit]/m;

	  hitdata.tx[nhit] =  (*gemHC)[idx]->GetPos().getY()/m;
	  hitdata.ty[nhit] = -(*gemHC)[idx]->GetPos().getX()/m;
	  hitdata.txp[nhit] =  (*gemHC)[idx]->GetYp();
	  hitdata.typ[nhit] = -(*gemHC)[idx]->GetXp();

//	  printf("GEM HIT %d (%f) %f %f\n", (*gemHC)[idx]->GetGEMID(), lz[nhit]/cm, lx[nhit]/cm, ly[nhit]/cm );
	  nhit++;
      }

      if( nhit >= 3 ){

	  // Perform fitting

	  CLHEP::HepMatrix mymat(nhit*2, 4);
	  CLHEP::HepMatrix sigmamat(nhit*2, nhit*2);
	  // Go x0,y0,x1,y1...

	  CLHEP::HepVector hitv(nhit*2);
	  for( i = 0; i < nhit; i++ ){
	      mymat[2*i][0] = 1.0;
	      mymat[2*i][1] = lz[i];
	      mymat[2*i][2] = 0.0;
	      mymat[2*i][3] = 0.0;

	      mymat[2*i+1][0] = 0.0;
	      mymat[2*i+1][1] = 0.0;
	      mymat[2*i+1][2] = 1.0;
	      mymat[2*i+1][3] = lz[i];

	      hitv[2*i]   = lx[i];
	      hitv[2*i+1] = ly[i];

	      sigmamat[2*i][2*i]     = 1.0/fGEMsigma[2*i];
	      sigmamat[2*i+1][2*i+1] = 1.0/fGEMsigma[2*i+1];
	  }

	  CLHEP::HepMatrix alpha = mymat.T()*sigmamat*sigmamat*mymat;
	  ierr = -10;
	  alpha.invert(ierr);

	  if( ierr == 0 ){
	      CLHEP::HepMatrix fitmat = alpha*mymat.T();

	      CLHEP::HepVector track = fitmat*sigmamat*sigmamat*hitv;

	      // Switch to "BigBite coordinates"
	      // Larger momentum is correlated to larger x
	      // larger angle is correlated with smaller y
	      trdata.x  = track[2]/m;
	      trdata.xp = track[3];
	      trdata.y  = -track[0]/m;
	      trdata.yp = -track[1];

	      trdata.tx  = ty/m;
	      trdata.txp = typ;
	      trdata.ty  = -tx/m;
	      trdata.typ = -txp;

	      trdata.gemtr = 1;

//	      printf("Reconstructed track = (%f, %f) (%f, %f)\n\n", trdata.x, trdata.y, trdata.xp, trdata.yp);
	      for( i = 0; i < nhit; i++ ){
		  double dx = track[0] + track[1]*lz[i] - lx[i];
		  double dy = track[2] + track[3]*lz[i] - ly[i];
//		  printf("Position deviations = %f um\n", sqrt(dx*dx+dy*dy)/um);

		  hitdata.dx[i] =  dy/m;
		  hitdata.dy[i] = -dx/m;
	      }
	      hitdata.ndata = nhit;
	  }
      }
  }

  fIO->SetTrackData(trdata);
  fIO->SetHitData(hitdata);
  fIO->FillTree();

  return;
}



