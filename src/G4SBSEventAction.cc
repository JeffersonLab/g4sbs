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

#include "TMatrix.h"
#include "TVector.h"

#include "G4SBSIO.hh"

#define MAXHIT 2000

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
	printf("Event %8d\r", ev->GetEventID());
	fflush(stdout);
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

  int i;

  tr_t trdata;
  cal_t caldata;

  caldata.bcndata = 0;
  caldata.hcndata = 0;

  if(bbcalHC) {
      if( bbcalHC->entries() > 0 ){
	  hasbb = true;

	  double xsum = 0.0;
	  double ysum = 0.0;
	  double esum = 0.0;

	  caldata.bcndata = bbcalHC->entries();
	  
	  if( bbcalHC->entries() > MAXHITDATA ){
	      G4cerr << "WARNING:  Number of hits exceeds array length - truncating" << G4endl;
	      caldata.bcndata = MAXHITDATA;
	  }


	  for( i = 0; i < caldata.bcndata; i++ ){
	      xsum += (*bbcalHC)[i]->GetPos().x()*(*bbcalHC)[i]->GetEdep();
	      ysum += (*bbcalHC)[i]->GetPos().y()*(*bbcalHC)[i]->GetEdep();
	      esum += (*bbcalHC)[i]->GetEdep();

	      caldata.bcx[i] = (*bbcalHC)[i]->GetPos().x()/cm;
	      caldata.bcy[i] = (*bbcalHC)[i]->GetPos().y()/cm;
	      caldata.bcz[i] = (*bbcalHC)[i]->GetPos().z()/cm;
	      caldata.bce[i] = (*bbcalHC)[i]->GetEdep()/GeV;

	      caldata.bcvx[i] = (*bbcalHC)[i]->GetVertex().x()/cm;
	      caldata.bcvy[i] = (*bbcalHC)[i]->GetVertex().y()/cm;
	      caldata.bcvz[i] = (*bbcalHC)[i]->GetVertex().z()/cm;

	      caldata.bcpid[i] = (*bbcalHC)[i]->GetPID();
	      caldata.bcmid[i] = (*bbcalHC)[i]->GetMID();
	      caldata.bctrid[i] = (*bbcalHC)[i]->GetTrID();
	  }
	  
	  //  Energy weighted sum
	  trdata.bcx = xsum/esum/cm; 
	  trdata.bcy = ysum/esum/cm;
      }
  }
  
  if(hcalHC) {
      if( hcalHC->entries() > 0 ){
	  hashcal = true;

	  double xsum = 0.0;
	  double ysum = 0.0;

	  double xlsum = 0.0;
	  double ylsum = 0.0;
	  double zlsum = 0.0;

	  double esum = 0.0;

	  caldata.hcndata = hcalHC->entries();

	  if( hcalHC->entries() > MAXHITDATA ){
	      G4cerr << "WARNING:  Number of hits exceeds array length - truncating" << G4endl;
	      caldata.hcndata = MAXHITDATA;
	  }

	  for( i = 0; i < caldata.hcndata; i++ ){
	      xsum += (*hcalHC)[i]->GetPos().x()*(*hcalHC)[i]->GetEdep();
	      ysum += (*hcalHC)[i]->GetPos().y()*(*hcalHC)[i]->GetEdep();
	      xlsum += (*hcalHC)[i]->GetLabPos().x()*(*hcalHC)[i]->GetEdep();
	      ylsum += (*hcalHC)[i]->GetLabPos().y()*(*hcalHC)[i]->GetEdep();
	      zlsum += (*hcalHC)[i]->GetLabPos().y()*(*hcalHC)[i]->GetEdep();
	      esum += (*hcalHC)[i]->GetEdep();

	      if( (*hcalHC)[i]->GetMID() == 0 ){
		  trdata.hct = (*hcalHC)[i]->GetTime()/ns + CLHEP::RandGauss::shoot(0.0, fevgen->GetToFres());
	      }

	      caldata.hcx[i] = (*hcalHC)[i]->GetPos().x()/cm;
	      caldata.hcy[i] = (*hcalHC)[i]->GetPos().y()/cm;
	      caldata.hcz[i] = (*hcalHC)[i]->GetPos().z()/cm;
	      caldata.hce[i] = (*hcalHC)[i]->GetEdep()/GeV;

	      caldata.hcvx[i] = (*hcalHC)[i]->GetVertex().x()/cm;
	      caldata.hcvy[i] = (*hcalHC)[i]->GetVertex().y()/cm;
	      caldata.hcvz[i] = (*hcalHC)[i]->GetVertex().z()/cm;

	      caldata.hcpid[i] = (*hcalHC)[i]->GetPID();
	      caldata.hcmid[i] = (*hcalHC)[i]->GetMID();
	      caldata.hctrid[i] = (*hcalHC)[i]->GetTrID();
	  }


      	  trdata.hcx = xsum/esum/cm;
	  trdata.hcy = ysum/esum/cm;
	  trdata.hclx = xlsum/esum/cm;
	  trdata.hcly = ylsum/esum/cm;
	  trdata.hclz = zlsum/esum/cm;

	  G4ThreeVector avglab = G4ThreeVector( trdata.hclx, trdata.hcly, trdata.hclz);


	  // Calculate expected time of flight
	  G4ThreeVector q3m = fevgen->GetBeamP()-fevgen->GetElectronP();
	  G4ThreeVector path = avglab-fevgen->GetV();
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
  /*
  if( !hasbb && !hashcal ){
      return;
  }
  */


  trdata.hcal = hashcal;
  trdata.bb   = hasbb;

  trdata.x = trdata.y = trdata.xp = trdata.yp = -1e9;
  trdata.gemtr = 0;
  int idx, nhit, gid;

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
	  hitdata.trid[nhit] = (*gemHC)[idx]->GetTrID();
	  hitdata.pid[nhit] = (*gemHC)[idx]->GetPID();
	  hitdata.mid[nhit] = (*gemHC)[idx]->GetMID();
	  hitdata.x[nhit] =  ly[nhit]/m;
	  hitdata.y[nhit] = -lx[nhit]/m;
	  hitdata.z[nhit] = lz[nhit]/m;
	  hitdata.p[nhit] = (*gemHC)[idx]->GetMom()/GeV;
	  hitdata.edep[nhit] = (*gemHC)[idx]->GetEdep()/GeV;

	  hitdata.vx[nhit] =  (*gemHC)[idx]->GetVertex().getX()/m;
	  hitdata.vy[nhit] =  (*gemHC)[idx]->GetVertex().getY()/m;
	  hitdata.vz[nhit] =  (*gemHC)[idx]->GetVertex().getZ()/m;

	  hitdata.tx[nhit] =  (*gemHC)[idx]->GetPos().getY()/m;
	  hitdata.ty[nhit] = -(*gemHC)[idx]->GetPos().getX()/m;
	  hitdata.txp[nhit] =  (*gemHC)[idx]->GetYp();
	  hitdata.typ[nhit] = -(*gemHC)[idx]->GetXp();

//	  printf("GEM HIT %d (%f) %f %f\n", (*gemHC)[idx]->GetGEMID(), lz[nhit]/cm, lx[nhit]/cm, ly[nhit]/cm );
	  nhit++;

	  if( nhit == MAXHITDATA ){
	      G4cerr << "WARNING:  Number of hits exceeds array length - truncating" << G4endl;
	      break;
	  }
      }

      if( nhit >= 3 ){

	  // Perform fitting

	  TMatrix mymat(nhit*2, 4);
	  TMatrix sigmamat(nhit*2, nhit*2);
	  // Go x0,y0,x1,y1...

	  TVector hitv(nhit*2);
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

	  TMatrix mytrans = mymat;

	  mytrans.T();

	  TMatrix alpha = mytrans*sigmamat*sigmamat*mymat;


	  if( alpha.Determinant() != 0.0 ){
	      alpha.Invert();

	      TMatrix fitmat = alpha*mytrans;

	      TVector track = fitmat*sigmamat*sigmamat*hitv;

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
  fIO->SetCalData(caldata);
  fIO->SetHitData(hitdata);
  fIO->FillTree();

  return;
}



