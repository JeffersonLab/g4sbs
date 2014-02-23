#include "TObjString.h"

#include <TFile.h>
#include <TTree.h>
#include <TH2F.h>
#include <TClonesArray.h>

#include "G4SBSGlobalField.hh"
#include "G4SBSRun.hh"
#include "G4SBSIO.hh"
#include <assert.h>

G4SBSIO::G4SBSIO(){
    fTree = NULL;
    InitializeTree();
    // Default filename
    strcpy(fFilename, "g4sbsout.root");
    fFile = NULL;

    gendata.Ebeam = 2.2;
    gendata.thbb = 40.0*deg;
    gendata.dbb = 1.5;
    gendata.thhcal = 39.4*deg;
    gendata.dhcal = 17.0;
}

G4SBSIO::~G4SBSIO(){
    if( fTree ){delete fTree;}
    fTree = NULL;

    if( fFile ){delete fFile;}
    fFile = NULL;
}


void G4SBSIO::InitializeTree(){
    if( fFile ){
	fFile->Close();
	delete fFile;
    }

    fFile = new TFile(fFilename, "RECREATE");

    if( fTree ){ delete fTree; }

    fTree = new TTree("T", "Geant4 SBS Simulation");
    fTree->Branch("ev", &evdata, "count/D:rate/D:solang/D:sigma/D:W2/D:xbj/D:Q2/D:th/D:ph/D:Aperp/D:Apar/D:vx/D:vy/D:vz/D:ep/D:np/D:epx/D:epy/D:epz/D:npx/D:npy/D:npz/D:nth/D:nph/D:pmperp/D:pmpar/D:pmparsm/D:z/D:phperp/D:phih/D:MX/D:nucl/I:fnucl/I:hadr/I");
    fTree->Branch("tr", &trdata, "x/D:y/D:xp/D:yp/D:tx/D:ty/D:txp/D:typ/D:hcal/I:bb/I:gemtr/I:hcx/D:hcy/D:bcx/D:bcy/D:hct/D:hctex/D:hclx/D:hcly/D:hclz/D:hcdang/D");
    fTree->Branch("gen", &gendata, "thbb/D:thhcal/D:dbb/D:dhcal/D:Ebeam/D");
    // Tedious, but we want dynamically scaled
    fTree->Branch("ht.ndata", &hitdata.ndata, "ht.ndata/I");
    fTree->Branch("ht.gid", &hitdata.gid, "ht.gid[ht.ndata]/I");
    fTree->Branch("ht.x", &hitdata.x, "ht.x[ht.ndata]/D");
    fTree->Branch("ht.y", &hitdata.y, "ht.y[ht.ndata]/D");
    fTree->Branch("ht.z", &hitdata.z, "ht.z[ht.ndata]/D");
    fTree->Branch("ht.vx", &hitdata.vx, "ht.vx[ht.ndata]/D");
    fTree->Branch("ht.vy", &hitdata.vy, "ht.vy[ht.ndata]/D");
    fTree->Branch("ht.vz", &hitdata.vz, "ht.vz[ht.ndata]/D");
    fTree->Branch("ht.dx", &hitdata.dx, "ht.dx[ht.ndata]/D");
    fTree->Branch("ht.dy", &hitdata.dy, "ht.dy[ht.ndata]/D");
    fTree->Branch("ht.p", &hitdata.p, "ht.p[ht.ndata]/D");

    fTree->Branch("ht.trid", &hitdata.trid, "ht.trid[ht.ndata]/I");
    fTree->Branch("ht.pid", &hitdata.pid, "ht.pid[ht.ndata]/I");
    fTree->Branch("ht.mid", &hitdata.mid, "ht.mid[ht.ndata]/I");
    fTree->Branch("ht.edep", &hitdata.edep, "ht.edep[ht.ndata]/D");

    fTree->Branch("ht.tx", &hitdata.tx, "ht.tx[ht.ndata]/D");
    fTree->Branch("ht.ty", &hitdata.ty, "ht.ty[ht.ndata]/D");
    fTree->Branch("ht.txp", &hitdata.txp, "ht.txp[ht.ndata]/D");
    fTree->Branch("ht.typ", &hitdata.typ, "ht.typ[ht.ndata]/D");

    fTree->Branch("hc.ndata", &caldata.hcndata, "hc.ndata/I");
    fTree->Branch("hc.x", &caldata.hcx, "hc.x[hc.ndata]/D");
    fTree->Branch("hc.y", &caldata.hcy, "hc.y[hc.ndata]/D");
    fTree->Branch("hc.z", &caldata.hcz, "hc.z[hc.ndata]/D");
    fTree->Branch("hc.e", &caldata.hce, "hc.e[hc.ndata]/D");
    fTree->Branch("hc.vx", &caldata.hcvx, "hc.vx[hc.ndata]/D");
    fTree->Branch("hc.vy", &caldata.hcvy, "hc.vy[hc.ndata]/D");
    fTree->Branch("hc.vz", &caldata.hcvz, "hc.vz[hc.ndata]/D");
    fTree->Branch("hc.trid", &caldata.hctrid, "hc.trid[hc.ndata]/I");
    fTree->Branch("hc.mid", &caldata.hcmid, "hc.mid[hc.ndata]/I");
    fTree->Branch("hc.pid", &caldata.hcpid, "hc.pid[hc.ndata]/I");

    fTree->Branch("bc.ndata", &caldata.bcndata, "bc.ndata/I");
    fTree->Branch("bc.x", &caldata.bcx, "bc.x[bc.ndata]/D");
    fTree->Branch("bc.y", &caldata.bcy, "bc.y[bc.ndata]/D");
    fTree->Branch("bc.z", &caldata.bcz, "bc.z[bc.ndata]/D");
    fTree->Branch("bc.e", &caldata.bce, "bc.e[bc.ndata]/D");
    fTree->Branch("bc.vx", &caldata.bcvx, "bc.vx[bc.ndata]/D");
    fTree->Branch("bc.vy", &caldata.bcvy, "bc.vy[bc.ndata]/D");
    fTree->Branch("bc.vz", &caldata.bcvz, "bc.vz[bc.ndata]/D");
    fTree->Branch("bc.trid", &caldata.bctrid, "bc.trid[bc.ndata]/I");
    fTree->Branch("bc.mid", &caldata.bcmid, "bc.mid[bc.ndata]/I");
    fTree->Branch("bc.pid", &caldata.bcpid, "bc.pid[bc.ndata]/I");

    //Declare RICH-related branches of the tree:
    //richdata are stored as STL vectors (basically dynamically sized arrays). Newer ROOT versions know how to handle this, older may not.
    fTree->Branch("RICH_nhits", &(richdata.nhits_RICH), "nhits_RICH/I");
    fTree->Branch("RICH_pmt", &(richdata.PMTnumber) );
    fTree->Branch("RICH_row", &(richdata.row) );
    fTree->Branch("RICH_col", &(richdata.col) );
    fTree->Branch("RICH_nphe", &(richdata.NumPhotoelectrons) );
    fTree->Branch("RICH_tavg", &(richdata.Time_avg) );
    fTree->Branch("RICH_trms", &(richdata.Time_rms) );
    fTree->Branch("RICH_mID", &(richdata.mTrackNo) );
    fTree->Branch("RICH_vol", &(richdata.volume_flag) );
    fTree->Branch("RICH_xhit", &(richdata.xhit) );
    fTree->Branch("RICH_yhit", &(richdata.yhit) );
    fTree->Branch("RICH_zhit", &(richdata.zhit) );
    fTree->Branch("RICH_pxhit", &(richdata.pxhit) );
    fTree->Branch("RICH_pyhit", &(richdata.pyhit) );
    fTree->Branch("RICH_pzhit", &(richdata.pzhit) );
    fTree->Branch("RICH_vxhit", &(richdata.pvx) );
    fTree->Branch("RICH_vyhit", &(richdata.pvy) );
    fTree->Branch("RICH_vzhit", &(richdata.pvz) );
    fTree->Branch("RICH_vpxhit", &(richdata.ppx) );
    fTree->Branch("RICH_vpyhit", &(richdata.ppy) );
    fTree->Branch("RICH_vpzhit", &(richdata.ppz) );
    fTree->Branch("RICH_ntracks", &(richdata.ntracks_RICH), "ntracks_RICH/I");
    fTree->Branch("RICH_mPID", &(richdata.mPID) );
    fTree->Branch("RICH_mvx", &(richdata.mvx) );
    fTree->Branch("RICH_mvy", &(richdata.mvy) );
    fTree->Branch("RICH_mvz", &(richdata.mvz) );
    fTree->Branch("RICH_mpx", &(richdata.mpx) );
    fTree->Branch("RICH_mpy", &(richdata.mpy) );
    fTree->Branch("RICH_mpz", &(richdata.mpz) );

    return;
}

void G4SBSIO::FillTree(){
    if( !fTree ){ 
	fprintf(stderr, "Error %s: %s line %d - Trying to fill non-existant tree\n", __PRETTY_FUNCTION__, __FILE__, __LINE__ );
	return; 
    }

    fTree->Fill();
}

void G4SBSIO::WriteTree(){
    assert( fFile );
    assert( fTree );
    if( !fFile->IsOpen() ){
	G4cerr << "ERROR: " << __FILE__ << " line " << __LINE__ << ": TFile not open" << G4endl;
	exit(1);
    }

    fFile->cd();
    fTree->Write("T", TObject::kOverwrite);

    G4SBSRun::GetRun()->GetData()->Write("run_data", TObject::kOverwrite);

    // Produce and write out field map graphics
    fGlobalField->DebugField();
    for( vector<TH2F *>::iterator it = fGlobalField->fFieldPlots.begin(); it!= fGlobalField->fFieldPlots.end(); it++ ){
	(*it)->Write((*it)->GetName(), TObject::kOverwrite );
    }

    fTree->ResetBranchAddresses();
    delete fTree;
    fTree = NULL;

    fFile->Close();
    delete fFile;
    fFile = NULL;

    return;
}




