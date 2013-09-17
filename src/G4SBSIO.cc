#include "G4SBSIO.hh"

#include <TFile.h>
#include <TTree.h>
#include <TClonesArray.h>

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
    delete fTree;
    fTree = NULL;
}


void G4SBSIO::InitializeTree(){
    if( fTree ){ delete fTree; }

    fTree = new TTree("T", "Geant4 SBS Simulation");
    fTree->Branch("ev", &evdata, "count/D:rate/D:solang/D:sigma/D:W2/D:xbj/D:Q2/D:th/D:ph/D:Aperp/D:Apar/D:vx/D:vy/D:vz/D:ep/D:np/D:epx/D:epy/D:epz/D:npx/D:npy/D:npz/D:nth/D:nph/D:pmperp/D:pmpar/D:pmparsm/D:nucl/I:fnucl/I");
    fTree->Branch("tr", &trdata, "x/D:y/D:xp/D:yp/D:tx/D:ty/D:txp/D:typ/D:hcal/I:bb/I:gemtr/I:hcx/D:hcy/D:bcx/D:bcy/D:hct/D:hctex/D:hclx/D:hcly/D:hclz/D:hcdang/D");
    fTree->Branch("gen", &gendata, "thbb/D:thhcal/D:dbb/D:dhcal/D:Ebeam/D");
    // Tedious, but we want dynamically scaled
    fTree->Branch("ht.ndata", &hitdata.ndata, "ht.ndata/I");
    fTree->Branch("ht.gid", &hitdata.gid, "ht.gid[ht.ndata]/I");
    fTree->Branch("ht.x", &hitdata.x, "ht.x[ht.ndata]/D");
    fTree->Branch("ht.y", &hitdata.y, "ht.y[ht.ndata]/D");
    fTree->Branch("ht.z", &hitdata.z, "ht.z[ht.ndata]/D");
    fTree->Branch("ht.dx", &hitdata.dx, "ht.dx[ht.ndata]/D");
    fTree->Branch("ht.dy", &hitdata.dy, "ht.dy[ht.ndata]/D");

    fTree->Branch("ht.tx", &hitdata.tx, "ht.tx[ht.ndata]/D");
    fTree->Branch("ht.ty", &hitdata.ty, "ht.ty[ht.ndata]/D");
    fTree->Branch("ht.txp", &hitdata.txp, "ht.txp[ht.ndata]/D");
    fTree->Branch("ht.typ", &hitdata.typ, "ht.typ[ht.ndata]/D");

    fTree->Branch("hc.ndata", &caldata.hcndata, "hc.ndata/I");
    fTree->Branch("hc.x", &caldata.hcx, "hc.x[hc.ndata]/D");
    fTree->Branch("hc.y", &caldata.hcy, "hc.y[hc.ndata]/D");
    fTree->Branch("hc.e", &caldata.hce, "hc.e[hc.ndata]/D");
    fTree->Branch("hc.trid", &caldata.hctrid, "hc.trid[hc.ndata]/D");
    fTree->Branch("hc.mid", &caldata.hcmid, "hc.mid[hc.ndata]/D");
    fTree->Branch("hc.pid", &caldata.hcpid, "hc.pid[hc.ndata]/D");

    fTree->Branch("bc.ndata", &caldata.bcndata, "bc.ndata/I");
    fTree->Branch("bc.x", &caldata.bcx, "bc.x[bc.ndata]/D");
    fTree->Branch("bc.y", &caldata.bcy, "bc.y[bc.ndata]/D");
    fTree->Branch("bc.e", &caldata.bce, "bc.e[hc.ndata]/D");
    fTree->Branch("bc.trid", &caldata.bctrid, "bc.trid[hc.ndata]/D");
    fTree->Branch("bc.mid", &caldata.bcmid, "bc.mid[hc.ndata]/D");
    fTree->Branch("bc.pid", &caldata.bcpid, "bc.pid[hc.ndata]/D");

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
    fFile = new TFile(fFilename, "RECREATE");
    fFile->cd();
    fTree->Write("T", TObject::kOverwrite);
    fFile->Close();
    delete fFile;
    fFile = NULL;

    return;
}




