void rootlogon(){
  gSystem->AddIncludePath(" -I/home/rachel/geant4/g4sbs_mtpc/g4sbs/include");
  gSystem->AddIncludePath(" -I/home/rachel/geant4/g4sbs_mtpc/g4sbs/include");
    FileStat_t buf;
    if( !gSystem->GetPathInfo("/home/rachel/geant4/g4sbs_mtpc/g4sbs/libg4sbsroot.so", buf) ){
	gSystem->Load("/home/rachel/geant4/g4sbs_mtpc/g4sbs/libg4sbsroot.so" ) ;
    }
    if( !gSystem->GetPathInfo("/home/rachel/geant4/g4sbs_mtpc/g4sbs/libg4sbsroot.dylib", buf) ){
	gSystem->Load("/home/rachel/geant4/g4sbs_mtpc/g4sbs/libg4sbsroot.dylib" ) ;
    }
}

