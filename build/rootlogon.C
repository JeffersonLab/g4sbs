void rootlogon(){
  gSystem->AddIncludePath(" -I/home/meriem/g4sbs/include");
  gSystem->AddIncludePath(" -I/home/meriem/g4sbs/include");
    FileStat_t buf;
    if( !gSystem->GetPathInfo("/home/meriem/g4sbs/libg4sbsroot.so", buf) ){
	gSystem->Load("/home/meriem/g4sbs/libg4sbsroot.so" ) ;
    }
    if( !gSystem->GetPathInfo("/home/meriem/g4sbs/libg4sbsroot.dylib", buf) ){
	gSystem->Load("/home/meriem/g4sbs/libg4sbsroot.dylib" ) ;
    }
}

