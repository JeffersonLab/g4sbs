void rootlogon(){
    FileStat_t buf;
    if( !gSystem->GetPathInfo("libg4sbsroot.so", buf) ){
	gSystem->Load("libg4sbsroot.so" ) ;
    }
    if( !gSystem->GetPathInfo("libg4sbsroot.dylib", buf) ){
	gSystem->Load("libg4sbsroot.dylib" ) ;
    }
}

