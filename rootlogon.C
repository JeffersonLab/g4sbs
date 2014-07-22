void rootlogon(){
  gSystem->AddIncludePath(" -I${PROJECT_SOURCE_DIR}/include");
  gSystem->AddIncludePath(" -I${PROJECT_BINARY_DIR}/include");
    FileStat_t buf;
    if( !gSystem->GetPathInfo("${PROJECT_BINARY_DIR}/libg4sbsroot.so", buf) ){
	gSystem->Load("${PROJECT_BINARY_DIR}/libg4sbsroot.so" ) ;
    }
    if( !gSystem->GetPathInfo("${PROJECT_BINARY_DIR}/libg4sbsroot.dylib", buf) ){
	gSystem->Load("${PROJECT_BINARY_DIR}/libg4sbsroot.dylib" ) ;
    }
}

