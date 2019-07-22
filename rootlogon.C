void rootlogon(){
  //Check for a pre-existing rootlogon.C file in the current working directory and/or
  //the current user's home directory:

  FILE *ftest = fopen( "./rootlogon.C","r");

  if( ftest != NULL ){
    gROOT->ProcessLine(".x ./rootlogon.C");
    fclose( ftest );
  }

  ftest = fopen("~/rootlogon.C","r");

  if( ftest != NULL ){
    gROOT->ProcessLine(".x ~/rootlogon.C");
    fclose( ftest );
  }
  
  //gSystem->AddIncludePath(" -I${PROJECT_SOURCE_DIR}/include");
  gSystem->AddIncludePath(" -I${CMAKE_INSTALL_PREFIX}/include");
  gSystem->AddIncludePath(" -I${PROJECT_BINARY_DIR}/include");
  
  FileStat_t buf;
  if( !gSystem->GetPathInfo("${PROJECT_BINARY_DIR}/libg4sbsroot.so", buf) ){
    gSystem->Load("${PROJECT_BINARY_DIR}/libg4sbsroot.so" ) ;
  }
  if( !gSystem->GetPathInfo("${PROJECT_BINARY_DIR}/libg4sbsroot.dylib", buf) ){
    gSystem->Load("${PROJECT_BINARY_DIR}/libg4sbsroot.dylib" ) ;
  }
}

