void rootlogon(){
  //Check for a pre-existing rootlogon.C file in the current working directory and/or
  //the current user's home directory:

  TString cwd = gSystem->GetWorkingDirectory();
  TString username = gSystem->Getenv("USER");
  TString homedir = gSystem->GetHomeDirectory(username.Data());

  TString fnametemp = homedir + "/rootlogon.C";

  if( homedir != cwd ){
  
    FILE *ftest = fopen( fnametemp.Data() ,"r");

  // if( ftest != NULL ){
  //   gROOT->ProcessLine(".x ./rootlogon.C");
  //   fclose( ftest );
  // }

  
  
  // ftest = fopen("~/rootlogon.C","r");
  //For the "run in build directory" approach, only search the home directory, because otherwise
  //this will create an infinite loop when 

    if( ftest != NULL ){
      cout << "found rootlogon.C in current user's home directory, executing..." << endl;

      TString cmd = ".x "+fnametemp;
      
      gROOT->ProcessLine(cmd.Data());
      fclose( ftest );
    }
  }
  //gSystem->AddIncludePath(" -I${PROJECT_SOURCE_DIR}/include");
  //gSystem->AddIncludePath(" -I${CMAKE_INSTALL_PREFIX}/include");
  gSystem->AddIncludePath(" -I${PROJECT_BINARY_DIR}/include");

  TString libname = "${PROJECT_BINARY_DIR}/libg4sbsroot.so";

  SysInfo_t sysinfo;
  gSystem->GetSysInfo( &sysinfo );

  TString OSname = sysinfo.fOS;

  if( OSname == "Darwin" ){
    libname.ReplaceAll(".so",".dylib");
  }
  
  FileStat_t buf;
  if( !gSystem->GetPathInfo(libname.Data(), buf) ){
    gSystem->Load( libname.Data() ) ;
  }
}

