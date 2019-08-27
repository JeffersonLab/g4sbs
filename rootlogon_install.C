void rootlogon(){
  //Check for a pre-existing rootlogon.C file in the current working directory and/or
  //the current user's home directory:

  TString cwd = gSystem->GetWorkingDirectory();
  TString username = gSystem->Getenv("USER");
  TString homedir = gSystem->GetHomeDirectory(username.Data());

  //check current working directory first:
  TString fnametemp = cwd + "/rootlogon.C";
  
  FILE *ftest = fopen( fnametemp.Data(),"r");

  if( ftest != NULL ){
    cout << "found rootlogon.C in current working directory, executing..." << endl;
    TString cmd = ".x " + fnametemp;
    gROOT->ProcessLine(cmd.Data());
    fclose( ftest );
  }

  if( homedir != cwd ){ //also check in user's home directory:
    fnametemp = homedir + "/rootlogon.C";
    
    ftest = fopen(fnametemp.Data(),"r");

    if( ftest != NULL ){
      cout << "found rootlogon.C in current user's home directory, executing..." << endl;
      TString cmd = ".x " + fnametemp;
      
      gROOT->ProcessLine(cmd.Data());
      fclose( ftest );
    }
  }
  
  //gSystem->AddIncludePath(" -I${PROJECT_SOURCE_DIR}/include");
  //gSystem->AddIncludePath(" -I${CMAKE_INSTALL_PREFIX}/include");
  gSystem->AddIncludePath(" -I${CMAKE_INSTALL_PREFIX}/include");

  TString libname = "${CMAKE_INSTALL_PREFIX}/lib/libg4sbsroot.so";

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

