void rootlogon(){
  gSystem->AddIncludePath(" -I/home/obrecht/Documents/GEn/gen_dev_rhel/g4sbs/include");
  gSystem->AddIncludePath(" -I/home/obrecht/Documents/GEn/gen_dev_rhel/g4sbs/include");
    FileStat_t buf;
    if( !gSystem->GetPathInfo("/home/obrecht/Documents/GEn/gen_dev_rhel/g4sbs/libg4sbsroot.so", buf) ){
	gSystem->Load("/home/obrecht/Documents/GEn/gen_dev_rhel/g4sbs/libg4sbsroot.so" ) ;
    }
    if( !gSystem->GetPathInfo("/home/obrecht/Documents/GEn/gen_dev_rhel/g4sbs/libg4sbsroot.dylib", buf) ){
	gSystem->Load("/home/obrecht/Documents/GEn/gen_dev_rhel/g4sbs/libg4sbsroot.dylib" ) ;
    }
}

