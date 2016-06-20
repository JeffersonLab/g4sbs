#ifndef __G4SBSTEXTFILE_HH
#define __G4SBSTEXTFILE_HH

#define __STRLEN 1024

#include "TObject.h"

class G4SBSTextFile : public TObject {
    public:
	 G4SBSTextFile();
	 G4SBSTextFile(const G4SBSTextFile &);
	 const G4SBSTextFile& operator=(const G4SBSTextFile &);
	 G4SBSTextFile(const char *);
	~G4SBSTextFile();

	 void copyFileIn(const char *);

	void Print(Option_t *) const;
	void Print() const;

	const char *GetFilename(){ return fFilename; }
	unsigned long long int GetBufferSize(){ return fBufferSize; }
	
	void Recreate(const char *fn = NULL, bool clobber = false);
	void RecreateInDir(const char *path, bool clobber = false);

	char* GetBuffer(){ return fBuffer;}

    private:
	int fFilenameSize;
	char *fFilename;

	unsigned long long int fBufferSize;
	char *fBuffer;

	const char *GetBaseFile(const char *fp = NULL);

	ClassDef(G4SBSTextFile, 1);
};

#endif//__G4SBSTEXTFILE_HH
