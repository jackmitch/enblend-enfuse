
#include "deghosting.h"

#include <vigra/imageinfo.hxx>

using namespace vigra;

namespace deghosting {
    
    void Deghosting::loadImages(std::vector<std::string>& newInputFiles) throw(NoImages, BadDimensions) {
        if (newInputFiles.empty())
            throw NoImages();
        const ImageImportInfo firstInfo = ImageImportInfo(newInputFiles[0].c_str());
        const int width = firstInfo.width();
        const int height = firstInfo.height();
        for (unsigned int i = 0; i< newInputFiles.size(); i++) {
            ImageImportInfo tmpInfo = ImageImportInfo(newInputFiles[i].c_str());
            if ((width != tmpInfo.width()) || (height != tmpInfo.height()))
                throw BadDimensions();
        }
        inputFiles = newInputFiles;
    }

    void Deghosting::setFlags(const uint16_t newFlags) {
        flags = newFlags;
    }

    void Deghosting::setDebugFlags(const uint16_t newFlags) {
        debugFlags = newFlags;
    }

    void Deghosting::setIterationNum(const int newIterations) {
        iterations = newIterations;
    }

    void Deghosting::setCameraResponse(EMoR newResponse) {
        response = newResponse;
    }

    void Deghosting::setVerbosity(int newVerbosity) {
        verbosity = newVerbosity;
    }

}
