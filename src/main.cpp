#include "common.h"
#include "meshio.h"
#include "hierarchy.h"
#include "timer.h"

template <typename DerivedW>
bool writeDMAT(
        const std::string file_name,
        const Eigen::MatrixBase<DerivedW> & W,
        const bool ascii) {
    FILE * fp = fopen(file_name.c_str(),"wb");
    if(fp == NULL)
    {
        fprintf(stderr,"IOError: writeDMAT() could not open %s...",file_name.c_str());
        return false;
    }
    if(ascii)
    {
        // first line contains number of rows and number of columns
        fprintf(fp,"%d %d\n",(int)W.cols(),(int)W.rows());
        // Loop over columns slowly
        for(int j = 0;j < W.cols();j++)
        {
            // loop over rows (down columns) quickly
            for(int i = 0;i < W.rows();i++)
            {
                fprintf(fp,"%0.17lg\n",(double)W(i,j));
            }
        }
    }else
    {
        // write header for ascii
        fprintf(fp,"0 0\n");
        // first line contains number of rows and number of columns
        fprintf(fp,"%d %d\n",(int)W.cols(),(int)W.rows());
        // reader assumes the binary part is double precision
        Eigen::MatrixXd Wd = W.template cast<double>();
        fwrite(Wd.data(),sizeof(double),Wd.size(),fp);
        //// Loop over columns slowly
        //for(int j = 0;j < W.cols();j++)
        //{
        //  // loop over rows (down columns) quickly
        //  for(int i = 0;i < W.rows();i++)
        //  {
        //    double d = (double)W(i,j);
        //    fwrite(&d,sizeof(double),1,fp);
        //  }
        //}
    }
    fclose(fp);
    return true;
}

void rosy_process(char *input, Float scale, int smooth_iter) {

    MultiResolutionHierarchy mRes;
    Timer<> timer;
    timer.beginStage("data pre-processing");
    mRes.load(input);

    mRes.build();

    mRes.setScale(scale);
    timer.endStage();

    timer.beginStage("rosy optimization");

    int mLevelIterations = 0;
    int mMaxIterations = smooth_iter;
    int mLevel = mRes.levels()-2;

    while (true) {
        mRes.smoothOrientationsTri(mLevel, true, true, true);
        mLevelIterations++;

        if (mLevelIterations >= mMaxIterations) {
            mLevelIterations = 0;
            if (mLevel == 0) {
                break;
            }

            mLevel--;
            mRes.prolongOrientations(mLevel);
        }
    }
    timer.endStage();
    writeDMAT("output.dmat",mRes.mQ[0], true);
}

int main(int argc, char **argv) {
	char batchInput[300] = "/Users/zhjiang/Workspace/robust_hex_dominant_meshing/Nefertiti.obj";
	char batchOutput[300] = "output.obj";
	Float scale = 3;
    uint32_t smooth_iter = 10;
    tbb::task_scheduler_init init(1);

    rosy_process(batchInput, scale, smooth_iter);

	return EXIT_SUCCESS;
}