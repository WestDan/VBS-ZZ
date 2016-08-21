#include "xAODRootAccess/Init.h"
#include "SampleHandler/SampleHandler.h"
#include "SampleHandler/ToolsDiscovery.h"
#include "SampleHandler/DiskListLocal.h"
#include "EventLoop/Job.h"
#include "EventLoop/DirectDriver.h"
#include "EventLoopGrid/PrunDriver.h"

#include <PATInterfaces/SystematicCode.h>
#include "PATInterfaces/CorrectionCode.h"
#include "PATInterfaces/SystematicsUtil.h"

#include "MyAnalysis/MyxAODAnalysis.h"

#include <iostream>     // std::cout
#include <fstream> 

using namespace std;

int main( int argc, char* argv[] ) {

   // Take the submit directory from the input if provided:
   std::string submitDir = "submitDir";
   if( argc > 1 ) submitDir = argv[ 1 ];


   ifstream dataFile(argv[2], ios::in); 
   if(!dataFile) 
     { cout << "ERROR: data file cannot be open!" << endl; return -1; }

   string dataset;
   getline(dataFile, dataset);
   
   //cout << dataset << endl;
   
   string data1;
   if(argc>=4)
     data1 = argv[3];
   
   // Set up the job for xAOD access:
   xAOD::Init().ignore();
   
   // Construct the samples to run on:
   SH::SampleHandler sh;
   SH::DiskListLocal list( dataset.c_str() );
   //SH::scanDir( sh, "/afs/cern.ch/atlas/project/PAT/xAODs/r5591/" );
   //SH::scanDir( sh, "/afs/cern.ch/user/c/cgeng/work/data/" );
   //SH::scanDir( sh, dataset.c_str() );
   //SH::scanDir( sh, list, "AOD.01512664._00000[1-2].pool.root.1" );
   //SH::scanDir( sh, list, "AOD.01512664._000001.pool.root.1" );
  
   if(argc>=4) { 
     SH::scanDir( sh, list, data1 );
   }
   else
     SH::scanDir( sh, list, "DAOD_HIGG2D1.06607565._000001.pool.root.1" );
   
   // Set the name of the input TTree. It's always "CollectionTree"
   // for xAOD files.
   sh.setMetaString( "nc_tree", "CollectionTree" );
   
   // Print what we found:
   sh.print();
   
   // Create an EventLoop job:
   EL::Job job;
   job.sampleHandler( sh );
   
   // Add our analysis to the job:
   MyxAODAnalysis* alg = new MyxAODAnalysis("physics");
 
   job.algsAdd( alg );
   
   // Run the job using the local/direct driver:
   EL::DirectDriver driver;
   driver.submit( job, submitDir );
   
   return 0;
}
