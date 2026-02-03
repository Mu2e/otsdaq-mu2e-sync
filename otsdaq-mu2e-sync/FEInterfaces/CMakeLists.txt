include(otsdaq::FEInterface)

cet_make_library(
  LIBRARY_NAME ROCMobileSyncInterface
  SOURCE ROCMobileSyncInterface_interface.cc
  LIBRARIES PUBLIC
  otsdaq-mu2e::ROCCore
  otsdaq-mu2e::ROCPolarFireCoreInterface
  )

cet_build_plugin(ROCMobileSyncInterface
  otsdaq::FEInterface
  LIBRARIES REG
  otsdaq-mu2e::ROCCore
  otsdaq-mu2e::ROCPolarFireCoreInterface
  otsdaq::Macros
  )
 
install_headers()
install_source()
