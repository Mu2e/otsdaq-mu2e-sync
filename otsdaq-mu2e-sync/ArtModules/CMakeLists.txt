
cet_build_plugin(SyncAna art::module LIBRARIES REG
art_root_io::TFileService_service
artdaq-core-mu2e::Data
canvas::canvas
otsdaq-mu2e::ArtModules
otsdaq::NetworkUtilities
Offline::Mu2eUtilities
Offline::RecoDataProducts
ROOT::Hist
ROOT::Tree
ROOT::Core
ROOT::RIO
ROOT::Gui
)

install_headers()
install_source()
install_fhicl(SUBDIRS fcl)
