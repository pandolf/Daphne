#ifndef DaphneCommon_h
#define DaphneCommon_h


#include "TPaveText.h"


class DaphneCommon {

 public:

  DaphneCommon();
  ~DaphneCommon();

  static void setStyle();

  static TPaveText* getLabelTop( float lumi );
  static TPaveText* getLabelTopSimulation( float lumi );
  static TPaveText* getLabelCMS( const std::string& text="CMS" );
  static TPaveText* getLabelTop( const std::string& text="CMS Preliminary, #sqrt{s} = 13 TeV" );
  static TPaveText* getLabelTopSimulation( const std::string& text="CMS Simulation, #sqrt{s} = 13 TeV" );

  static std::string getLumiText( float lumi );

  static std::string etaText( float eta );

  static std::vector<int> colors();

};

#endif
