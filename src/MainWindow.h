#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include "TimeFrequencyView.h"
#include "sbsms.h"
#include "../JuceLibraryCode/JuceHeader.h"

using namespace _sbsms_;

class ContentComponent;

class MainWindow 
: public DocumentWindow
{
public:
  MainWindow();
  ~MainWindow();
	void openFile();
  void getSaveFile();
	void openAudio(File &name);
	MenuBarModel *menuBar;
  void closeButtonPressed() override;  

	ContentComponent *contentComponent;
};


#endif
