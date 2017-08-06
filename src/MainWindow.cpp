#include "MainWindow.h"
#include "../JuceLibraryCode/JuceHeader.h"
#include "import.h"
#include "sbsms.h"
#include "convert.h"
#include "pcm.h"

class MainMenuModel : public MenuBarModel 
{
public:
  MainMenuModel(MainWindow *main) : main(main) {}

  MainWindow *main;

  StringArray getMenuBarNames()  {
    const char* const names[] = {"File"};
    return StringArray ((const char**) names,1);
  }

  PopupMenu getMenuForIndex(int menuIndex, const String &menuName) {
    PopupMenu menu;
    switch(menuIndex) {
    case 0:
      menu.addItem(1,"Open...");
      menu.addItem(2,"Open Recent...");
      break;
    }
    return menu;
  }

  void menuItemSelected(int menuItemID, int menuIndex) {
    switch(menuIndex) {
    case 0:
      switch(menuIndex) {
      case 0:
        main->openFile();
        break;
      case 1:
        break;
      }
      break;
    }
  }
};


class ContentComponent : public Component 
{
public:
	ContentComponent() {
    tfView = NULL;
  }
    
  ~ContentComponent() {
    delete tfView;
  }

  void resized() {
    Rectangle<int> r(getLocalBounds());
    if(tfView) tfView->setBounds(r);
  }
	
	void setSBSMS(SBSMS *sbsms) {
		if(tfView) {
      delete tfView;
    }
    tfView = new TimeFrequencyView(sbsms);
		addAndMakeVisible(tfView);
    resized();
	}

	TimeFrequencyView *tfView;	
};

MainWindow :: MainWindow() : DocumentWindow("SBSMSer",Colours::lightgrey,allButtons,true)
{
	setUsingNativeTitleBar (true);
  setResizable (true, false);
  setResizeLimits (400, 400, 10000, 10000);
  menuBar = new MainMenuModel(this);
#ifdef __APPLE__
  MenuBarModel::setMacMainMenu(menuBar,NULL,"Open Recent...");
#else
  setMenuBar(menuBar);
#endif
  setBounds ((int) (0.05f * getParentWidth()),
             (int) (0.05f * getParentHeight()),
             jmax (850, (int) (0.9f * getParentWidth())),
             jmax (600, (int) (0.9f * getParentHeight())));

	contentComponent = new ContentComponent();
  setContentOwned (contentComponent, false);
	setVisible (true);

  StringArray args = JUCEApplicationBase::getCommandLineParameterArray();
  if(args.size() >= 1) {
    File file(args[0]);
    openAudio(file);
    if(args.size() >= 5) {
      contentComponent->tfView->setView(atoi(args[1].getCharPointer()),atoi(args[2].getCharPointer()),atof(args[3].getCharPointer()),atof(args[4].getCharPointer()));
    }
  }
}

MainWindow :: ~MainWindow()
{
#ifdef __APPLE__
  MenuBarModel::setMacMainMenu(NULL,NULL);
#else
  setMenuBar(NULL);
#endif
  delete menuBar;
  clearContentComponent();
}

void MainWindow :: openAudio(File &filename)
{  
  const char *in = filename.getFullPathName().getCharPointer();
  //sbsms_convert(in,"bla.wav",true,true,progressCB,NULL,0.5,0.5,1.0,1.0,1.0);
  SBSMS *sbsms = createCachedSBSMS(in);
  contentComponent->setSBSMS(sbsms);
  addKeyListener(contentComponent->tfView);
}

void MainWindow :: openFile()
{
	File lastDirectory = File::getCurrentWorkingDirectory().getFullPathName();
  FileChooser myChooser ("Open file...",
                         lastDirectory,
                         "*.wav;*.aif;*.aiff;*.mp3");
if(myChooser.browseForFileToOpen()) {
  File file = myChooser.getResult();
  if (file.existsAsFile()) {
    if(file.hasFileExtension(".sbsms")) {


    }   
  }
 }
}

void MainWindow::closeButtonPressed()
{
  JUCEApplication::getInstance()->systemRequestedQuit();
}
