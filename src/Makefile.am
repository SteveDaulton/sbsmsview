bin_PROGRAMS = SBSMSView
EXTRA_programs = SBSMSView

FRAMEWORKS = \
	 -framework Accelerate \
	 -framework AudioToolbox \
	 -framework Carbon \
	 -framework Cocoa \
	 -framework CoreAudio \
	 -framework CoreMIDI \
	 -framework DiscRecording \
	 -framework IOKit \
	 -framework OpenGL \
	 -framework QTKit \
	 -framework QuartzCore \
	 -framework QuickTime \
	 -framework WebKit \
	 -framework Accelerate \
	 -framework AudioToolbox \
	 -framework Carbon \
	 -framework Cocoa \
	 -framework CoreAudio \
	 -framework CoreMIDI \
	 -framework DiscRecording \
	 -framework IOKit \
	 -framework OpenGL \
	 -framework QTKit \
	 -framework QuartzCore \
	 -framework QuickTime \
	 -framework WebKit

SBSMSView_LDADD = @SNDFILE_LIBS@ @MAD_LIBS@ @AVCODEC_LIBS@ @AVFORMAT_LIBS@ @AVUTIL_LIBS@
SBSMSView_LDFLAGS = $(os_ldflags)  $(FRAMEWORKS)
SBSMSView_CXXFLAGS = -I../include -I/usr/local/include -I../JuceLibraryCode -I../JuceLibraryCode/modules @SBSMS_CFLAGS@ @AVFORMAT_CFLAGS@ @AVCODEC_CFLAGS@ @AVUTIL_CFLAGS@ -DDEBUG -DPNG_NO_CONFIG_H -g -D__STDC_CONSTANT_MACROS -D__STDC_FORMAT_MACROS -ObjC++ -std=c++11

gl_matrix_SOURCES = gl-matrix/gl-matrix.h \
	gl-matrix/gl-matrix.h \
	gl-matrix/mat3.c \
	gl-matrix/mat4.c \
	gl-matrix/quat.c \
	gl-matrix/vec3.c \
	gl-matrix/str.c

juce_SOURCES = \
	../JuceLibraryCode/modules/juce_audio_basics/juce_audio_basics.cpp \
	../JuceLibraryCode/modules/juce_audio_devices/juce_audio_devices.cpp \
	../JuceLibraryCode/modules/juce_audio_formats/juce_audio_formats.cpp \
	../JuceLibraryCode/modules/juce_audio_processors/juce_audio_processors.cpp \
	../JuceLibraryCode/modules/juce_core/juce_core.cpp \
	../JuceLibraryCode/modules/juce_cryptography/juce_cryptography.cpp \
	../JuceLibraryCode/modules/juce_data_structures/juce_data_structures.cpp \
	../JuceLibraryCode/modules/juce_events/juce_events.cpp \
	../JuceLibraryCode/modules/juce_graphics/juce_graphics.cpp \
	../JuceLibraryCode/modules/juce_gui_basics/juce_gui_basics.cpp \
	../JuceLibraryCode/modules/juce_gui_extra/juce_gui_extra.cpp \
	../JuceLibraryCode/modules/juce_opengl/juce_opengl.cpp \
	../JuceLibraryCode/modules/juce_video/juce_video.cpp

SBSMSView_SOURCES = \
	sms.cpp \
	track.cpp \
	trackpoint.cpp \
	resample.cpp \
	sbsms.cpp sbsms.h\
	fft.cpp \
	grain.cpp \
	buffer.cpp \
	dBTable.cpp \
	slide.cpp \
	dBTable.h \
	sincCoeffs.h \
	synthTable.cpp synthTable.h \
	sse.h \
	subband.cpp \
	buffer.h \
	fft.h \
	utils.h \
	trackpoint.h \
	grain.h \
	real.h \
	../include/sbsms.h \
	sms.h \
	subband.h \
	track.h \
	convert.cpp convert.h\
	pcm.cpp pcm.h \
	mp3.cpp mp3.h \
	ffmpeg.cpp ffmpeg.h \
	mp3tech.cpp mp3tech.h \
	import.cpp import.h \
	audiobuffer.cpp audiobuffer.h \
	Main.cpp \
	MainWindow.cpp MainWindow.h \
	TimeFrequencyView.cpp TimeFrequencyView.h \
	Audio.cpp Audio.h \
	$(gl_matrix_SOURCES) \
	$(juce_SOURCES)

all: $(top_builddir)/approot/SBSMSView.app $(top_builddir)/mac/Info.plist $(top_builddir)/mac/PkgInfo

$(top_builddir)/approot/SBSMSView.app: SBSMSView
	mkdir -p $(top_builddir)/approot/SBSMSView.app/Contents/MacOS
	mkdir -p $(top_builddir)/approot/SBSMSView.app/Contents/Resources
	cp -f SBSMSView $(top_builddir)/approot/SBSMSView.app/Contents/MacOS/
	cp -f $(top_builddir)/mac/Info.plist $(top_builddir)/approot/SBSMSView.app/Contents/
	cp -f $(top_builddir)/mac/PkgInfo $(top_builddir)/approot/SBSMSView.app/Contents/

