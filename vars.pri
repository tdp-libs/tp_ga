TARGET = tp_ga
TEMPLATE = lib

DEFINES += TP_GA_LIBRARY

SOURCES += src/Globals.cpp
HEADERS += inc/tp_ga/Globals.h

SOURCES += src/RefineArray.cpp
HEADERS += inc/tp_ga/RefineArray.h

SOURCES += src/RANSAC.cpp
HEADERS += inc/tp_ga/RANSAC.h

SOURCES += src/RANSACRefineArray.cpp
HEADERS += inc/tp_ga/RANSACRefineArray.h
