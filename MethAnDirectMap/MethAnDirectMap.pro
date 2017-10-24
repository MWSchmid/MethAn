#-------------------------------------------------
#
# Project created by QtCreator 2014-08-18T13:01:52
#
#-------------------------------------------------

#QT       += core
QT       += core gui xml widgets
#QT       -= gui

TARGET = MethAnDirectMap
CONFIG   += console
CONFIG   -= app_bundle

TEMPLATE = app

HEADERS += ../../Rcount/source_V2/p502_SOURCE/dataStructure/databaseitem.h \
    ../../Rcount/source_V2/p502_SOURCE/dataStructure/databasereader.h \
    ../../Rcount/source_V2/p502_SOURCE/dataStructure/databasewriter.h \
    ../../Rcount/source_V2/p502_SOURCE/dataStructure/database.h \
    #../../Rcount/source_V2/p502_SOURCE/dataStructure/databasemodel.h \
    ../../Rcount/source_V2/p502_SOURCE/dataStructure/temporaryresults.h \
    ../../Rcount/source_V2/p502_SOURCE/dataStructure/mappingtreeitem.h \

SOURCES += ../../Rcount/source_V2/p502_SOURCE/dataStructure/databaseitem.cpp \
    ../../Rcount/source_V2/p502_SOURCE/dataStructure/databasereader.cpp \
    ../../Rcount/source_V2/p502_SOURCE/dataStructure/databasewriter.cpp \
    ../../Rcount/source_V2/p502_SOURCE/dataStructure/database.cpp \
    #../../Rcount/source_V2/p502_SOURCE/dataStructure/databasemodel.cpp \
    ../../Rcount/source_V2/p502_SOURCE/dataStructure/temporaryresults.cpp \
    main.cpp

INCLUDEPATH+= ../../Rcount/source_V2/zlib-1.2.8 \
              ../../Rcount/source_V2/

LIBS += -lz

QMAKE_CXXFLAGS += -D_GLIBCXX_USE_CXX11_ABI=0
