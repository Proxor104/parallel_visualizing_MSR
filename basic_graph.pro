QT += widgets

HEADERS       = headers.h
SOURCES       = main.cpp \
                classes.cpp\
                in_out.cpp\
                proc_time.cpp\
                window.cpp\
                calculations.cpp
                

QMAKE_CXXFLAGS += -std=c++0x -pthread
QMAKE_CXXFLAGS += -O3 -mfpmath=sse -fstack-protector-all -g -W -Wall -Wextra -Wunused -Wcast-align -pedantic -pedantic-errors -Wfloat-equal -Wpointer-arith -Wformat-security -Wmissing-format-attribute -Wformat=1 -Wwrite-strings -Wcast-align -Wno-long-long -Woverloaded-virtual -Wnon-virtual-dtor -Wcast-qual -Wno-suggest-attribute=format

LIBS += -pthread
