CC = g++
CFLAGS  = -g -Wall `root-config  --libs --cflags`
TARGET = CCM
CCM_LIB = CheckCCM.cpp CCM.cpp Cross_correlation.cpp
USER_SOURCE = main_example.cpp

all: CCM

CCM: $(USER_SOURCE)
	$(CC) $(CFLAGS) $(CCM_LIB) $(USER_SOURCE) -o $(TARGET)

clean:
	$(RM) $(TARGET)

