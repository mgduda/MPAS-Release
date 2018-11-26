PWD=$(shell pwd)
EXE_NAME=dummy_model
NAMELIST_SUFFIX=dummy
override CPPFLAGS += -DCORE_DUMMY

report_builds:
	@echo "CORE=dummy"
