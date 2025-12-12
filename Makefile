CONTAINER ?= $(container)

MAKEFILE_DIR := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))
build_directory := $(MAKEFILE_DIR)/build
install_directory := $(MAKEFILE_DIR)/install

.PHONY: all setupenv install build clean

all: setupenv install build

# Root prefixes
INSTALL_PREFIX := $(install_directory)
LOCAL_PREFIX   := $(install_directory)

INSTALL_BIN_DIR     := $(INSTALL_PREFIX)/bin
INSTALL_LIB_DIR     := $(INSTALL_PREFIX)/lib
INSTALL_INCLUDE_DIR := $(INSTALL_PREFIX)/include

LOCAL_BIN_DIR       := $(LOCAL_PREFIX)/bin
LOCAL_LIB_DIR       := $(LOCAL_PREFIX)/lib
LOCAL_INCLUDE_DIR   := $(LOCAL_PREFIX)/include


# (Optional) Source files – currently unused

BIN_SRCS      := $(wildcard $(INSTALL_BIN_DIR)/*)
LIB_SRCS      := $(wildcard $(INSTALL_LIB_DIR)/*)
INCLUDE_SRCS  := $(wildcard $(INSTALL_INCLUDE_DIR)/*)

BIN_LINKS     := $(patsubst $(INSTALL_BIN_DIR)/%,$(LOCAL_BIN_DIR)/%,$(BIN_SRCS))
LIB_LINKS     := $(patsubst $(INSTALL_LIB_DIR)/%,$(LOCAL_LIB_DIR)/%,$(LIB_SRCS))
INCLUDE_LINKS := $(patsubst $(INSTALL_INCLUDE_DIR)/%,$(LOCAL_INCLUDE_DIR)/%,$(INCLUDE_SRCS))


build:
	cmake . -B$(build_directory) -DCMAKE_INSTALL_PREFIX=$(install_directory)
	cmake --build $(build_directory)

install: build
	cmake --install $(build_directory)

setupenv: install
	@echo "Setting up environment under $(MAKEFILE_DIR)..."
	. ./end.sh

clean:
	rm -rf $(build_directory) $(install_directory) end.sh
