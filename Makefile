CONTAINER ?= $(container)

build_directory := build
install_directory := $(CURDIR)/install

.PHONY: build install all clean install_extra link_dirs

ifeq ($(CONTAINER),podman)
# Root prefixes
INSTALL_PREFIX := $(install_directory)
LOCAL_PREFIX   := /ratpac-setup/local

INSTALL_BIN_DIR     := $(INSTALL_PREFIX)/bin
INSTALL_LIB_DIR     := $(INSTALL_PREFIX)/lib
INSTALL_INCLUDE_DIR := $(INSTALL_PREFIX)/include

LOCAL_BIN_DIR       := $(LOCAL_PREFIX)/bin
LOCAL_LIB_DIR       := $(LOCAL_PREFIX)/lib
LOCAL_INCLUDE_DIR   := $(LOCAL_PREFIX)/include

# (Optional) Source files � currently unused
BIN_SRCS      := $(wildcard $(INSTALL_BIN_DIR)/*)
LIB_SRCS      := $(wildcard $(INSTALL_LIB_DIR)/*)
INCLUDE_SRCS  := $(wildcard $(INSTALL_INCLUDE_DIR)/*)

BIN_LINKS     := $(patsubst $(INSTALL_BIN_DIR)/%,$(LOCAL_BIN_DIR)/%,$(BIN_SRCS))
LIB_LINKS     := $(patsubst $(INSTALL_LIB_DIR)/%,$(LOCAL_LIB_DIR)/%,$(LIB_SRCS))
INCLUDE_LINKS := $(patsubst $(INSTALL_INCLUDE_DIR)/%,$(LOCAL_INCLUDE_DIR)/%,$(INCLUDE_SRCS))

install_extra: install link_dirs
	ln -sfn $(INSTALL_PREFIX)/share/end $(LOCAL_PREFIX)/share/end

link_dirs: $(BIN_LINKS) $(LIB_LINKS) $(INCLUDE_LINKS)

# Generic rule: link from install ? local
$(LOCAL_PREFIX)/%: $(INSTALL_PREFIX)/%
	ln -sfn $(abspath $<) $@
else
install_extra:
	@true
endif

build:
	cmake . -B$(build_directory) -DCMAKE_INSTALL_PREFIX=$(install_directory)
	cmake --build $(build_directory)

install: build
	cmake --install $(build_directory)

all: build install install_extra
	ln -sf ../../cformat.sh ./.git/hooks/pre-commit

clean:
	rm -rf $(build_directory) $(install_directory) end.sh

