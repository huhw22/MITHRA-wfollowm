SHELL = /bin/sh
COMP = mpic++

CPPFLAGS +=
CFLAGS += -std=c++11
CFLAGS += -O3
LDFLAGS +=
LDLIBS +=

# 如果系统装了 pkg-config 和 hdf5，这两行会自动生效
CPPFLAGS += $(shell pkg-config --cflags hdf5 2>/dev/null)
LDLIBS   += $(shell pkg-config --libs hdf5 2>/dev/null)

PREFIX ?= .
LIB_DIR = $(PREFIX)/lib
INC_DIR = $(PREFIX)/include/mithra/

EXEC = prj/MITHRA
SRC_DIR = src
SRCS := $(wildcard $(SRC_DIR)/*.cpp)
HDRS := $(wildcard $(SRC_DIR)/*.h)
OBJ_DIR = obj
OBJS := $(SRCS:$(SRC_DIR)/%.cpp=$(OBJ_DIR)/%.o)
LIB = $(OBJ_DIR)/libmithra.a

.PHONY: all clean debug install

all: $(EXEC) $(LIB)

debug: CFLAGS += -g
debug: all

install: all install-incs install-lib

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	@mkdir -p $(@D)
	$(COMP) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

$(EXEC): $(OBJS)
	$(COMP) $^ -o $@ $(LDFLAGS) $(LDLIBS)

$(LIB): $(OBJS)
	@mkdir -p $(@D)
	$(AR) rcs $@ $^

install-incs:
	@mkdir -p $(INC_DIR)
	install $(HDRS) $(INC_DIR)

install-lib:
	@mkdir -p $(LIB_DIR)
	install $(LIB) $(LIB_DIR)

clean:
	$(RM) $(EXEC) $(OBJS) $(LIB)
	rmdir $(OBJ_DIR) 2>/dev/null || true