CC = mpicc.mpich
AR = ar rcs
CFLAGS = -O2 -Wall -Wextra
PREFIX = $(HOME)/.local
LIB_NAME = libvptcomm.a
SRC_DIR = src
INCLUDE_DIR = include
TEST_DIR = tests
BUILD_DIR = build

SRCS = $(wildcard $(SRC_DIR)/*.c)
HEADERS = $(wildcard $(INCLUDE_DIR)/*.h)
OBJS = $(patsubst $(SRC_DIR)/%.c, $(BUILD_DIR)/%.o, $(SRCS))
TESTS = $(patsubst $(TEST_DIR)/%.c, $(BUILD_DIR)/%, $(wildcard $(TEST_DIR)/*.c))

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.c | $(BUILD_DIR)
	$(CC) $(CFLAGS) -c $< -I./$(INCLUDE_DIR) -o $@

$(LIB_NAME): $(OBJS)
	$(AR) $@ $^

$(BUILD_DIR)/%: $(TEST_DIR)/%.c $(LIB_NAME)
	$(CC) $(CFLAGS) $< -o $@ -I./$(INCLUDE_DIR) -L. -lvptcomm -lm

all: $(LIB_NAME) $(TESTS)

test: all
	for test in $(TESTS); do mpirun.mpich -np 8 $$test; done

install: $(LIB_NAME)
	mkdir -p $(PREFIX)/lib $(PREFIX)/include/vptcomm
	cp $(LIB_NAME) $(PREFIX)/lib/
	cp $(HEADERS) $(PREFIX)/include/vptcomm/

uninstall:
	rm -f $(PREFIX)/lib/$(LIB_NAME)
	rm -rf $(PREFIX)/include/vptcomm

clean:
	rm -rf $(BUILD_DIR) $(LIB_NAME)

.PHONY: all test install uninstall clean

