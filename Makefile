FLAGS = -Wall -g
TARGET = ./bin/main.exe
OBJECTS := $(patsubst ./src/%.cpp, ./bin/%.o, $(wildcard ./src/*.cpp))

all : $(TARGET)

$(TARGET) : $(OBJECTS)
	g++ $(FLAGS) -o $(TARGET) $(OBJECTS)

./bin/%.o : ./src/%.cpp
	g++ $(FLAGS) -c $^ -o $@

clean :
	rm -f $(OBJECTS) $(TARGET)

run : $(TARGET)
	$(TARGET)