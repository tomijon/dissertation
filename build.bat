@echo off

if not exist build/ mkdir build/

if not exist build/shaders/ mkdir build/shaders/

xcopy "src/shaders" "build/shaders" /E /I /Y

g++ -std=c++11 -o build/Application.exe src/*.cpp src/*.c -I./include -L./lib -lglfw3 -lopengl32 -lgdi32 
