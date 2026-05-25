cmake_minimum_required(VERSION 3.15)

project(GEOstar)

set(CMAKE_CXX_STANDARD 20)

find_package(glfw3 CONFIG REQUIRED)
find_package(glad CONFIG REQUIRED)
find_package(glm CONFIG REQUIRED)

set(SHARED_SOURCES
        Tensor.cpp
        Solver.cpp
        Eval.cpp

        imgui/imgui.cpp
        imgui/imgui_draw.cpp
        imgui/imgui_tables.cpp
        imgui/imgui_widgets.cpp
        imgui/imgui_demo.cpp

        imgui/backends/imgui_impl_glfw.cpp
        imgui/backends/imgui_impl_opengl3.cpp

        implot/implot.cpp
        implot/implot_items.cpp
)

add_executable(GEOstarRK4
        mainRK4.cpp
        ${SHARED_SOURCES}
)

target_include_directories(GEOstarRK4 PRIVATE
        imgui
        imgui/backends
        implot
)

target_link_libraries(GEOstarRK4 PRIVATE
        glfw
        glad::glad
        glm::glm
        opengl32
)