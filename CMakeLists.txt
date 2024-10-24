cmake_minimum_required(VERSION 3.16)
project(MyProject VERSION 0.1 LANGUAGES CXX)

set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_BUILD_TYPE Release)

find_package(QT NAMES Qt6 Qt5 REQUIRED COMPONENTS Widgets)
find_package(Qt${QT_VERSION_MAJOR} REQUIRED COMPONENTS Widgets)

# 原有的可执行文件
add_executable(Matrix Matrix.cpp)
add_executable(test_LpVariable ./test/test_LpVariable.cpp LpVariable.cpp)
add_executable(test_LpProblem ./test/test_LpProblem.cpp LpProblem.cpp LpVariable.cpp)
add_executable(test_SimplexAlgorithm ./test/test_SimplexAlgorithm.cpp SimplexAlgorithm.cpp LpProblem.cpp LpVariable.cpp)
add_executable(LP1 ./LP_test/LP1.cpp SimplexAlgorithm.cpp LpProblem.cpp LpVariable.cpp)
add_executable(LP2 ./LP_test/LP2.cpp SimplexAlgorithm.cpp LpProblem.cpp LpVariable.cpp)
add_executable(LP3 ./LP_test/LP3.cpp SimplexAlgorithm.cpp LpProblem.cpp LpVariable.cpp)
add_executable(LP4 ./LP_test/LP4.cpp SimplexAlgorithm.cpp LpProblem.cpp LpVariable.cpp)
add_executable(LP5 ./LP_test/LP5.cpp SimplexAlgorithm.cpp LpProblem.cpp LpVariable.cpp)

# Qt GUI 部分
set(PROJECT_SOURCES
        main.cpp
        mainwindow.cpp
        mainwindow.h
        mainwindow.ui
        SimplexAlgorithm.cpp
        LpProblem.cpp
        LpVariable.cpp
)

if(${QT_VERSION_MAJOR} GREATER_EQUAL 6)
    qt_add_executable(main
            MANUAL_FINALIZATION
            ${PROJECT_SOURCES}
    )
else()
    add_executable(main ${PROJECT_SOURCES})
endif()

target_link_libraries(main PRIVATE Qt${QT_VERSION_MAJOR}::Widgets)

# 设置属性
set_target_properties(main PROPERTIES
        MACOSX_BUNDLE TRUE
        WIN32_EXECUTABLE TRUE
)

include(GNUInstallDirs)
install(TARGETS main
        BUNDLE DESTINATION .
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

if(QT_VERSION_MAJOR EQUAL 6)
    qt_finalize_executable(main)
endif()
