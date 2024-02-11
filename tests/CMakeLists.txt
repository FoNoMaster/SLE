find_package(GTest REQUIRED)

enable_testing()

add_executable(run_test tests.cpp ../src/tridiagonal_solver.cpp)
target_link_libraries(run_test GTest::gtest_main)

add_test(Alltest run_test)
