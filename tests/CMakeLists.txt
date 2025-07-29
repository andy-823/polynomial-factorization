Include(FetchContent)

FetchContent_Declare(
  Catch2
  GIT_REPOSITORY https://github.com/catchorg/Catch2.git
  GIT_TAG        v3.5.2
)

FetchContent_MakeAvailable(Catch2)

add_executable(test_galois_field test_galois_field.cpp)
target_link_libraries(test_galois_field PRIVATE factorization)
target_link_libraries(test_galois_field PRIVATE Catch2::Catch2WithMain)

add_executable(test_galois_field_element test_galois_field_element.cpp)
target_link_libraries(test_galois_field_element PRIVATE factorization)
target_link_libraries(test_galois_field_element PRIVATE Catch2::Catch2WithMain)

add_executable(test_polynomial test_polynomial.cpp)
target_link_libraries(test_polynomial PRIVATE factorization)
target_link_libraries(test_polynomial PRIVATE Catch2::Catch2WithMain)

add_executable(test_solver test_solver.cpp)
target_link_libraries(test_solver PRIVATE factorization)
target_link_libraries(test_solver PRIVATE Catch2::Catch2WithMain)