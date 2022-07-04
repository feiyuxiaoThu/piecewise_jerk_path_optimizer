cmake_minimum_required(VERSION 3.20)
project(QP C)

set(CMAKE_C_STANDARD 99)

include_directories(include)
set(SRC src/amd_1.c  src/amd_2.c src/amd_aat.c src/amd_control.c src/amd_defaults.c src/amd_dump.c src/amd_global.c src/amd_info.c
        src/amd_order.c src/amd_post_tree.c src/amd_postorder.c src/amd_preprocess.c src/amd_valid.c src/Auxilary.c src/ldl.c src/matrix_gen_utils.c
        src/qpSWIFT.c src/timer.c)

link_libraries(m) # add for math.h
#LINK_LIBRARIES(m)
add_executable(QP runqp.c ${SRC})
