CMAKE_MINIMUM_REQUIRED(VERSION 3.4)

SET(EXECUTABLE_OUTPUT_PATH ../../cgal_test/bin)

## ----- change cgal to our result format -----

SET(APP_NAME EfficientRANSAC)

ADD_EXECUTABLE( ${APP_NAME}
        main_cgal2ourres.cpp
        )

TARGET_LINK_LIBRARIES( ${APP_NAME}
        ${CGAL_LIB}
        )

TARGET_COMPILE_OPTIONS( ${APP_NAME} PRIVATE -std=c++11 )


## ----- write cgal points -----

SET(APP_NAME WritePoints)

ADD_EXECUTABLE( ${APP_NAME}
        main_cgal_write_points.cpp
        )

TARGET_LINK_LIBRARIES( ${APP_NAME}
        ${CGAL_LIB}
        )

TARGET_COMPILE_OPTIONS( ${APP_NAME} PRIVATE -std=c++11 )
