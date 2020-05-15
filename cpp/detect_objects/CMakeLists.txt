CMAKE_MINIMUM_REQUIRED(VERSION 3.4)

SET(EXECUTABLE_OUTPUT_PATH ../../detect_objects/bin)

## ----- test object detector with different DetectorSettings -----

SET(APP_NAME DetectObjects)

ADD_EXECUTABLE( ${APP_NAME}
        main_voting_sample.cpp
        )

TARGET_LINK_LIBRARIES( ${APP_NAME}
        tinyply
        pcd_lib
        )
        
TARGET_COMPILE_OPTIONS( ${APP_NAME} PRIVATE -std=c++11 )

## ----- detect and refine planes with "optimal" settings -----

SET(APP_NAME PlaneDetectRefine)

ADD_EXECUTABLE( ${APP_NAME}
        main_plane_dr.cpp
        )

TARGET_LINK_LIBRARIES( ${APP_NAME}
        ${CERES_LIBRARIES}
        tinyply
        pcd_lib
        )
        
TARGET_COMPILE_OPTIONS( ${APP_NAME} PRIVATE -std=c++11 )

## ----- detect and refine spheres with "optimal" settings -----

SET(APP_NAME SphereDetectRefine)

ADD_EXECUTABLE( ${APP_NAME}
        main_sphere_dr.cpp
        )

TARGET_LINK_LIBRARIES( ${APP_NAME}
        ${CERES_LIBRARIES}
        tinyply
        pcd_lib
        )
        
TARGET_COMPILE_OPTIONS( ${APP_NAME} PRIVATE -std=c++11 )

## ----- detect and refine cylinders with "optimal" settings -----

SET(APP_NAME CylDetectRefine)

ADD_EXECUTABLE( ${APP_NAME}
        main_cylinder_dr.cpp
        )

TARGET_LINK_LIBRARIES( ${APP_NAME}
        ${CERES_LIBRARIES}
        tinyply
        pcd_lib
        )
        
TARGET_COMPILE_OPTIONS( ${APP_NAME} PRIVATE -std=c++11 )

## ----- detect and refine cones with "optimal" settings -----

SET(APP_NAME ConeDetectRefine)

ADD_EXECUTABLE( ${APP_NAME}
        main_cone_dr.cpp
        )

TARGET_LINK_LIBRARIES( ${APP_NAME}
        ${CERES_LIBRARIES}
        tinyply
        pcd_lib
        )
        
TARGET_COMPILE_OPTIONS( ${APP_NAME} PRIVATE -std=c++11 )

## ----- refine plane parameters -----

SET(APP_NAME OptPlane)

ADD_EXECUTABLE( ${APP_NAME}
        main_refine_planes.cpp
        )

TARGET_LINK_LIBRARIES( ${APP_NAME}
        ${CERES_LIBRARIES}
        tinyply
        pcd_lib
        )
        
TARGET_COMPILE_OPTIONS( ${APP_NAME} PRIVATE -std=c++11 )

## ----- refine sphere parameters -----

SET(APP_NAME OptSphere)

ADD_EXECUTABLE( ${APP_NAME}
        main_refine_spheres.cpp
        )

TARGET_LINK_LIBRARIES( ${APP_NAME}
        ${CERES_LIBRARIES}
        tinyply
        pcd_lib
        )
        
TARGET_COMPILE_OPTIONS( ${APP_NAME} PRIVATE -std=c++11 )

## ----- refine cylinder parameters -----

SET(APP_NAME OptCylinder)

ADD_EXECUTABLE( ${APP_NAME}
        main_refine_cylinders.cpp
        )

TARGET_LINK_LIBRARIES( ${APP_NAME}
        ${CERES_LIBRARIES}
        tinyply
        pcd_lib
        )
        
TARGET_COMPILE_OPTIONS( ${APP_NAME} PRIVATE -std=c++11 )
