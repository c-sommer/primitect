CMAKE_MINIMUM_REQUIRED(VERSION 3.4)

SET(EXECUTABLE_OUTPUT_PATH ../../track_object/bin)


## ----- track sphere in video sequence (our method) -----

SET(APP_NAME TrackSphere)

ADD_EXECUTABLE( ${APP_NAME}
        main_track_sphere.cpp
        )

TARGET_LINK_LIBRARIES( ${APP_NAME}
        tinyply
        pcd_lib
        )
        
TARGET_COMPILE_OPTIONS( ${APP_NAME} PRIVATE -std=c++11 )

# ## ----- track sphere in video sequence (Efficient RANSAC))-----

# SET(APP_NAME TrackSphereRansac)

# ADD_EXECUTABLE( ${APP_NAME}
#         main_track_sphere_ransac.cpp
#         )

# TARGET_LINK_LIBRARIES( ${APP_NAME}
#         ${CGAL_LIB}
#         )
        
# TARGET_COMPILE_OPTIONS( ${APP_NAME} PRIVATE -std=c++11 )


## ----- track cylinder in video sequence (our method) -----

SET(APP_NAME TrackCylinder)

ADD_EXECUTABLE( ${APP_NAME}
        main_track_cylinder.cpp
        )

TARGET_LINK_LIBRARIES( ${APP_NAME}
        tinyply
        pcd_lib
        )
        
TARGET_COMPILE_OPTIONS( ${APP_NAME} PRIVATE -std=c++11 )

# ## ----- track cylinder in video sequence (Efficient RANSAC))-----

# SET(APP_NAME TrackCylinderRansac)

# ADD_EXECUTABLE( ${APP_NAME}
#         main_track_cylinder_ransac.cpp
#         )

# TARGET_LINK_LIBRARIES( ${APP_NAME}
#         ${CGAL_LIB}
#         )
        
# TARGET_COMPILE_OPTIONS( ${APP_NAME} PRIVATE -std=c++11 )


## ----- track cone in video sequence (our method) -----

SET(APP_NAME TrackCone)

ADD_EXECUTABLE( ${APP_NAME}
        main_track_cone.cpp
        )

TARGET_LINK_LIBRARIES( ${APP_NAME}
        tinyply
        pcd_lib
        )
        
TARGET_COMPILE_OPTIONS( ${APP_NAME} PRIVATE -std=c++11 )

# ## ----- track cone in video sequence (Efficient RANSAC))-----

# SET(APP_NAME TrackConeRansac)

# ADD_EXECUTABLE( ${APP_NAME}
#         main_track_cone_ransac.cpp
#         )

# TARGET_LINK_LIBRARIES( ${APP_NAME}
#         ${CGAL_LIB}
#         )
        
# TARGET_COMPILE_OPTIONS( ${APP_NAME} PRIVATE -std=c++11 )

