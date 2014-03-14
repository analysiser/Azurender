add_executable(p3 main.cpp raytracer.cpp Photon.cpp)
target_link_libraries(p3 application math scene tinyxml ${SDL_LIBRARY}
                      ${PNG_LIBRARIES} ${OPENGL_LIBRARIES} ${GLUT_LIBRARIES}
                      ${GLEW_LIBRARIES})

if(APPLE)
    target_link_libraries(p3 SDLmain)
endif()

install(TARGETS p3 DESTINATION ${PROJECT_SOURCE_DIR}/..)
