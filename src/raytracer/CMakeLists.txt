add_executable(raytracer main.cpp  raytracer.cpp raytracer.hpp Photon.cpp Photon.hpp azReflection.cpp azReflection.hpp montecarlo.cpp montecarlo.hpp  Utils.h constants.h)
target_link_libraries(raytracer application math scene tinyxml ${SDL_LIBRARY}
                      ${PNG_LIBRARIES} ${OPENGL_LIBRARIES} ${GLUT_LIBRARIES}
                      ${GLEW_LIBRARIES})

if(APPLE)
    target_link_libraries(raytracer SDLmain)
endif()

install(TARGETS raytracer DESTINATION ${PROJECT_SOURCE_DIR}/..)
