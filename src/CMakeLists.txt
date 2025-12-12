add_library( END SHARED
  End.cc
  GeoEndFactory.cc
  #DichroiconArrayFactory.cc
  HitmanProc.cc
  NtupleProc.cc
  LaserballGenerator.cc
  MyDAQProc.cc
  ForcedTriggerProc.cc
)

target_include_directories( END PUBLIC
  $<BUILD_INTERFACE:${END_INCLUDE_DIR}>
  $<INSTALL_INTERFACE:include>)

file(COPY include/ DESTINATION ${END_INCLUDE_DIR})

add_executable(end end.cpp)
target_link_libraries(end PUBLIC 
  END ${ROOT_LIBRARIES} ${Geant4_LIBRARIES} ${Ratpac_LIBRARIES} ${NLOPT_LIBRARIES})

set_target_properties(END end
  PROPERTIES
  ARCHIVE_OUTPUT_DIRECTORY ${END_LIBRARY_DIR}
  LIBRARY_OUTPUT_DIRECTORY ${END_LIBRARY_DIR})

install(TARGETS END end 
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  RUNTIME DESTINATION bin)

install(DIRECTORY ${END_INCLUDE_DIR}
  DESTINATION ${CMAKE_INSTALL_PREFIX})
