set(SOURCE_LIST
    BoundaryHeadEvent.cpp
    BoundaryTailEvent.cpp
    InjectionEvent.cpp
    ChangeInputMixtureEvent.cpp
    MergeBifurcationEvent.cpp
    MergeChannelEvent.cpp
    TimeStepEvent.cpp
    StopSimulationEvent.cpp
)

set(HEADER_LIST
    BoundaryHeadEvent.h
    BoundaryTailEvent.h
    Event.h
    InjectionEvent.h
    ChangeInputMixtureEvent.h
    MergeBifurcationEvent.h
    MergeChannelEvent.h
    TimeStepEvent.h
    StopSimulationEvent.h
)

target_sources(${TARGET_NAME} PRIVATE ${SOURCE_LIST} ${HEADER_LIST})
target_include_directories(${TARGET_NAME} PRIVATE ${CMAKE_CURRENT_LIST_DIR})