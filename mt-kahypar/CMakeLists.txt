add_subdirectory(datastructures)
add_subdirectory(io)
add_subdirectory(partition)
add_subdirectory(utils)

foreach(modtarget IN LISTS TARGETS_WANTING_ALL_SOURCES)
	target_precompile_headers(${modtarget} PUBLIC definitions.h)
	target_precompile_headers(${modtarget} PUBLIC partition/context.h)
endforeach()
